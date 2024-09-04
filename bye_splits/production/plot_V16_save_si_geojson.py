import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from collections import defaultdict
import os
import numpy as np
import json

def rotate_point(x, y, angle):
    """Rotate a point (x, y) by angle degrees around the origin."""
    angle_rad = np.radians(angle)
    new_x = x * np.cos(angle_rad) - y * np.sin(angle_rad)
    new_y = x * np.sin(angle_rad) + y * np.cos(angle_rad)
    return new_x, new_y

def transform_U_V(U, V):
    """Transform U and V values as per the new rule."""
    new_U = -V
    new_V = U - V
    return new_U, new_V

def parse_and_plot_geometry(filename_geom, filename_rings):
    # Read the hexagon geometry file
    with open(filename_geom, 'r') as file:
        lines = file.readlines()
    
    # Group polygons by layer and cassette
    layers = defaultdict(lambda: defaultdict(list))
    
    # Define colors for each cassette
    colors = [
        'lightblue', 'darkturquoise', 'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'orange',
        'purple', 'brown', 'pink', 'gray'
    ]
    
    # Process each line in the geometry file
    for line in lines:
        parts = line.split()
        layer = int(parts[0])  # layer number
        center_x = float(parts[2]) / 10  # x coordinate of the center
        center_y = float(parts[3]) / 10  # y coordinate of the center
        N = int(parts[5])  # number of vertices
        cassette = parts[-9]  # cassette number
        U = int(parts[-2])  # U identifier (assuming U is an integer)
        V = int(parts[-1])  # V identifier (assuming V is an integer)
        
        # Extract vertices
        vertices = []
        for i in range(N):
            x = float(parts[6 + 2 * i]) / 10
            y = float(parts[6 + 2 * i + 1]) / 10
            
            # Rotate vertices for specified layers
            if layer in [28, 30, 32]:
                x, y = rotate_point(x, y, 30)
            
            vertices.append((x, y))
        
        # Rotate center for specified layers
        if layer in [28, 30, 32]:
            center_x, center_y = rotate_point(center_x, center_y, 30)
        
        # Store vertices and center info in the corresponding layer and cassette
        layers[layer][cassette].append((vertices, (center_x, center_y), U, V, N))
    
    # Read the rings file
    rings_data = defaultdict(list)
    with open(filename_rings, 'r') as file:
        lines = file.readlines()
    
    # Conditions for red lines
    red_line_conditions = {
        (34, 37): [25],
        (38, 43): [17, 25, 33],
        (44, 47): [5, 17, 25, 33]
    }

    # Process each line in the rings file, skipping headers and handling transitions
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            continue
        parts = line.split()
        layer = int(parts[0])
        second_column = int(parts[1])
        inner_radius = float(parts[10]) / 10  # Convert mm to cm
        outer_radius = float(parts[11]) / 10  # Convert mm to cm
        rings_data[layer].append((second_column, inner_radius, outer_radius))
    
    # Create a directory to save the plots
    output_dir = 'plots_geomV16_detailed_all'
    os.makedirs(output_dir, exist_ok=True)
    
    geojson_data = {
        "type": "FeatureCollection",
        "features": []
    }
    
    # Plot each layer separately
    for layer, cassettes in layers.items():
        plt.figure()
        
        # Plot circles for the layer
        if layer in rings_data:
            most_inner_radius = float('inf')
            most_outer_radius = float('-inf')
            for second_column, inner_radius, outer_radius in rings_data[layer]:
                color = 'grey'
                linewidth = 0.1
                for (start_layer, end_layer), columns in red_line_conditions.items():
                    if start_layer <= layer <= end_layer and second_column in columns:
                        color = 'red'
                        linewidth = 0.5
                        break
                circle_inner = plt.Circle((0, 0), inner_radius, color='grey', linewidth=0.1, fill=False)
                circle_outer = plt.Circle((0, 0), outer_radius, color=color, linewidth=linewidth, fill=False)
                plt.gca().add_patch(circle_inner)
                plt.gca().add_patch(circle_outer)
                if inner_radius < most_inner_radius:
                    most_inner_radius = inner_radius
                if outer_radius > most_outer_radius:
                    most_outer_radius = outer_radius
        
        # Draw lines every 10 degrees for layers >= 34
        if layer >= 34:
            for angle in range(0, 360, 10):
                angle_rad = np.radians(angle)
                x_inner = most_inner_radius * np.cos(angle_rad)
                y_inner = most_inner_radius * np.sin(angle_rad)
                x_outer = most_outer_radius * np.cos(angle_rad)
                y_outer = most_outer_radius * np.sin(angle_rad)
                plt.plot([x_inner, x_outer], [y_inner, y_outer], color='black', linestyle='-', linewidth=0.5)
        
        for cassette, figures in cassettes.items():
            color = colors[int(cassette) % len(colors)]
            for vertices, center, U, V, N in figures:
                polygon = Polygon(vertices)
                x, y = polygon.exterior.xy
                plt.plot(x, y, label=f'Layer {layer}', color=color)
                plt.plot(center[0], center[1], 'o', color=color, markersize=1)  # Plot the center
                plt.text(center[0], center[1], f'{U},{V}', fontsize=4, ha='center', color='black')  # Add U, V identifiers

                # Add feature to geojson data
                feature = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [vertices + [vertices[0]]]
                    },
                    "properties": {
                        "Layer": layer,
                        "N": N,
                        "wu": U,
                        "wv": V
                    }
                }
                geojson_data["features"].append(feature)

                # First rotation by 120 degrees counterclockwise and transform U and V values
                rotated_vertices_1 = [rotate_point(vx, vy, 120) for vx, vy in vertices]
                rotated_center_1 = rotate_point(center[0], center[1], 120)
                new_U_1, new_V_1 = transform_U_V(U, V)

                if layer >= 34 and layer % 2 == 0:
                    new_U_1 = new_U_1 + 1
                    new_V_1 = new_V_1 + 1 

                if layer >= 35 and layer % 2 == 1:
                    new_U_1 = new_U_1 - 1
                    new_V_1 = new_V_1 - 1 

                else:
                    new_U_1 = new_U_1
                    new_V_1 = new_V_1 
                
                # Plot the first rotated module
                rotated_polygon_1 = Polygon(rotated_vertices_1)
                rx1, ry1 = rotated_polygon_1.exterior.xy
                plt.plot(rx1, ry1, label=f'Layer {layer} (Rotated 1)', color=color)
                plt.plot(rotated_center_1[0], rotated_center_1[1], 'o', color=color, markersize=1)  # Plot the rotated center
                plt.text(rotated_center_1[0], rotated_center_1[1], f'{new_U_1},{new_V_1}', fontsize=4, ha='center', color='black')  # Add transformed U, V identifiers

                # Add feature to geojson data
                feature_rotated_1 = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [rotated_vertices_1 + [rotated_vertices_1[0]]]
                    },
                    "properties": {
                        "Layer": layer,
                        "N": N,
                        "wu": new_U_1,
                        "wv": new_V_1
                    }
                }
                geojson_data["features"].append(feature_rotated_1)

                # Second rotation by another 120 degrees counterclockwise and transform U and V values
                rotated_vertices_2 = [rotate_point(vx, vy, 120) for vx, vy in rotated_vertices_1]
                rotated_center_2 = rotate_point(rotated_center_1[0], rotated_center_1[1], 120)
                new_U_2, new_V_2 = transform_U_V(new_U_1, new_V_1)

                if layer >= 34 and layer % 2 == 0:
                    new_U_2 = new_U_2+1
                    new_V_2 = new_V_2+1

                if layer >= 35 and layer % 2 == 1:
                    new_U_2 = new_U_2-1
                    new_V_2 = new_V_2-1

                else:
                    new_U_2 = new_U_2
                    new_V_2 = new_V_2
            
                # Plot the second rotated module
                rotated_polygon_2 = Polygon(rotated_vertices_2)
                rx2, ry2 = rotated_polygon_2.exterior.xy
                plt.plot(rx2, ry2, label=f'Layer {layer} (Rotated 2)', color=color)
                plt.plot(rotated_center_2[0], rotated_center_2[1], 'o', color=color, markersize=1)  # Plot the rotated center
                plt.text(rotated_center_2[0], rotated_center_2[1], f'{new_U_2},{new_V_2}', fontsize=4, ha='center', color='black')  # Add transformed U, V identifiers

                # Add feature to geojson data
                feature_rotated_2 = {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [rotated_vertices_2 + [rotated_vertices_2[0]]]
                    },
                    "properties": {
                        "Layer": layer,
                        "N": N,
                        "wu": new_U_2,
                        "wv": new_V_2
                    }
                }
                geojson_data["features"].append(feature_rotated_2)
        
        # Customize plot
        plt.xlabel('X coordinate [cm]')
        plt.ylabel('Y coordinate [cm]')
        plt.title(f'Complete Geom- Layer {layer}')
        plt.grid(True)
        plt.axis('equal')
        
        # Save the plot
        #plt.savefig(os.path.join(output_dir, f'layer_{layer}.png'), dpi=500)
        plt.close()

    # Save the geojson data to a file
    with open('hexagon_geometry_V16.json', 'w') as json_file:
        json.dump(geojson_data, json_file, indent=4)

def main():
    filename_geom = 'geomdetailed10052021.txt'
    filename_rings = 'tiles_posts_pattern_spaces-scenario13.txt'
    parse_and_plot_geometry(filename_geom, filename_rings)

if __name__ == "__main__":
    main()