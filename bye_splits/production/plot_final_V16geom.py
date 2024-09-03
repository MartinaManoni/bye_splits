import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import json
import os

def plot_geojson_figures(geojson_file, eta_type):
    with open(geojson_file, 'r') as file:
        geojson_data = json.load(file)

    if eta_type == 'GT0':
        geojson_data = invert_coordinates(geojson_data)

    layers = {}

    # Group features by layer
    for feature in geojson_data['features']:
        layer = feature['properties']['Layer']
        if layer not in layers:
            layers[layer] = []
        layers[layer].append(feature)

    # Create directory to save plots
    output_dir = f'plots_geojson_ETA_{eta_type}'
    os.makedirs(output_dir, exist_ok=True)

    # Plot each layer separately
    for layer, features in layers.items():
        plt.figure()

        for feature in features:
            vertices = feature['geometry']['coordinates'][0]
            U = feature['properties']['wu']
            V = feature['properties']['wv']

            polygon = Polygon(vertices)
            x, y = polygon.exterior.xy
            plt.plot(x, y, label=f'Layer {layer}', color='black')

            # Calculate centroid of the polygon
            centroid_x, centroid_y = polygon.centroid.xy
            centroid_x = centroid_x[0]
            centroid_y = centroid_y[0]

            plt.text(centroid_x, centroid_y, f'{U},{V}', fontsize=5, ha='center', va='center', color='red')  # Add U, V identifiers at centroid

        # Customize plot
        plt.xlabel('X coordinate [cm]')
        plt.ylabel('Y coordinate [cm]')
        plt.title(f'Layer {layer} - Silicon - Eta {eta_type}')
        plt.grid(True)
        plt.axis('equal')
        #plt.legend()

        # Save the plot
        plt.savefig(os.path.join(output_dir, f'layer_{layer}.png'), dpi=500)
        plt.close()

def invert_coordinates(geojson_data):
    inverted_features = []

    for feature in geojson_data['features']:
        inverted_feature = feature.copy()
        inverted_coordinates = []

        for coord_pair in feature['geometry']['coordinates'][0]:
            x, y = coord_pair
            inverted_coordinates.append([-x, y])
        
        inverted_feature['geometry']['coordinates'] = [inverted_coordinates]
        inverted_features.append(inverted_feature)

    inverted_geojson = geojson_data.copy()
    inverted_geojson['features'] = inverted_features

    return inverted_geojson

def main():
    geojson_file = 'hexagon_geometry_V16.json'
    
    # Ask the user which plot to generate: EtaLT0 or EtaGT0
    eta_type = input("Enter the plot type ('LT0' for EtaLT0, 'GT0' for EtaGT0): ").strip().upper()
    
    if eta_type not in ['LT0', 'GT0']:
        print("Invalid input. Please enter 'LT0' or 'GT0'.")
        return

    plot_geojson_figures(geojson_file, eta_type)

if __name__ == "__main__":
    main()
