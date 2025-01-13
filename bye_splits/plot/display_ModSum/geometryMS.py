_all_ = [ ]

import os
from pathlib import Path
import sys

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import json
import h5py
import re
import plotly.express.colors as px
import random
import logging
import time
import multiprocessing
from shapely.strtree import STRtree


log = logging.getLogger(__name__)

#from bye_splits.plot.display_plotly import yaml, np, pd, go, dbc
from utils import params, common 
from data_handle.data_process import *
from data_handle.geometry import GeometryData
import pandas as pd
import numpy as np
import yaml
import processingMS
import plotMS
import algosMS
import resolutionMS
import matplotlib
from scipy.spatial.distance import cdist
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon, Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as matPoly
from matplotlib.collections import PatchCollection
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

class Geometry():
    def __init__(self):
        self.ds_geom = GeometryData(reprocess=False, logger=log, library='plotly').provide()

    #TOWERS BINS 

    def create_arc(self, x_start, y_start, x_stop, y_stop, center):
        # Calculate radius based on the distance from the center to the starting point
        radius = np.linalg.norm([x_start - center[0], y_start - center[1]])

        # Calculate the start and end angles
        start_angle = np.arctan2(y_start - center[1], x_start - center[0])
        end_angle = np.arctan2(y_stop - center[1], x_stop - center[0])

        # Ensure that the end angle is greater than the start angle
        if end_angle < start_angle:
            print("end_angle < start_angle")
            end_angle += 2 * np.pi

        # Generate angles for the arc
        theta = np.linspace(start_angle, end_angle, 10)

        # Calculate points on the arc
        x = center[0] + radius * np.cos(theta)
        y = center[1] + radius * np.sin(theta)

        # Combine the starting point, arc, and end point into a single list of points
        arc_points = np.column_stack([x, y])
        arc_line = geom.LineString(arc_points)
        return arc_line

    def create_bin_polygon(self, bin_x, bin_y, bin_eta, bin_phi, center):
        center_x = np.mean(bin_x)
        center_y = np.mean(bin_y)
        arcs = []
        #lines = []

        #find index of points with same eta but different phi
        idx_pairs = []
        idx_pairs_no = []
        for i in range(4):
            for j in range(i + 1, 4):
                if bin_eta[i] == bin_eta[j] and bin_phi[i] != bin_phi[j]:
                    idx_pairs.append((i, j))
                if bin_eta[i] != bin_eta[j] and bin_phi[i] == bin_phi[j]:
                    idx_pairs_no.append((i, j))

        for start_idx, end_idx in idx_pairs:
            arc = self.create_arc(bin_x[start_idx], bin_y[start_idx], bin_x[end_idx], bin_y[end_idx], center)
            arcs.append(arc)

        #for start_idx, end_idx in idx_pairs_no:
           # straight_lines = geom.LineString([[bin_x[start_idx], bin_y[start_idx]], [bin_x[end_idx], bin_y[end_idx]]])
           # lines.append(straight_lines)  

        # Combine arcs and lines to form the polygon
        points = []
        for arc in arcs:
            points.extend(arc.coords[:])  

        # Arrange the points in clockwise order around the center
        # (Note the role reversal: the "y-coordinate" is the first function parameter, the "x-coordinate" is the second.)
        points_clockwise = sorted(points, key=lambda p: np.arctan2(p[1] - center_y, p[0] - center_x))

        # Remove duplicate points while preserving order
        unique_points = []
        prev_point = None
        for point in points_clockwise:
            if point != prev_point:
                unique_points.append(point)
                prev_point = point

        # Ensure that the polygon is closed
        if unique_points[0] != unique_points[-1]:
            unique_points.append(unique_points[0])

        polygon = geom.Polygon(unique_points)
        #self.plot_polygon(polygon)
        return polygon
    
    def plot_polygon(self,polygon):
        x, y = polygon.exterior.xy
        plt.plot(x, y)
        plt.fill(x, y, alpha=0.5)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Polygon')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.grid(True)
        plt.show()

    def plot_polygons(self, polygons):
        plt.figure()
        for polygon in polygons[:40]:
            x, y = polygon.exterior.xy
            plt.plot(x, y)
            plt.fill(x, y, alpha=0.5)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Polygons')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.grid(True)
        plt.show()

    def save_bin_geo(self, df_bins, output_file, output_file_vertex):
            print("saving geometry bins to geojson")
            # Lists to store GeoJSON features
            features_bin_poly_with_arcs = []
            features_vertices = []

            # Iterate over each bin in df_bins
            for _, row in df_bins.iterrows():
                layer = row['Layer']
                z = row['Z_value']
                center = np.array([0., 0.])
                
                # Iterate over bins in the current layer
                for bin_info in row['Bins']:
                    # Extract bin vertices
                    x_vertices = bin_info['X_Vertices']
                    y_vertices = bin_info['Y_Vertices']
                    eta_vertices = bin_info['Eta_Vertices']
                    phi_vertices = bin_info['Phi_Vertices']

                    eta_vertices_list = eta_vertices.tolist() if isinstance(eta_vertices, np.ndarray) else eta_vertices
                    phi_vertices_list = phi_vertices.tolist() if isinstance(phi_vertices, np.ndarray) else phi_vertices

                    # Create a tower bin with arcs using Shapely
                    bin_polygon = self.create_bin_polygon(x_vertices, y_vertices, eta_vertices, phi_vertices, center)

                    # Convert Shapely polygon to GeoJSON feature
                    feature_bin_poly_with_arcs = {
                        'type': 'Bin',
                        'geometry': mapping(bin_polygon),  # Convert polygon to GeoJSON geometry
                        'properties': {
                            'Layer': layer,
                            'Z_value':z,
                            'Eta_vertices': eta_vertices_list,
                            'Phi_vertices': phi_vertices_list,
                        }
                    }
                    features_bin_poly_with_arcs.append(feature_bin_poly_with_arcs)  # Add feature to the list


                    # Create a tower bin with with only the four vertices using Shapely
                    vertices_polygon = Polygon(zip(x_vertices, y_vertices))


                    # Convert Shapely polygon to GeoJSON feature with vertices as geometry
                    feature_vertices = {
                        'type': 'Feature',
                        'geometry': mapping(vertices_polygon),  # Convert polygon to GeoJSON geometry
                        'properties': {
                            'Layer': layer,
                            'Z_value': z,
                            'Eta_vertices': eta_vertices_list,
                            'Phi_vertices': phi_vertices_list,
                        }
                    }
                    features_vertices.append(feature_vertices)  # Add feature to the list

            # Create GeoJSON FeatureCollection
            feature_collection = {
                'type': 'FeatureCollection',
                'features': features_bin_poly_with_arcs
            }
            # Write GeoJSON data to file
            with open(output_file, 'w') as f:
                json.dump(feature_collection, f, indent=4)

            print("GeoJSON with bin poly with arcs saved to:", output_file)

            # Create GeoJSON FeatureCollection for polygon features
            feature_collection_vertices = {
                'type': 'FeatureCollection',
                'features': features_vertices
            }

            # Write GeoJSON data with bin_polygon geometry to file
            with open(output_file_vertex, 'w') as f_polygon:
                json.dump(feature_collection_vertices, f_polygon, indent=4)

            print("GeoJSON with bin poly, only vertices, saved to:", output_file_vertex)

    #V11 SILICON GEOMETRY
    def save_bin_hex(self, output_file):
        features = []  # List to store GeoJSON features
        existing_properties = set()

        # Three different x/y shifts for CE-E (subdet 1), CE-H (subdet 2) for even and odd layers
        diff_x_subdet1 = -1.387
        diff_y_subdet1 = -0.601

        diff_x_subdet2_even = -1.387
        diff_y_subdet2_even = -0.745

        diff_x_subdet2_odd = -1.387
        diff_y_subdet2_odd = -0.457

        # Iterate over each bin in self.ds_geom['si']
        for _, row in self.ds_geom['si'].iterrows():
            layer = row['ts_layer']
            z = row['z']
            wu = row['ts_wu']
            wv = row['ts_wv']

            hex_x, hex_y = row['hex_x'], row['hex_y']

            if layer <= 28:
                hex_x_shifted = [x + diff_x_subdet1 for x in hex_x]
                hex_y_shifted = [y + diff_y_subdet1 for y in hex_y]
            elif layer % 2 == 0:
                hex_x_shifted = [x + diff_x_subdet2_even for x in hex_x]
                hex_y_shifted = [y + diff_y_subdet2_even for y in hex_y]
            else:
                hex_x_shifted = [x + diff_x_subdet2_odd for x in hex_x]
                hex_y_shifted = [y + diff_y_subdet2_odd for y in hex_y]

            hex_polygon = Polygon([(x, y) for x, y in zip(hex_x_shifted, hex_y_shifted)])
            
            feature = {
                'type': 'Feature',
                'geometry': mapping(hex_polygon),  # Convert polygon to GeoJSON geometry
                'properties': {
                    'Layer': layer,
                    'z': z,
                    'wu': wu,
                    'wv': wv
                }
            }

            feature_properties = (layer, z, wu, wv)  # Tuple representing the properties of the feature
            if feature_properties not in existing_properties:  # Check if properties already exist
                existing_properties.add(feature_properties)  # Add properties to the set
                features.append(feature)  # Add feature to the list

        # Create GeoJSON FeatureCollection
        feature_collection = {
            'type': 'FeatureCollection',
            'features': features
        }
        
        print("saving hexagons geometry in geojson file")
        # Write GeoJSON data to file
        with open(output_file, 'w') as f:
            json.dump(feature_collection, f, indent=4)

    # SILICON V16 GEOMETRY
    def create_si_mod_geometry_V16(self):
        """
        Implementing the V16 geometry for the silicon part (hexagonal modules which include partials)
        """
        with open('/grid_mnt/vol_home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/silicon_geometry_V16_eta_GT0.json', 'r') as file:
            geojson_data = file.read()

        geojson = json.loads(geojson_data)

        # Initialize lists to hold the extracted data
        layers, waferus, wafervs, hex_xs, hex_ys = [], [], [], [], []

        # Extract data from GeoJSON features
        for feature in geojson['features']:
            properties = feature['properties']
            geometry = feature['geometry']['coordinates'][0]

            layers.append(properties['Layer'])
            waferus.append(properties['wu'])
            wafervs.append(properties['wv'])

            hex_x = [coord[0] for coord in geometry]
            hex_y = [coord[1] for coord in geometry]

            hex_xs.append(hex_x)
            hex_ys.append(hex_y)

        # Create a DataFrame from the extracted data
        df_geo = pd.DataFrame({
            'ts_layer': layers,
            'ts_wu': waferus,
            'ts_wv': wafervs,
            'hex_x': hex_xs,
            'hex_y': hex_ys
        })

        return df_geo



    # SCINTILLATOR GEOMETRY

    def implement_scint_modules_id(self, df):
        """
        Implementing the ieta (ts_ieta) and iphi (ts_iphi) identifier for scintillator modules based on CMSSW
        https://github.com/jbsauvan/cmssw/commit/a8a2c2c9a0cfe2745a6a523e2b3818e2124b80f0#diff-78500392da6bd14745d9ab56db0d0b182e62ba8dffae5616bd7afa71539b216dR980
        """
        hSc_tcs_per_module_phi = 4
        hSc_back_layers_split = 8
        hSc_front_layers_split = 12
        hSc_layer_for_split = 40

        df['ts_iphi'] = (df['tc_cv'] - 1) // hSc_tcs_per_module_phi + 1 #Phi index 1-12

        split = hSc_front_layers_split
        if df['tc_layer'].iloc[0] > hSc_layer_for_split:
            split = hSc_back_layers_split

        df['ts_ieta'] = df.apply(lambda row: 0 if row['tc_cu'] <= split else 1, axis=1)

        return df

    def order_vertices_clockwise(self, vertices):
        centroid = Polygon(vertices).centroid
        centroid_x, centroid_y = centroid.x, centroid.y
        # Convert vertices list to NumPy array
        vertices_array = np.array(vertices)
        # Calculate angles with respect to centroid
        angles = np.arctan2(vertices_array[:, 1] - centroid_y, vertices_array[:, 0] - centroid_x)
        # Sort vertices based on angles
        sorted_indices = np.argsort(angles)
        ordered_vertices = vertices_array[sorted_indices]

        return ordered_vertices.tolist()

    def find_corresponding_y(self, x, all_coords):
        closest_y = None
        min_distance = float('inf')
        for coord in all_coords:
            for x_val, y_val in coord:
                if x_val == x:
                    distance = abs(y_val)
                    if distance < min_distance:
                        min_distance = distance
                        closest_y = y_val
        return closest_y

    def find_corresponding_x(self, y, all_coords):
        closest_x = None
        min_distance = float('inf')
        for coord in all_coords:
            for x_val, y_val in coord:
                if y_val == y:
                    distance = abs(x_val)
                    if distance < min_distance:
                        min_distance = distance
                        closest_x = x_val
        return closest_x

    def save_scint_mod_geo(self, merged_df, output_file):
        print("Saving scintillator module geometries to GeoJSON")
        features_scint_mod_poly = []

        # Iterate over each row in the merged DataFrame
        for _, row in merged_df.iterrows():
            ts_ieta = row['ts_ieta']
            ts_iphi = row['ts_iphi']
            tc_layer = row['tc_layer']
            ordered_vertices_list = row['vertices_clockwise']
            #print(ordered_vertices_list)

            # Convert string representation of vertices to list of floats
            ordered_vertices = [list(map(float, vertex)) for vertex in ordered_vertices_list]

            # Create a Shapely Polygon from the ordered vertices
            polygon = Polygon(ordered_vertices)

            # Convert Shapely polygon to GeoJSON feature
            feature_scint_mod_poly = {
                    'type': 'Feature',
                    'geometry': mapping(polygon),  # Convert polygon to GeoJSON geometry
                    'properties': {
                        'Layer': tc_layer,
                        'ieta': ts_ieta,
                        'iphi': ts_iphi,
                    }
                }
            features_scint_mod_poly.append(feature_scint_mod_poly)  # Add feature to the list

        # Create GeoJSON FeatureCollection
        feature_collection = {
            'type': 'FeatureCollection',
            'features': features_scint_mod_poly
        }

        # Write GeoJSON data to file
        with open(output_file, 'w') as f:
            json.dump(feature_collection, f, indent=4)

        print("GeoJSON with scintillator module geometries saved to:", output_file)


    def create_scint_mod_geometry(self, df, save_geojson=False):
        """
        Function that creates the Modules Scintillator geometry and outputs a DataFrame
        with hex_x and hex_y coordinates representing the ordered vertices.

        If save_geojson = True, the geometry is saved into a geojson file.
        """
        # Group DataFrame by ts_ieta, ts_iphi, and tc_layer
        grouped = df.groupby(['ts_ieta', 'ts_iphi', 'tc_layer'])
        merged_geometries = []

        for (ts_ieta, ts_iphi, tc_layer), group in grouped:
            # Create a list to store individual diamond polygons
            diamond_polygons = []
            # Iterate over rows in the group and create Shapely Polygon objects
            for _, row in group.iterrows():
                # Define vertices of the diamond polygon
                vertices = [
                    (row['diamond_x'][0], row['diamond_y'][0]),
                    (row['diamond_x'][1], row['diamond_y'][1]),
                    (row['diamond_x'][2], row['diamond_y'][2]),
                    (row['diamond_x'][3], row['diamond_y'][3]),
                ]
                # Create a Polygon object
                polygon = Polygon(vertices)
                # Add polygon to the list
                diamond_polygons.append(polygon)

            # Extract coordinates of all sub-polygons directly from the diamond polygons
            all_coords = [list(polygon.exterior.coords) for polygon in diamond_polygons]

            # Flatten the list of coordinates
            flat_coords = [coord for sublist in all_coords for coord in sublist]

            # Find min and max x and y coordinates for all sub-polygons
            x_coords, y_coords = zip(*flat_coords)
            x_min, x_max = min(x_coords), max(x_coords)
            y_min, y_max = min(y_coords), max(y_coords)

            # Find the corresponding y-coordinate or x-coordinate
            Y_x_min = self.find_corresponding_y(x_min, all_coords)
            Y_x_max = self.find_corresponding_y(x_max, all_coords)
            X_y_min = self.find_corresponding_x(y_min, all_coords)
            X_y_max = self.find_corresponding_x(y_max, all_coords)

            # Order the vertices in a clockwise manner
            ordered_vertices = self.order_vertices_clockwise([(x_min, Y_x_min), (x_max, Y_x_max), (X_y_max, y_max), (X_y_min, y_min)])

            # Check for repeated vertices (if needed)
            if len(set(map(tuple, ordered_vertices))) < len(ordered_vertices):
                print(f"Repeated vertices found in Layer {tc_layer}: {ordered_vertices}")

            # Separate the ordered vertices into hex_x and hex_y
            # The 'hex' naming is mantained to be consistent with silicon geometry in order to merge events from both sicintillator and silicon geom,
            # even if here we are not conidering hexagons but quadrangles.
            hex_x, hex_y = zip(*ordered_vertices)

            # Append the result to the merged geometries list
            merged_geometries.append({
                'ts_ieta': ts_ieta,
                'ts_iphi': ts_iphi,
                'tc_layer': tc_layer,
                'hex_x': list(hex_x),  # Store ordered x-coordinates
                'hex_y': list(hex_y)   # Store ordered y-coordinates
            })

        # Create a DataFrame from the merged geometries
        merged_df = pd.DataFrame(merged_geometries)

        # Optionally save the geojson
        if save_geojson:
            self.save_scint_mod_geo(merged_df, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo_FINAL.geojson')

        return merged_df