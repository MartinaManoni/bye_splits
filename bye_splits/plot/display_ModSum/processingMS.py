_all_ = [ ]

import os
from pathlib import Path
import sys

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

#from dash import dcc, html
import json
#import dash_daq as daq
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
import matplotlib
from scipy.spatial.distance import cdist
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon, Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as matPoly
from matplotlib.collections import PatchCollection
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
#matplotlib.use("TkAgg")

class Processing():
    def __init__(self):
        with open(params.CfgPath, "r") as afile:
            self.cfg = yaml.safe_load(afile)

        self.ds_geom = GeometryData(reprocess=False, logger=log, library='plotly').provide()
        self.filename = None
        self.list_events = None

    def random_event(self, f):
        return random.choice(self.list_events)

    def filestore(self):
        self.filename = params.LocalStorage + "/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"
        #"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/new_OKAY/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"

    def get_data_new(self, event=None):
        # Check if the HDF5 file exists
        self.filestore()
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"HDF5 file '{self.filename}' not found.")
        # Read the list of events from the HDF5 file
        with h5py.File(self.filename, 'r') as file:
            self.list_events = list(set([key.split('_')[1] for key in file.keys()]))
        print(f"List of events in HDF5 file: {self.list_events}")   
        print(f"Number of events in HDF5 file:",  len(self.list_events))  

        # If no event is specified, choose a random one
        if event is None:
            event = random.choice(self.list_events)
            print(f"Selected random event: {event}")
            processed_event_df = self.get_event(event)

        elif event == '-1':
            events_to_process = self.list_events
            print(f"Processing all events: {events_to_process}")
            processed_event_df = self.get_allevents(events_to_process)

        else:
            print(f"Selected event: {event}")
            # If the specified event is not in the list, raise an exception
            if str(event) not in self.list_events:
                raise ValueError(f"Event {event} not found in the HDF5 file.")
            processed_event_df = self.get_event(event)
        return processed_event_df

    def get_event(self, event):
        with h5py.File(self.filename, 'r') as file:
            for key in file.keys():
                if event in key and 'ts' in key:
                    dataset = file[key]
                    column_names = [str(col) for col in dataset.attrs['columns']]
                    # Convert the dataset to a Pandas DataFrame with proper column names
                    df_ts = pd.DataFrame(dataset[:], columns=column_names)
                    df_ts['event'] = int(event)
                    silicon_df_proc= self.process_event(df_ts)

                    return silicon_df_proc   
    
    def get_allevents(self, events_to_process):
        all_ts_event_dfs = []
        all_tc_event_dfs = []

        with h5py.File(self.filename, 'r') as file:
            for event in events_to_process:
                for key in file.keys():
                    if event in key:
                        if 'ts' in key:
                            dataset = file[key]
                            column_names = [str(col) for col in dataset.attrs['columns']]
                            df_ts = pd.DataFrame(dataset[:], columns=column_names)
                            df_ts['event'] = int(event)  # Add a new column for the event number
                            all_ts_event_dfs.append(df_ts)  # Append each ts event data frame to the list

                        #FIXME - Temporary processing tc data
                        elif 'tc' in key:
                            dataset = file[key]
                            column_names = [str(col) for col in dataset.attrs['columns']]
                            df_tc = pd.DataFrame(dataset[:], columns=column_names)
                            df_tc['event'] = int(event)  # Add a new column for the event number
                            all_tc_event_dfs.append(df_tc)  # Append each tc event data frame to the list
        
        combined_ts_df = pd.concat(all_ts_event_dfs, ignore_index=True)
        combined_tc_df = pd.concat(all_tc_event_dfs, ignore_index=True)
        print("tc data dataframe", combined_tc_df.columns)
        
        # Process the combined ts DataFrame
        processed_combined_ts_df = self.process_event(combined_ts_df) #FIXME - Temporary adding tc data (combined_ts_df,combined_tc_df)
        return processed_combined_ts_df
    
    def read_hdf5_structure(self,file_path, group_name='df', indent=0):

        def _print_group_structure(group, indent):
            for name in group:
                print(" " * indent + name)
                if isinstance(group[name], h5py.Group):
                    _print_group_structure(group[name], indent + 2)

        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                _print_group_structure(group, indent)

    def read_all_block0_values(self,file_path, group_name='df'):
        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                if 'block0_values' in group:
                    block0_values = group['block0_values']
                    print("Values in 'block0_values':")
                    for row in block0_values:
                        print(row)
                else:
                    print("Error: 'block0_values' not found in group {}.".format(group_name))
            else:
                print("Error: Group '{}' not found in HDF5 file.".format(group_name))
    
    def get_genpart_data(self, file_path, block0_value, group_name='df'):
        block1_data = []
        events = []
        print(block0_value)

        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                if 'block1_items' in group and 'block1_values' in group:
                    block0_values = group['block0_values'][:, 0]
                    print(" block0_values - events", block0_values)
                    block1_items = [item.decode() for item in group['block1_items']]
                    block1_values = group['block1_values']
                    if block0_value == '-1':
                        for idx, value in enumerate(block0_values):
                            block1_data.append(block1_values[idx])
                            events.append(value)
                    elif int(block0_value) in block0_values:
                        idx = block0_values.tolist().index(int(block0_value))
                        block1_data.append(block1_values[idx])
                        events.append(block0_value)
                    else:
                        print(f"Block0 value '{block0_value}' not found.")

                    # Create DataFrame and include event numbers as a column
                    df = pd.DataFrame(block1_data, columns=block1_items)
                    df['event'] = events

                    return df
                else:
                    print("Error: 'block0_values', 'block1_items', or 'block1_values' not found in group {}.".format(group_name))
            else:
                print("Error: Group '{}' not found in HDF5 file.".format(group_name))


    def store_event(self, path, data):
        if isinstance(data, dict):
            for key, item in data.items():
                if isinstance(item, pd.DataFrame):
                    item.to_hdf(self.filename, path + str(key))
                else:
                    raise ValueError('Cannot save %s type'%type(item))
        else:
            data.to_hdf(self.filename, path)

    def process_event(self, df_ts):
        #(self, df_ts, df_tc)
        print("Process events ...")
        print("Filled data ts columns", df_ts.columns)
        #print("DATA tf ts pt", df_ts['ts_pt'].sum())
        print("DATA tf ts pt", df_ts['ts_pt'].sum())

        ts_keep = {'waferu'       : 'ts_wu',
                   'waferv'       : 'ts_wv',
                   'layer'        : 'ts_layer'}
        
        sci_update = {'triggercellieta' : 'tc_cu',
                      'triggercelliphi' : 'tc_cv',
                      'layer'           : 'tc_layer',
                      'waferu'          : 'tc_wu',
                      'waferv'          : 'tc_wv'}
        
        self.ds_geom['si']  = self.ds_geom['si'].rename(columns=ts_keep)
        self.ds_geom['sci'] = self.ds_geom['sci'].rename(columns=sci_update)

        #SILICON
        print("GEOMETRY silicon columns", self.ds_geom['si'].columns)
        #print("GEO scint columns", self.ds_geom['sci'].columns)
        #mask= self.ds_geom['sci']['tc_layer']== 42
        #print("GEO scint columns", self.ds_geom['sci'][mask])
        
        ds_si_geo = self.ds_geom['si'].drop_duplicates(subset=['ts_layer', 'ts_wu', 'ts_wv'], keep='first')

        #plotMS.plot_layers(df_ts, ds_new)
        #plotMS.plot_layers_sci(df_tc, self.ds_geom['sci'])

        # Merge between data and geometry SILICON
        silicon_df = pd.merge(left=df_ts, right=ds_si_geo, how='inner',
                                      on=['ts_layer', 'ts_wu', 'ts_wv'])
        silicon_df = silicon_df.drop(['triggercellu','triggercellv','waferorient', 'waferpart','diamond_x', 'diamond_y'], axis=1)

        # Shifting hexagons vertices based on difference between wx_center/wy_center (byesplit) and ts_x/ts_y (CMSSW)
        shifted_hex_df = self.shift_hex_values(silicon_df, self.ds_geom['si'], df_ts)

        #SCINTILLATOR
        #Adding scintillator modules to scintillator geometry - reproducing what is done at the moment in CMSSW:
        #https://github.com/hgc-tpg/cmssw/blob/hgc-tpg-devel-CMSSW_14_0_0_pre1/L1Trigger/L1THGCal/plugins/geometries/HGCalTriggerGeometryV9Imp3.cc#L989

        self.add_scint_modules_var(self.ds_geom['sci'])
        #scint_mod_geom = self.create_scint_mod_geometry(self.ds_geom['sci'])
        #plotMS.plot_scint_tiles(self.ds_geom['sci'])
        #print("Scintillator geometry", scint_mod_geom.columns)

        '''scintillator_df = pd.merge(left=df_tc, right=self.ds_geom['sci'], how='inner', #FIXME
                                           on=['tc_layer', 'tc_wu', 'tc_cu', 'tc_cv'])'''

        #print("shifted_hex_df")
        #print(shifted_hex_df[['hex_x', 'hex_y']])
        #print("\nsilicon_df")
        #print(silicon_df[['hex_x', 'hex_y']])

        return shifted_hex_df

    def add_scint_modules_var(self, df):
        #Adding to the dataframe ieta (ts_ieta) and iphi (ts_iphi) identifier for scintillator modules
        hSc_tcs_per_module_phi = 4
        hSc_back_layers_split = 8
        hSc_front_layers_split = 12
        hSc_layer_for_split = 40

        df['ts_iphi'] = (df['tc_cv'] - 1) // hSc_tcs_per_module_phi

        split = hSc_front_layers_split
        if df['tc_layer'].iloc[0] > hSc_layer_for_split:
            split = hSc_back_layers_split

        df['ts_ieta'] = df.apply(lambda row: 0 if row['tc_cu'] <= split else 1, axis=1)

    def create_scint_mod_geometry(self, df):
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

            # Find the corresponding y-coordinate or x-coordinate #FIXME
            Y_x_min = self.find_corresponding_y(x_min, all_coords)
            Y_x_max = self.find_corresponding_y(x_max, all_coords)
            X_y_min = self.find_corresponding_x(y_min, all_coords)
            X_y_max = self.find_corresponding_x(y_max, all_coords)

            ordered_vertices = self.order_vertices_clockwise([(x_min, Y_x_min), (x_max, Y_x_max), (X_y_max, y_max), (X_y_min, y_min)])

            # Check for repeated vertices
            #if len(set(map(tuple, ordered_vertices))) < len(ordered_vertices):
                #print(f"Repeated vertices found in Layer {tc_layer}: {ordered_vertices}")

            # Append the result to the new DataFrame
            merged_geometries.append({'ts_ieta': ts_ieta, 'ts_iphi': ts_iphi, 'tc_layer': tc_layer,
                                  'geometry': diamond_polygons, 'vertices_clockwise': ordered_vertices})

        merged_df = pd.DataFrame(merged_geometries)

        self.save_scint_mod_geo(merged_df, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/scint_modules_geo.geojson')

        return merged_df

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


    def shift_hex_values(self, silicon_df_proc, df_geom, df_ts):
        '''shifting hexagons vertices based on difference between wx_center/wy_center (byesplit) and ts_x/ts_y (CMSSW),
           in order to maintain as center of the hexagon the orginal ts_x/ts_y. 
           It returns the dataframe with the applied shifts'''
        
        # Three different x/y shifts for CE-E (subdet 1), CE-H (subdet 2) for even and odd layers
        diff_x_subdet1 = -1.387
        diff_y_subdet1 = -0.601

        diff_x_subdet2_even = -1.387
        diff_y_subdet2_even = -0.745

        diff_x_subdet2_odd = -1.387
        diff_y_subdet2_odd = -0.457

        # Make a copy of the DataFrame
        shifted_df = silicon_df_proc.copy()
        
        # Iterate over each row
        for idx, row in shifted_df.iterrows():
            # Apply shift for the hex_x and hex_y columns
            if row['ts_layer'] <= 28:
                shifted_df.at[idx, 'hex_x'] = [v + diff_x_subdet1 for v in row['hex_x']]
                shifted_df.at[idx, 'hex_y'] = [v + diff_y_subdet1 for v in row['hex_y']]
                shifted_df.at[idx, 'wx_center'] = row['wx_center'] + diff_x_subdet1
                shifted_df.at[idx, 'wy_center'] = row['wy_center'] + diff_y_subdet1

            if (row['ts_layer'] > 28 and row['ts_layer'] %2 ==0):
                shifted_df.at[idx, 'hex_x'] = [v + diff_x_subdet2_even for v in row['hex_x']]
                shifted_df.at[idx, 'hex_y'] = [v + diff_y_subdet2_even for v in row['hex_y']]
                shifted_df.at[idx, 'wx_center'] = row['wx_center'] + diff_x_subdet2_even
                shifted_df.at[idx, 'wy_center'] = row['wy_center'] + diff_x_subdet2_even
            if (row['ts_layer'] > 28 and row['ts_layer'] %2 !=0):
                shifted_df.at[idx, 'hex_x'] = [v + diff_x_subdet2_odd for v in row['hex_x']]
                shifted_df.at[idx, 'hex_y'] = [v + diff_y_subdet2_odd for v in row['hex_y']]
                shifted_df.at[idx, 'wx_center'] = row['wx_center'] + diff_x_subdet2_odd
                shifted_df.at[idx, 'wy_center'] = row['wy_center'] + diff_x_subdet2_odd

        # Plot shifted modules for chosen layer
        #layer_number= 15
        #plotMS.plot_shifted_modules(silicon_df_proc, shifted_df, df_geom, df_ts, layer_number)
        
        return shifted_df    
    
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

        #print("unique points", unique_points)    

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
    

    def process_layer(self, layer_name, data_layer, bins_layer, cart2sph):
        # Create spatial index for bins
        bin_polygons = [Polygon(bin_feature['geometry']['coordinates'][0]) for bin_feature in bins_layer]
        #print("bin_polygons", bin_polygons)
        str_tree = STRtree(bin_polygons)

        results = []

        for hex_row in data_layer.itertuples():
            hex_x = hex_row.hex_x
            hex_y = hex_row.hex_y

            ts_x = hex_row.ts_x  # CMSSW correct modules position
            ts_y = hex_row.ts_y
            ts_z = round(hex_row.ts_z, 2)

            wx_center = hex_row.wx_center  # byesplit modules position
            wy_center = hex_row.wy_center

            ts_pt= hex_row.ts_pt

            hex_polygon = Polygon(zip(hex_x, hex_y))
            #print("hex_polygon", hex_polygon)

            hex_centroid = hex_polygon.centroid
            hex_centroid_x = hex_centroid.x
            hex_centroid_y = hex_centroid.y

            hex_eta_centroid, hex_phi_centroid = cart2sph(hex_centroid_x, hex_centroid_y, z=ts_z)

            bins_overlapping = []

            # Query spatial index for potential overlapping bins
            possible_bins = str_tree.query(hex_polygon)
            possible_bins_polygons = str_tree.geometries.take(possible_bins)
            #print("possible bins", possible_bins, possible_bins_polygons)

            for bin_polygon in possible_bins_polygons:
                overlap_area = hex_polygon.intersection(bin_polygon).area
                percentage_overlap = overlap_area / hex_polygon.area if overlap_area > 0 else 0

                if percentage_overlap > 0:
                    bin_index = bin_polygons.index(bin_polygon)
                    bin_feature = bins_layer[bin_index]
                    
                    x_vertices = np.array([point[0] for point in bin_feature['geometry']['coordinates'][0]])
                    y_vertices = np.array([point[1] for point in bin_feature['geometry']['coordinates'][0]])
                    eta_vertices = np.array(bin_feature['properties'].get('Eta_vertices'))
                    phi_vertices = np.array(bin_feature['properties'].get('Phi_vertices'))
                    centroid_eta = np.mean(eta_vertices)
                    centroid_phi = np.mean(phi_vertices)

                    bins_overlapping.append({
                        'x_vertices': x_vertices,
                        'y_vertices': y_vertices,
                        'eta_vertices': eta_vertices,
                        'phi_vertices': phi_vertices,
                        'centroid_eta': centroid_eta,
                        'centroid_phi': centroid_phi,
                        'percentage_overlap': percentage_overlap
                    })

            results.append({
                'hex_x': hex_x,
                'hex_y': hex_y,
                'hex_x_centroid': hex_centroid_x,
                'hex_y_centroid': hex_centroid_y,
                'hex_eta_centroid': hex_eta_centroid,
                'hex_phi_centroid': hex_phi_centroid,
                'ts_pt': ts_pt,
                'layer': layer_name,
                'bins_overlapping': bins_overlapping
            })
        return results

    #ALGORITHMS

    def baseline_by_event(self, df_hexagon_info, subdet):
        print("Baseline OKAY algorithm ...")
        start_time = time.time()
        
        # Initialize a list to store the results for all events
        all_event_rows = []

        # Get unique event indices
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        # Iterate over each event
        for event in unique_events:
            print(f"Processing event {event}...")
            
            # Extract the DataFrame for the current event
            df_event = df_hexagon_info.loc[event]

            # Initialize a dictionary to store the summed pt for each bin by layer for the current event
            bin_pt_by_layer = {}

            # Get unique layer indices for the current event
            unique_layers = df_event.index.get_level_values('layer').unique()

            # Filter layers if required
            if subdet == 'CEE':
                print("CEE subdet ...")
                unique_layers = unique_layers[unique_layers < 29]
            elif subdet == 'CEH':
                print("CEH subdet ...")
                unique_layers = unique_layers[unique_layers >= 29]
            else:
                # No layer selection
                print("CEE + CEH ...")
                unique_layers = unique_layers

            # Iterate over each layer for the current event
            for layer_idx in unique_layers:
                layer_df = df_event.loc[df_event.index.get_level_values('layer') == layer_idx]
                #print("layer_df", layer_df )

                bin_pt = {}  # Initialize pt for bins in the current layer

                # Iterate over each hexagon in the layer
                for hex_row in layer_df.itertuples():
                    if not hex_row.bins_overlapping: # Skip hexagons with no overlapping bins
                        continue

                    # Prepare array of hexagon eta/phi centroids
                    # Prepare array of hexagon eta/phi centroids
                    hex_centroid = np.array([hex_row.hex_eta_centroid, hex_row.hex_phi_centroid])

                    # Prepare array of bin eta/phi centroids and pt
                    bin_centroids = []
                    bin_pts = []

                    for bin_info in hex_row.bins_overlapping:
                        bin_eta_centroid = bin_info['centroid_eta']
                        bin_phi_centroid = bin_info['centroid_phi']
                        bin_centroid = np.array([bin_eta_centroid, bin_phi_centroid])
                        bin_centroids.append(bin_centroid)

                        # Directly use the pt of the hexagon for the nearest bin
                        #bin_pts.append(hex_row.ts_pt)
                        bin_pts.append(hex_row.ts_pt)

                    bin_centroids = np.array(bin_centroids)

                    # Compute pairwise distances between hexagon and bins
                    distances = np.linalg.norm(bin_centroids - hex_centroid, axis=1)

                    # Find the nearest bin for the current hexagon
                    nearest_bin_idx = np.argmin(distances)

                    # Assign pt of hexagons to the nearest bins
                    nearest_bin_info = hex_row.bins_overlapping[nearest_bin_idx]
                    bin_key = (tuple(nearest_bin_info['eta_vertices']), tuple(nearest_bin_info['phi_vertices']))

                    if bin_key not in bin_pt:
                        bin_pt[bin_key] = 0

                    bin_pt[bin_key] += bin_pts[nearest_bin_idx]

                bin_pt_by_layer[layer_idx] = bin_pt

            # Convert the bin_pt_by_layer dictionary into a list of dictionaries
            event_rows = []

            for layer_idx, pt_dict in bin_pt_by_layer.items():
                for bin_key, pt in pt_dict.items():
                    row = {
                        'event': event,
                        'layer': layer_idx,
                        'eta_vertices': bin_key[0],
                        'phi_vertices': bin_key[1],
                        'pt': pt
                    }
                    event_rows.append(row)

            all_event_rows.extend(event_rows)

        end_time = time.time()
        print("Execution time loop - baseline", end_time - start_time)

        # Convert the list of dictionaries to a DataFrame
        df_baseline = pd.DataFrame(all_event_rows, columns=['event', 'layer', 'eta_vertices', 'phi_vertices', 'pt'])

        # Group by eta_vertices and phi_vertices and sum the pt
        df_baseline = df_baseline.groupby(['event', 'eta_vertices', 'phi_vertices']).agg({'pt': 'sum'}).reset_index()

        #print("df_baseline", df_baseline)
        return df_baseline


    def area_overlap_by_event(self, df_hexagon_info, subdet):
        print("Area overlap algorithm ...")
        # Initialize a list to store the results for all events
        all_event_rows = []
        # Get unique event indices
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        start_time = time.time()

        for event in unique_events:
            print(f"Processing event {event}...")
            df_event = df_hexagon_info.loc[event]
            bin_pt_by_layer = {}
            unique_layers = df_event.index.get_level_values('layer').unique()

            if subdet == 'CEE':
                print("CEE subdet ...")
                unique_layers = unique_layers[unique_layers < 29]
            elif subdet == 'CEH':
                print("CEH subdet ...")
                unique_layers = unique_layers[unique_layers >= 29]
            else:
                print("CEE + CEH subdet ...")
                unique_layers = unique_layers

            for layer_idx in unique_layers:
                layer_df = df_event.loc[df_event.index.get_level_values('layer') == layer_idx]

                bin_pt = {}  # Initialize pt for bins in the current layer

                for hex_row in layer_df.itertuples():
                    if not hex_row.bins_overlapping: # Skip hexagons with no overlapping bins
                        continue

                    for bin_info in hex_row.bins_overlapping:
                        bin_info['pt'] = hex_row.ts_pt * bin_info['percentage_overlap']
                        bin_key = (tuple(bin_info['eta_vertices']), tuple(bin_info['phi_vertices']))

                        if bin_key not in bin_pt:
                            bin_pt[bin_key] = 0

                        bin_pt[bin_key] += bin_info['pt']

                    bin_pt_by_layer[layer_idx] = bin_pt

            # Convert the bin_pt_by_layer dictionary into a list of dictionaries
            event_rows = []

            for layer_idx, pt_dict in bin_pt_by_layer.items():
                for bin_key, pt in pt_dict.items():
                    row = {
                        'event': event,
                        'layer': layer_idx,
                        'eta_vertices': bin_key[0],
                        'phi_vertices': bin_key[1],
                        'pt': pt
                    }
                    event_rows.append(row)

            all_event_rows.extend(event_rows)

        end_time = time.time()
        print("Execution time loop - area overlap by event", end_time - start_time)

        # Convert the list of dictionaries to a DataFrame
        df_2 = pd.DataFrame(all_event_rows, columns=['event', 'layer', 'eta_vertices', 'phi_vertices', 'pt'])

        # Group by eta_vertices and phi_vertices and sum the pt
        df_2 = df_2.groupby(['event', 'eta_vertices', 'phi_vertices']).agg({'pt': 'sum'}).reset_index()
        return df_2

    def area_overlap_8Towers_by_event(self, df_hexagon_info, subdet):
        print("Algorithm 1/8s towers ...")
        start_time = time.time()

        # Initialize a list to store the results for all events
        all_event_results = []

        # Get unique event indices
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        # Iterate over each event
        for event in unique_events:
            print(f"Processing event {event}...")

            # Extract the DataFrame for the current event
            df_event = df_hexagon_info.loc[event]

            # Get unique layer indices for the current event
            unique_layers = df_event.index.get_level_values('layer').unique()

            if subdet == 'CEE':
                print("CEE subdet ...")
                unique_layers = unique_layers[unique_layers < 29]
            elif subdet == 'CEH':
                print("CEH subdet ...")
                unique_layers = unique_layers[unique_layers >= 29]
            else:
                print("CEE + CEH subdet ...")

            # Iterate over unique layer indices
            for layer_idx in unique_layers:
                # Filter DataFrame for the current layer
                layer_df = df_event[df_event.index.get_level_values('layer') == layer_idx]

                # Iterate over each row of the filtered DataFrame for the current layer
                for _, row in layer_df.iterrows():
                    # Initialize pt values for all bins to 0
                    for bin_info in row['bins_overlapping']:
                        bin_info['pt'] = 0.0
                        bin_info['event'] = event  # Ensure event information is added
                        bin_info['layer'] = layer_idx  # Ensure layer information is added

                    total_pt = row['ts_pt']
                    #print("Total pt hexagon:", total_pt)

                    # Filter out bins with area overlap greater than 0
                    filtered_bins = [bin_info for bin_info in row['bins_overlapping'] if bin_info['percentage_overlap'] > 0]
                    num_bins = len(filtered_bins)

                    # Sort bins based on percentage overlap in descending order
                    sorted_bins = sorted(filtered_bins, key=lambda x: x['percentage_overlap'], reverse=True)

                    # Calculate pt distribution
                    remaining_pt = total_pt
                    top_bins = sorted_bins[:8] if num_bins > 8 else sorted_bins
                    for bin_info in top_bins:
                        percent = round(bin_info['percentage_overlap'] * 8)
                        percent = max(1, percent)  # Ensure at least 1
                        pt_fraction = percent / 8
                        pt_assigned = min(remaining_pt, pt_fraction * total_pt)
                        bin_info['pt'] = pt_assigned
                        remaining_pt -= pt_assigned
                        if remaining_pt <= 0:
                            #print("pt finished")
                            break  # Stop assigning pt if all available pt is exhausted

                    # Append results for the current layer to the event_results list
                    all_event_results.extend(row['bins_overlapping'])

        #print("All events results:", all_event_results)
        # Calculate the execution time
        end_time = time.time()
        print("Execution time loop - area overlap by event:", end_time - start_time)

        # Extract pt information for each bin and construct DataFrame
        flattened_bins = []
        for bin_info in all_event_results:
            eta_vertices = bin_info['eta_vertices']
            phi_vertices = bin_info['phi_vertices']
            pt = bin_info['pt']
            event = bin_info['event']
            layer = bin_info['layer']
            flattened_bins.append({'event': event, 'layer': layer, 'eta_vertices': tuple(eta_vertices), 'phi_vertices': tuple(phi_vertices), 'pt': pt})

        # Convert the list of dictionaries to a DataFrame
        flattened_bins_df = pd.DataFrame(flattened_bins)

        # Group by event, eta_vertices, and phi_vertices and sum the pt
        df_over_final = flattened_bins_df.groupby(['event', 'eta_vertices', 'phi_vertices']).agg({'pt': 'sum'}).reset_index()

        total_pt = df_over_final['pt'].sum()
        #print("Total pt sum:", total_pt)

        # Print intermediate results
        #print("Flattened bins:", flattened_bins_df)
        #print("Final bins:", df_over_final)

        print("Execution time loop - area 8towers:", end_time - start_time)

        return df_over_final

    def sph2cart(self, eta, phi, z=322.):
        ''' Useful conversion: Spherical coordinates to cartesian coordinates (x, y)  '''
        theta = 2*np.arctan(np.exp(-eta))
        r = z / np.cos(theta)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        return x, y
    
    def cart2sph(self, x, y, z=322.):
        ''' Useful conversion: Cartesian coordinates to spherical coordinates (eta, phi) '''
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r)
        eta = -np.log(np.tan(theta / 2))
        phi = np.arctan2(y, x)
        return eta, phi

    def find_closest_vertices(self, center_eta, center_phi, eta_values, phi_values):
        eta_mesh, phi_mesh = np.meshgrid(eta_values, phi_values)
        points = np.column_stack((eta_mesh.flatten(), phi_mesh.flatten()))
        distances = cdist([[center_eta, center_phi]], points)
        closest_indices = np.argsort(distances.flatten())[:4]
        closest_vertices = points[closest_indices]

        #sort vertices to be consecutive 
        angles = np.arctan2(closest_vertices[:, 1] - center_phi, closest_vertices[:, 0] - center_eta)
        sorted_indices = np.argsort(angles)
        closest_vertices = closest_vertices[sorted_indices]

        return closest_vertices
    
    
    def create_bin_df_new(self, kw):
            
        phi_values = np.linspace(kw['MinPhi'], kw['MaxPhi'], kw['NbinsPhi'] + 1)
        eta_values = np.linspace(kw['EtaMin'], kw['EtaMax'], kw['NbinsEta'] + 1)

        bins_info = []

        # Iterate through each layer
        for layer in range(1, 51):  # Loop through all layers from 1 to 50
            # Filter dataframe for the current layer based on the conditions
            if (layer <= 29 and layer % 2 != 0) or (layer >= 30 and layer <= 50):
                # Filter dataframe for the current layer
                layer_df = self.ds_geom['si'][self.ds_geom['si']['ts_layer'] == layer].copy()
                
                layer_df['z_appr'] = layer_df['z']
                # Extract unique 'z' values for the current layer
                unique_z_values = layer_df['z_appr'].unique()[:1]   #select only the first entry cause there are multiple z values for approximation reasons
                
                layer_bins = []

                # Calculate the bin widths for eta and phi
                d_phi = (kw['MaxPhi'] - kw['MinPhi']) / kw['NbinsPhi']
                d_eta = (kw['EtaMax'] - kw['EtaMin']) / kw['NbinsEta']

                for phi in phi_values[:-1]:  # Exclude the last value to avoid duplicating the center
                    for eta in eta_values[:-1]:  # Exclude the last value to avoid duplicating the center
                        for z_value in unique_z_values:
                            # Calculate the center of the rectangular bin
                            center_phi = phi + d_phi / 2
                            center_eta = eta + d_eta / 2

                            # Find the 4 closest vertices to the center
                            vertices = self.find_closest_vertices(center_eta, center_phi, eta_values, phi_values)

                            # Convert spherical coordinates to cartesian for vertices and center
                            x_vertices, y_vertices = zip(*[self.sph2cart(v[0], v[1], z=z_value) for v in vertices])

                            layer_bins.append({
                                'X_Vertices': x_vertices, 
                                'Y_Vertices': y_vertices,
                                'Eta_Vertices': vertices[:, 0], 
                                'Phi_Vertices': vertices[:, 1]
                            })
                
                bins_info.append({
                    'Layer': layer,
                    'Z_value': unique_z_values.tolist(),
                    'Bins': layer_bins
                })

        df_bins = pd.DataFrame(bins_info)

        '''for layer in range(1, 4):  # Iterate over the first two layers
            if (layer % 2 != 0):
                bins = df_bins.loc[df_bins['Layer'] == layer, 'Bins'].iloc[0][:5]  # Select the first 5 entries of 'Bins'
                print(f'Layer {layer} - First 5 entries of Bins:')
                for bin_info in bins:
                    print(bin_info)
                print('\n')'''

        return df_bins
    

    #SAVING

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

    def read_geojson_files(self, bins_geojson_file, hexagons_geojson_file, scint_geojson_file):
        with open(bins_geojson_file, 'r') as f:
            bins_data = json.load(f)
        with open(hexagons_geojson_file, 'r') as f:
            hexagons_data = json.load(f)
        with open(scint_geojson_file, 'r') as f:
            scint_data = json.load(f)
        return bins_data, hexagons_data, scint_data

    def update_bins(self, baseline_df, bins_df):
        # Filter the bins from the first layer
        first_layer_bins = bins_df[bins_df['Layer'] == 1]
        all_bins = []

        for _, row in first_layer_bins.iterrows():
            for bin_info in row['Bins']:
                eta_vertices = tuple(bin_info['Eta_Vertices'])
                phi_vertices = tuple(bin_info['Phi_Vertices'])
                all_bins.append((eta_vertices, phi_vertices))

        all_bins_df = pd.DataFrame(all_bins, columns=['eta_vertices', 'phi_vertices'])

        # Merge baseline_df with all_bins_df to find bins that are missing in baseline_df
        merged_df = all_bins_df.merge(baseline_df, how='left', on=['eta_vertices', 'phi_vertices'])

        # Fill NaN pt values with 0 for bins that are not in baseline_df
        merged_df['pt'] = merged_df['pt'].fillna(0)

        # No need to group by layer since we're only interested in the first layer bins
        final_df = merged_df

        #print("Final df_baseline", final_df)
        return final_df

    def apply_update_to_each_event(self, df, bins_df):
        all_event_dfs = []

        unique_events = df['event'].unique()

        for event in unique_events:
            #print(f"Processing event {event}...")
            
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            #Add bins that have 0 overlap with the hexagons with pt 0
            updated_event_df = self.update_bins(event_df, bins_df)

            # Add the event column back
            updated_event_df['event'] = event
            
            all_event_dfs.append(updated_event_df)

        # Concatenate all the event DataFrames together
        final_df = pd.concat(all_event_dfs, ignore_index=True)
        #print("SUPER FINAL", final_df)

        # Group by eta_vertices and phi_vertices and sum the pt
        summed_df = final_df.groupby(['eta_vertices', 'phi_vertices']).agg({'pt': 'sum'}).reset_index()
        #print("SUMMED FINAL", summed_df)

        return final_df, summed_df

    def save_eta_phi_differences(self, results_df, filename):
        eta_diffs = results_df['eta_diff']
        phi_diffs = results_df['phi_diff']
        
        with open(filename, 'w') as file:
            file.write("eta_diff,phi_diff\n")
            for eta_diff, phi_diff in zip(eta_diffs, phi_diffs):
                file.write(f"{eta_diff},{phi_diff}\n")

    #SPLITTING
            
    def ModSumToTowers(self, kw, data, subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, hdf5_filename, data_gen, bins_df):
        print('Mod sum to towers')

        print("data", data.columns)
        print("ts_en", "gen_pt")
        print(data['ts_pt'].sum(), data_gen['gen_pt'])
        #print("data pt TS", data['ts_pt'].sum())
        #print("data gen pt", data_gen['gen_en'])


        #overlap_data = self.read_hdf5_file(hdf5_filename)
        #print("OVERLAP DATA INPUT ", overlap_data)
        #print("OVERLAP DATA INPUT columns ", overlap_data.columns)

        hexagon_info_df = self.eval_hex_bin_overlap_OK(data, bin_geojson_filename)
        #print("QUIIIIII, hexagon_info_df", hexagon_info_df)
        #print("QUIIIIII, hexagon_info_df COLUMNS", hexagon_info_df.columns)

        if algo == 'baseline':
            df_algo = self.baseline_by_event(hexagon_info_df, subdet)
            
        elif algo == 'area_overlap':
            df_algo = self.area_overlap_by_event(hexagon_info_df, subdet)

        elif algo == '8towers':
            df_algo = self.area_overlap_8Towers_by_event(hexagon_info_df, subdet)
        
        else:
            raise ValueError("Invalid algorithm specified. Choose 'baseline', 'area_overlap' or '8towers'.")

        df , df_sum = self.apply_update_to_each_event(df_algo, bins_df)

        print("Eta/Phi resolution evaluation")
        df_resolution = self.eval_eta_phi_photon_resolution(df, data_gen, window_size=8, subwindow_size=3)
        self.save_eta_phi_differences(df_resolution, f'{algo}_{particle}_{event}_{subdet}_eta_phi_resolution_hadd_123_9x9_OK.txt')
        plotMS.plot_energy_ratio_histogram()
        #plotMS.plot_eta_phi_resolution(df_resolution, algo, event, particle, subdet)
        plotMS.plot_towers_eta_phi_grid(df_sum, data_gen, algo, event, particle, subdet)

    def compute_overlap(self, hex_polygon, hex_centroid, hex_properties, bins_layer):
        hex_centroid_x, hex_centroid_y = hex_centroid.x, hex_centroid.y
        hex_eta_centroid, hex_phi_centroid = self.cart2sph(hex_centroid_x, hex_centroid_y, hex_properties['ts_z'])

        # Create an R-tree spatial index for fast nearest-neighbor queries
        #bins_tree = STRtree([bin_feature['geometry'] for bin_feature in bins_layer])

        # Find the 35 closest bins using the spatial index
        #closest_bins = bins_tree.nearest(hex_centroid, 35)

        # Compute distances to all bin centroids and sort them
        bin_distances = [
            (bin_feature, hex_centroid.distance(bin_feature['geometry'].centroid))
            for bin_feature in bins_layer
        ]
        bin_distances.sort(key=lambda x: x[1])  # Sort by distance
        closest_bins = [bin_dist[0] for bin_dist in bin_distances[:35]] #takes total of 1.8 sec per event 

        hex_info = []
        for bin_feature in closest_bins:
            bin_polygon = bin_feature['geometry']
            overlap_area = hex_polygon.intersection(bin_polygon).area
            if overlap_area > 0:
                percentage_overlap = overlap_area / hex_polygon.area
                x_vertices = np.array([point[0] for point in bin_polygon.exterior.coords])
                y_vertices = np.array([point[1] for point in bin_polygon.exterior.coords])
                eta_vertices = np.array(bin_feature['properties'].get('Eta_vertices'))
                phi_vertices = np.array(bin_feature['properties'].get('Phi_vertices'))
                centroid_eta = np.mean(eta_vertices)
                centroid_phi = np.mean(phi_vertices)
                hex_info.append({
                    'x_vertices': x_vertices,
                    'y_vertices': y_vertices,
                    'eta_vertices': eta_vertices,
                    'phi_vertices': phi_vertices,
                    'centroid_eta': centroid_eta,
                    'centroid_phi': centroid_phi,
                    'percentage_overlap': percentage_overlap
                })

        return {
            'hex_x': hex_properties['hex_x'],
            'hex_y': hex_properties['hex_y'],
            'hex_x_centroid': hex_centroid_x,
            'hex_y_centroid': hex_centroid_y,
            'hex_eta_centroid': hex_eta_centroid,
            'hex_phi_centroid': hex_phi_centroid,
            'ts_pt': hex_properties['ts_pt'],
            'bins_overlapping': hex_info
        }


    def eval_hex_bin_overlap_OK(self, data, df_bin):
        print("Evaluating overlap between hexagons and towers bins layer by layer")


        with open(df_bin) as f:
            bin_geojson = json.load(f)

        bin_features = bin_geojson['features']
        layer_names = set(bin_feature['properties']['Layer'] for bin_feature in bin_features)
        event_groups = data.groupby('event')

        # Precompute bin geometries
        bins_by_layer = {
            layer_name: [
                {
                    'geometry': Polygon(bin_feature['geometry']['coordinates'][0]),
                    'properties': bin_feature['properties']
                } for bin_feature in bin_features if bin_feature['properties']['Layer'] == layer_name
            ] for layer_name in layer_names
        }

        # Precompute hexagon geometries
        hex_geometries = {
            (event, layer_name): [
                {
                    'hex_properties': {
                        'hex_x': hex_row.hex_x,
                        'hex_y': hex_row.hex_y,
                        'ts_z': round(hex_row.ts_z, 2),
                        'ts_pt': hex_row.ts_pt
                    },
                    'hex_polygon': Polygon(zip(hex_row.hex_x, hex_row.hex_y)),
                    'hex_centroid': Polygon(zip(hex_row.hex_x, hex_row.hex_y)).centroid
                } for hex_row in event_data[event_data['ts_layer'] == layer_name].itertuples()
            ] for event, event_data in event_groups for layer_name in layer_names
        }

        hierarchical_data = []

        start_time = time.time()
        for event, event_data in event_groups:
            event_info = {'event': event, 'layers': []}
            
            for layer_name in layer_names:
                layer_data = {'layer': layer_name, 'hexagons': []}
                bins_layer = bins_by_layer[layer_name]
                hex_info = [
                    self.compute_overlap(hex_geometry['hex_polygon'], hex_geometry['hex_centroid'], hex_geometry['hex_properties'], bins_layer)
                    for hex_geometry in hex_geometries[(event, layer_name)]
                ]
                
                layer_data['hexagons'] = hex_info
                event_info['layers'].append(layer_data)
            
            hierarchical_data.append(event_info)

        end_time = time.time()
        print("LOOP OVERLAP", end_time - start_time)

        hexagon_info = []
        max_overlapping_bins = 0
        for event_data in hierarchical_data:
            event = event_data['event']
            for layer in event_data['layers']:
                layer_idx = int(layer['layer'])  # Ensure the layer index is an integer
                for hexagon in layer['hexagons']:
                    hex_info = {
                        'event': event,
                        'layer': layer_idx,
                        'hex_x': hexagon['hex_x'],
                        'hex_y': hexagon['hex_y'],
                        'hex_x_centroid': hexagon['hex_x_centroid'],
                        'hex_y_centroid': hexagon['hex_y_centroid'],
                        'hex_eta_centroid': hexagon['hex_eta_centroid'],
                        'hex_phi_centroid': hexagon['hex_phi_centroid'],
                        'ts_pt': hexagon['ts_pt'],
                        'bins_overlapping': []
                    }

                    for bin_info in hexagon['bins_overlapping']:
                        bin_data = {
                            'x_vertices': bin_info['x_vertices'],
                            'y_vertices': bin_info['y_vertices'],
                            'eta_vertices': bin_info['eta_vertices'],
                            'phi_vertices': bin_info['phi_vertices'],
                            'centroid_eta': bin_info['centroid_eta'],
                            'centroid_phi': bin_info['centroid_phi'],
                            'percentage_overlap': bin_info['percentage_overlap']
                        }
                        hex_info['bins_overlapping'].append(bin_data)
                    
                    hexagon_info.append(hex_info)

                    # Update the maximum number of overlapping bins
                    #num_overlapping_bins = len(hexagon['bins_overlapping'])
                    #if num_overlapping_bins > max_overlapping_bins:
                        #max_overlapping_bins = num_overlapping_bins

        #print("max number of overlapping bins", max_overlapping_bins)
        df_hexagon_info = pd.DataFrame(hexagon_info)
        df_hexagon_info.set_index(['event', 'layer'], inplace=True)  # Now, 'event' and 'layer' are part of the index, not regular columns. This enables hierarchical organization

        return df_hexagon_info

    def find_particle_bin_and_evaluate_windows_2(self, baseline_df, genpart_df, window_size=8, subwindow_size=3):
        particle_eta = genpart_df['gen_eta'].iloc[0]
        particle_phi = genpart_df['gen_phi'].iloc[0]
        event_number = genpart_df['event'].iloc[0]
        print("gen part dataframe",  genpart_df)
        #print("gen part dataframe",  genpart_df.columns)

        # Find the bin where the generated particle is located
        baseline_df = baseline_df.copy()
        baseline_df['eta_center'] = baseline_df['eta_vertices'].apply(lambda x: np.mean(x))
        baseline_df['phi_center'] = baseline_df['phi_vertices'].apply(lambda x: np.mean(x))

        # Find the bin with the minimum distance to the particle's eta and phi
        particle_bin_idx = ((baseline_df['eta_center'] - particle_eta).abs() + (baseline_df['phi_center'] - particle_phi).abs()).idxmin()
        particle_bin = baseline_df.loc[particle_bin_idx]

        #print("particle bin", particle_bin)

        # Get the minimum vertex of the particle's bin
        min_eta_vertex = np.min(particle_bin['eta_vertices'])
        min_phi_vertex = np.min(particle_bin['phi_vertices'])

        # Define the window size in eta and phi
        bin_eta_size = np.max(particle_bin['eta_vertices']) - min_eta_vertex
        bin_phi_size = np.max(particle_bin['phi_vertices']) - min_phi_vertex

        # Calculate eta_min and phi_min based on the minimum vertex of the particle's bin
        eta_min = min_eta_vertex - (window_size // 2) * bin_eta_size
        eta_max = min_eta_vertex + (1+ (window_size // 2)) * bin_eta_size
        phi_min = min_phi_vertex - (window_size // 2) * bin_phi_size
        phi_max = min_phi_vertex + (1+(window_size // 2)) * bin_phi_size

        # Ensure eta and phi bounds are within the specified range
        eta_min = max(eta_min, 1.305)
        eta_max = min(eta_max, 3.045)
        phi_min = max(phi_min, -3.141593)
        phi_max = min(phi_max, 3.141593)

        window_bins = baseline_df[
            (baseline_df['eta_center'] >= eta_min) & (baseline_df['eta_center'] <= eta_max) &
            (baseline_df['phi_center'] >= phi_min) & (baseline_df['phi_center'] <= phi_max)
        ]

        # Create an empty list to store the energies for each 3x3 subwindow
        subwindow_energies = []

        # Step size for subwindow iteration
        step_eta = bin_eta_size
        step_phi = bin_phi_size

        max_pt = -1
        best_subwindow = None
        best_subwindow_data = {}

        # Iterate over each possible 3x3 subwindow within the 12x12 window
        for i in range(window_size - subwindow_size + 2):
            for j in range(window_size - subwindow_size + 2):
                eta_start = eta_min + i * step_eta
                eta_end = eta_start + subwindow_size * step_eta
                phi_start = phi_min + j * step_phi
                phi_end = phi_start + subwindow_size * step_phi

                eta_start = max(eta_start, 1.305)
                eta_end = min(eta_end, 3.045)
                phi_start = max(phi_start, -3.141593)
                phi_end = min(phi_end, 3.141593)

                subwindow = window_bins[
                    (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] <= eta_end) &
                    (window_bins['phi_center'] >= phi_start) & (window_bins['phi_center'] <= phi_end)
                ]
                #print("subwindow", subwindow.columns)
                #print("subwindow", subwindow.columns)
                if len(subwindow) == subwindow_size * subwindow_size:
                    total_pt = subwindow['pt'].sum()
                    subwindow_energies.append(total_pt)

                    particle_eta_1 = genpart_df['gen_eta'].iloc[0]
                    particle_phi_1 = genpart_df['gen_phi'].iloc[0]

                    #plotMS.plot_window_with_subwindows(window_bins, eta_min, eta_max, phi_min, phi_max, eta_start, eta_end, phi_start, phi_end, particle_eta_1, particle_phi_1)

                    if total_pt > max_pt:
                        max_pt = total_pt
                        best_subwindow = subwindow
                        best_subwindow_data = {
                        'subwindow': best_subwindow,
                        'total_pt': max_pt
                    }

        #print("subwindow data", best_subwindow_data)
        if best_subwindow_data:
            best_subwindow = best_subwindow_data['subwindow']
            best_subwindow_pt = best_subwindow_data['total_pt']

            # Calculate the weighted eta position
            if best_subwindow_data['total_pt'] < 0.001:
                # If total pt is very small, use the center of the central bin
                central_bin_index = 4  # Fifth row (assuming zero-indexed)
                weighted_eta = best_subwindow.iloc[central_bin_index]['eta_center']
                weighted_phi = best_subwindow.iloc[central_bin_index]['phi_center']
            else:
                # Otherwise, calculate weighted eta and phi positions
                total_pt = best_subwindow['pt'].sum()
                weighted_eta = (best_subwindow['pt'] * best_subwindow['eta_center']).sum() / total_pt
                weighted_phi = (best_subwindow['pt'] * best_subwindow['phi_center']).sum() / total_pt

            prova_phi = best_subwindow.iloc[4]['phi_center']
            prova_eta = best_subwindow.iloc[4]['eta_center']

            #print("prova eta e phi", prova_phi, prova_eta)

        # Calculate the difference with the actual eta of the generated particle
        eta_diff = weighted_eta - particle_eta
        phi_diff = weighted_phi - particle_phi

        # Check for infinite or NaN values
        if not np.isfinite(eta_diff) or not np.isfinite(phi_diff):
            print("Warning: Infinite or NaN values detected in eta_diff or phi_diff. Handling edge case.")
            # Handle edge case (for example, set to zero or NaN)
            eta_diff = 0.0
            phi_diff = 0.0  # or phi_diff = np.nan

        # Save required data to files --> devo lavorare tutto in pt e non in pt altriemnti non riesco a confrontare con la gen part
        genpart_pt = genpart_df['gen_pt'].iloc[0]
        #best_subwindow_pt = max_pt
        pt_ratio = best_subwindow_pt / genpart_pt if genpart_pt != 0 else np.nan
        with open('genpart_pt.txt', 'a') as f1, open('best_subwindow_pt_overlap.txt', 'a') as f2, open('pt_ratio_overlap.txt', 'a') as f3:
            f1.write(f"{genpart_pt}\n")
            f2.write(f"{best_subwindow_pt}\n")
            f3.write(f"{pt_ratio}\n")

        if pt_ratio > 1.3:
            print("event", event_number)
            print("genpart_pt", genpart_pt)
            print("best_subwindow_pt", best_subwindow_pt)
            print("pt_ratio", pt_ratio)

        events_gt_2_5 = []
        events_gt_2_5_and_eta_diff_lt_minus_0_2 = []
        # Save eta_diff and phi_diff to different files based on gen_eta condition
        if particle_eta < 2:
            with open('diffs_lt_2.txt', 'a') as f:
                f.write(f"{eta_diff},{phi_diff}\n")
        elif particle_eta > 2.5:
            with open('diffs_gt_2_5.txt', 'a') as f:
                f.write(f"{eta_diff},{phi_diff}\n")
                events_gt_2_5.append(event_number)
            if eta_diff < -0.2:
                events_gt_2_5_and_eta_diff_lt_minus_0_2.append(event_number)
                print(f"Event info: {event_number}, eta_diff: {eta_diff}, phi_diff: {phi_diff}")

        #Save event numbers to files
        with open('events_gt_2_5.txt', 'a') as f:
            for event in events_gt_2_5:
                f.write(f"{event}\n")

        with open('events_gt_2_5_and_eta_diff_lt_minus_0_2.txt', 'a') as f:
            for event in events_gt_2_5_and_eta_diff_lt_minus_0_2:
                f.write(f"{event}\n")

        return subwindow_energies, eta_diff, phi_diff

    def eval_eta_phi_photon_resolution(self, df, genpart_df, window_size=8, subwindow_size=3):
        all_results = []
        #print("df IN RESOLUTION", df.columns) #contains pt columns

        #print("GEN PART COLUMNS", genpart_df.columns)
        #print(genpart_df['event'])

        unique_events = df['event'].unique()
        #print("len(unique_events)", len(unique_events))

        for event in unique_events:
            #print(f"Processing event {event}...")
            
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            #print("event_df", event_df, "genpart_df['event']", genpart_df['gen_eta'])

            if len(unique_events) ==1: 
                event_genpart_df = genpart_df
            else: 
                event_genpart_df = genpart_df[genpart_df['event'] == event]

            #print("event_genpart_df", event_genpart_df)

            # Apply the find_particle_bin_and_evaluate_windows function
            subwindow_energies, eta_diff, phi_diff = self.find_particle_bin_and_evaluate_windows_2(event_df, event_genpart_df, window_size, subwindow_size)

            # Store the results
            result = {
                'event': event,
                'subwindow_energies': subwindow_energies,
                'eta_diff': eta_diff,
                'phi_diff': phi_diff
            }
            all_results.append(result)

        return pd.DataFrame(all_results)