_all_ = [ ]

import os
from pathlib import Path
import sys

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

from dash import dcc, html
import json
#import dash_daq as daq
import h5py
import re
import plotly.express.colors as px
import random
import logging
import time
import multiprocessing

log = logging.getLogger(__name__)

from bye_splits.plot.display_plotly import yaml, np, pd, go, dbc
from utils import params, common 
from data_handle.data_process import *
from data_handle.geometry import GeometryData
import pandas as pd
import processingMS
import plotMS
from scipy.spatial.distance import cdist
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon, Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as matPoly
from matplotlib.collections import PatchCollection



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
        processed_combined_ts_df = self.process_event(combined_ts_df,combined_tc_df) #FIXME - Temporary adding tc data
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
        print(block0_value)

        with h5py.File(file_path, 'r') as f:
            if group_name in f:
                group = f[group_name]
                if 'block1_items' in group and 'block1_values' in group:
                    block0_values = group['block0_values'][:, 0]
                    #print(" block0_values - events", block0_values)
                    block1_items = [item.decode() for item in group['block1_items']]
                    block1_values = group['block1_values']
                    if block0_value == '-1':
                        for idx, value in enumerate(block0_values):
                            block1_data.append(block1_values[idx])
                    elif int(block0_value) in block0_values:
                        idx = block0_values.tolist().index(int(block0_value))
                        block1_data.append(block1_values[idx])
                    else:
                        print(f"Block0 value '{block0_value}' not found.")

                    df = pd.DataFrame(block1_data, columns=block1_items)
                    #print("gen part data: ", df)
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

    def process_event(self, df_ts, df_tc):
        print("Process events ...")
        print("Filled data ts columns", df_ts.columns)
        print("DATA tc columns", df_tc.columns)

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
        scint_mod_geom = self.create_scint_mod_geometry(self.ds_geom['sci'])
        #plotMS.plot_scint_modules(geom_poly_scint_develop)
        #plotMS.plot_scint_tiles(self.ds_geom['sci'])
        print("Scintillator geometry", scint_mod_geom.columns)

        scintillator_df = pd.merge(left=df_tc, right=self.ds_geom['sci'], how='inner', #FIXME
                                           on=['tc_layer', 'tc_wu', 'tc_cu', 'tc_cv'])

        print("shifted_hex_df")
        print(shifted_hex_df[['hex_x', 'hex_y']])
        print("\nsilicon_df")
        print(silicon_df[['hex_x', 'hex_y']])

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

            for polygon_coords in all_coords:
                for i in range(len(polygon_coords)):
                    x, y = polygon_coords[i]
                    polygon_coords[i] = (round(x, 3), round(y, 3))

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


    def eval_hex_bin_overlap(self, data, df_bin, hdf5_filename):
        num_cores = multiprocessing.cpu_count()
        #print("Number of CPU cores:", num_cores)

        print("Evaluating overlap between hexagons and towers bins layer by layer")
        hexagon_info = [] #list of dictionaries 
        
        # Load GeoJSON file containing bin polygons
        with open(df_bin) as f:
            bin_geojson = json.load(f)
        
        # Extract bin features
        bin_features = bin_geojson['features']
        
        # Get unique layer names from bin features
        layer_names = set(bin_feature['properties']['Layer'] for bin_feature in bin_features)

        start_time = time.time()
        # Iterate over each layer
        for layer_name in layer_names:
            # Filter hexagons and bins belonging to the current layer
            data_layer = data[data['ts_layer'] == layer_name]
            bins_layer = [bin_feature for bin_feature in bin_features if bin_feature['properties']['Layer'] == layer_name]
            
            # Iterate over each hexagon in the current layer
            for hex_row in data_layer.itertuples():
                hex_x = hex_row.hex_x
                hex_y = hex_row.hex_y

                ts_x = hex_row.ts_x #CMSSW correct modules position   #FIXME
                ts_y = hex_row.ts_y
                ts_z = round(hex_row.ts_z, 2)

                wx_center = hex_row.wx_center #byesplit modules position
                wy_center = hex_row.wy_center
        
                ts_mipPt = hex_row.ts_mipPt
                
                hex_polygon = Polygon(zip(hex_x, hex_y))

                hex_centroid = hex_polygon.centroid
                hex_centroid_x = hex_centroid.x
                hex_centroid_y = hex_centroid.y

                hex_eta_centroid, hex_phi_centroid = self.cart2sph(hex_centroid_x, hex_centroid_y, z=ts_z)

                # List to store information about bins overlapping with this hexagon
                bins_overlapping = []
                
                # Iterate over each bin feature in the current layer
                for bin_feature in bins_layer:
                    bin_polygon = Polygon(bin_feature['geometry']['coordinates'][0])
                    
                    # Calculate overlap area between hexagon and bin
                    overlap_area = hex_polygon.intersection(bin_polygon).area
                  
                    #Execution time: 0.00003 seconds
                    percentage_overlap = overlap_area / hex_polygon.area if overlap_area > 0 else 0
                    
                    # If there is an overlap, store the information
                    if percentage_overlap > 0:
                        # Extract bin vertices
                        '''x_vertices = [point[0] for point in bin_feature['geometry']['coordinates'][0]]
                        y_vertices = [point[1] for point in bin_feature['geometry']['coordinates'][0]]

                        eta_vertices = bin_feature['properties'].get('Eta_vertices')  # Extract eta_vertices if available
                        phi_vertices = bin_feature['properties'].get('Phi_vertices')  # Extract phi_vertices if available
                        centroid_eta = sum(eta_vertices) / len(eta_vertices)
                        centroid_phi = sum(phi_vertices) / len(phi_vertices)'''


                        x_vertices = np.array([point[0] for point in bin_feature['geometry']['coordinates'][0]])
                        y_vertices = np.array([point[1] for point in bin_feature['geometry']['coordinates'][0]])
                        eta_vertices = np.array(bin_feature['properties'].get('Eta_vertices'))
                        phi_vertices = np.array(bin_feature['properties'].get('Phi_vertices'))
                        centroid_eta = np.mean(eta_vertices)
                        centroid_phi = np.mean(phi_vertices)

                        #bin_properties = bin_feature.get('properties', {})
                        
                        bins_overlapping.append({
                            'x_vertices': x_vertices,
                            'y_vertices': y_vertices,
                            'eta_vertices': eta_vertices,
                            'phi_vertices': phi_vertices,
                            'centroid_eta': centroid_eta,
                            'centroid_phi': centroid_phi,
                            'percentage_overlap': percentage_overlap
                        })
                
                # Append information about this hexagon to hexagon_info
                hexagon_info.append({
                    'hex_x': hex_x,
                    'hex_y': hex_y,
                    'hex_x_centroid': hex_centroid_x, #FIXME for now taking the centroid but should distinguish between the CMSSW position and byesplit position of modules
                    'hex_y_centroid': hex_centroid_y,
                    'hex_eta_centroid': hex_eta_centroid,
                    'hex_phi_centroid': hex_phi_centroid,
                    'ts_mipPt': ts_mipPt,
                    'layer': layer_name,
                    'bins_overlapping': bins_overlapping
                })

        end_time = time.time()
        execution_time = end_time - start_time
        print("Execution time loop: {:.5f} seconds".format(execution_time))        
        
        start_time1 = time.time()
        with h5py.File(hdf5_filename, 'w') as hf:
            for layer_idx, hex_info in enumerate(hexagon_info):
                layer_group = hf.require_group(f'layer_{hex_info["layer"]}')
                hex_group = layer_group.create_group(f'hex_{layer_idx}')
                hex_group.create_dataset('hex_x', data=hex_info['hex_x'])
                hex_group.create_dataset('hex_y', data=hex_info['hex_y'])
                hex_group.create_dataset('hex_x_centroid', data=hex_info['hex_x_centroid'])
                hex_group.create_dataset('hex_y_centroid', data=hex_info['hex_y_centroid'])
                hex_group.create_dataset('hex_eta_centroid', data=hex_info['hex_eta_centroid'])
                hex_group.create_dataset('hex_phi_centroid', data=hex_info['hex_phi_centroid'])
                hex_group.create_dataset('ts_mipPt', data=hex_info['ts_mipPt'])
                for idx, bin_info in enumerate(hex_info['bins_overlapping']):
                    bin_group = hex_group.create_group(f'bin_{idx}')
                    bin_group.create_dataset('x_vertices', data=bin_info['x_vertices'])
                    bin_group.create_dataset('y_vertices', data=bin_info['y_vertices'])
                    bin_group.create_dataset('eta_vertices', data=bin_info['eta_vertices'])
                    bin_group.create_dataset('phi_vertices', data=bin_info['phi_vertices'])
                    bin_group.create_dataset('centroid_eta', data=bin_info['centroid_eta'])
                    bin_group.create_dataset('centroid_phi', data=bin_info['centroid_phi'])
                    bin_group.create_dataset('percentage_overlap', data=bin_info['percentage_overlap'])
        
        end_time1 = time.time()
        execution_time1 = end_time1 - start_time1
        print("Execution time save hdf5: {:.5f} seconds".format(execution_time1)) 
        
        
        return hexagon_info
        
    def read_hdf5_file(self, hdf5_filename):
        print('reading hdf5 file')
        start_time1 = time.time()
        hexagon_info = []
        with h5py.File(hdf5_filename, 'r') as hf:
            for layer_key in hf.keys():
                layer_group = hf[layer_key]
                layer_idx = int(layer_key.split('_')[1])  # Extract layer index from key

                for hex_key in layer_group.keys():
                    hex_group = layer_group[hex_key]

                    hex_info = {
                        'layer': layer_idx,
                        'hex_x': hex_group['hex_x'][:],
                        'hex_y': hex_group['hex_y'][:],
                        'hex_x_centroid': hex_group['hex_x_centroid'][()],
                        'hex_y_centroid': hex_group['hex_y_centroid'][()],
                        'hex_eta_centroid': hex_group['hex_eta_centroid'][()],
                        'hex_phi_centroid': hex_group['hex_phi_centroid'][()],
                        'ts_mipPt': hex_group['ts_mipPt'][()],
                        'bins_overlapping': []
                    }

                    for bin_key in hex_group.keys():
                        if bin_key.startswith('bin_'):
                            bin_group = hex_group[bin_key]
                            bin_info = {
                                'x_vertices': bin_group['x_vertices'][:],
                                'y_vertices': bin_group['y_vertices'][:],
                                'eta_vertices': bin_group['eta_vertices'][:] if 'eta_vertices' in bin_group else None,
                                'phi_vertices': bin_group['phi_vertices'][:] if 'phi_vertices' in bin_group else None,

                                'centroid_eta': bin_group['centroid_eta'][()] if 'centroid_eta' in bin_group else None,
                                'centroid_phi': bin_group['centroid_phi'][()] if 'centroid_phi' in bin_group else None,
                                'percentage_overlap': bin_group['percentage_overlap'][()]
                            }
                            hex_info['bins_overlapping'].append(bin_info)

                    hexagon_info.append(hex_info)

        df_hexagon_info = pd.DataFrame(hexagon_info)

        end_time1 = time.time()
        execution_time1 = end_time1 - start_time1
        print("Execution time read hdf5: {:.5f} seconds".format(execution_time1)) 
        return df_hexagon_info
    
    #ALGORITHMS

    def baseline(self, df_hexagon_info, subdet):
        print("Baseline algorithm ...")
        start_time = time.time()
        
        # Initialize a dictionary to store the summed mipPt for each bin by layer
        bin_mipPt_by_layer = {}

        # Convert DataFrame to numpy array for faster computation
        hexagon_info_array = df_hexagon_info.to_numpy()

        # Get unique layer indices
        unique_layers = np.unique(hexagon_info_array[:, 0])

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

        # Iterate over each layer
        for layer_idx in unique_layers:
            layer_df = hexagon_info_array[hexagon_info_array[:, 0] == layer_idx]

            bin_mipPt = {}  # Initialize mipPt for bins in the current layer

            # Iterate over each hexagon in the layer
            for hex_row in layer_df:
                if not hex_row[8]:  # Skip hexagons with no overlapping bins
                    continue

                # Prepare array of hexagon eta/phi centroids
                hex_centroid = np.array([hex_row[5], hex_row[6]])

                # Prepare array of bin eta/phi centroids and mipPt
                bin_centroids = []
                bin_mipPts = []

                for bin_info in hex_row[8]:
                    bin_eta_centroid = bin_info['centroid_eta']
                    #print("bin_eta_centroid",bin_eta_centroid)
                    bin_phi_centroid = bin_info['centroid_phi']
                    bin_centroid = np.array([bin_eta_centroid, bin_phi_centroid])
                    bin_centroids.append(bin_centroid)

                    # Directly use the mipPt of the hexagon for the nearest bin
                    bin_mipPts.append(hex_row[7])

                bin_centroids = np.array(bin_centroids)

                # Compute pairwise distances between hexagon and bins
                distances = np.linalg.norm(bin_centroids - hex_centroid, axis=1)

                # Find the nearest bin for the current hexagon
                nearest_bin_idx = np.argmin(distances)

                # Assign mipPt of hexagons to the nearest bins
                nearest_bin_info = hex_row[8][nearest_bin_idx]
                bin_key = (tuple(nearest_bin_info['eta_vertices']), tuple(nearest_bin_info['phi_vertices']))

                if bin_key not in bin_mipPt:
                    bin_mipPt[bin_key] = 0

                bin_mipPt[bin_key] += bin_mipPts[nearest_bin_idx]

            bin_mipPt_by_layer[layer_idx] = bin_mipPt

        # Convert the bin_mipPt_by_layer dictionary into a list of dictionaries
        rows = []

        for layer_idx, mipPt_dict in bin_mipPt_by_layer.items():
            for bin_key, mipPt in mipPt_dict.items():
                row = {
                    'layer': layer_idx,
                    'eta_vertices': bin_key[0],
                    'phi_vertices': bin_key[1],
                    'mipPt': mipPt
                }
                rows.append(row)

        end_time = time.time()
        print("Execution time loop - baseline", end_time-start_time)

        # Convert the list of dictionaries to a numpy array
        df_baseline = np.array([(row['layer'], row['eta_vertices'], row['phi_vertices'], row['mipPt']) for row in rows],
                            dtype=[('layer', int), ('eta_vertices', object), ('phi_vertices', object), ('mipPt', float)])

        df_baseline = pd.DataFrame(df_baseline, columns=['layer', 'eta_vertices', 'phi_vertices', 'mipPt'])
        df_baseline = df_baseline.groupby(['eta_vertices', 'phi_vertices']).agg({'mipPt': 'sum'}).reset_index()

        print("df_baseline", df_baseline)
        return df_baseline
    

    def area_overlap(self, df_overlap, subdet):
        print("Area overlap algorithm ...")

        start_time = time.time()

        if subdet == 'CEE':
            print("CEE subdet ...")
            df_overlap = df_overlap[df_overlap['layer'] < 29]
        elif subdet == 'CEH':
            print("CEH subdet ...")
            df_overlap = df_overlap[df_overlap['layer'] >= 29]
        else:
            print("CEE + CEH subdet ...")
            df_overlap = df_overlap


        # Group DataFrame by 'layer' column
        grouped = df_overlap.groupby('layer')
        
        # Iterate over each layer group
        for layer_idx, layer_df in grouped:
            #print("layer idx", layer_idx)
            # Iterate over each row of the layer DataFrame
            for index, row in layer_df.iterrows():
                # Calculate the number of bins associated with the hexagon
                num_bins = len(row['bins_overlapping'])
                #print("num bins", num_bins)

                # Update the mipPt value for each bin associated with the hexagon
                for bin_info in row['bins_overlapping']:
                    bin_info['mipPt'] = row['ts_mipPt'] * bin_info['percentage_overlap']
                    #print("hex mip, %, bin mip", row['ts_mipPt']," ", bin_info['percentage_overlap'] ," ", bin_info['mipPt'] )
        
        #print("grouped df", grouped)
        #print("DF OVERLAP ALGO", df_overlap.columns)
        #print("DF OVERLAP ALGO", df_overlap)

        flattened_bins = []

        # Iterate over each row of the DataFrame
        for index, row in df_overlap.iterrows():
            # Iterate over each bin in the 'bins_overlapping' column
            for bin_info in row['bins_overlapping']:
                # Extract eta and phi vertices
                eta_vertices = bin_info['eta_vertices']
                phi_vertices = bin_info['phi_vertices']
                mipPt = bin_info['mipPt']
                
                # Append to flattened_bins list
                flattened_bins.append({'eta_vertices': eta_vertices, 'phi_vertices': phi_vertices, 'mipPt': mipPt})

        # Convert the list of dictionaries to a DataFrame
        flattened_bins_df = pd.DataFrame(flattened_bins)
        # Convert arrays to tuples

        flattened_bins_df['eta_vertices'] = flattened_bins_df['eta_vertices'].apply(tuple)
        flattened_bins_df['phi_vertices'] = flattened_bins_df['phi_vertices'].apply(tuple)

        df_over_final = flattened_bins_df.groupby(['eta_vertices', 'phi_vertices']).agg({'mipPt': 'sum'}).reset_index()

        print("flattened bins", flattened_bins_df)   
        print("final bins", df_over_final)  

        total_mipPt = df_over_final['mipPt'].sum()
        print("Total mipPt sum:", total_mipPt)

        end_time = time.time()
        print("Execution time loop - area overlap not opt", end_time-start_time)

        return df_over_final

    def area_overlap_8Towers(self, df_hexagon_info, subdet):
        print("Algoirthm 1/8s towers ...")
        start_time = time.time()

        if subdet == 'CEE':
            print("CEE subdet ...")
            df_hexagon_info = df_hexagon_info[df_hexagon_info['layer'] < 29]
        elif subdet == 'CEH':
            print("CEH subdet ...")
            df_hexagon_info = df_hexagon_info[df_hexagon_info['layer'] >= 29]
        else:
            print("CEE + CEH subdet ...")

        # Iterate over unique layer indices
        for layer_idx in df_hexagon_info['layer'].unique():
            # Filter DataFrame for the current layer
            layer_df = df_hexagon_info[df_hexagon_info['layer'] == layer_idx]

            # Iterate over each row of the filtered DataFrame for the current layer
            for index, row in layer_df.iterrows():
                # Initialize mipPt values for all bins to 0
                for bin_info in row['bins_overlapping']:
                    bin_info['mipPt'] = 0.

                total_energy = row['ts_mipPt']
                print("total energy hexagon", total_energy)

                # Calculate the number of bins associated with the hexagon
                # Filter out bins with area overlap greater than 0
                filtered_bins = [bin_info for bin_info in row['bins_overlapping'] if bin_info['percentage_overlap'] > 0]
                num_bins = len(filtered_bins)

                # Sort bins based on percentage overlap in descending order
                sorted_bins = sorted(filtered_bins, key=lambda x: x['percentage_overlap'], reverse=True)

                # Check if the number of bins is less than 8
                if num_bins <= 8:
                    print("number of bins <= 8:", num_bins )
                    remaining_energy = total_energy
                    for bin_info in sorted_bins:
                        percent = round(bin_info['percentage_overlap'] * 8) 
                        if percent == 0:
                            percent = 1
                        energy_fraction = percent/ 8
                        energy_assigned = min(remaining_energy, energy_fraction * total_energy)
                        bin_info['mipPt'] = energy_assigned
                        remaining_energy -= energy_assigned
                        if remaining_energy <= 0:
                            print("energy finished")
                            break  # Stop assigning energy if all available energy is exhausted
                else:
                    print("number of bins > 8:", num_bins)
                    # Select the top 8 bins with the highest percentage overlap
                    top_bins = sorted_bins[:8]
                    remaining_energy = total_energy
                    for bin_info in top_bins:
                        percent = round(bin_info['percentage_overlap'] * 8) 
                        if percent == 0:
                            percent = 1
                        energy_fraction = percent/ 8
                        energy_assigned = min(remaining_energy, energy_fraction * total_energy)
                        bin_info['mipPt'] = energy_assigned
                        remaining_energy -= energy_assigned
                        if remaining_energy <= 0:
                            print("energy finished")
                            break  # Stop assigning energy if all available energy is exhausted
        flattened_bins = []

        # Iterate over each row of the DataFrame
        for index, row in df_hexagon_info.iterrows():
            # Iterate over each bin in the 'bins_overlapping' column
            for bin_info in row['bins_overlapping']:
                # Extract eta and phi vertices
                eta_vertices = bin_info['eta_vertices']
                phi_vertices = bin_info['phi_vertices']
                mipPt = bin_info['mipPt']
                
                # Append to flattened_bins list
                flattened_bins.append({'eta_vertices': eta_vertices, 'phi_vertices': phi_vertices, 'mipPt': mipPt})

        # Convert the list of dictionaries to a DataFrame
        flattened_bins_df = pd.DataFrame(flattened_bins)
        # Convert arrays to tuples

        flattened_bins_df['eta_vertices'] = flattened_bins_df['eta_vertices'].apply(tuple)
        flattened_bins_df['phi_vertices'] = flattened_bins_df['phi_vertices'].apply(tuple)

        df_over_final = flattened_bins_df.groupby(['eta_vertices', 'phi_vertices']).agg({'mipPt': 'sum'}).reset_index()

        total_mipPt = df_over_final['mipPt'].sum()
        print("Total mipPt sum:", total_mipPt)

        print("flattened bins", flattened_bins_df)   
        print("final bins", df_over_final)  

        end_time = time.time()
        print("Execution time loop - area 8towers", end_time-start_time)
                   
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

    #SPLITTING
            
    def ModSumToTowers(self, kw, data, subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, hdf5_filename, data_gen):
        print('Mod sum to towers') 

        #plotMS.plot_bins_from_geojson(bin_geojson_filename, 'plot_layers')
        #plotMS.plot_bins_and_hexagons_from_geojson(bin_geojson_filename, hex_geojson_filename, 'plot_layers')

        #plotMS.plot_hexagons_from_geojson('/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/hexagons_CMSSW.geojson', 'plot_layers')
        
        #plotMS.plot_hex_bin_overlap_save(hdf5_filename, '/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/plot_layers/plot_overlap/')

        overlap_data = self.read_hdf5_file(hdf5_filename)

        if algo == 'baseline':
            df = self.baseline(overlap_data, subdet)
        
        elif algo == 'area_overlap':
            df = self.area_overlap(overlap_data, subdet)
            #df= self.area_overlap_opt(overlap_data)

        elif algo == '8towers':
            df = self.area_overlap_8Towers(overlap_data, subdet)
        
        else:
            raise ValueError("Invalid algorithm specified. Choose 'baseline', 'area_overlap' or '8towers'.")

        plotMS.plot_towers_eta_phi_grid(df, data_gen, algo, event, particle, subdet)
