_all_ = [ ]

import os
from pathlib import Path
import sys

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

from dash import dcc, html
#import dash_daq as daq
import h5py
import re
import plotly.express.colors as px
import random
import logging
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
from shapely.geometry import Polygon

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
        print(f"number of events in HDF5 file:",  len(self.list_events))   
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
                        elif 'tc' in key:
                            dataset = file[key]
                            column_names = [str(col) for col in dataset.attrs['columns']]
                            df_tc = pd.DataFrame(dataset[:], columns=column_names)
                            df_tc['event'] = int(event)  # Add a new column for the event number
                            all_tc_event_dfs.append(df_tc)  # Append each tc event data frame to the list
        
        combined_ts_df = pd.concat(all_ts_event_dfs, ignore_index=True)
        #combined_tc_df = pd.concat(all_tc_event_dfs, ignore_index=True)
        
        # Process the combined ts DataFrame
        processed_combined_ts_df = self.process_event(combined_ts_df)
        return processed_combined_ts_df
    
    
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
        print("Process events")
        print("DATA ts columns", df_ts.columns)
        #print("DATA tc columns", df_tc.columns)

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

        print("GEO silicon columns", self.ds_geom['si'].columns)
        print("GEO scint columns", self.ds_geom['sci'].columns)
        mask= self.ds_geom['sci']['tc_layer']== 42
        print("GEO scint columns", self.ds_geom['sci'][mask])
        
        ds_new = self.ds_geom['si'].drop_duplicates(subset=['ts_layer', 'ts_wu', 'ts_wv'], keep='first')

        #scintillator_df = pd.merge(left=df_tc, right=self.ds_geom['sci'], how='inner',
                                           #on=['tc_layer', 'tc_wu', 'tc_cu', 'tc_cv'])
        #plotMS.plot_layers(df_ts, ds_new)
        #plotMS.plot_layers_sci(df_tc, self.ds_geom['sci'])


        # Merge between data and geometry
        silicon_df = pd.merge(left=df_ts, right=ds_new, how='inner',
                                      on=['ts_layer', 'ts_wu', 'ts_wv'])
        silicon_df = silicon_df.drop(['triggercellu','triggercellv','waferorient', 'waferpart','diamond_x', 'diamond_y'], axis=1)
                
        #mask_ok = (silicon_df['ts_layer'] >= 28) & (silicon_df['ts_layer'] %2== 0) & (silicon_df['ts_wu'] == 0) & (silicon_df['ts_wv'] == 6)
        #print("Rows for layers >= 28:")
        #print(silicon_df[mask_ok])

        # Shifting hexagons vertices based on difference between wx_center/wy_center (byesplit) and ts_x/ts_y (CMSSW)
        #shifted_hex_df = self.shift_hex_values(silicon_df, self.ds_geom['si'], df_ts)
        processed_df = self.project_data(silicon_df) #TO MODIFY
        return processed_df
            
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
            if (row['ts_layer'] > 28 and row['ts_layer'] %2 ==0):
                shifted_df.at[idx, 'hex_x'] = [v + diff_x_subdet2_even for v in row['hex_x']]
                shifted_df.at[idx, 'hex_y'] = [v + diff_y_subdet2_even for v in row['hex_y']]   
            if (row['ts_layer'] > 28 and row['ts_layer'] %2 !=0):
                shifted_df.at[idx, 'hex_x'] = [v + diff_x_subdet2_odd for v in row['hex_x']]
                shifted_df.at[idx, 'hex_y'] = [v + diff_y_subdet2_odd for v in row['hex_y']]      

        # Plot shifted modules for chosen layer
        #layer_number= 44
        #plotMS.plot_shifted_modules(silicon_df_proc, shifted_df, df_geom, df_ts, layer_number)
        
        return shifted_df        

    '''def project_data(self, shifted_hex_df):
        # Group by hex_x and hex_y, then aggregate the ts_mipPt values by sum
        #print("ORIGINAL DATA MIP PT", shifted_hex_df['ts_mipPt'])
        shifted_hex_df['hex_x'] = shifted_hex_df['hex_x'].apply(tuple)
        shifted_hex_df['hex_y'] = shifted_hex_df['hex_y'].apply(tuple)
        projected_df = shifted_hex_df.groupby(['hex_x', 'hex_y', 'ts_x', 'ts_y', 'ts_wu',
       'ts_wv']).agg({'ts_mipPt': 'sum'}).reset_index()
        #print("PROJECTED HEXAGONS", projected_df['ts_mipPt'])
        return projected_df'''
    
    def project_data(self, shifted_hex_df):
        # Define conditions for each group
        conditions = [
            (shifted_hex_df['ts_layer'] <= 28),
            (shifted_hex_df['ts_layer'] > 28) & (shifted_hex_df['ts_layer'] % 2 != 0),
            (shifted_hex_df['ts_layer'] > 28) & (shifted_hex_df['ts_layer'] % 2 == 0)
        ]
        
        # Define corresponding group labels
        groups = ['CEE', 'CEH_odd', 'CEH_even']
        
        # Apply conditions and group labels
        shifted_hex_df['group'] = np.select(conditions, groups)
        
        # Apply tuple conversion to hex_x and hex_y
        shifted_hex_df['hex_x'] = shifted_hex_df['hex_x'].apply(tuple)
        shifted_hex_df['hex_y'] = shifted_hex_df['hex_y'].apply(tuple)

        # Round ts_x and ts_y to the third decimal number
        shifted_hex_df['ts_x'] = shifted_hex_df['ts_x'].round(3)
        shifted_hex_df['ts_y'] = shifted_hex_df['ts_y'].round(3)
        
        # Group by hex_x, hex_y, ts_layers, ts_x, ts_y, ts_wu, ts_wv, and group
        # Aggregate the ts_mipPt values by sum
        final_projected_df = shifted_hex_df.groupby(['hex_x', 'hex_y', 'wx_center', 'wy_center', 'ts_wu', 'ts_wv', 'group']).agg({'ts_mipPt': 'sum'}).reset_index()
        #use 'ts_x', 'ts_y' for CMSSW geom
       
        print("columns projected df", final_projected_df.columns)
        print("projected df", final_projected_df)

        #plotMS.plot_group_positions(final_projected_df)

        # Find rows within the same group that have the same ts_wu and ts_wv
        for group_name, group_df in final_projected_df.groupby('group'):
            duplicated_rows = group_df[group_df.duplicated(subset=['ts_wu', 'ts_wv'], keep=False)]
            if not duplicated_rows.empty:
                print(f"Duplicate rows in group {group_name}:")
                print(duplicated_rows)

        return final_projected_df
    
    
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


    def eval_hex_bin_overlap(self, data, df_bin):
        '''function that evaluates the area overlap between each hexagon and the eta/phi bins and stores the information in a dataframe'''
        #funzione che funziona sia che io gli dia il data masked che tutti i dati per cEE/CEH even e odd perche analizza un esagono per volta 
        print("EVAL_HEX_BIN_OVERLAP")
        hexagon_info = []
        polygons = []

        # Iterate over each hexagon in df
        for hex_row in data.itertuples():
            hex_x, hex_y = hex_row.hex_x, hex_row.hex_y
            hex_polygon = Polygon([(x, y) for x, y in zip(hex_x, hex_y)])

            hex_centroid = hex_polygon.centroid
            hex_x_centroid = hex_centroid.x
            hex_y_centroid = hex_centroid.y

            hex_eta_centroid, hex_phi_centroid = self.cart2sph(hex_x_centroid, hex_y_centroid, z=322.)
         
            # List to store information about bins overlapping with this hexagon
            bins_overlapping = []
            center = np.array([0., 0.])

            # Iterate over each bin in df_bin
            for bin_row in df_bin.itertuples():
                bin_x, bin_y = bin_row.bin_vertex_X, bin_row.bin_vertex_Y
                bin_eta, bin_phi = bin_row.eta_vertices, bin_row.phi_vertices
                #print("binx", bin_x)

                #create polygon with arc of circumference
                bin_polygon = self.create_bin_polygon(bin_x, bin_y, bin_eta, bin_phi, center)
                polygons.append(bin_polygon)
                #print("Polygon ARCO DI CIRCONFERENZA:", polygon)

                #bin_polygon = Polygon([(x, y) for x, y in zip(bin_x, bin_y)]) #polygon with straight lines

                # Calculate centroid of the bin
                bin_centroid = bin_polygon.centroid
                bin_centroid_x = bin_centroid.x
                bin_centroid_y = bin_centroid.y

                # Calculate the overlap area between the bin and hexagon
                overlap_area = hex_polygon.intersection(bin_polygon).area   #WORKING BUT SLOW--> OPTIMIZE
                percentage_overlap = overlap_area / hex_polygon.area if overlap_area > 0 else 0

                # If there is an overlap, store the information
             
                bins_overlapping.append({
                        'centroid_X': bin_centroid_x,
                        'centroid_Y': bin_centroid_y,
                        'percentage_overlap': percentage_overlap,
                        'eta_vertices': bin_row.eta_vertices,
                        'phi_vertices': bin_row.phi_vertices
                    })

            # Append the information about this hexagon to the hexagon_info list
            hexagon_info.append({
                'hex_x': hex_row.hex_x,
                'hex_y': hex_row.hex_y,
                'hex_x_centroid': hex_x_centroid,
                'hex_y_centroid': hex_y_centroid,
                'hex_eta_centroid': hex_eta_centroid,
                'hex_phi_centroid': hex_phi_centroid,
                'ts_x': hex_row.ts_x,  # Add ts_x column
                'ts_y': hex_row.ts_y,  # Add ts_y column
                'ts_mipPt': hex_row.ts_mipPt,  # Add ts_mipPt column
                'bins_overlapping': bins_overlapping
            })

        # Create a DataFrame from hexagon_info
        df_overlap = pd.DataFrame(hexagon_info)

        # Print hexagon coordinates, mipPt, and number of overlapping bins
        for idx, row in df_overlap.iterrows():
            overlapping_bins = [bin_info for bin_info in row['bins_overlapping'] if bin_info['percentage_overlap'] > 0]
            num_overlapping_bins = len(overlapping_bins)
            print(f"Hexagon {idx}:")
            print(f"  Coordinates: ({row['hex_x_centroid']}, {row['hex_y_centroid']})")
            print(f"  Coordinates eta/phi: ({row['hex_eta_centroid']}, {row['hex_phi_centroid']})")
            print(f"  MipPt: {row['ts_mipPt']}")
            print(f"  Number of overlapping bins: {num_overlapping_bins}")

            '''print("  Bins overlapping:")
            for bin_info in row['bins_overlapping']:
                if bin_info['percentage_overlap'] > 0:
                    print(f"    Centroid: ({bin_info['centroid_X']}, {bin_info['centroid_Y']})")
                    print(f"    Percentage overlap: {bin_info['percentage_overlap']}")'''
            
        #self.plot_polygons(polygons)   
        return df_overlap


    def ModSumToEightTowers(self, df_overlap): 
        print("new algo")
        '''assign the module energy to maximum 8 towers depending on eta values, 
        if 2.5<hex_eta_centroid<3--> 8 towers, 
        if 2<hex_eta_centroid<2.5--> 5 towers, 
        if 1.5<hex_eta_centroid<2--> 3 towers,  '''
        # Initialize a dictionary to store the ts_mipPt values for each bin
        mip_pt_for_bins = {}
        hexagon_coordinates = df_overlap[['ts_x', 'ts_y']]
        first_hexagon = df_overlap.iloc[0]
        unique_bin_coordinates = []
    
        # Iterate over the bins overlapping with the first hexagon
        for entry in first_hexagon['bins_overlapping']:
            bin_coordinate2 = (entry['centroid_X'], entry['centroid_Y'])
            # Check if the bin coordinate is already in the list
            if bin_coordinate2 not in unique_bin_coordinates:
                # If not, add it to the list
                unique_bin_coordinates.append(bin_coordinate2)

        bin_coordinates = np.array(unique_bin_coordinates)

        # Compute distances using cdist
        distances = cdist(hexagon_coordinates, bin_coordinates)

        #print("Distances:")
        #print(distances)
        #print(distances.shape)
        #print()
      
        '''for row_idx, row in enumerate(distances):
            # Find the minimum distance in the current row
            min_distance_index = np.argmin(row)
            min_distance = np.min(row)
            # Get the centroid of the bin with the minimum distance
            min_distance_bin_centroid = bin_coordinates[min_distance_index]
            # Get the ts_x and ts_y coordinates of the hexagon
            hexagon_ts_x = hexagon_coordinates.iloc[row_idx]['ts_x']
            hexagon_ts_y = hexagon_coordinates.iloc[row_idx]['ts_y']'''

        # Iterate over each row in df_overlap
        for idx, row in enumerate(df_overlap.itertuples()):
            ts_mipPt = row.ts_mipPt
            hex_x_centroid = row.ts_x
            hex_y_centroid = row.ts_y

            if 2.5 <= row.hex_eta_centroid <= 3:
                print("eta 8 bins ", row.hex_eta_centroid)
                num_bins = min(8, len(row.bins_overlapping))  # set the number of bins
                if num_bins <8:
                    print("for 2.5 <= eta <= 3 superimposed bins are lower than 8, num b in is: ", num_bins) # Ensure we have at least 8 bins

                print("ts_mipPt", ts_mipPt)
                energy_per_tower = ts_mipPt/num_bins
                print("Energy per tower", energy_per_tower)

                nearest_bin_indices = np.argsort(distances[idx])[:num_bins]

                # Assign 1/8 of the energy to each of the 8 nearest bins
                for nearest_idx in nearest_bin_indices:
                    nearest_bin_info = row.bins_overlapping[nearest_idx]

                    # Use centroid_X and centroid_Y as keys in mip_pt_for_bins dictionary
                    bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

                    # If the bin centroid is not already in the dictionary, initialize it with 0
                    if bin_key not in mip_pt_for_bins:
                        mip_pt_for_bins[bin_key] = 0

                    # Add energy_per_bin value to the corresponding bin centroid
                    mip_pt_for_bins[bin_key] += energy_per_tower

                    nearest_bin_info['mipPt'] = mip_pt_for_bins[bin_key]
                    
                print("bins energy distribution", mip_pt_for_bins)
               
            if 2 <= row.hex_eta_centroid < 2.5:
                print(" eta 5 bins:", row.hex_eta_centroid)
                num_bins = min(5, len(row.bins_overlapping))  # Ensure we have at most 5 bins
                if num_bins < 5:
                    print("For 2 <= eta < 2.5, superimposed bins are lower than 5, num bins is:", num_bins)

                # Calculate energy per tower
                print("ts_mipPt", ts_mipPt)
                energy_per_tower = ts_mipPt/8.
                print("Energy per tower:", energy_per_tower)

                # Sort the distances and get the indices of the nearest bins
                nearest_bin_indices = np.argsort(distances[idx])[:num_bins]

                # Assign energy to the nearest bins
                for i, nearest_idx in enumerate(nearest_bin_indices):
                    nearest_bin_info = row.bins_overlapping[nearest_idx]

                    # Use centroid_X and centroid_Y as keys in mip_pt_for_bins dictionary
                    bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

                    # Determine energy portion based on the index
                    if i < 3:
                        energy_portion = 2 * energy_per_tower
                    else:
                        energy_portion = 1 * energy_per_tower

                    # If the bin centroid is not already in the dictionary, initialize it with 0
                    if bin_key not in mip_pt_for_bins:
                        mip_pt_for_bins[bin_key] = 0

                    # Add energy portion to the corresponding bin centroid
                    print("energy portion", energy_portion)
                    mip_pt_for_bins[bin_key] += energy_portion

                    nearest_bin_info['mipPt'] = mip_pt_for_bins[bin_key]

                print("Bins energy distribution:", mip_pt_for_bins) 

            elif 1.5 <= row.hex_eta_centroid < 2:
                print("eta 3 bins ", row.hex_eta_centroid)
                num_bins = min(3, len(row.bins_overlapping))  # Ensure we have at most 3 bins
                if num_bins <3:
                    print("for 2.5 <= eta <= 3 superimposed bins are lower than 8, num b in is: ", num_bins) # Ensure we have at most 8 bins
                
                # Calculate energy per tower
                print("ts_mipPt", ts_mipPt)    
                energy_per_tower = ts_mipPt/8.
                print("Energy per tower:", energy_per_tower)

                # Sort the distances and get the indices of the nearest bins
                nearest_bin_indices = np.argsort(distances[idx])[:num_bins]

                # Assign energy to the nearest bins
                for i, nearest_idx in enumerate(nearest_bin_indices):
                    nearest_bin_info = row.bins_overlapping[nearest_idx]

                    # Use centroid_X and centroid_Y as keys in mip_pt_for_bins dictionary
                    bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

                    # Determine energy portion based on the index
                    if i < 2:
                        energy_portion = 3 * energy_per_tower
                    else:
                        energy_portion = 2 * energy_per_tower

                    # If the bin centroid is not already in the dictionary, initialize it with 0
                    if bin_key not in mip_pt_for_bins:
                        mip_pt_for_bins[bin_key] = 0

                    # Add energy portion to the corresponding bin centroid
                    print("energy portion", energy_portion)
                    mip_pt_for_bins[bin_key] += energy_portion

                    nearest_bin_info['mipPt'] = mip_pt_for_bins[bin_key]

                print("Bins energy distribution:", mip_pt_for_bins) 
        
        for row in df_overlap['bins_overlapping']:
            for bin_info in row:
                if 'mipPt' not in bin_info:
                    bin_info['mipPt'] = 0    

        for row in df_overlap.itertuples():
                    print("Bins Overlapping:")
                    print("  ts_mipPt:", row.ts_mipPt)
                    total_mipPt = 0
                    for bin_info in row.bins_overlapping:
                        if bin_info['mipPt']>0:
                            print("  Bin Centroid X:", bin_info['centroid_X'])
                            print("  Bin Centroid Y:", bin_info['centroid_Y'])
                            print("  mipPt:", bin_info['mipPt'])
                            total_mipPt += bin_info['mipPt']

                    print("Total mipPt for Hexagon:", total_mipPt)

        return df_overlap

    def assign_mip_pt_to_bins_centroid(self, df_overlap):
        # Initialize a dictionary to store the ts_mipPt values for each bin
        mip_pt_for_bins = {}
        hexagon_coordinates = df_overlap[['ts_x', 'ts_y']]
        unique_bin_coordinates = []

        first_hexagon = df_overlap.iloc[0]
    
        # Iterate over the bins overlapping with the first hexagon
        for entry in first_hexagon['bins_overlapping']:
            bin_coordinate = (entry['centroid_X'], entry['centroid_Y'])
            # Check if the bin coordinate is already in the list
            if bin_coordinate not in unique_bin_coordinates:
                # If not, add it to the list
                unique_bin_coordinates.append(bin_coordinate)

        bin_coordinates = np.array(unique_bin_coordinates)
        # Compute distances between each hexagon and each bin
        distances = cdist(hexagon_coordinates, bin_coordinates)
        # Find the index of the nearest bin centroid for each hexagon
        nearest_bin_indices = np.argmin(distances, axis=1)
        # Set a unique set of keys to identify bins
        unique_bin_keys = set()

        # Iterate over each row in df_overlap
        for idx, row in enumerate(df_overlap.itertuples()):
            # Extract ts_mipPt value for the current hexagon
            ts_mipPt = row.ts_mipPt

            # Get the nearest bin index and info for the current hexagon
            nearest_idx = nearest_bin_indices[idx]
            nearest_bin_info = row.bins_overlapping[nearest_idx]
            
            # Use tuple directly for bin key
            bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

            # Check if bin_key is already in mip_pt_for_bins or not
            if bin_key not in unique_bin_keys:
                # If not, add it to the set of unique bin keys and initialize it in mip_pt_for_bins
                unique_bin_keys.add(bin_key)
                mip_pt_for_bins[bin_key] = 0
                
            # Add ts_mipPt value to the corresponding bin centroid
            mip_pt_for_bins[bin_key] += ts_mipPt

        print("Set of unique bins with associated mipPt", mip_pt_for_bins)

        # Iterate over each row in df_overlap to update 'mipPt' value
        for idx, row in enumerate(df_overlap.itertuples()):
            nearest_idx = nearest_bin_indices[idx]
            nearest_bin_info = row.bins_overlapping[nearest_idx]
            bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

            # If the bin key has been encountered before, skip updating 'mipPt'
            if bin_key in unique_bin_keys:
                # Update 'bin_mipPt' value in df_overlap with ts_mipPt from mip_pt_for_bins
                nearest_bin_info['mipPt'] = mip_pt_for_bins.get(bin_key, 0)
                unique_bin_keys.remove(bin_key)
        
        for row in df_overlap['bins_overlapping']:
            for bin_info in row:
                if 'mipPt' not in bin_info:
                    bin_info['mipPt'] = 0

        '''for row in df_overlap.itertuples():
                print("Bins Overlapping:")
                total_mipPt = 0
                for bin_info in row.bins_overlapping:
                    if bin_info['mipPt']>0:
                        print("  Bin Centroid X:", bin_info['centroid_X'])
                        print("  Bin Centroid Y:", bin_info['centroid_Y'])
                        print("  mipPt:", bin_info['mipPt'])
                        total_mipPt += bin_info['mipPt']
                if total_mipPt>0:        
                    print("Total mipPt for Hexagon:", total_mipPt)'''

        return df_overlap
    

    def assign_mipPt_to_bins_area_overlap(self, df_overlap): 
        # Iterate over each row of the DataFrame
        for index, row in df_overlap.iterrows():
            # Calculate the number of bins associated with the hexagon
            num_bins = len(row['bins_overlapping'])

            # Update the mipPt value for each bin associated with the hexagon
            for bin_info in row['bins_overlapping']:
                #print(bin_info)
                bin_info['mipPt'] = row['ts_mipPt'] * bin_info['percentage_overlap']

        return df_overlap
    
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

    def create_bin_df(self, kw):
        
        phi_values = np.linspace(kw['MinPhi'], kw['MaxPhi'], kw['NbinsPhi'] + 1)
        eta_values = np.linspace(kw['EtaMin'], kw['EtaMax'], kw['NbinsEta'] + 1)
        
        grid_bin = []
        # Calculate the bin widths for eta and phi
        d_phi = (kw['MaxPhi'] - kw['MinPhi']) / kw['NbinsPhi']
        d_eta = (kw['EtaMax'] - kw['EtaMin']) / kw['NbinsEta']

        for phi in phi_values[:-1]:  # Exclude the last value to avoid duplicating the center
            for eta in eta_values[:-1]:  # Exclude the last value to avoid duplicating the center
                # Calculate the center of the rectangular bin
                center_phi = phi + d_phi / 2
                center_eta = eta + d_eta / 2

                # Find the 4 closest vertices to the center
                vertices = self.find_closest_vertices(center_eta, center_phi, eta_values, phi_values)

                #for v in vertices:
                    #print("v0", v[0]) #all x
                    #print("v1", v[1]) #all y 

                x_vertices, y_vertices = zip(*[self.sph2cart(v[0], v[1], z=322.) for v in vertices])

                #print("x_vertices", x_vertices, "y_vertices", y_vertices)
 
                x_center, y_center = self.sph2cart(center_eta, center_phi, z=322.)

                # Create a Shapely polygon for the bin
                bin_polygon = Polygon(zip(x_vertices, y_vertices))
                
                # Ensure the polygon is valid
                if not bin_polygon.is_valid:
                    bin_polygon = bin_polygon.buffer(0)
                
                # Calculate the centroid using Shapley method
                bin_centroid = bin_polygon.centroid
                
                grid_bin.append({
                    'center_X': x_center, 'center_Y': y_center,
                    'centroid_X': bin_centroid.x, 'centroid_Y': bin_centroid.y,
                    'bin_vertex_X': x_vertices, 'bin_vertex_Y': y_vertices,
                    'eta_vertices': vertices[:, 0], 'phi_vertices': vertices[:, 1]
                })
                
        df_bin= pd.DataFrame(grid_bin)
        print("BIN_DF", df_bin.columns)

        return df_bin 


    def create_grid_df(self,kw): #only used for plotting
        # Implement the logic for creating the grid DataFrame based on the parameters in kw
        # For demonstration purposes, a simple example is provided here
        phi_values = np.linspace(kw['MinPhi'], kw['MaxPhi'], kw['NbinsPhi'] + 1)
        eta_values = np.linspace(kw['EtaMin'], kw['EtaMax'], kw['NbinsEta'] + 1)
        
        grid_data = []

        for phi in phi_values:
            for eta in eta_values:
                x, y = self.sph2cart(eta,phi,z=322.)
                grid_data.append({'X': x, 'Y': y, 'Phi': phi, 'Eta': eta})

    
        df_grid = pd.DataFrame(grid_data)
        return df_grid      


    # MODULE SPLITTING

    def ModSumToTowers(self, kw, data, df_bin, algo, subdet, event, particle):    
        print('mod sum to towers')                                                                                                               
        df_overlap = self.eval_hex_bin_overlap(data, df_bin) 
        #plotMS.plot_single_event_proj(kw, data, df_bin, "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/20x72_marco_geo_mod.png")

        if algo == 'baseline':
            print('baseline algorithm')
            df_centroid = self.assign_mip_pt_to_bins_centroid(df_overlap)
            plotMS.plot_grid_area_overlap(df_centroid, algo, subdet, event, particle)

        elif algo == 'area_overlap':
            print('area_overlap algorithm')
            df_area_overlap = self.assign_mipPt_to_bins_area_overlap(df_overlap) 
            plotMS.plot_grid_area_overlap(df_area_overlap, algo, subdet, event, particle) 

        elif algo == 'algo3':
            print('algo3')
            df = self.ModSumToEightTowers(df_overlap)
            plotMS.plot_grid_area_overlap(df, algo, subdet, event, particle)     

        else:
            raise ValueError("Invalid algorithm specified. Choose 'baseline', 'area_overlap' or algo3.")