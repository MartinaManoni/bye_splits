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
from shapely.geometry import Polygon
import pandas as pd
import processingMS
import plotMS
from scipy.spatial.distance import cdist
from shapely.geometry import Polygon

class Processing():
    def __init__(self):
        with open(params.CfgPath, "r") as afile:
            self.cfg = yaml.safe_load(afile)

        self.ds_geom = GeometryData(reprocess=True, logger=log, library='plotly').provide()
        #print("qui", self.ds_geom["si"].columns)
        self.filename = None
        self.list_events = None


    def random_event(self, f):
        return random.choice(self.list_events)


    def filestore(self):
        self.filename = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/new_algos3/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"


    def get_data(self, event='4773'):
        #Check if the HDF5 file exists
        self.filestore() 
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"HDF5 file '{self.filename}' not found.")

        # Read the list of events from the HDF5 file
        with h5py.File(self.filename, 'r') as file:
            for key in file.keys():
                #print("key", key)
                self.list_events = [key.split('_')[1] for key in file.keys()]
        print(f"List of events in HDF5 file: {self.list_events}")    
        # Determine the event to use (either specified or a random one)
        event = event or self.random_event(self.filename)
        print(f"Selected event: {event}")
        # If the specified event is not in the list, raise an exception
        if str(event) not in self.list_events:
            raise ValueError(f"Event {event} not found in the HDF5 file.")

        # Get the data for the specified event
        processed_event_df = self.get_event(event)
        return processed_event_df, event
    

    def get_data_new(self, event=None):
        # Check if the HDF5 file exists
        self.filestore()
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"HDF5 file '{self.filename}' not found.")

        # Read the list of events from the HDF5 file
        with h5py.File(self.filename, 'r') as file:
            self.list_events = list(set([key.split('_')[1] for key in file.keys()]))
        print(f"List of events in HDF5 file: {self.list_events}")    
        
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
        all_event_dfs = []
        print("get all events")
        
        for event in events_to_process:
            with h5py.File(self.filename, 'r') as file:
                for key in file.keys():
                    if event in key and 'ts' in key:
                        dataset = file[key]
                        column_names = [str(col) for col in dataset.attrs['columns']]
                        df_ts = pd.DataFrame(dataset[:], columns=column_names)
                        df_ts['event'] = int(event)  # Add a new column for the event number
                        all_event_dfs.append(df_ts)  # Append each event data frame to the list

        # Concatenate all event data frames into a single data frame
        combined_df = pd.concat(all_event_dfs, ignore_index=True)
        
        # Process the combined DataFrame
        processed_combined_df = self.process_event(combined_df)
        
        return processed_combined_df
                
    
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
        print("process events")
        ts_keep = {'waferu'       : 'ts_wu',
                   'waferv'       : 'ts_wv',
                   'layer'        : 'ts_layer'}
        
        self.ds_geom['si']  = self.ds_geom['si'].rename(columns=ts_keep)
       
        ds_new = self.ds_geom['si'].drop_duplicates(subset=['ts_layer', 'ts_wu', 'ts_wv'], keep='last')
        silicon_df = pd.merge(left=df_ts, right=ds_new, how='inner',
                                      on=['ts_layer', 'ts_wu', 'ts_wv'])
        silicon_df = silicon_df.drop(['triggercellu','triggercellv','waferorient', 'waferpart','diamond_x', 'diamond_y'], axis=1)

        shifted_hex_df = processingMS.Processing().shift_hex_values_for_all_entries(silicon_df)
        print("shifted_hex_df columns", shifted_hex_df.columns)

        #print("SHIFTED DATA", shifted_hex_df)
        processed_df = processingMS.Processing().project_data(shifted_hex_df)

        #print("rows ts", df_ts.shape[0])
        #print("ts columns second ", silicon_df.columns)
        #print("rows GEOM", silicon_df.shape[0])
        #silicon_df = silicon_df.drop(['waferorient', 'waferpart'], axis=1)
        #print("columns ",silicon_df.columns)
        #print("columns ",silicon_df['ts_x'], silicon_df['ts_y'], silicon_df['wx_center'], silicon_df['wy_center'], silicon_df['hex_x'][0], silicon_df['hex_y'][0] )
        #print("columns ",silicon_df['y'])
        print("ZERO")
        return processed_df #this dataframe is the data + geometry dataframe, the hexagons are shifted and projected only on one plane 
    

    def eval_and_print_geo_diff(self,silicon_df_proc):
        # Calculate differences and print them
        differences = silicon_df_proc[['ts_x', 'wx_center', 'ts_y', 'wy_center']].apply(
            lambda row: (row['ts_x'] - row['wx_center'], row['ts_y'] - row['wy_center']), axis=1
        )
        
        for idx, (diff_x, diff_y) in enumerate(differences):
            print(f"Difference for entry {idx + 1}:")
            print(f"  ts_x - wx_center: {diff_x}")
            print(f"  ts_y - wy_center: {diff_y}")
            
            
    def shift_hex_values_for_all_entries(self, silicon_df_proc):
        diff_x = silicon_df_proc.at[0, 'ts_x'] - silicon_df_proc.at[0, 'wx_center']
        diff_y = silicon_df_proc.at[0, 'ts_y'] - silicon_df_proc.at[0, 'wy_center']
        
        # Make a copy of the DataFrame
        modified_df = silicon_df_proc.copy()
        
        # Iterate over each row
        for idx, row in modified_df.iterrows():
            # Modify the hex_x and hex_y columns
            modified_df.at[idx, 'hex_x'] = [v + diff_x for v in row['hex_x']]
            modified_df.at[idx, 'hex_y'] = [v + diff_y for v in row['hex_y']]

        '''print("Hex values for all entries after shifting:")
        for idx, row in modified_df.iterrows():
            print(f"For entry {idx}:")
            print(f"  hex_x: {row['hex_x']}")
            print(f"  hex_y: {row['hex_y']}")'''
        return modified_df        


    def project_data(self, shifted_hex_df):
        # Group by hex_x and hex_y, then aggregate the ts_mipPt values by sum
        #print("ORIGINAL DATA MIP PT", shifted_hex_df['ts_mipPt'])
        shifted_hex_df['hex_x'] = shifted_hex_df['hex_x'].apply(tuple)
        shifted_hex_df['hex_y'] = shifted_hex_df['hex_y'].apply(tuple)
        projected_df = shifted_hex_df.groupby(['hex_x', 'hex_y', 'ts_x', 'ts_y', 'ts_wu',
       'ts_wv']).agg({'ts_mipPt': 'sum'}).reset_index()
        #print("PROJECTED HEXAGONS", projected_df['ts_mipPt'])
        return projected_df
    

    def eval_hex_bin_overlap(self, data, df_bin):
        hexagon_info = []

        # Iterate over each hexagon in aggregated_df
        for hex_row in data.itertuples():
            hex_x, hex_y = hex_row.hex_x, hex_row.hex_y
            hex_polygon = Polygon([(x, y) for x, y in zip(hex_x, hex_y)])

            # List to store information about bins overlapping with this hexagon
            bins_overlapping = []

            # Iterate over each bin in df_bin
            for bin_row in df_bin.itertuples():
                bin_x, bin_y = bin_row.bin_vertex_X, bin_row.bin_vertex_Y
                bin_polygon = Polygon([(x, y) for x, y in zip(bin_x, bin_y)])

                   # Calculate centroid of the bin
                bin_centroid = bin_polygon.centroid
                bin_centroid_x = bin_centroid.x
                bin_centroid_y = bin_centroid.y

                # Calculate the overlap area between the bin and hexagon
                overlap_area = hex_polygon.intersection(bin_polygon).area
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
                'ts_x': hex_row.ts_x,  # Add ts_x column
                'ts_y': hex_row.ts_y,  # Add ts_y column
                'ts_mipPt': hex_row.ts_mipPt,  # Add ts_mipPt column
                'bins_overlapping': bins_overlapping
            })

        # Create a DataFrame from hexagon_info
        df_overlap = pd.DataFrame(hexagon_info)
        #print("df_overlap col:", df_overlap.columns)
        #print("df_overlap bins overapping columns", df_overlap['bins_overlapping'])
        for col in df_overlap['bins_overlapping'].iloc[0][0].keys():
            print("Column bin overalp:", col)

        hexagon_index = 0
        hexagon_info = df_overlap.loc[hexagon_index]
        hex_x = hexagon_info['hex_x']
        hex_y = hexagon_info['hex_y']
        bins_overlapping = hexagon_info['bins_overlapping']

        #print("Hexagon coordinates:", hex_x, hex_y)
        for bin_info in bins_overlapping:
            centroid_x = bin_info['centroid_X']
            centroid_y = bin_info['centroid_Y']
            percentage_overlap = bin_info['percentage_overlap']
            #print("Bin centroid:", centroid_x, centroid_y)
            #print("Percentage overlap:", percentage_overlap)
        return df_overlap

    def assign_mip_pt_to_bins_centroid(self, df_overlap):
        # Initialize a dictionary to store the ts_mipPt values for each bin centroid
        mip_pt_for_bins = {}

        # Calculate distances between each (ts_x, ts_y) from df_overlap and each (centroid_X, centroid_Y) from bins_overlapping
        distances = cdist(df_overlap[['ts_x', 'ts_y']], 
                        np.array([[entry['centroid_X'], entry['centroid_Y']] for sublist in df_overlap['bins_overlapping'] for entry in sublist]))

        # Find the index of the nearest bin centroid for each hexagon
        nearest_bin_indices = np.argmin(distances, axis=1)

        # Iterate over each row in df_overlap
        for idx, row in enumerate(df_overlap.itertuples()):
            # Extract ts_mipPt value for the current hexagon
            ts_mipPt = row.ts_mipPt

            # Get the nearest bin index for the current hexagon
            nearest_idx = nearest_bin_indices[idx]

            # Get the nearest bin info
            nearest_bin_info = row.bins_overlapping[nearest_idx]
            
            # Use centroid_X and centroid_Y as keys in mip_pt_for_bins dictionary
            bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])

            # If the bin centroid is not already in the dictionary, initialize it with 0
            if bin_key not in mip_pt_for_bins:
                mip_pt_for_bins[bin_key] = 0
                
            # Add ts_mipPt value to the corresponding bin centroid
            mip_pt_for_bins[bin_key] += ts_mipPt

        # Update df_overlap with the mip_pt_for_bins information
        for idx, row in enumerate(df_overlap.itertuples()):
            nearest_idx = nearest_bin_indices[idx]
            nearest_bin_info = row.bins_overlapping[nearest_idx]
            bin_key = (nearest_bin_info['centroid_X'], nearest_bin_info['centroid_Y'])
            
            # Update 'bin_mipPt' value in df_overlap with ts_mipPt from mip_pt_for_bins
            nearest_bin_info['mipPt'] = mip_pt_for_bins.get(bin_key, 0)
        # Assign mipPt = 0 to bins for which mipPt is not assigned
            
        for row in df_overlap['bins_overlapping']:
            for bin_info in row:
                if 'mipPt' not in bin_info:
                    bin_info['mipPt'] = 0

            # Iterate over each row in df_overlap
        for idx, row in df_overlap.iterrows():
            print(f"Row {idx}:")
            # Iterate over each bin overlapping with the hexagon
            for bin_info in row['bins_overlapping']:
                #print("centroid_X:", bin_info['centroid_X'])
                #print("centroid_Y:", bin_info['centroid_Y'])
                #print("percentage_overlap:", bin_info['percentage_overlap'])
                #print("eta_vertices:", bin_info['eta_vertices'])
                #print("phi_vertices:", bin_info['phi_vertices'])
                # Check if 'bin_mipPt' exists before accessing
                #if bin_info['mipPt']>0:
                    print("bin_mipPt:", bin_info['mipPt'])
                    #print("bin_mipPt:", bin_info['centroid_X'])
                    #print("bin_mipPt:", bin_info['centroid_Y'])
            print("\n")

        return df_overlap

    def assign_mipPt_to_bins_area_overlap(self, df_overlap):  #in questo modo assegno a tutti i bin la stessa energia (energia totale del modulo divisio il numero di bin)
        #merged_df = pd.merge(hexagon_info_df, data, on=['hex_x', 'hex_y'], how='left')
        #print("Merged DataFrame COLUMNS:", merged_df.columns)
        #print("Merged DataFrame:", merged_df)

        # Iterate over each row of the DataFrame
        for index, row in df_overlap.iterrows():
            # Calculate the number of bins associated with the hexagon
            num_bins = len(row['bins_overlapping'])
            #print("number of bins", num_bins)

            # Calculate the mipPt for each bin associated with the hexagon
            #mipPt_per_bin = row['ts_mipPt'] / num_bins
            #print("mipPt per bin", mipPt_per_bin )

            # Update the mipPt value for each bin associated with the hexagon
            for bin_info in row['bins_overlapping']:
                #print(bin_info)
                bin_info['mipPt'] = row['ts_mipPt'] * bin_info['percentage_overlap']

        
        print("df_overlap AREA PERC:", df_overlap.columns)
        print("df_overlap bins overapping columns AREA PERC", df_overlap['bins_overlapping'])
        for col in df_overlap['bins_overlapping'].iloc[0][0].keys():
            print("Column bin overalp AREA PERC:", col)
        #for index, row in df_overlap.iterrows():
            #num_bins = len(row['bins_overlapping'])
            #print("number of bins", num_bins)

            #for bin_info in row['bins_overlapping']:
                #print(bin_info)

        return df_overlap
    

    #BINNING

    def sph2cart(self, eta, phi, z=322.):
        ''' Useful conversion '''
        theta = 2*np.arctan(np.exp(-eta))
        r = z / np.cos(theta)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        return x, y

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


    #MODULE SPLITTING
    #method that splits the ModSums in towers according to different algorithms and plots the event (grid + modules)

    def ModSumToTowers(self, kw, data, df_bin, algo='baseline'):                                                                                                                   
        df_overlap = self.eval_hex_bin_overlap(data,df_bin) 
        plotMS.plot_single_event_proj(kw, data, df_bin, "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/today_prova.png")

        if algo == 'baseline':
            print('baseline algorithm')
            df_centroid = self.assign_mip_pt_to_bins_centroid(df_overlap)
            plotMS.plot_grid_area_overlap(df_centroid)

        elif algo == 'area_overlap':
            print('area_overlap algorithm')
            df_area_overlap = self.assign_mipPt_to_bins_area_overlap(df_overlap) 
            plotMS.plot_grid_area_overlap(df_area_overlap) 

        else:
            raise ValueError("Invalid algorithm specified. Choose 'baseline' or 'area_overlap'.")