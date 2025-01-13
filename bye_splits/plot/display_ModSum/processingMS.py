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
import matplotlib.patches as PolygonPatch
import random
import logging
import time
import multiprocessing
from shapely.strtree import STRtree

import matplotlib.pyplot as plt
import pandas as pd


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
import geometryMS
import matplotlib
from scipy.spatial.distance import cdist
from shapely import geometry as geom
from shapely.geometry import Polygon, mapping, shape, MultiPolygon, Point
from shapely.ops import unary_union
from matplotlib.patches import Polygon as matPoly
from matplotlib.collections import PatchCollection
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
import uproot
import random

class Processing():
    def __init__(self):
        with open(params.CfgPath, "r") as afile:
            self.cfg = yaml.safe_load(afile)

        self.ds_geom = GeometryData(reprocess=False, logger=log, library='plotly').provide()
        self.filename = None
        self.list_events = None
        self.algorithms = algosMS.Algorithms()
        self.resolution = resolutionMS.Resolution()
        self.geometry= geometryMS.Geometry()


    def get_gen_particles(self, root_file_path, particle, n=None, event=None):
        print("Read root file and create DF for generated particles")

        # Open the ROOT file and navigate to the TTree
        with uproot.open(root_file_path) as file:
            tree = file["l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple"]

            if particle == 'jets':
                # Extract the necessary branches and rename them
                branches = {
                    "good_genjet_eta": "gen_eta",
                    "good_genjet_phi": "gen_phi",
                    "good_genjet_energy": "gen_en",
                    "good_genjet_pt": "gen_pt",
                    "event": "event",
                    "jet_flag":"jet_flag",
                }

            else:
                # Extract the necessary branches and rename them
                branches = {
                    "good_genpart_exeta": "gen_eta",
                    "good_genpart_exphi": "gen_phi",
                    "good_genpart_energy": "gen_en",
                    "good_genpart_pt": "gen_pt",
                    "event": "event"
                }

            # Load branches into a dictionary of NumPy arrays
            arrays = tree.arrays(list(branches.keys()), library="np")

            # Flatten arrays and create a DataFrame
            data = {new_name: np.concatenate(arrays[old_name]) for old_name, new_name in branches.items() if old_name != "event"}
            #data["event"] = np.repeat(arrays["event"], [len(arr) for arr in arrays["good_genpart_exeta"]])

            #df = pd.DataFrame(data)
            #baseline_selections = (df['gen_eta'] > 1.7) & (df['gen_eta'] < 2.8)
            if particle == 'jets':
                data["event"] = np.repeat(arrays["event"], [len(arr) for arr in arrays["good_genjet_eta"]])
                df = pd.DataFrame(data)
                baseline_selections = (df['gen_eta'] > 1.7) & (df['gen_eta'] < 2.8) & (df['jet_flag']>=0)
            else:
                data["event"] = np.repeat(arrays["event"], [len(arr) for arr in arrays["good_genpart_exeta"]])
                df = pd.DataFrame(data)
                baseline_selections = (df['gen_eta'] > 1.7) & (df['gen_eta'] < 2.8)

            df = df[baseline_selections]

            # Ensure necessary columns are present in the dataframe
            '''if 'gen_eta' in df.columns and 'gen_phi' in df.columns and 'gen_pt' in df.columns:
                # Plot gen_eta distribution
                plt.figure(figsize=(12, 4))

                plt.subplot(1, 3, 1)
                plt.hist(df['gen_eta'], bins=30, color='#4682B4', alpha=0.7)
                plt.xlabel('gen_eta')
                plt.ylabel('Count')
                plt.title('Distribution of gen_eta')

                # Plot gen_phi distribution
                plt.subplot(1, 3, 2)
                plt.hist(df['gen_phi'], bins=30, color='#4682B4', alpha=0.7)
                plt.xlabel('gen_phi')
                plt.ylabel('Count')
                plt.title('Distribution of gen_phi')

                # Plot gen_phi distribution
                plt.subplot(1, 3, 3)
                print("LEN PT", len(df['gen_pt']))
                plt.hist(df['gen_pt'], bins=30, color='#4682B4', alpha=0.7)
                plt.xlabel('gen_pt')
                plt.ylabel('Count')
                plt.title('Distribution of gen_pt')
                #plt.yscale('log')

                plt.tight_layout()
                plt.show()
            else:
                print("Columns 'gen_eta' and 'gen_phi' are required for plotting.")

            #print("ciao", df)'''

            # Get unique events
            unique_events = df['event'].unique()
            print("unique_events", unique_events)
            print("unique_events", len(unique_events))

            # Event selection logic
            if event is None:
                # Randomly choose one event if none is specified
                selected_event = random.choice(unique_events)
                print(f"Randomly selected event: {selected_event}")
                selected_events = [selected_event]
            elif event == '-1':
                # Include all events if event == '-1'
                selected_events = unique_events
                if n is not None:
                    # Limit to the first n events if specified
                    selected_events = unique_events[:n]
                    print(f"Selected first {n} events.")
            else:
                # Process specific event
                event = int(event)  # Ensure the event ID is an integer
                if event not in unique_events:
                    raise ValueError(f"Event {event} not found in the ROOT file.")
                print(f"Selected event: {event}")
                selected_events = [event]

            print("selected_events qui", selected_events)
            # Filter the DataFrame to only include the selected events
            filtered_df = df[df['event'].isin(selected_events)]

            filtered_df = filtered_df.reset_index(drop=True)

            # Return the filtered DataFrame and list of selected events
            return filtered_df, selected_events

    def read_root_and_create_dataframe(self, root_file_path, subdet, selected_events):
            print("Read root file and create DF")
            # Open the ROOT file and navigate to the TTree
            with uproot.open(root_file_path) as file:
                tree = file["l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple"]

                # Extract all branches with their new names
                ts_branches = {
                    "good_ts_layer": "ts_layer",
                    "good_ts_mipPt": "ts_mipPt",
                    "good_ts_pt": "ts_pt",
                    "good_ts_waferu": "ts_wu",
                    "good_ts_waferv": "ts_wv",
                    "good_ts_x": "ts_x",
                    "good_ts_y": "ts_y",
                    "good_ts_z": "ts_z",
                    "good_ts_eta": "ts_eta",
                    "good_ts_phi": "ts_phi",
                    "good_ts_energy": "ts_energy",
                    "good_ts_subdet": "ts_subdet",
                    "event": "event"
                }

                # Load branches into a dictionary of NumPy arrays
                ts_arrays = tree.arrays(list(ts_branches.keys()), library="np")

                # Process ts_* data : Flatten arrays and create a DataFrame
                ts_data = {new_name: np.concatenate(ts_arrays[old_name]) for old_name, new_name in ts_branches.items() if old_name != "event"}
                ts_data["event"] = np.repeat(ts_arrays["event"], [len(arr) for arr in ts_arrays["good_ts_layer"]])
                ts_df = pd.DataFrame(ts_data)

                # Check if events should be filtered
                if selected_events is not None:
                    # If events are specified, filter the DataFrame to only include those events
                    ts_df = ts_df[ts_df['event'].isin(selected_events)]
                else:
                    #df = df #[df['event']==493403]
                    #print("EVENTS NEUTRINOS", df['event'])
                    unique_events_count = ts_df['event'].nunique()
                    print("Number of unique events:", unique_events_count)

                # Apply baseline selections for ts_df (you can modify this selection as needed)
                baseline_selections_ts = (ts_df['ts_z'] > 0) & (ts_df['ts_mipPt'] > 0.5)  # Example: Ensure positive z position
                ts_df = ts_df[baseline_selections_ts]

                # Reset the index after filtering
                ts_df = ts_df.reset_index(drop=True)

                # Process the events if needed (you can adapt this step)
                ts_df = self.process_event_V16(ts_df, subdet)

                # Return both DataFrames (ts_df and tc_df)
                return ts_df


    

    def random_event(self, f):
        return random.choice(self.list_events)
    



    def filestore(self,geom=None, particle=None):
        if geom =='V11':
            self.filename = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/DoublePhotonsPU0_hadd_123_energy/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"
        elif geom =='V16' and particle == 'photons':
            self.filename = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/SinglePhotonPU0V16/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"
        elif geom =='V16' and particle == 'pions':
            print("Considering pions data")
            self.filename = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/SinglePionPU0V16/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"
        else:
            self.filename = None  # Or raise an error if geom is required
            print("Unrecognized geom value. Please provide 'V11' or 'V16'.")
            #"/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/data/new_OKAY/fill_prova_SEL_all_REG_Si_SW_1_SK_default_CA_min_distance_NEV_100.hdf5"

    def get_data_new(self, event=None, n=None, geom=None, subdet=None, particle=None):
        # Check if the HDF5 file exists
        self.filestore(geom, particle)
        if not os.path.exists(self.filename):
            raise FileNotFoundError(f"HDF5 file '{self.filename}' not found.")
        # Read the list of events from the HDF5 file
        with h5py.File(self.filename, 'r') as file:
            self.list_events = list(set([key.split('_')[1] for key in file.keys()]))
        #print(f"List of events in HDF5 file: {self.list_events}")   
        print(f"Number of events in HDF5 file:",  len(self.list_events)) 

        # If no event is specified, choose a random one
        if event is None:
            events_to_process = random.choice(self.list_events)
            print(f"Selected random event: {events_to_process}")
            processed_event_df = self.get_event(events_to_process,geom, subdet)

        elif event == '-1':
            events_to_process = self.list_events
            if n is not None:
                print("Number of selected events (n)", n)
                events_to_process = events_to_process[:n]
            #print(f"Processing events: {events_to_process}")
            processed_event_df = self.get_allevents(events_to_process, geom, subdet)

        else:
            print(f"Selected event: {event}")
            # If the specified event is not in the list, raise an exception
            if str(event) not in self.list_events:
                raise ValueError(f"Event {event} not found in the HDF5 file.")
            processed_event_df = self.get_event(event,geom, subdet)
            events_to_process = event
        return processed_event_df, events_to_process

    def get_event(self, event, geom, subdet):
        with h5py.File(self.filename, 'r') as file:
            for key in file.keys():
                if event in key and 'ts' in key:
                    dataset = file[key]
                    column_names = [str(col) for col in dataset.attrs['columns']]
                    print("column_names", column_names)
                    # Convert the dataset to a Pandas DataFrame with proper column names
                    df_ts = pd.DataFrame(dataset[:], columns=column_names)
                    df_ts['event'] = int(event)
                    #print("df_ts STOP", df_ts.columns)
                    #print("df_ts LAYERS", sorted(df_ts['ts_layer'].unique()) )
                    #filtered_df_sub3 = df_ts[df_ts['ts_subdet'] == 3]
                    #unique_layers_subdet_3 = sorted(filtered_df_sub3['ts_layer'].unique())
                    #print("Unique Layers for  == 3:", unique_layers_subdet_3)

                    #print("df_ts ETA", df_ts['ts_eta'].unique())
                    #print("df_ts values", df_ts)
                    print("DF_TS_TODAY", df_ts)
                    print("DF_TS_TODAY", df_ts.columns)
                    if geom == 'V11':
                        silicon_df= self.process_event(df_ts, subdet)
                        print("silicon_df_proc V11", silicon_df.columns)
                    if geom == 'V16':
                        silicon_df= self.process_event_V16(df_ts, subdet)
                        print("silicon_df_V16", silicon_df.columns)
                        print("silicon_df_V16", silicon_df)

                    return silicon_df
    
    def get_allevents(self, events_to_process, geom, subdet):
        all_ts_event_dfs = []
        all_tc_event_dfs = []

        event_key_map = {}

        with h5py.File(self.filename, 'r') as file:
            for event in events_to_process:
                for key in file.keys():
                    if event in key:

                        # Track keys associated with the current event
                        if event not in event_key_map:
                            event_key_map[event] = []
                        event_key_map[event].append(key)

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

        #combined_tc_df = pd.concat(all_tc_event_dfs, ignore_index=True)
        #print("tc data dataframe", combined_tc_df.columns)
        
        # Process the combined ts DataFrame
        if geom == 'V11':
            processed_combined_ts_df = self.process_event(combined_ts_df, subdet) #FIXME - Temporary adding tc data (combined_ts_df,combined_tc_df)
        if geom == 'V16': 
            processed_combined_ts_df = self.process_event_V16(combined_ts_df, subdet)
            lenght = len(processed_combined_ts_df)
            print("lenght df",lenght)
            print("output DF V16 process event", processed_combined_ts_df.columns)
        return processed_combined_ts_df
    
    def get_genpart_data(self, file_path, block0_value, events_to_process, n, group_name='df'):
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
                        if n is not None:
                            # Convert events_to_process to integers for comparison
                            matching_indices = [i for i, v in enumerate(block0_values) if str(v) in events_to_process]
                            print("block0_values inside", block0_values)
                            for idx in matching_indices:
                                block1_data.append(block1_values[idx])
                                events.append(block0_values[idx])
                            '''for idx, value in enumerate(events_to_process):
                                block1_data.append(block1_values[idx])
                                events.append(value)'''
                        else:
                            # Otherwise, select all values
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
                    print("GEN PART DF", df)
                    return df, df['event']
                else:
                    print("Error: 'block0_values', 'block1_items', or 'block1_values' not found in group {}.".format(group_name))
            else:
                print("Error: Group '{}' not found in HDF5 file.".format(group_name))

    def process_event(self, df_ts):
        """
        Merging V11 data with geometry, separately for silicon and scintillator
        """
        print("Processing events V11...")

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

        #SILICON V11

        ds_si_geo = self.ds_geom['si'].drop_duplicates(subset=['ts_layer', 'ts_wu', 'ts_wv'], keep='first')

        #plotMS.plot_layers(df_ts, ds_new)
        #plotMS.plot_layers_sci(df_tc, self.ds_geom['sci'])

        # Merge dataframes
        silicon_df = pd.merge(left=df_ts, right=ds_si_geo, how='inner',
                                      on=['ts_layer', 'ts_wu', 'ts_wv'])
        silicon_df = silicon_df.drop(['triggercellu','triggercellv','waferorient', 'waferpart','diamond_x', 'diamond_y'], axis=1) #anche qui tagliare sul subdet!!

        # Shifting hexagons vertices based on difference between wx_center/wy_center (byesplit) and ts_x/ts_y (CMSSW)
        shifted_hex_df = self.shift_hex_values(silicon_df, self.ds_geom['si'], df_ts)

        #SCINTILLATOR V11

        self.geometry.implement_scint_modules_id(self.ds_geom['sci'])
        scint_mod_geom = self.geometry.create_scint_mod_geometry(self.ds_geom['sci'], save_geojson = False)
        #plotMS.plot_scint_tiles(self.ds_geom['sci'])

        sci_merge = {'ts_ieta' : 'ts_eta',
                      'ts_iphi' : 'ts_phi',
                      'tc_layer' : 'ts_layer'}

        scint_mod_geom = scint_mod_geom.rename(columns=sci_merge)

         # Merge dataframes
        '''scintillator_df = pd.merge(left=df_ts, right=scint_mod_geom, how='inner', #FIXME ma qui devo scelgiere solo IL 3° SUBDET dei dati che non HO!! devo aggiungerlo!!
                                           on=['ts_layer', 'ts_eta', 'ts_phi'])'''

        return shifted_hex_df

    def process_event_V16(self, df_ts, subdet):
        """
        Merging V16 data with geometry, separately for silicon and scintillator
        """
        print("**Processing V16 geometry ...")
        print("df_ts QUII", df_ts.columns)

        if subdet == 1 or subdet ==2:
            print(f"**Processing SUBDET {subdet}")

            #SILICON V16
            df_geo_si_V16 = self.geometry.create_si_mod_geometry_V16()
            df_si_merged = pd.merge(left=df_ts, right=df_geo_si_V16, how='inner', on=['ts_layer', 'ts_wu', 'ts_wv'])
            df_si_merged_subdet = df_si_merged[df_si_merged['ts_subdet'] == subdet]
            print("df_si_merged_subdet", df_si_merged_subdet)

            # Collect all unique ts_layers
            unique_ts_layers = sorted(df_si_merged_subdet['ts_layer'].unique())

            # Print the unique sequence of ts_layers
            print(f"Unique ts_layers for subdet {subdet}: {unique_ts_layers}")
            #plotMS.plot_layer_hexagons(df_si_merged, 11, x,y)
            print(f"df_si_merged_subdet", df_si_merged_subdet.columns)
            return df_si_merged_subdet

        elif subdet == 3:
            print(f"**Processing subdet {subdet}: CEH --> only scintillator part")

            #SCINTILLATOR V16
            sci_update = {'triggercellieta' : 'tc_cu',
                        'triggercelliphi' : 'tc_cv',
                        'layer'           : 'tc_layer',
                        'waferu'          : 'tc_wu',
                        'waferv'          : 'tc_wv'}

            self.ds_geom['sci'] = self.ds_geom['sci'].rename(columns=sci_update)
            df = self.geometry.implement_scint_modules_id(self.ds_geom['sci'])
            scint_mod_geom = self.geometry.create_scint_mod_geometry(df, save_geojson = False)

            #print("scint_mod_geom", scint_mod_geom)
            #print("scint_mod_geom_2", scint_mod_geom_2)
            #plotMS.plot_scint_tiles(self.ds_geom['sci'])

            sci_merge = {'ts_ieta' : 'ts_wu',
                        'ts_iphi' : 'ts_wv',
                        'tc_layer' : 'ts_layer'}

            scint_mod_geom = scint_mod_geom.rename(columns=sci_merge)
            df_ts_filtered = df_ts[df_ts['ts_subdet'] == 3] #filtering to select only SUBDET 3 corresponding to the scintillator part

            # Merge dataframes
            scintillator_df = pd.merge(left=df_ts_filtered, right=scint_mod_geom, how='inner',
                                            on=['ts_layer', 'ts_wu', 'ts_wv'])
            return scintillator_df

        elif subdet == 4:
            print("Processing combined subdet 2 and 3 (all CEH)")

            # Get the DataFrames for subdet 2 and 3
            df_subdet_2 = self.process_event_V16(df_ts, subdet=2)
            df_subdet_3 = self.process_event_V16(df_ts, subdet=3)

            # Combine the two DataFrames
            df_combined = pd.concat([df_subdet_2, df_subdet_3], ignore_index=True)

            #print("df_combined COL", df_combined.columns)
            #print("df_combined ALL", df_combined)

            return df_combined

        elif subdet == 5:

            print("Processing combined subdet 1, 2 and 3 (CEE+ CEH)")

            # Get the DataFrames for subdet 2 and 3
            df_subdet_1 = self.process_event_V16(df_ts, subdet=1)
            df_subdet_2 = self.process_event_V16(df_ts, subdet=2)
            df_subdet_3 = self.process_event_V16(df_ts, subdet=3)

            # Combine the two DataFrames
            df_combined = pd.concat([df_subdet_1, df_subdet_2, df_subdet_3], ignore_index=True)

            #print("df_combined COL", df_combined.columns)
            #print("df_combined ALL", df_combined)

            return df_combined

        else:
            print(f"Subdet {subdet} is not recognized.")
            return None

    def shift_hex_values(self, silicon_df_proc, df_geom, df_ts):
        """Shifting hexagons vertices based on difference between wx_center/wy_center (byesplit) and ts_x/ts_y (CMSSW),
           in order to maintain as center of the hexagon the orginal ts_x/ts_y. 
           It returns the dataframe with the applied shifts"""
        
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
    
    
    def create_and_save_tower_bins(self, kw, df_data, geom ):
        print("using df data to calculate bins", df_data.columns)
        if geom=='V16':
            max_layer = 47
        else:
            max_layer = 50
            
        phi_values = np.linspace(kw['MinPhi'], kw['MaxPhi'], kw['NbinsPhi'] + 1)
        eta_values = np.linspace(kw['EtaMin'], kw['EtaMax'], kw['NbinsEta'] + 1)

        bins_info = []

        # Iterate through each layer
        for layer in range(1, max_layer + 1):  # Loop through all layers from 1 to 50
            # Filter dataframe for the current layer based on the conditions
            if (layer <= 26 and layer % 2 != 0) or (layer >= 27 and layer <= max_layer):
                # Filter dataframe for the current layer
                layer_df = df_data[df_data['ts_layer'] == layer].copy()
                print("layer_df", layer_df.columns)

                unique_z_values = layer_df['ts_z'].unique()[:1]   #select only the first entry cause there are multiple z values for approximation reasons
                unique_layer_values = layer_df['ts_layer'].unique()[:1] 
                print("LAYER:", unique_layer_values, " Z:", unique_z_values)
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

        self.geometry.save_bin_geo(df_bins, f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_with_arcs_correct.geojson', f'/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/geojson/bins_only_vertices_correct.geojson')

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

    def apply_update_to_each_event(self, df, geojson_bins):
        print("apply update to each event")
        #print("bin df", bins_df)
        #print("bin df", bins_df.columns)

        # Load the GeoJSON file
        with open(geojson_bins) as f:  # Replace 'geojson_bins.json' with the actual file path
            bin_geojson = json.load(f)

        # Initialize lists to hold extracted data
        bins_info = []

        # Iterate over each feature in the GeoJSON
        for feature in bin_geojson['features']:
            layer = feature['properties']['Layer']
            z_values = feature['properties']['Z_value']

            # Extract geometry coordinates (X and Y vertices)
            coordinates = feature['geometry']['coordinates'][0]
            x_vertices = [coord[0] for coord in coordinates]
            y_vertices = [coord[1] for coord in coordinates]

            # Extract Eta and Phi vertices from properties
            eta_vertices = feature['properties']['Eta_vertices']
            phi_vertices = feature['properties']['Phi_vertices']

            # Prepare the layer bins data
            layer_bins = [{
                'X_Vertices': x_vertices,
                'Y_Vertices': y_vertices,
                'Eta_Vertices': eta_vertices,
                'Phi_Vertices': phi_vertices
            }]

            # Append to bins_info
            bins_info.append({
                'Layer': layer,
                'Z_value': z_values,
                'Bins': layer_bins
            })

        # Create DataFrame from the bins_info list
        df_bins = pd.DataFrame(bins_info)
        #print("df_ORA",df_bins )

        # Group the DataFrame by 'Layer' and combine 'Bins' by concatenating lists
        grouped_bins = df_bins.groupby('Layer').agg({
            'Z_value': 'first',  # Retain the first Z_value for each layer
            'Bins': lambda x: sum(x, [])  # Concatenate lists of dictionaries in the 'Bins' column
        }).reset_index()

        #print("df_ORA",grouped_bins["Bins"] )





        #df bins qua da non uare, leggere goejson e ricreare stesso df
        all_event_dfs = []

        unique_events = df['event'].unique()

        for event in unique_events:
            #print(f"Processing event {event}...")
            
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            #Add bins that have 0 overlap with the hexagons with pt 0
            updated_event_df = self.update_bins(event_df, grouped_bins)

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
            
    def ModSumToTowers(self, kw, data, subdet, event, particle, algo, bin_geojson_filename, hex_geojson_filename, data_gen, geom):
        print('Mod sum to towers')

        #print("data", data.columns)
        #print("data gen", data_gen.columns)

        #print("ts_en", "gen_pt")
        #print(data['ts_pt'].sum(), data_gen['gen_pt'])

        #self.check_pt(data, data_gen)
        #print("data pt TS", data['ts_pt'].sum())
        #print("data gen pt", data_gen['gen_en'])


        #overlap_data = self.read_hdf5_file(hdf5_filename)
        #print("DATA OGGI ", data)
        #print("OVERLAP DATA INPUT columns ", overlap_data.columns)

        #hexagon_info_df = self.eval_hex_bin_overlap_scint(data, bin_geojson_filename, geom)

        #hexagon_info_df = self.eval_hex_bin_overlap_OK(data, bin_geojson_filename, geom)
        filename_precomputed_silicon = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/hex_bin_precomputed_overlaps_corrected.json"
        filename_precomputed_scint = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/tiles_bin_precomputed_overlaps_corrected.json"
        hexagon_info_df= self.eval_hex_bin_overlap_with_precomputed_jsons(data,filename_precomputed_silicon, filename_precomputed_scint ,geom)
        #print("eval overlap dataframe", hexagon_info_df.columns)
        print(hexagon_info_df.index.get_level_values('event').unique())

        #plotMS.plot_hex_bins(hexagon_info_df)

        #print("QUIIIIII, hexagon_info_df COLUMNS", hexagon_info_df.columns)

        if algo == 'baseline':
            df_algo = self.algorithms.baseline_by_event(hexagon_info_df, subdet)
            
        elif algo == 'area_overlap':
            df_algo = self.algorithms.area_overlap_by_event(hexagon_info_df, subdet)

        elif algo == '8towers':
            df_algo = self.algorithms.area_overlap_8Towers_by_event(hexagon_info_df, subdet)
        
        elif algo == '16towers':
            df_algo = self.algorithms.area_overlap_16Towers_by_event(hexagon_info_df, subdet)

        else:
            raise ValueError("Invalid algorithm specified. Choose 'baseline', 'area_overlap' or '8towers'.")

        df , df_sum = self.apply_update_to_each_event(df_algo, bin_geojson_filename)

        #print("Eta/Phi resolution evaluation...")
        #print("df", df)
        #print("df columns", df.columns)
        #print("data_gen", data_gen)
        #print("data_gen columns", data_gen.columns)

        if particle == "pions" or particle == "jets":
            results_df, jets = self.resolution.perform_clustering_antikt_matched(df, data_gen, f'{algo}_{particle}_{event}_{subdet}_PIONS_results.txt')
            #print(results_df)

        elif particle == "photons":
            results_df = self.resolution.eval_eta_phi_photon_resolution(df, data_gen, algo, subdet, window_size=12, subwindow_size=9)
            self.resolution.save_eta_phi_differences(results_df, f'{algo}_{particle}_{event}_{subdet}_eta_phi_resolution_12w_9sub.txt')

        elif particle == "neutrinos":
            print("Sei arrivato!---antikt for neutrinos")
            results_df, jet_counts = self.resolution.perform_clustering_antikt(df,f'{algo}_{particle}_{event}_{subdet}_results_2Ntuples.txt' )
            #print("RESULTS NEUTRINO",results_df)
            #print("COUNTS",jet_counts)


        #plotMS.plot_energy_ratio_histogram()
        #plotMS.plot_eta_phi_resolution(df_resolution, algo, event, particle, subdet)
        #print("data_gen", data_gen.columns)
        #print("df_sum", df_sum.columns)
        #print("event", event)
        #plotMS.plot_towers_eta_phi_grid(df_sum, data_gen, algo, event, particle, subdet, results_df)
        #plotMS.plot_towers_xy_grid(df_sum, data_gen, algo, event, particle, subdet)

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
        closest_bins = [bin_dist[0] for bin_dist in bin_distances[:35]] #takes total of 1.8 sec per event #35
        
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

    def eval_hex_bin_overlap_with_precomputed_jsons(self, data, precomputed_json_file, scintillator_json_file, geom):
        """
        Evaluate overlaps between hexagons and bins using precomputed JSON data.

        :param data: DataFrame containing hexagon information.
        :param silicon_json_file: Path to the JSON file with precomputed overlaps for subdet 1 and 2.
        :param scintillator_json_file: Path to the JSON file with precomputed overlaps for subdet 3.
        :return: A DataFrame with hexagon and bin overlap information.
        """
        print("Loading precomputed overlap data...")

        # Load the two precomputed JSON files
        with open(precomputed_json_file) as f:
            precomputed_data = json.load(f)

        with open(scintillator_json_file) as f:
            scintillator_data = json.load(f)

        # Index both precomputed data files
        print("Indexing precomputed data...")
        precomputed_index = {
            (entry['hex_layer'], entry['hex_wu'], entry['hex_wv']): entry
            for entry in precomputed_data
        }

        scintillator_index = {
            (entry['hex_layer'], entry['hex_wu'], entry['hex_wv']): entry
            for entry in scintillator_data
        }

        print("Processing hexagon data...")
        hexagon_info = []

        for hex_row in data.itertuples():
            # Determine which index to use based on the subdet
            if hex_row.ts_subdet in [1, 2]:
                hex_key = (hex_row.ts_layer, hex_row.ts_wu, hex_row.ts_wv)
                precomputed_entry = precomputed_index.get(hex_key)
            elif hex_row.ts_subdet == 3:
                hex_key = (hex_row.ts_layer, hex_row.ts_wu, hex_row.ts_wv)
                precomputed_entry = scintillator_index.get(hex_key)
            else:
                continue  # Skip if subdet is not in the expected range

            if precomputed_entry is not None:  # If a match is found
                hex_info = {
                    'event': hex_row.event,
                    'layer': hex_row.ts_layer,
                    'hex_x': hex_row.hex_x,
                    'hex_y': hex_row.hex_y,
                    'hex_eta_centroid': precomputed_entry['hex_eta_centroid'],
                    'hex_phi_centroid': precomputed_entry['hex_phi_centroid'],
                    'ts_pt': hex_row.ts_pt,
                    'bins_overlapping': []
                }

                # Add overlapping bin details
                for bin_overlap in precomputed_entry['overlapping_bins']:
                    bin_info = {
                        'percentage_overlap': bin_overlap['percentage_overlap'],
                        'bin_layer': bin_overlap['bin_layer'],
                        'centroid_eta': bin_overlap['bin_eta_centroid'],
                        'centroid_phi': bin_overlap['bin_phi_centroid'],
                        'eta_vertices': bin_overlap['bin_properties']['Eta_vertices'],
                        'phi_vertices': bin_overlap['bin_properties']['Phi_vertices']
                    }
                    hex_info['bins_overlapping'].append(bin_info)

                hexagon_info.append(hex_info)

        # Convert hexagon_info into a structured DataFrame
        print("Converting to DataFrame...")
        df_hexagon_info = pd.DataFrame(hexagon_info)
        df_hexagon_info.set_index(['event', 'layer'], inplace=True)

        return df_hexagon_info