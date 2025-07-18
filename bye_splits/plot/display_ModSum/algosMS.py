_all_ = [ ]

import os
import sys
import logging
import time
import pandas as pd
import numpy as np
import json
from data_handle.data_process import *

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)
log = logging.getLogger(__name__)

class Algorithms():
    def __init__(self):
        pass

    def baseline_by_event(self, df_hexagon_info, subdet):
        print("Baseline algorithm ...")
        start_time = time.time()
        
        # Initialize a list to store the results for all events
        all_event_rows = []

        # Get unique event indices
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        # Iterate over each event
        for event in unique_events:
            #print(f"Processing event {event}...")
            
            # Extract the DataFrame for the current event
            df_event = df_hexagon_info.loc[event]

            # Initialize a dictionary to store the summed pt for each bin by layer for the current event
            bin_pt_by_layer = {}

            # Get unique layer indices for the current event
            unique_layers = df_event.index.get_level_values('layer').unique()

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
            #print(f"Processing event {event}...")
            df_event = df_hexagon_info.loc[event]
            bin_pt_by_layer = {}
            unique_layers = df_event.index.get_level_values('layer').unique()

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


    def splitted_MS_Towers_by_event(self, df_hexagon_info, num_towers = 8):
        print(f"Splitting Mod Sums into {num_towers} towers ...")
        start_time = time.time()

        # Initialize a list to store the results for all events
        all_event_results = []

        # Get unique event indices
        unique_events = df_hexagon_info.index.get_level_values('event').unique()

        # Iterate over each event
        for event in unique_events:

            # Extract the DataFrame for the current event
            df_event = df_hexagon_info.loc[event]

            # Get unique layer indices for the current event
            unique_layers = df_event.index.get_level_values('layer').unique()

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
                    top_bins = sorted_bins[:num_towers] if num_bins > num_towers  else sorted_bins
                    for bin_info in top_bins:
                        percent = round(bin_info['percentage_overlap'] * num_towers )
                        percent = max(1, percent)  # Ensure at least 1
                        pt_fraction = percent / num_towers 
                        pt_assigned = min(remaining_pt, pt_fraction * total_pt)
                        bin_info['pt'] = pt_assigned
                        remaining_pt -= pt_assigned
                        if remaining_pt <= 0:
                            break  # Stop assigning pt if all available pt is exhausted

                    # Append results for the current layer to the event_results list
                    all_event_results.extend(row['bins_overlapping'])

        # Calculate the execution time
        end_time = time.time()
        print("Execution time loop:", end_time - start_time)

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

        print("Execution time loop:", end_time - start_time)

        return df_over_final

    # Function to load bin geometry from a JSON file
    def load_bins_from_json(bins_json_path):
        print("load bins from json")
        with open(bins_json_path, "r") as f:
            bins_data = json.load(f)

        # Parse bin information
        bins = []
        for feature in bins_data["features"]:
            eta_vertices = feature["properties"]["Eta_vertices"]
            phi_vertices = feature["properties"]["Phi_vertices"]
            centroid_eta = np.mean(eta_vertices)
            centroid_phi = np.mean(phi_vertices)
            layer = feature["properties"]["Layer"]
            bins.append({
                "eta_vertices": eta_vertices,
                "phi_vertices": phi_vertices,
                "centroid_eta": centroid_eta,
                "centroid_phi": centroid_phi,
                "layer": layer
            })

        return pd.DataFrame(bins)

    # Function to assign points to bins layer by layer
    def assign_points_to_closest_bin_layer_by_layer(df, bins_df):
        print("assign_points_to_closest_bin_layer_by_layer")
        results = []

        # Process by event and tc_layer
        for (event, layer), group in df.groupby(["event", "tc_layer"]):
            #print(event, layer)
            # For each point in this group (layer), find the closest bin
            for _, point in group.iterrows():
                # Calculate the distance to each bin's centroid
                distances = np.sqrt(
                    (bins_df["centroid_eta"] - point["tc_eta"]) ** 2 +
                    (bins_df["centroid_phi"] - point["tc_phi"]) ** 2
                )

                # Find the closest bin (the one with the minimum distance)
                closest_bin_idx = distances.idxmin()
                closest_bin = bins_df.loc[closest_bin_idx]

                # Record result for this point, assign energy of the point to the closest bin
                results.append({
                    "event": event,
                    "layer": layer,
                    "eta_vertices": closest_bin["eta_vertices"],
                    "phi_vertices": closest_bin["phi_vertices"],
                    "pt": point["tc_pt"]  # Assign energy of the point to the closest bin
                })

        return pd.DataFrame(results)

    # Function to group by event, eta_vertices, and phi_vertices, and sum the pt values
    def group_and_sum_pts(final_df):
        print("group_and_sum_pts")
        # Group by event, eta_vertices, and phi_vertices and sum the pt values
        # Convert eta_vertices and phi_vertices from lists to tuples (hashable)
        final_df['eta_vertices'] = final_df['eta_vertices'].apply(tuple)
        final_df['phi_vertices'] = final_df['phi_vertices'].apply(tuple)

        # Group by event, eta_vertices_tuple, and phi_vertices_tuple and sum the pt values
        df_baseline = final_df.groupby(['event', 'eta_vertices', 'phi_vertices']).agg({'pt': 'sum'}).reset_index()

        return df_baseline


    # Main function that integrates all steps
    def process_and_assign_points_to_bins(self,df_STCs, bins_json_path):
        # Load bins geometry from JSON
        bins_df = Algorithms.load_bins_from_json(bins_json_path)

        # Assign points to bins layer by layer
        final_df = Algorithms.assign_points_to_closest_bin_layer_by_layer(df_STCs, bins_df)

        # Group by event, eta_vertices, and phi_vertices, and sum the pt values
        df_baseline = Algorithms.group_and_sum_pts(final_df)

        # Return the final DataFrame with summed pt values
        return df_baseline
