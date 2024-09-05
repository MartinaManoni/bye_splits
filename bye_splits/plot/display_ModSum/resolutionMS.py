_all_ = [ ]

import os
import sys
import plotly.express.colors as px
import logging
from data_handle.data_process import *
import pandas as pd
import numpy as np

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)
log = logging.getLogger(__name__)

class Resolution():
    def __init__(self):
        pass

    def find_particle_bin_and_evaluate_windows(self, baseline_df, genpart_df, window_size=8, subwindow_size=5):
        particle_eta = genpart_df['gen_eta'].iloc[0]
        particle_phi = genpart_df['gen_phi'].iloc[0]
        event_number = genpart_df['event'].iloc[0]
        print( particle_eta,particle_phi, event_number)
        print("gen part dataframe",  genpart_df)
        print("baseline_df dataframe",  baseline_df)

        # Find the bin where the generated particle is located
        baseline_df = baseline_df.copy()
        baseline_df['eta_center'] = baseline_df['eta_vertices'].apply(lambda x: np.mean(x))
        baseline_df['phi_center'] = baseline_df['phi_vertices'].apply(lambda x: np.mean(x))

        # Find the bin with the minimum distance to the particle's eta and phi
        particle_bin_idx = ((baseline_df['eta_center'] - particle_eta).abs() + (baseline_df['phi_center'] - particle_phi).abs()).idxmin()
        particle_bin = baseline_df.loc[particle_bin_idx]

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
    
    def eval_eta_phi_photon_resolution(self, df, genpart_df, window_size=8, subwindow_size=5):
        all_results = []

        unique_events = df['event'].unique()
        print("unique_events", unique_events)

        for event in unique_events:
            print("event", event)
            
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            event_genpart_df = genpart_df[genpart_df['event'] == event]
            print("event_genpart_df", event_genpart_df)

            # Apply the find_particle_bin_and_evaluate_windows function
            subwindow_energies, eta_diff, phi_diff = self.find_particle_bin_and_evaluate_windows(event_df, event_genpart_df, window_size, subwindow_size)

            # Store the results
            result = {
                'event': event,
                'subwindow_energies': subwindow_energies,
                'eta_diff': eta_diff,
                'phi_diff': phi_diff
            }
            all_results.append(result)

        return pd.DataFrame(all_results)
    
    def check_pt(self, data, data_gen):
        print("check pt")
        # Initialize a list to store ratios
        ratios = []

        # Group data by events
        data_grouped = data.groupby('event')
        data_gen_grouped = data_gen.groupby('event')

        for event_id in data_grouped.groups.keys():
            event_data = data_grouped.get_group(event_id)
            event_gen_data = data_gen_grouped.get_group(event_id)

            # Ensure there is only one gen entry per event for simplicity
            if len(event_gen_data) == 1:
                gen_pt = event_gen_data['gen_pt'].values[0]
                event_data_ts_pt_sum = event_data['ts_pt'].sum()

                ratio = event_data_ts_pt_sum / gen_pt
                ratios.append(ratio)

                if ratio > 3:
                    print(f"Event {event_id} has ratio > 3")
                    print(f"ts_pt: {event_data_ts_pt_sum}, gen_pt: {gen_pt}")
                    print(f"gen_phi: {event_gen_data['gen_phi'].values[0]}, gen_eta: {event_gen_data['gen_eta'].values[0]}")

        # Plot histogram of ratios
        bins = np.linspace(0, 10, 11)  # 10 bins between 0 and 10
        plt.hist(ratios, bins=bins, edgecolor='black')
        plt.xlabel('ts_pt / gen_pt')
        plt.ylabel('Frequency')
        plt.title('Histogram of ts_pt/gen_pt Ratios per Event')
        plt.xlim(0, 10)  # Set x-axis limit between 0 and 10
        plt.show()
    
    def save_eta_phi_differences(self, results_df, filename):
        eta_diffs = results_df['eta_diff']
        phi_diffs = results_df['phi_diff']
        
        with open(filename, 'w') as file:
            file.write("eta_diff,phi_diff\n")
            for eta_diff, phi_diff in zip(eta_diffs, phi_diffs):
                file.write(f"{eta_diff},{phi_diff}\n")