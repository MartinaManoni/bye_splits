_all_ = [ ]

import os
import sys
import plotly.express.colors as px
import logging
from data_handle.data_process import *
import pandas as pd
import numpy as np
import plotMS
import fastjet

parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)
log = logging.getLogger(__name__)

class Resolution():
    def __init__(self):
        pass

    def find_particle_bin_and_evaluate_windows(self, baseline_df, genpart_df, algo, subdet, window_size=12, subwindow_size=9):
        particle_eta = genpart_df['gen_eta'].iloc[0]
        particle_phi = genpart_df['gen_phi'].iloc[0]
        event_number = genpart_df['event'].iloc[0]

        window_bins_part1 = None
        window_bins_part2 = None    
        #print( particle_eta,particle_phi, event_number)
        #print("gen part dataframe",  genpart_df)
        #print("baseline_df dataframe",  baseline_df)

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

        # Wrap-around for phi_min and phi_max
        phi_min = (phi_min + np.pi) % (2 * np.pi) - np.pi
        phi_max = (phi_max + np.pi) % (2 * np.pi) - np.pi



        if phi_min > phi_max:
            # Part 1: From phi_min to +pi (wrap around)
            window_bins_part1 = baseline_df[
                (baseline_df['eta_center'] >= eta_min) & (baseline_df['eta_center'] <= eta_max) &
                (baseline_df['phi_center'] >= phi_min) & (baseline_df['phi_center'] <= np.pi)
            ]
            # Part 2: From -pi to phi_max (wrap around)
            window_bins_part2 = baseline_df[
                (baseline_df['eta_center'] >= eta_min) & (baseline_df['eta_center'] <= eta_max) &
                (baseline_df['phi_center'] >= -np.pi) & (baseline_df['phi_center'] <= phi_max)
            ]
            window_bins = pd.concat([window_bins_part1, window_bins_part2], ignore_index=True)

        else:

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

                # Wrap-around logic for phi_start and phi_end
                phi_start = (phi_start + np.pi) % (2 * np.pi) - np.pi
                phi_end = (phi_end + np.pi) % (2 * np.pi) - np.pi

                # Select subwindow, handling the wrap-around case
                if phi_start > phi_end:
                    # Wrap-around case: Split the selection into two parts
                    subwindow_part1 = window_bins[
                        (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] <= eta_end) &
                        (window_bins['phi_center'] >= phi_start) & (window_bins['phi_center'] <= np.pi)
                    ]
                    subwindow_part2 = window_bins[
                        (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] <= eta_end) &
                        (window_bins['phi_center'] >= -np.pi) & (window_bins['phi_center'] <= phi_end)
                    ]
                    subwindow = pd.concat([subwindow_part1, subwindow_part2], ignore_index=True)

                else:
                    # Normal case
                    subwindow = window_bins[
                        (window_bins['eta_center'] >= eta_start) & (window_bins['eta_center'] <= eta_end) &
                        (window_bins['phi_center'] >= phi_start) & (window_bins['phi_center'] <= phi_end)
                    ]


                if len(subwindow) == subwindow_size * subwindow_size:
                    total_pt = subwindow['pt'].sum()
                    subwindow_energies.append(total_pt)

                    particle_eta_1 = genpart_df['gen_eta'].iloc[0]
                    particle_phi_1 = genpart_df['gen_phi'].iloc[0]

                    # Inside the loop or function body
                    #if window_bins_part1 is not None and window_bins_part2 is not None:
                        #plotMS.plot_window_with_wraparound(window_bin=None, window_bins_part1=window_bins_part1, window_bins_part2=window_bins_part2, eta_min= eta_min,eta_max= eta_max, phi_min= phi_min, phi_max=phi_max, particle_eta_1=particle_eta_1, particle_phi_1=particle_phi_1, eta_start=eta_start, eta_end=eta_end, phi_start=phi_start, phi_end=phi_end)
                    
                    #elif window_bins is not None:
                        #plotMS.plot_window_with_wraparound(window_bins= window_bins, window_bins_part1=None, window_bins_part2=None, eta_min= eta_min,eta_max= eta_max, phi_min= phi_min, phi_max=phi_max, particle_eta_1=particle_eta_1, particle_phi_1=particle_phi_1, eta_start=eta_start, eta_end=eta_end, phi_start=phi_start, phi_end=phi_end)
                    
                    if total_pt > max_pt:
                        max_pt = total_pt
                        best_subwindow = subwindow
                        best_subwindow_data = {
                        'subwindow': best_subwindow,
                        'total_pt': max_pt
                    }

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

                # Handle phi wrap-around using complex numbers

                phi_values = best_subwindow['phi_center']
                weights = best_subwindow['pt']
                # Convert phi values to complex numbers: e^(i*phi)
                complex_phi = np.exp(1j * phi_values)
                # Calculate the weighted mean in the complex plane
                weighted_complex_mean = (weights * complex_phi).sum() / weights.sum()
                # Convert the weighted mean back to an angle
                weighted_phi = np.angle(weighted_complex_mean)

        # Calculate the difference
        eta_diff = weighted_eta - particle_eta
        phi_diff_0 = weighted_phi - particle_phi
        phi_diff = (phi_diff_0 + np.pi) % (2 * np.pi) - np.pi

        # Check for infinite or NaN values
        if not np.isfinite(eta_diff) or not np.isfinite(phi_diff):
            print("Warning: Infinite or NaN values detected in eta_diff or phi_diff. Handling edge case.")
            # Handle edge case (for example, set to zero or NaN)
            eta_diff = 0.0
            phi_diff = 0.0  # or phi_diff = np.nan

        # Save required data to files
        genpart_pt = genpart_df['gen_pt'].iloc[0]
        pt_ratio = best_subwindow_pt / genpart_pt if genpart_pt != 0 else np.nan

        # Check if the file is empty or doesn't exist
        if not os.path.exists(f'pt_scale_&_res_{algo}_{subdet}.txt') or os.stat(f'pt_scale_&_res_{algo}_{subdet}.txt').st_size == 0:
            with open(f'pt_scale_&_res_{algo}_{subdet}.txt', 'a') as f4:
                # Write the header
                f4.write("best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, weighted_eta,weighted_phi, event_number \n")

        with open(f'genpart_pt_{algo}_{subdet}.txt', 'a') as f1, open(f'best_subwindow_pt_overlap_{algo}_{subdet}.txt', 'a') as f2, open(f'pt_ratio_overlap_{algo}_{subdet}.txt', 'a') as f3, open(f'pt_scale_&_res_{algo}_{subdet}.txt', 'a') as f4:
            f1.write(f"{genpart_pt}\n")
            f2.write(f"{best_subwindow_pt}\n")
            f3.write(f"{pt_ratio}\n")
            f4.write(f"{best_subwindow_pt}, {genpart_pt}, {pt_ratio}, {eta_diff}, {phi_diff}, {particle_eta}, {particle_phi}, {weighted_eta}, {weighted_phi}, {event_number}\n")

        return subwindow_energies, eta_diff, phi_diff
    
    def eval_eta_phi_photon_resolution(self, df, genpart_df, algo, subdet, window_size=12, subwindow_size=9):
        all_results = []
        print("eval_eta_phi_photon_resolution")
        print("df", df)
        print("genpart_df", genpart_df)

        unique_events = df['event'].unique()
        print("unique_events", unique_events)

        for event in unique_events:
            #print("event", event)
            
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            #print("genpart_df", genpart_df)
            genpart_df['event'] = genpart_df['event'].astype(int)
            event_genpart_df = genpart_df[genpart_df['event'] == event]
            #print("event_genpart_df", event_genpart_df)

            #print(df['event'].dtype)
            #print(genpart_df['event'].dtype)


            # Apply the find_particle_bin_and_evaluate_windows function
            subwindow_energies, eta_diff, phi_diff = self.find_particle_bin_and_evaluate_windows(event_df, event_genpart_df, algo, subdet, window_size, subwindow_size)

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

    def perform_clustering_antikt_matched(self, df, genpart_df, ouptput_txt):
        """
        This function clusters particle data using the Anti-kt algorithm and matches reconstructed jets to generated particles.

        Key steps:
        1. Construct pseudo-jets from particle data and perform clustering (radius 0.4).
        2. Match reconstructed jets to generated particles within a DeltaR < 0.1.
        3. Calculate jet resolutions (eta, phi, pT) and record unmatched particles.
        4. Save results to a txt file.

        Returns:
        - A DataFrame of the matched jets.
        - A list of reconstructed jets for each event.
        """

        #print("df", df)
        #print("genpart_df", genpart_df)
        all_results = []
        all_jets = []

        # Get unique events from the data
        unique_events = df['event'].unique()
        print("unique_events:", unique_events)

        # Ensure the 'event' column in genpart_df is of integer type
        genpart_df['event'] = genpart_df['event'].astype(int)

        # Iterate over each unique event
        for event in unique_events:
            # Extract the DataFrame for the current event from both dataframes
            event_df = df[df['event'] == event]
            event_genpart_df = genpart_df[genpart_df['event'] == event]

            # Prepare the PseudoJet list for this event
            pseudojet_data = []
            particle_info = []
            for index, row in event_df.iterrows():
                pt = row['pt']
                if pt <= 0:
                    continue  # Skip PseudoJet creation if pt is <= 0

                # Compute eta_center and phi_center using the mean of vertices
                eta_center = np.mean(row['eta_vertices'])
                phi_center = np.mean(row['phi_vertices'])

                # Convert (pt, eta, phi) to (px, py, pz, E)
                px = pt * np.cos(phi_center)
                py = pt * np.sin(phi_center)
                pz = pt * np.sinh(eta_center)
                energy = (px**2 + py**2 + pz**2)**0.5

                # Create a PseudoJet and associate particle information
                pseudojet = fastjet.PseudoJet(px, py, pz, energy)
                pseudojet_data.append(pseudojet)
                particle_info.append({'eta_vertices': row['eta_vertices'], 'phi_vertices': row['phi_vertices']})

            # Perform clustering using the Anti-kt algorithm
            jet_definition = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
            cluster_sequence = fastjet.ClusterSequence(pseudojet_data, jet_definition)
            inclusive_jets = cluster_sequence.inclusive_jets()

            # Store the jets for this event
            all_jets.append(inclusive_jets)

            # Evaluate the resolution between reconstructed jets and generated particles
            for _, gen_row in event_genpart_df.iterrows():
                gen_eta = gen_row['gen_eta']
                gen_phi = gen_row['gen_phi']
                gen_pt = gen_row['gen_pt']

                # Find the closest jet within a radius of 0.1 in eta/phi
                matched_jet = None
                min_distance = 0.1
                for jet in inclusive_jets:
                    jet_eta = jet.eta()
                    jet_phi = jet.phi()
                    jet_phi = (jet_phi + np.pi) % (2 * np.pi) - np.pi
                    jet_pt = jet.pt()

                    # Compute the distance in eta/phi
                    delta_eta = jet_eta - gen_eta
                    delta_phi = (jet_phi - gen_phi + np.pi) % (2 * np.pi) - np.pi
                    distance = np.sqrt(delta_eta**2 + delta_phi**2)

                    if distance < min_distance:
                        matched_jet = jet
                        min_distance = distance

                # Handle the case where no jet is found within the threshold
                if matched_jet is None:
                     all_results.append({
                        'event': event,
                        'gen_eta': gen_eta,
                        'gen_phi': gen_phi,
                        'gen_pt': gen_pt,
                        'reco_eta': jet_eta,
                        'reco_phi': jet_phi,
                        'reco_pt': jet_pt,
                        'eta_diff': None,
                        'phi_diff': None,
                        'pt_ratio': None,
                        'matched': False  # Add a flag for unmatched particles
                    })

                # If a matching jet is found, compute resolutions and store results
                if matched_jet is not None:
                    jet_eta = matched_jet.eta()
                    jet_phi = matched_jet.phi()
                    jet_phi = (jet_phi + np.pi) % (2 * np.pi) - np.pi
                    jet_pt = matched_jet.pt()
                    jet_E = matched_jet.E()

                    # Calculate pt_ratio (reconstructed jet pt / generated particle pt)
                    pt_ratio = jet_pt / gen_pt if gen_pt != 0 else None  # Avoid division by zero

                    # Gather eta and phi vertices of the particles contributing to the jet
                    constituents = matched_jet.constituents()
                    constituent_eta_vertices = []
                    constituent_phi_vertices = []

                    #print("particle_info", particle_info)
                    for constituent in constituents:
                        idx = constituent.user_index()  # Retrieve original index
                        constituent_eta_vertices.append(particle_info[idx]['eta_vertices'])
                        constituent_phi_vertices.append(particle_info[idx]['phi_vertices'])

                    # Compute the difference in eta and phi (resolution)
                    eta_resolution = jet_eta - gen_eta
                    phi_resolution_no_wrap = jet_phi - gen_phi
                    phi_resolution = (phi_resolution_no_wrap + np.pi) % (2 * np.pi) - np.pi

                    # Collect results in a dictionary
                    all_results.append({
                        'event': event,
                        'gen_eta': gen_eta,
                        'gen_phi': gen_phi,
                        'gen_pt': gen_pt,
                        'reco_eta': jet_eta,
                        'reco_phi': jet_phi,
                        'reco_pt': jet_pt,
                        'eta_diff': eta_resolution,
                        'phi_diff': phi_resolution,
                        'pt_ratio': pt_ratio,
                        'matched': True  # Add a flag for matched particles
                    })

        # Convert all results to a DataFrame for easier analysis
        results_df = pd.DataFrame(all_results)
        #print("results_df", results_df)

        # Save the results to a txt file
        results_df.to_csv(ouptput_txt, sep=',', index=False)
        return results_df, all_jets

    def perform_clustering_antikt(self, df, output_txt):
        print("Performing antikt clustering for Minimum Bias sample")
        #print("Input DataFrame:", df.columns)
        all_results = []

        # Get unique events from the data
        unique_events = df['event'].unique()
        #print("Unique events:", unique_events)

        # Iterate over each unique event
        for event in unique_events:
            # Extract the DataFrame for the current event
            event_df = df[df['event'] == event]

            # Prepare the PseudoJet list for this event
            pseudojet_data = []
            for index, row in event_df.iterrows():
                pt = row['pt']
                if pt <= 0:
                    continue  # Skip PseudoJet creation if pt is <= 1 (selecting TTs only above 1GeV)

                # Compute eta_center and phi_center using the mean of vertices
                eta_center = np.mean(row['eta_vertices'])
                phi_center = np.mean(row['phi_vertices'])

                # Convert (pt, eta, phi) to (px, py, pz, E)
                px = pt * np.cos(phi_center)
                py = pt * np.sin(phi_center)
                pz = pt * np.sinh(eta_center)
                energy = (px**2 + py**2 + pz**2)**0.5

                # Create a PseudoJet
                pseudojet = fastjet.PseudoJet(px, py, pz, energy)
                pseudojet_data.append(pseudojet)

            # Perform clustering using the Anti-kt algorithm
            jet_definition = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
            cluster_sequence = fastjet.ClusterSequence(pseudojet_data, jet_definition)
            inclusive_jets = cluster_sequence.inclusive_jets()


            # Store the jets for this event
            for jet in inclusive_jets:
                all_results.append({
                    'event': event,
                    'reco_eta': jet.eta(),
                    'reco_phi':  (jet.phi() + np.pi) % (2 * np.pi) - np.pi,
                    'reco_pt': jet.pt(),
                })

        # Convert all results to a DataFrame for easier analysis
        results_df = pd.DataFrame(all_results)
        #print("Results DataFrame:", results_df)

        # Save the results to a txt file
        results_df.to_csv(output_txt, sep=',', index=False)

        # Count and print the number of jets reconstructed per event
        jet_counts = results_df.groupby('event').size()
        #print("Jet counts per event:", jet_counts)

        return results_df, jet_counts