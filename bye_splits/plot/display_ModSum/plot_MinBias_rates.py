import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep

def compute_and_plot_jet_based_rates(file_path):
    # Load the txt file, treating it like a CSV
    data = pd.read_csv(file_path, delimiter=',')  # Adjust delimiter if necessary

    # Calculate the total number of events and reconstructed jets
    total_events = len(data['event'].unique())
    total_jets = len(data)
    average_jets_per_event = total_jets / total_events
    
    print(f"Total events: {total_events}")
    print(f"Total reconstructed jets: {total_jets}")
    print(f"Average number of reconstructed jets per event: {average_jets_per_event:.2f}")
    
    # Define the pT cut thresholds
    pt_thresholds = np.arange(1, 201, 5)
    total_events = len(data)
    
    # Initialize arrays for rates and errors
    rates = []
    y_errors = []
    x_errors = []

    # Scaling factor
    scaling_factor = 2340 * 11.245
    
    # Compute rates and errors
    for i, pt_cut in enumerate(pt_thresholds[:-1]):  # Exclude the last bin for x-errors
        entries_with_cut = len(data[data['reco_pt'] >= pt_cut])
        rate = entries_with_cut / total_events
        actual_rate = rate * scaling_factor  # Apply scaling factor
        y_error = np.sqrt(rate * scaling_factor)

        rates.append(actual_rate)
        
        y_errors.append(y_error)
        
        # Calculate horizontal error (half-width of the bin)
        x_error = (pt_thresholds[i + 1] - pt_cut) / 2
        x_errors.append(x_error)
    
    print("rates", rates)
    print("rates err", y_errors)
    print("x err", len(x_errors))
    print("pt_thresholds", len(pt_thresholds))

    # Plot the rates with error bars
    plt.style.use(mplhep.style.CMS)
    mplhep.cms.label('Private work', data=True, rlabel='Minimum bias')
    plt.errorbar(
        pt_thresholds[:-1],  # Exclude last point for x-axis
        rates,          # Exclude last point for y-axis (rates)
        yerr=y_errors,  # Vertical error
        xerr=x_errors,       # Horizontal error
        fmt='o',             # Marker style
        label='Jets-based rate',
        linestyle='-',
    )
    plt.xlabel(r'$p_{T}^{jets}$ [GeV]')
    plt.ylabel('Rate [kHz]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.yscale('log')  # Log scale for clarity if needed
    plt.savefig('rate_jet_based_area_overlap_TT.pdf')
    plt.savefig('rate_jet_based_area_overlap_TT.png')
    plt.close()

def compute_and_plot_event_based_rates(file_path, scaling_factor=2340 * 11.245): #2340 * 11.245
    """Compute and plot rate using the number of events with at least one jet passing the threshold."""
    # Load the txt file, treating it like a CSV
    data = pd.read_csv(file_path, delimiter=',')  # Adjust delimiter if necessary

    # Calculate the total number of events
    total_events = len(data['event'].unique())
    
    # Define the pT cut thresholds
    pt_thresholds = np.arange(1, 201, 5)
    
    # Initialize arrays for rates and errors
    rates = []
    y_errors = []
    x_errors = []
    
    # Compute rates and errors
    for i, pt_cut in enumerate(pt_thresholds[:-1]):  # Exclude the last bin for x-errors
        events_with_cut = len(data[data['reco_pt'] >= pt_cut]['event'].unique())
        print(total_events)
        rate = events_with_cut / total_events
        actual_rate = rate * scaling_factor
        y_error = np.sqrt(rate * scaling_factor)
        
        rates.append(actual_rate)
        y_errors.append(y_error)
        
        # Horizontal error (half-width of the bin)
        x_error = (pt_thresholds[i + 1] - pt_cut) / 2
        x_errors.append(x_error)

    print(rates)
    print(y_errors)
    
    # Plot the rates with error bars
    plt.style.use(mplhep.style.CMS)
    mplhep.cms.label('Private work', data=True, rlabel='Minimum bias')
    
    plt.errorbar(
        pt_thresholds[:-1],  # Exclude last point for x-axis
        rates,               # Event-based rates
        yerr=y_errors,       # Vertical error
        xerr=x_errors,       # Horizontal error
        fmt='s',             # Marker style
        label='Event-based rate',
        linestyle='-',
    )
    
    plt.xlabel(r'$p_{T}^{jets}$ [GeV]')
    plt.ylabel('Rate [kHz]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.yscale('log')  # Log scale for clarity if needed
    plt.savefig('rate_event_based_area_overlap_TT.pdf')
    plt.savefig('rate_event_based_area_overlap_TT.png')
    plt.close()



# Execute the function
compute_and_plot_jet_based_rates('area_overlap_neutrinos_-1_5_results_newgeom.txt')
compute_and_plot_event_based_rates('area_overlap_neutrinos_-1_5_results_newgeom.txt')
