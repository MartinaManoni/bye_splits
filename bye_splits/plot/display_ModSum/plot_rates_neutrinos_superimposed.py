import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep

# Hardcoded labels for up to 4 files
'''file_labels = [
    "Baseline: TT > 1 GeV",
    "Split in 1/16s: TT > 1 GeV"
]'''

file_labels = [
    "No cut",
    "TT > 1 GeV",
    "TT > 2 GeV",
    "TT > 3 GeV"
]

# Hardcoded marker styles and colors
markers = ['o', 's', 'v', '^']
colors = ['C0', 'C1', 'C2', 'C3']

SCALING_FACTOR = 2340 * 11.245

def compute_jet_based_rates(data):
    total_events = len(data['event'].unique())
    pt_thresholds = np.arange(1, 201, 5)
    rates, y_errors, x_errors = [], [], []

    for i, pt_cut in enumerate(pt_thresholds[:-1]):
        entries_with_cut = len(data[data['reco_pt'] >= pt_cut])
        rate = entries_with_cut / total_events
        actual_rate = rate * SCALING_FACTOR
        y_error = np.sqrt(rate * SCALING_FACTOR)
        rates.append(actual_rate)
        y_errors.append(y_error)
        x_errors.append((pt_thresholds[i + 1] - pt_cut) / 2)

    return pt_thresholds[:-1], rates, y_errors, x_errors

def compute_event_based_rates(data):
    total_events = len(data['event'].unique())
    pt_thresholds = np.arange(1, 201, 5)
    rates, y_errors, x_errors = [], [], []

    for i, pt_cut in enumerate(pt_thresholds[:-1]):
        events_with_cut = len(data[data['reco_pt'] >= pt_cut]['event'].unique())
        rate = events_with_cut / total_events
        actual_rate = rate * SCALING_FACTOR
        y_error = np.sqrt(rate * SCALING_FACTOR)
        rates.append(actual_rate)
        y_errors.append(y_error)
        x_errors.append((pt_thresholds[i + 1] - pt_cut) / 2)

    return pt_thresholds[:-1], rates, y_errors, x_errors

def plot_superimposed_rates(file_paths, particle_label, plot_type='jet'):
    """
    Plot rates from multiple files on the same plot with dynamic CMS rlabel.
    - particle_label: string to include in CMS rlabel (e.g., 'Muon', 'Electron')
    """
    plt.style.use(mplhep.style.CMS)
    mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'{particle_label} PU200', fontsize=15)
    
    for idx, file_path in enumerate(file_paths[:4]):
        data = pd.read_csv(file_path, delimiter=',')
        label = file_labels[idx]
        marker = markers[idx]
        color = colors[idx]

        if plot_type == 'jet':
            x, y, yerr, xerr = compute_jet_based_rates(data)
        elif plot_type == 'event':
            x, y, yerr, xerr = compute_event_based_rates(data)
        else:
            raise ValueError("plot_type must be 'jet' or 'event'")

        plt.errorbar(
            x, y, yerr=yerr, xerr=xerr,
            fmt=marker, color=color, label=label, linestyle='-'
        )
    
    plt.xlabel(r'$p_{T}^{jets}$ [GeV]')
    plt.ylabel('Rate [kHz]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.yscale('log')
    
    outname = f'rate_{plot_type}_based_superimposed_16towers_STCS.pdf'
    plt.savefig(outname)
    plt.savefig(outname.replace('.pdf', '.png'))
    plt.close()


# Example usage
files = [
    '16towers_neutrinos_-1_5_results_2Ntuples.txt',
    '16towers_neutrinos_-1_5_results_2Ntuples_TT_1GeV.txt',
    '16towers_neutrinos_-1_5_results_2Ntuples_TT_2GeV.txt',
    '16towers_neutrinos_-1_5_results_2Ntuples_TT_3GeV.txt'
]

particle = "Neutrinos"  # Example particle name

plot_superimposed_rates(files, particle_label=particle, plot_type='jet')
plot_superimposed_rates(files, particle_label=particle, plot_type='event')


