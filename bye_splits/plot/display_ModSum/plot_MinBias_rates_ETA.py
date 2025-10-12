import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep

def compute_event_based_rates(data, pt_thresholds, scaling_factor=2340 * 11.245):
    total_events = len(data['event'].unique())
    rates, y_errors, x_errors = [], [], []

    for i, pt_cut in enumerate(pt_thresholds[:-1]):
        events_with_cut = len(data[data['reco_pt'] >= pt_cut]['event'].unique())
        rate = events_with_cut / total_events
        actual_rate = rate * scaling_factor
        y_error = np.sqrt(rate * scaling_factor)

        rates.append(actual_rate)
        y_errors.append(y_error)
        x_error = (pt_thresholds[i + 1] - pt_cut) / 2
        x_errors.append(x_error)

    return rates, y_errors, x_errors

def compare_event_based_rates(file1, file2, scaling_factor=2340 * 11.245):
    data1 = pd.read_csv(file1, delimiter=',')
    data2 = pd.read_csv(file2, delimiter=',')
    pt_thresholds = np.arange(1, 201, 5)

    # Print sums for both files
    for label, data in zip(["File1", "File2"], [data1, data2]):
        total_events = len(data['event'].unique())
        sum_reco_pt_all = data['reco_pt'].sum()
        sum_reco_pt_eta_lt2 = data[data['reco_eta'] < 2]['reco_pt'].sum()
        sum_reco_pt_eta_ge2 = data[data['reco_eta'] >= 2]['reco_pt'].sum()
        print(f"--- {label} ---")
        print(f"Total events: {total_events}")
        print(f"Sum of reco_pt (all): {sum_reco_pt_all:.2f}")
        print(f"Sum of reco_pt (reco_eta < 2): {sum_reco_pt_eta_lt2:.2f}")
        print(f"Sum of reco_pt (reco_eta >= 2): {sum_reco_pt_eta_ge2:.2f}")

    # Function to plot comparison for a subset
    def plot_comparison(subset1, subset2, label, suffix):
        rates1, yerr1, xerr1 = compute_event_based_rates(subset1, pt_thresholds, scaling_factor)
        rates2, yerr2, xerr2 = compute_event_based_rates(subset2, pt_thresholds, scaling_factor)

        plt.style.use(mplhep.style.CMS)
        mplhep.cms.label('Simulation Preliminary', data=True, rlabel='Minimum bias')

        plt.errorbar(pt_thresholds[:-1], rates1, yerr=yerr1, xerr=xerr1,
                     fmt='o', label=f'{label} (No module splitting)', linestyle='-')
        plt.errorbar(pt_thresholds[:-1], rates2, yerr=yerr2, xerr=xerr2,
                     fmt='s', label=f'{label} (1/16 module splitting)', linestyle='--')

        plt.xlabel(r'$p_{T}^{jets}$ [GeV]')
        plt.ylabel('Rate [kHz]')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig(f'compare_event_based_{suffix}3GeV.pdf')
        plt.savefig(f'compare_event_based_{suffix}3GeV.png')
        plt.close()

    # Compare for reco_eta < 2
    plot_comparison(data1[data1['reco_eta'] < 2], data2[data2['reco_eta'] < 2],
                    'Event-based rate (reco_eta < 2)', 'eta_lt2')

    # Compare for reco_eta >= 2
    plot_comparison(data1[data1['reco_eta'] >= 2], data2[data2['reco_eta'] >= 2],
                    'Event-based rate (reco_eta >= 2)', 'eta_ge2')


# Example usage
compare_event_based_rates(
    'baseline_neutrinos_-1_5_results_2Ntuples_TT_3GeV.txt',
    '16towers_neutrinos_-1_5_results_2Ntuples_TT_3GeV.txt'
)


