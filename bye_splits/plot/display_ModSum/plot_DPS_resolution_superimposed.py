import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal*")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import mplhep
import argparse
import os
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

mplhep.style.use("CMS")

plt.rcParams.update({
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 40,   # Bigger text
    "legend.framealpha": 0.5 # Semi-transparent background (optional)
})

# -------------------------
# Argument parser
# -------------------------
parser = argparse.ArgumentParser(description="Generate analysis plots with custom settings.")
parser.add_argument("--files", type=str, nargs='+', required=True,
                    help="List of input files (1 to 4).")
parser.add_argument("--algo", type=str, default="base_16_PIONS_SUPERIMPOSED")
parser.add_argument("--subdet", type=str, default="CEE_Mod_CEH_STC")
parser.add_argument("--events", type=str, default="4k")
parser.add_argument("--particle", type=str, default="Pions")
args = parser.parse_args()

# -------------------------
# Setup
# -------------------------
output_dir = f'plots_{args.particle}_{args.algo}_{args.subdet}'
os.makedirs(output_dir, exist_ok=True)

colors = ['#4682B4', '#E74C3C', '#2ECC71', '#FFA500']
assert len(args.files) <= 4, "You can only input up to 4 files."

# Hardcoded labels (must match number of input files)
file_labels = ["No module splitting", "1/16 module splitting", "TT > 2 GeV", "TT > 3 GeV"][:len(args.files)]

use_fit = False
use_eff_rms = True

# -------------------------
# Functions
# -------------------------
def effrms(df, c=0.68):
    """Compute half-width of the shortest interval
    containing a fraction 'c' of items in a 1D array or DataFrame.
    """
    out = {}
    for col in df:
        x = df[col].values
        x = np.sort(x, kind="mergesort")
        m = int(c * len(x)) + 1
        out[col] = [np.min(x[m:] - x[:-m]) / 2.0]
    return pd.DataFrame(out).iloc[0]

def load_file(filename):
    data = []
    with open(filename, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            values = line.strip().split(',')
            if len(values) >= 11:
                event = int(values[0])
                gen_eta = float(values[1]) if values[1] else np.nan
                gen_phi = float(values[2]) if values[2] else np.nan
                gen_pt = float(values[3]) if values[3] else np.nan
                reco_eta = float(values[4]) if values[4] else np.nan
                reco_phi = float(values[5]) if values[5] else np.nan
                reco_pt = float(values[6]) if values[6] else np.nan
                eta_diff = float(values[7]) if values[7] else np.nan
                phi_diff = float(values[8]) if values[8] else np.nan
                pt_ratio = float(values[9]) if values[9] else np.nan
                matched = 1.0 if values[10] == 'True' else 0.0
                non_matched = 1.0 if values[10] == 'False' else 0.0
                data.append((event, gen_eta, gen_phi, gen_pt,
                             reco_eta, reco_phi, reco_pt,
                             eta_diff, phi_diff, pt_ratio,
                             matched, non_matched))
    data = np.array(data)
    return {
        "gen_pt": data[:, 3],
        "pt_ratio": data[:, 9],
        "eta_diff": data[:, 7],
        "phi_diff": data[:, 8],
        "matched": data[:, 10].astype(bool),
        "non_matched": data[:, 11].astype(bool)
    }

# -------------------------
# Load datasets
# -------------------------
print("Load datasets...")
datasets = [load_file(f) for f in args.files]

# -------------------------
# 1) PT ratio distributions overlay with mean ± effRMS in legend
# -------------------------
plt.figure(figsize=(12, 8))

for i, data in enumerate(datasets):
    pt_ratios = data["pt_ratio"][data["matched"]]
    
    mean = np.mean(pt_ratios)
    pt_ratios_df = pd.DataFrame({'pt_ratio': pt_ratios})
    eff_rms_val = effrms(pt_ratios_df)['pt_ratio']
    print("eff_rms_val", eff_rms_val)
    mplhep.style.use("CMS")
    plt.hist(pt_ratios, bins=50, density=True, alpha=0.5, color=colors[i])
    plt.plot([], [], color=colors[i], alpha=0.5,
             label = f'{file_labels[i]}: \n$\\mu={mean:.2f}, \\sigma_{{\\rm eff}}={eff_rms_val:.2f}$')



plt.text(x=0.75, y=0.5,  # x>1 moves outside the right side
         s="Pions PU=0", 
         transform=plt.gca().transAxes,  # use axes coordinates
         fontsize=19,
         verticalalignment='center')
plt.legend(fontsize=18)
plt.xlabel(r'$p_{T}^{reco} / p_{T}^{gen} $', fontsize=25)
plt.ylabel('a.u', fontsize=25)
plt.xlim(0.25, 1.50)
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'', fontsize=25)
plt.grid(True)
plt.savefig(f'{output_dir}/pt_ratio_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/pt_ratio_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()

from scipy.optimize import curve_fit

# -------------------------
# 2) Resolution plot (σ/μ vs pT) + per-dataset superimposed histograms
# -------------------------
bin_edges = np.arange(0, 201, 20)

# Folder for combined histograms
hist_dir = os.path.join(output_dir, "hists")
os.makedirs(hist_dir, exist_ok=True)

# Precompute global x-range across all datasets
all_ratios = np.concatenate([data["pt_ratio"][data["matched"]] for data in datasets])
x_min, x_max = np.min(all_ratios), np.max(all_ratios)

# Define consistent bin edges for all histograms
n_hist_bins = 80
hist_bins = np.linspace(x_min, x_max, n_hist_bins + 1)

# Store σ/μ curves
sigma_mu_all = []
err_sigma_mu_all = []
bin_centers_all = []

# Define ranges to split (full, 0–100, 100–200)
range_splits = {
    "all_bins": (0, 200),
    "bins_0_100": (0, 100),
    "bins_100_200": (100, 200),
}

for i, data in enumerate(datasets):
    gen_pt = data["gen_pt"][data["matched"]]
    pt_ratios = data["pt_ratio"][data["matched"]]

    sigma_mu_values = []
    err_sigma_mu_values = []
    bin_centers = []

    # --- Loop over range_splits ---
    for range_label, (pt_min, pt_max) in range_splits.items():
        plt.figure(figsize=(8, 6))

        for j in range(len(bin_edges) - 1):
            if bin_edges[j] < pt_min or bin_edges[j+1] > pt_max:
                continue  # skip bins outside chosen range

            bin_mask = (gen_pt >= bin_edges[j]) & (gen_pt < bin_edges[j + 1])
            pt_ratios_in_bin = pt_ratios[bin_mask]

            if len(pt_ratios_in_bin) == 0:
                continue

            # Mean and effRMS
            pt_ratios_df = pd.DataFrame({'pt_ratio': pt_ratios_in_bin})
            eff_rms_val = effrms(pt_ratios_df)['pt_ratio']
            mu = np.mean(pt_ratios_in_bin)
            std = eff_rms_val

            sigma_mu = std / mu if mu != 0 else 0
            #err_sigma_mu = std / (np.sqrt(2 * len(pt_ratios_in_bin) - 2) * mu) if mu != 0 and len(pt_ratios_in_bin) > 1 else 0

            N = len(pt_ratios_in_bin)
            if N > 1 and mu != 0:
                # Uncertainty on std and mean
                err_sigma = std / np.sqrt(2*(N-1))     
                err_mu    = std / np.sqrt(N)           
                err_ratio = np.sqrt((err_sigma / mu)**2 + ((std * err_mu) / mu**2)**2)
            else:
                err_ratio = 0

            sigma_mu_values.append(sigma_mu)
            err_sigma_mu_values.append(err_ratio)
            print("err_sigma_mu_values", err_sigma_mu_values)
            bin_centers.append((bin_edges[j] + bin_edges[j+1]) / 2)

            # Histogram for this bin on same figure
            plt.hist(pt_ratios_in_bin, bins=hist_bins, histtype='step',
                     linewidth=1, density=True,
                     label = f"{bin_edges[j]}–{bin_edges[j+1]} GeV: $\\mu={mu:.3f}, \\sigma_{{\\rm eff}}={std:.3f}, \\sigma_{{\\rm eff}}/\\mu={std/mu:.3f}$")

        # Styling for this dataset + range
        plt.xlim(x_min, x_max)
        plt.title(f"{file_labels[i]}: pT ratio distributions ({pt_min}–{pt_max} GeV)")
        plt.xlabel("pT ratio")
        plt.ylabel("Counts")
        plt.grid(True)
        plt.legend(fontsize=18, framealpha=0.3)

        # Save figure
        safe_label = file_labels[i].replace(" ", "_").replace("/", "_")
        plt.savefig(f"{hist_dir}/{range_label}_{safe_label}.png")
        plt.close()

    sigma_mu_all.append(sigma_mu_values)
    err_sigma_mu_all.append(err_sigma_mu_values)
    bin_centers_all.append(bin_centers)

# --- Define the resolution fit function (new version) ---
def res_func(E, a, b, c):
    """Jet resolution parametrization: sigma/mu = sqrt((a/√E)^2 + (b/E)^2 + c^2)"""
    return np.sqrt((a / np.sqrt(E))**2 + (b / E)**2 + c**2)

def res_func_simplified(E, a, c):
    """Jet resolution parametrization: sigma/mu = sqrt((a/√E)^2 + (b/E)^2 + c^2)"""
    return np.sqrt((a / np.sqrt(E))**2 + c**2)
    
# --- Resolution overlay plot with ratio panel ---
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)


popt_list = []
perr_list = []

# --- Top panel: original overlay ---
ax0 = fig.add_subplot(gs[0])
for i in range(len(datasets)):
    x = np.array(bin_centers_all[i])
    y = np.array(sigma_mu_all[i])
    yerr = np.array(err_sigma_mu_all[i])

    # Fit the resolution function
    try:
        #popt, pcov = curve_fit(
            #res_func, x, y, sigma=yerr, absolute_sigma=True, p0=[1.0, 1.0, 0.01]
        #)

        '''popt, pcov = curve_fit(
            res_func, x, y, sigma=yerr, absolute_sigma=True,
            p0=[1.0, 0.0, 0.01],  # starting values
            bounds=([0, 0, 0], [np.inf, 0.5, 1.0])  # constrain parameters
        )'''

        popt, pcov = curve_fit(
            res_func_simplified, x, y, sigma=yerr, absolute_sigma=True,
            p0=[1.0, 0.01],  # starting values
            bounds=([0, 0], [np.inf, 1.0])  # constrain parameters
        )

        y_fit_vals = res_func_simplified(x, *popt)
        chi2 = np.sum(((y - y_fit_vals) / yerr)**2)
        dof = len(y) - len(popt)
        chi2_text = f"$\\chi^2$/dof = {chi2:.1f}/{dof}"
        print("chi2_text", chi2_text)

        # Extract fit parameters with uncertainties
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [0, 0]
        #param_text = (f"a={popt[0]:.2f}±{perr[0]:.2f}, "
                    #f"c={popt[1]:.3f}±{perr[1]:.3f}")


        # --- Print fit results with errors ---
        print(f"Fit results for {file_labels[i]}:")
        print(f"  a = {popt[0]:.4f} ± {perr[0]:.4f}")
        print(f"  c = {popt[1]:.4f} ± {perr[1]:.4f}")
        print(f"  chi2/dof = {chi2:.2f}/{dof}")

        # store them for later comparison
        popt_list.append(popt)
        perr_list.append(perr)
                

        # only do this if you have at least 2 datasets
        if len(popt_list) >= 2:
            a1, c1 = popt_list[0][0], popt_list[0][1]
            err_a1, err_c1 = perr_list[0][0], perr_list[0][1]

            a2, c2 = popt_list[1][0], popt_list[1][1]
            err_a2, err_c2 = perr_list[1][0], perr_list[1][1]

            # --- compare c ---
            delta_c = abs(c1 - c2)
            combined_err_c = np.sqrt(err_c1**2 + err_c2**2)
            significance_c = delta_c / combined_err_c
            print(f"Difference in c is {delta_c:.4f}, which is {significance_c:.2f}σ")

            # --- compare a ---
            delta_a = abs(a1 - a2)
            combined_err_a = np.sqrt(err_a1**2 + err_a2**2)
            significance_a = delta_a / combined_err_a
            print(f"Difference in a is {delta_a:.4f}, which is {significance_a:.2f}σ")
        # Make a string showing the fitted resolution function explicitly in terms of pT
        param_text = (r"$\sigma/p_T = \sqrt{"
                    f"({popt[0]:.2f} / \sqrt{{p_T}})^2 + "
                    f"{popt[1]:.3f}^2"
                    r"}$")

        print("Fit params:", param_text)

    except RuntimeError:
        print(f"Warning: Fit failed for dataset {file_labels[i]}")
        popt = [0, 0, 0]

    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = res_func_simplified(x_fit, *popt)

    mplhep.style.use("CMS")
    ax0.errorbar(x, y, yerr=yerr, fmt='o', color=colors[i], label=f"{file_labels[i]}", capsize=5)
    ax0.plot(x_fit, y_fit, '-', color=colors[i], alpha=0.6, label=f"{file_labels[i]}: {param_text}")

ax0.text(0.75, 0.5, s="Pions PU=0", transform=ax0.transAxes, fontsize=19, verticalalignment='center')
ax0.set_xticks(np.arange(0, 225, 25))
ax0.set_ylabel(r'$\sigma_{eff}/\mu$', fontsize=25)
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'', fontsize=25)
ax0.legend(fontsize=16)
ax0.grid(True)

# --- Bottom panel: ratio of first two datasets using existing points ---
ax1 = fig.add_subplot(gs[1], sharex=ax0)

if len(datasets) >= 2:
    # Use only overlapping bins
    n_bins = min(len(bin_centers_all[0]), len(bin_centers_all[1]))
    x_common = bin_centers_all[0][:n_bins]  # use first dataset's bin centers
    R1 = np.array(sigma_mu_all[0][:n_bins])
    R2 = np.array(sigma_mu_all[1][:n_bins])
    ratio = R1 / R2

    # Compute error bars on the ratio
    err1 = np.array(err_sigma_mu_all[0][:n_bins])
    err2 = np.array(err_sigma_mu_all[1][:n_bins])
    ratio_err = ratio * np.sqrt((err1 / R1)**2 + (err2 / R2)**2)

    # Plot points with error bars
    ax1.set_xticks(np.arange(0, 225, 25))
    mplhep.style.use("CMS")
    ax1.errorbar(x_common, ratio, yerr=ratio_err, fmt='o', color='#4682B4', markersize=5, capsize=3)
    ax1.axhline(1.0, color='#E74C3C', linestyle='--', linewidth=1)
    ax1.set_ylabel(f'Ratio', fontsize=25)
    ax1.set_xlabel(r'$p_{T}^{gen} [GeV]$', fontsize=25)
    ax1.grid(True)

# Remove x tick labels on top panel
plt.setp(ax0.get_xticklabels(), visible=False)

# Save figure
plt.savefig(f'{output_dir}/resolution_overlay_with_ratio_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/resolution_overlay_with_ratio_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()




# -------------------------
# 3) Matched percentage vs pT
# -------------------------
pt_bins = np.linspace(0, 200, 21)
plt.figure(figsize=(10, 6))

for i, data in enumerate(datasets):
    gen_pt = data["gen_pt"]
    mask = data["matched"]

    total_counts, _ = np.histogram(gen_pt, bins=pt_bins)
    matched_counts, _ = np.histogram(gen_pt[mask], bins=pt_bins)

    percentage_matched = (matched_counts / total_counts) * 100
    percentage_matched = np.nan_to_num(percentage_matched)

    lo_err, up_err = [], []
    for k in range(len(pt_bins) - 1):
        n = total_counts[k]
        s = matched_counts[k]
        if n == 0:
            lo_err.append(0)
            up_err.append(0)
            continue
        result = stats.binomtest(s, n, p=s/n)
        ci = result.proportion_ci(confidence_level=0.95)
        lo_err.append(percentage_matched[k] - ci.low * 100)
        up_err.append(ci.high * 100 - percentage_matched[k])

    pt_bin_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
    pt_bin_widths = (pt_bins[1:] - pt_bins[:-1]) / 2

    plt.errorbar(pt_bin_centers, percentage_matched,
                 xerr=pt_bin_widths, yerr=[lo_err, up_err],
                 fmt='o-', color=colors[i],
                 label=file_labels[i],
                 capsize=3, alpha=0.8)

plt.text(x=0.75, y=0.5,  # x>1 moves outside the right side
         s="Pions PU=0", 
         transform=plt.gca().transAxes,  # use axes coordinates
         fontsize=19,
         verticalalignment='center')
plt.xlabel(r'$p_{T}^{gen} [GeV]$', fontsize=25)
plt.ylabel('Percentage of Matched Pions (%)', fontsize=25)
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'', fontsize=25)
plt.legend(fontsize=18)
plt.grid(True)
plt.savefig(f'{output_dir}/matched_percentage_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/matched_percentage_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()

# -------------------------
# 4) Gen vs Reco pT distributions
# -------------------------
plt.figure(figsize=(12, 8))

# Overlay gen_pt distributions
for i, data in enumerate(datasets):
    gen_pt = data["gen_pt"][data["matched"]]  # only matched for fair comparison
    plt.hist(gen_pt, bins=50, histtype='step', linewidth=2,
             color=colors[i], label=f'{file_labels[i]} Gen pT')

plt.xlabel(r'$p_{T}^{gen}$ [GeV]', fontsize=25)
plt.ylabel('Counts (a.u)', fontsize=25)
plt.title('Generator-level pT Distributions', fontsize=22)
plt.legend(fontsize=18)
plt.grid(True)
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'', fontsize=25)
plt.savefig(f'{output_dir}/gen_pt_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/gen_pt_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()


plt.figure(figsize=(12, 8))

# Overlay reco_pt distributions
for i, data in enumerate(datasets):
    reco_pt = data["pt_ratio"][data["matched"]] * data["gen_pt"][data["matched"]]
    plt.hist(reco_pt, bins=50, histtype='step', linewidth=2,
             color=colors[i], label=f'{file_labels[i]} Reco pT')

plt.xlabel(r'$p_{T}^{reco}$ [GeV]', fontsize=20)
plt.ylabel('Counts (a.u)', fontsize=20)
plt.title('Reconstructed pT Distributions', fontsize=22)
plt.legend(fontsize=18)
plt.grid(True)
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'', fontsize=25)
plt.savefig(f'{output_dir}/reco_pt_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/reco_pt_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()
