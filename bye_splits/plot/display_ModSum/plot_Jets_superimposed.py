import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal*")
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import mplhep
import argparse
import os
# -------------------------
# Argument parser
# -------------------------
parser = argparse.ArgumentParser(description="Generate analysis plots with custom settings.")
parser.add_argument("--files", type=str, nargs='+', required=True,
                    help="List of input files (1 to 4).")
parser.add_argument("--algo", type=str, default="16towers_200PUantikt02_superimposed")
parser.add_argument("--subdet", type=str, default="CEE_Mod_CEH_STC")
parser.add_argument("--events", type=str, default="4k")
parser.add_argument("--particle", type=str, default="JETS")
args = parser.parse_args()

# -------------------------
# Setup
# -------------------------
output_dir = f'plots_{args.particle}_{args.algo}_{args.subdet}'
os.makedirs(output_dir, exist_ok=True)

colors = ['#4682B4', '#E74C3C', '#2ECC71', '#FFA500']
assert len(args.files) <= 4, "You can only input up to 3 files."

# Hardcoded labels (must match number of input files)
file_labels = ["No TT cut", "TT > 1 GeV", "TT > 2 GeV", "TT > 3 GeV"][:len(args.files)]

use_fit = False
use_eff_rms = True

# -------------------------
# Functions
# -------------------------
def effrms(resp_bin, c=0.68):
    resp_bin = np.sort(resp_bin, kind="mergesort")
    m = int(c * len(resp_bin)) + 1
    min_index = np.argmin(resp_bin[m:] - resp_bin[:-m])
    return resp_bin[min_index:min_index + m]

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
datasets = [load_file(f) for f in args.files]

# -------------------------
# 1) PT ratio distributions overlay with mean ± std in legend
# -------------------------
plt.figure(figsize=(12, 8))

for i, data in enumerate(datasets):
    pt_ratios = data["pt_ratio"][data["matched"]]
    
    # Compute mean and std
    mean = np.mean(pt_ratios)
    std = np.std(pt_ratios)
    
    # Plot histogram
    plt.hist(pt_ratios, bins=50, density=True, alpha=0.5,
             color=colors[i])
    
    # Use an invisible plot to include mean ± std in legend
    plt.plot([], [], color=colors[i], alpha=0.5,
             label=f'{file_labels[i]}: μ={mean:.2f}, σ={std:.2f}')

plt.xlabel(r'$p_{T}^{reco} / p_{T}^{gen} $')
plt.ylabel('Density')
plt.xlim(0, 3)  # limit x-axis to 0-3
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'{args.particle} PU200', fontsize=15)
plt.legend()
plt.grid(True)
plt.savefig(f'{output_dir}/pt_ratio_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/pt_ratio_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()

# -------------------------
# 2) Resolution plot (σ/μ vs pT)
# -------------------------
bin_edges = np.arange(0, 200, 20)
plt.figure(figsize=(10, 6))

for i, data in enumerate(datasets):
    gen_pt = data["gen_pt"][data["matched"]]
    pt_ratios = data["pt_ratio"][data["matched"]]

    sigma_mu_values = []
    err_sigma_mu_values = []

    for j in range(len(bin_edges) - 1):
        bin_mask = (gen_pt >= bin_edges[j]) & (gen_pt < bin_edges[j + 1])
        pt_ratios_in_bin = pt_ratios[bin_mask]

        if use_fit and len(pt_ratios_in_bin) > 0:
            mu, std = norm.fit(pt_ratios_in_bin)
        elif use_eff_rms and len(pt_ratios_in_bin) > 1:
            eff_rms_vals = effrms(pt_ratios_in_bin)
            mu = np.mean(eff_rms_vals)
            std = np.std(eff_rms_vals)
        elif len(pt_ratios_in_bin) > 0:
            mu = np.mean(pt_ratios_in_bin)
            std = np.std(pt_ratios_in_bin)
        else:
            mu, std = 0, 0

        sigma_mu = std / mu if mu != 0 else 0
        err_sigma_mu = std / (np.sqrt(2 * len(pt_ratios_in_bin) - 2) * mu) if mu != 0 and len(pt_ratios_in_bin) > 1 else 0
        sigma_mu_values.append(sigma_mu)
        err_sigma_mu_values.append(err_sigma_mu)

    plt.errorbar(bin_edges[:-1] + 10, sigma_mu_values,
                 yerr=err_sigma_mu_values,
                 fmt='o-', color=colors[i],
                 label=file_labels[i], capsize=5)

plt.xticks(np.arange(0, 220, 50))
plt.xlabel(r'$p_{T}^{gen} [GeV]$')
plt.ylabel(r'$\sigma/\mu$')
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'{args.particle} PU200', fontsize=15)
plt.legend()
plt.grid(True)
plt.savefig(f'{output_dir}/resolution_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/resolution_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
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

plt.xlabel(r'$p_{T}^{gen} [GeV]$')
plt.ylabel('Percentage of Matched Jets (%)')
mplhep.cms.label('Simulation Preliminary', data=True, rlabel=f'{args.particle} PU200', fontsize=15)
plt.legend()
plt.grid(True)
plt.savefig(f'{output_dir}/matched_percentage_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/matched_percentage_overlay_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()

