import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
import mplhep
import argparse
import os



# Argument parser for customization
parser = argparse.ArgumentParser(description="Generate analysis plots with custom settings.")
parser.add_argument("--algo", type=str, default="8towers", help="Algorithm used (e.g., baseline, area_overlap..).")
parser.add_argument("--subdet", type=str, default="CEE_CEH", help="Subdetector type (e.g., 1,2,3 ...")
parser.add_argument("--events", type=str, default="4k", help="Event description (e.g., 1k, 2k...")
parser.add_argument("--particle", type=str, default="Jets", help="Event description (e.g., 1k, 2k...")
args = parser.parse_args()

# Create output directory dynamically
output_dir = f'plots_{args.particle}_{args.algo}_{args.subdet}'
os.makedirs(output_dir, exist_ok=True)

# Define the effrms function
def effrms(resp_bin, c=0.68):
    """ Compute half-width of the shortest interval (min std)
    containing a fraction 'c' of items """
    resp_bin = np.sort(resp_bin, kind="mergesort")
    m = int(c * len(resp_bin)) + 1
    min_index = np.argmin(resp_bin[m:] - resp_bin[:-m])
    return resp_bin[min_index:min_index + m]

# Option to choose whether to use the Gaussian fit results or raw mean and std from the data
use_fit = False  # Set to False to use the raw mean and std from the data
use_eff_rms = True  # Set to True to use the eff_rms method

# Load the data from the file and read the header
data = []
with open('8towers_jets_-1_5_PIONS_results.txt', 'r') as f: #matched_filtered_reordered_results.txt
    # Read the header line
    header = f.readline().strip().split(',')
    print("Header:", header)  # Optional: Print the header to check if it's correct

    # Read the rest of the data
    for line in f:
        values = line.strip().split(',')
        #print(len(values))
        if len(values) >=11:
            # Convert the relevant columns to floats
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
            matched = 1.0 if values[10] == 'True' else 0.0  # Convert 'True'/'False' to numeric values
            non_matched = 1.0 if values[10] == 'False' else 0.0
            data.append((event, gen_eta, gen_phi,gen_pt,reco_eta,reco_phi,reco_pt,eta_diff,phi_diff,pt_ratio, matched, non_matched))
        else:
            print("line boh")
# Convert the data into a numpy array for easier processing
data = np.array(data)

# Extract the genpart_pt and pt_ratio columns
events = data[:, 0]
print("events", events)
gen_etas = data[:, 1]
gen_phis = data[:, 2]
gen_pts = data[:, 3]
print("gen_pts", gen_pts)
reco_etas = data[:, 4]
reco_phis = data[:, 5]
reco_pts = data[:, 6]
eta_diffs = data[:, 7]
phi_diffs = data[:, 8]
pt_ratios = data[:, 9]
matched = data[:, 10].astype(bool)
non_matched = data[:, 11].astype(bool)


matched= matched

print("pt_ratios", pt_ratios)


def plot_scatter(x, y, xlabel, ylabel, title, filename, xlims=None, ylims=None):
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, alpha=0.6, s=10, c='blue', label='All Events')
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    if xlims:
        plt.xlim(xlims)
    if ylims:
        plt.ylim(ylims)
    plt.grid(alpha=0.5)
    mplhep.cms.label('Private work', data=True, rlabel=f'{args.particle} PU0')
    plt.legend()
    plt.savefig(filename)
    plt.close()


# Gen vs Reco Eta (All Events)
plot_scatter(
    gen_etas, reco_etas,
    r'$\eta^{gen}$', r'$\eta^{reco}$',
    'Gen vs Reco Eta (All Events)',
    f'{output_dir}/gen_vs_reco_eta_all_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

plot_scatter(
    gen_etas[matched], reco_etas[matched],
    r'$\eta^{gen}$', r'$\eta^{reco}$',
    'Gen vs Reco Eta (Matched)',
    f'{output_dir}/gen_vs_reco_eta_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

# Gen vs Reco Eta (Non-Matched)
plot_scatter(
    gen_etas[non_matched], reco_etas[non_matched],
    r'$\eta^{gen}$', r'$\eta^{reco}$',
    'Gen vs Reco Eta (Non-Matched)',
    f'{output_dir}/gen_vs_reco_eta_non_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

# Gen vs Reco Phi (All Events)
plot_scatter(
    gen_phis, reco_phis,
    r'$\phi^{gen}$', r'$\phi^{reco}$',
    'Gen vs Reco Phi (All Events)',
    f'{output_dir}/gen_vs_reco_phi_all_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

# Gen vs Reco Phi (Matched)
plot_scatter(
    gen_phis[matched], reco_phis[matched],
    r'$\phi^{gen}$', r'$\phi^{reco}$',
    'Gen vs Reco Phi (Matched)',
    f'{output_dir}/gen_vs_reco_phi_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

# Gen vs Reco Phi (Non-Matched)
plot_scatter(
    gen_phis[non_matched], reco_phis[non_matched],
    r'$\phi^{gen}$', r'$\phi^{reco}$',
    'Gen vs Reco Phi (Non-Matched)',
    f'{output_dir}/gen_vs_reco_phi_non_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png'
)

# Gen vs Reco Pt (All Events)
plot_scatter(
    gen_pts, reco_pts,
    r'$p_T^{gen}$ [GeV]', r'$p_T^{reco}$ [GeV]',
    'Gen vs Reco Pt (All Events)',
    f'{output_dir}/gen_vs_reco_pt_all_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png',
    xlims=(0, 200), ylims=(0, 200)
)

# Gen vs Reco Pt (Matched)
plot_scatter(
    gen_pts[matched], reco_pts[matched],
    r'$p_T^{gen}$ [GeV]', r'$p_T^{reco}$ [GeV]',
    'Gen vs Reco Pt (Matched)',
    f'{output_dir}/gen_vs_reco_pt_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png',
    xlims=(0, 200), ylims=(0, 200)
)

# Gen vs Reco Pt (Non-Matched)
plot_scatter(
    gen_pts[non_matched], reco_pts[non_matched],
    r'$p_T^{gen}$ [GeV]', r'$p_T^{reco}$ [GeV]',
    'Gen vs Reco Pt (Non-Matched)',
    f'{output_dir}/gen_vs_reco_pt_non_matched_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png',
    xlims=(0, 200), ylims=(0, 200)
)


mask = matched
events = events[mask]
gen_etas = gen_etas[mask]
gen_phis = gen_phis[mask]
reco_pts = reco_pts[mask]
pt_ratios = pt_ratios[mask]

print("pt_ratios masked",pt_ratios )




plt.figure(figsize=(10, 6))

# Plot the histogram
plt.hist(gen_pts[mask], bins=30, color='#4682B4', alpha=0.7, edgecolor='black', density=True)

# Add labels and title
plt.xlabel(r'$p_{T}^{gen} \, [GeV]$', fontsize=14)
plt.ylabel('Normalized Density', fontsize=14)
plt.title('Distribution of Gen-Level Transverse Momentum ($p_{T}^{gen}$)', fontsize=16)

# Add grid and style
plt.grid(alpha=0.5)
mplhep.cms.label('Private work', data=True, rlabel=f'{args.particle} PU0')

# Save and display the plot
#plt.savefig('gen_pt_distribution.png')
#plt.savefig('gen_pt_distribution.pdf')
plt.close()






# Define the number of bins for pt
bin_edges = np.arange(0, 200, 20)
#bin_edges = np.arange(0, 220, 20)
print("bin_edges", bin_edges)

# Prepare to store results for each bin
sigma_mu_values = []
err_sigma_mu_values = []

# Loop over each bin and process the data
for i in range(len(bin_edges) - 1):
    # Find the events that fall into the current bin
    bin_mask = (gen_pts[mask] >= bin_edges[i]) & (gen_pts[mask] < bin_edges[i + 1])

    # Get the pt_ratios for the events in the current bin
    pt_ratios_in_bin = pt_ratios[bin_mask]

    # Calculate mu and sigma based on the chosen method
    if use_fit:
        print("Using Histogram Fitted values ...")
        # Use the Gaussian fit
        mu, std = norm.fit(pt_ratios_in_bin)
    elif use_eff_rms:
        print("Using eff rms ...")
        # Use the eff_rms method
        eff_rms_vals = effrms(pt_ratios_in_bin) if len(pt_ratios_in_bin) > 1 else [0]
        mu = np.mean(eff_rms_vals) if len(eff_rms_vals) > 1 else 0
        std = np.std(eff_rms_vals) if len(eff_rms_vals) > 1 else 0
    else:
        print("Using Histogram mean and sigma values ...")
        # Use the raw mean and standard deviation from the data
        mu = np.mean(pt_ratios_in_bin) if len(pt_ratios_in_bin) > 0 else 0
        std = np.std(pt_ratios_in_bin) if len(pt_ratios_in_bin) > 0 else 0

    # Calculate sigma/mu ratio
    sigma_mu = std / mu if mu != 0 else 0
    err_sigma_mu = std / (np.sqrt(2 * len(pt_ratios_in_bin) - 2) * mu) if mu != 0 and len(pt_ratios_in_bin) > 1 else 0
    sigma_mu_values.append(sigma_mu)
    err_sigma_mu_values.append(err_sigma_mu)

    # Plot the histogram and, optionally, the Gaussian fit
    plt.hist(pt_ratios_in_bin, bins=30, density=True, alpha=0.6, color='g')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)

    if use_fit:
        # Plot the Gaussian fit
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=2)

    plt.title(f'{"Fit" if use_fit else "Raw"} for pt bin {i + 1}: ({bin_edges[i]:.2f}, {bin_edges[i + 1]:.2f})')
    plt.xlabel('pt_ratio')
    plt.ylabel('Density')
    plt.close()

#################
## RESOLUTION ##
################
    
# Plot the sigma/mu ratio for each bin with error bars
    
plt.style.use(mplhep.style.CMS)
plt.errorbar(
    #range(1, len(sigma_mu_values) + 1),  # Bin numbers (x-axis)
    bin_edges[:-1] + 10,
    sigma_mu_values,  # Sigma/mu values (y-axis)
    yerr=err_sigma_mu_values,  # Error bars
    fmt='o',  # Marker style
    #linestyle='-',  # Line style
    #color='b',  # Line color
    color='#4682B4', # Marker and line color (light blue)
    ecolor='#4682B4',  
    capsize=5  # Length of the error bar caps
)

plt.xticks(np.arange(0, 220, 50))
plt.xlabel(r'$p_{T}^{gen} [GeV]$')
plt.ylabel(r'$\sigma/\mu$')
mplhep.cms.label('Private work', data=True, rlabel=f'{args.particle} PU0')
plt.grid(True)
plt.savefig(f'{output_dir}/scale_plot_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.savefig(f'{output_dir}/scale_plot_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.close()

#################
##    SCALE   ##
################

# Calculate mean and standard deviation of pt_ratios
mean_pt_ratio = np.mean(pt_ratios)
std_pt_ratio = np.std(pt_ratios)


# Now, plot the pt_ratio for all events without binning
plt.figure(figsize=(12, 8))

# Plot the pt_ratio distribution for all events
plt.hist(pt_ratios, bins=50, density=True, alpha=0.6, color='#4682B4')

# Add vertical lines for mean and +/- sigma
plt.axvline(mean_pt_ratio, color='red', linestyle='--', linewidth=1, label=f'Mean: {mean_pt_ratio:.2f}')
plt.axvline(mean_pt_ratio - std_pt_ratio, color='green', linestyle='--', linewidth=1, label=f'-1σ: {mean_pt_ratio - std_pt_ratio:.2f}')
plt.axvline(mean_pt_ratio + std_pt_ratio, color='green', linestyle='--', linewidth=1, label=f'+1σ: {mean_pt_ratio + std_pt_ratio:.2f}')

# Set plot labels and title
plt.xlabel(r'$p_{T}^{reco} / p_{T}^{gen} $')
plt.ylabel('Density')
mplhep.cms.label('Private work', data=True, rlabel=f'{args.particle} PU0')
plt.legend()
plt.grid(True)

# Save the new plot as a PDF
plt.savefig(f'{output_dir}/pt_ratio_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.savefig(f'{output_dir}/pt_ratio_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.close()


# Apply filters to select only matched data
mask = matched# Example filters
filtered_eta_diffs = eta_diffs[mask]
filtered_phi_diffs = phi_diffs[mask]

# Combine filtered_eta_diffs and filtered_phi_diffs into a single array
eta_phi_diff_data = np.column_stack((filtered_eta_diffs, filtered_phi_diffs))

# Save to a new text file
output_file = f'{output_dir}/filtered_eta_phi_diffs_{args.particle}_{args.algo}_{args.subdet}_{args.events}.txt'
np.savetxt(output_file, eta_phi_diff_data, fmt='%.8f', delimiter=',', header='eta_diff,phi_diff', comments='')

print(f"Filtered eta_diff and phi_diff data saved to {output_file}")


#########################################
# Percentage of matched pions per PT bins Binomail counting
#########################################

# Define pt bins (same as used for resolution plot)
pt_bins = np.linspace(0, 200, 21)  # Example: 20 bins between 0 and 200 GeV

# Compute total and matched counts per pt bin
total_counts, _ = np.histogram(gen_pts, bins=pt_bins)
matched_counts, _ = np.histogram(gen_pts[mask], bins=pt_bins)

# Calculate percentage of matched jets
percentage_matched = (matched_counts / total_counts) * 100
print("matched_counts", matched_counts, " |total_counts", total_counts, " |percentage_matched", percentage_matched )


# Replace NaN and inf values (e.g., if a bin has 0 total events)
percentage_matched = np.nan_to_num(percentage_matched)

# Compute error bars based on binomial statistics
lo_err = []
up_err = []

for i in range(len(pt_bins) - 1):
    k = matched_counts[i]  # Number of successes (matched)
    n = total_counts[i]    # Total trials (total)
    
    if n == 0:  # If no events in the bin, set errors to 0
        lo_err.append(0)
        up_err.append(0)
        continue
    
    # Binomial test to get the confidence interval for the proportion
    result = stats.binomtest(k, n, p=k/n)
    ci = result.proportion_ci(confidence_level=0.95)  # 95% confidence interval
    lo_err.append(percentage_matched[i] - ci.low * 100)  # Lower error
    up_err.append(ci.high * 100 - percentage_matched[i])  # Upper error

# Midpoints of pt bins for plotting
pt_bin_centers = (pt_bins[:-1] + pt_bins[1:]) / 2

# Create the plot
plt.figure(figsize=(10, 6))

# Compute the bin width for horizontal error bars
pt_bin_widths = (pt_bins[1:] - pt_bins[:-1]) / 2 

# Plot the percentage with error bars
plt.errorbar(pt_bin_centers, percentage_matched, 
             xerr=pt_bin_widths,
             yerr=[lo_err, up_err], 
             fmt='o', color='#4682B4', label=f'Matched {args.particle} (%)', capsize=3, alpha=0.7)

# Add labels and title
plt.xlabel(r'$p_{T}^{gen} \, [GeV]$', fontsize=14)
plt.ylabel('Percentage of Matched Jets (%)', fontsize=14)

# Add grid and style
plt.grid(alpha=0.5)
mplhep.cms.label('Private work', data=True, rlabel=f'{args.particle} PU0')

# Add legend
plt.legend()

# Save and display the plot
plt.savefig(f'{output_dir}/matched_percentage_pt_bins_with_errors_{args.particle}_{args.algo}_{args.subdet}_{args.events}.png')
plt.savefig(f'{output_dir}/matched_percentage_pt_bins_with_errors_{args.particle}_{args.algo}_{args.subdet}_{args.events}.pdf')
plt.close()