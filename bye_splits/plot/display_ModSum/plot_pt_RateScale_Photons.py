import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import mplhep
import os
import argparse

# Argument parser for customization
parser = argparse.ArgumentParser(description="Generate analysis plots with custom settings.")
parser.add_argument("--algo", type=str, default="16towers", help="Algorithm used (e.g., 16t, DNN).")
parser.add_argument("--subdet", type=str, default="CEE_CEH", help="Subdetector type (e.g., Jets, Calorimeter).")
parser.add_argument("--events", type=str, default="5k", help="Events")
args = parser.parse_args()

# Directory to save all plots
output_dir = f"plots_PHOTONS_NEW_SUBWIND_{args.algo}_{args.subdet}_{args.events}"
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
with open('pt_scale_&_res_8towers_5.txt', 'r') as f:
    # Read the header line
    header = f.readline().strip().split(', ')
    print("Header:", header)  # Optional: Print the header to check if it's correct

    # Read the rest of the data
    for line in f:
        values = line.strip().split(', ')
        print(len(values))
        if len(values) >= 9:
            # Convert the relevant columns to floats
            best_subwindow_pt = float(values[0])
            genpart_pt = float(values[1])
            pt_ratio = float(values[2])
            eta_diff = float(values[3])
            phi_diff = float(values[4])
            particle_eta = float(values[5])
            particle_phi = float(values[6])
            reco_eta = float(values[7])
            reco_phi = float(values[8])
            event = float(values[9])
            data.append((best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, reco_eta, reco_phi, event))
        else:
            print("line boh")
# Convert the data into a numpy array for easier processing
data = np.array(data)

# Filter out events where pt_ratio > 2 or pt_ratio < -2
#pt_ratios = data[:, 2]
#mask = (pt_ratios < 2)
#data = data[mask]


# Extract the genpart_pt and pt_ratio columns
genpart_pts = data[:, 1]
pt_ratios = data[:, 2]
eta_diffs = data[:, 3]
phi_diffs = data[:, 4]
particle_etas = data[:, 5]
particle_phis = data[:, 6]
reco_etas = data[:, 7]
reco_phis = data[:, 8]


# Filter out events where pt_ratio > 2
'''mask = (pt_ratios >= 0.75) & (pt_ratios <= 2)

# Apply the mask to extract the corresponding eta_diff and phi_diff
filtered_eta_diffs = eta_diffs[mask]
filtered_phi_diffs = phi_diffs[mask]

# Write the filtered data to a new file
with open('eta_phi_diff_filtered_pt_ratio_ab2_bl075.txt', 'w') as output_file:
    # Write the header
    output_file.write("eta_diff, phi_diff\n")
    
    # Write each filtered eta_diff and phi_diff pair
    for eta, phi in zip(filtered_eta_diffs, filtered_phi_diffs):
        output_file.write(f"{eta}, {phi}\n")

print("Filtered data written to eta_phi_diff_filtered_pt_ratio_ab2_bl075.5.txt")'''


# Define the number of bins for pt
#num_bins = 10
#bin_edges = np.linspace(np.min(genpart_pts), np.max(genpart_pts), num_bins + 1)
bin_edges = np.arange(0, 220, 20)

# Prepare to store results for each bin
sigma_mu_values = []
err_sigma_mu_values = []

# Loop over each bin and process the data
for i in range(len(bin_edges) - 1):
    # Find the events that fall into the current bin
    bin_mask = (genpart_pts >= bin_edges[i]) & (genpart_pts < bin_edges[i + 1])

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
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.grid(True)
plt.savefig(os.path.join(output_dir,f'scale_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'scale_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()

#################
##    SCALE   ##
################
print("Values where pt_ratio < 0.5:")
for best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, reco_eta, reco_phi, event in data:
    if pt_ratio < 0.5:
        #print(f"best_subwindow_pt: {best_subwindow_pt}, genpart_pt: {genpart_pt}, pt_ratio: {pt_ratio}, eta_diff :{eta_diff}, phi_diff :{phi_diff}, part eta:{particle_eta}, part phi: {particle_phi}  ")
        print(f"best_subwindow_pt: {best_subwindow_pt}, genpart_pt: {genpart_pt}, pt_ratio: {pt_ratio}, eta_diff :{eta_diff}, phi_diff :{phi_diff}, part eta:{particle_eta}, part phi: {particle_phi},reco_eta: {reco_eta}, reco_phi: {reco_phi}, event: {event} ")

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
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.legend()
plt.grid(True)

# Save the new plot as a PDF
plt.savefig(os.path.join(output_dir,f'pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()

#################################
## OUTLIERS ETA DISTRIBUTION   ##
#################################

#Now, plot the particle_eta distribution for pt_ratio < 0.5
plt.figure(figsize=(12, 8))

# Filter the data for pt_ratio < 0.5
#pt_ratio_mask = pt_ratios < 0.5
filtered_particle_etas = particle_etas #[pt_ratio_mask]

# Plot the particle_eta distribution for pt_ratio < 0.5
plt.hist(filtered_particle_etas, bins=50, density=False, alpha=0.6, color='#4682B4')

# Set plot labels and title
plt.xlabel(r'$\eta$') #(for $p_{T}^{reco} / p_{T}^{gen} < 0.5$)
plt.ylabel('Counts')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.legend()
plt.grid(True)

# Save the particle_eta distribution plot as a PDF
plt.savefig(os.path.join(output_dir,f'particle_eta_distribution_ALL_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'particle_eta_distribution_ALL_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()


## OUTLIERS PHI DISTRIBUTION ##
#################################

# Now, plot the particle_phi distribution for pt_ratio < 0.5
plt.figure(figsize=(12, 8))

# Filter the data for pt_ratio < 0.5
#pt_ratio_mask = pt_ratios < 0.5
filtered_particle_phis = particle_phis #[pt_ratio_mask]

# Plot the particle_phi distribution for pt_ratio < 0.5
plt.hist(filtered_particle_phis, bins=50, density=False, alpha=0.6, color='#4682B4')

# Set plot labels and title
plt.xlabel(r'$\phi$ ') #(for $p_{T}^{reco} / p_{T}^{gen} < 0.5$)
plt.ylabel('Counts')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.legend()
plt.grid(True)

# Save the particle_phi distribution plot as a PDF
plt.savefig(os.path.join(output_dir,f'particle_phi_distribution_ALL_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'particle_phi_distribution_ALL_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()



## OUTLIERS PT DISTRIBUTION ##
#################################

# Now, plot the particle_phi distribution for pt_ratio < 0.5
plt.figure(figsize=(12, 8))

plt.hist(genpart_pts, bins=50, density=False, alpha=0.6, color='#4682B4')

# Set plot labels and title
plt.xlabel(r'$p_{T}$') #(for $p_{T}^{reco} / p_{T}^{gen} < 0.5$)
plt.ylabel('Counts')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.legend()
plt.grid(True)

# Save the particle_phi distribution plot as a PDF
plt.savefig(os.path.join(output_dir,f'particle_pT_distribution_pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'particle_pT_distribution_pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()


##############################################
## CORRELATION PLOT: ETA vs PHI             ##
##############################################

# Create a scatter plot for filtered particle_eta vs particle_phi
plt.figure(figsize=(12, 8))

# Plot a scatter plot
plt.scatter(particle_etas, particle_phis, alpha=0.6, color='#4682B4', s=10, edgecolor='k')

# Set plot labels and title
plt.xlabel(r'$\eta$ ') #(for $p_{T}^{reco} / p_{T}^{gen} < 0.5$)
plt.ylabel(r'$\phi$ ') #(for $p_{T}^{reco} / p_{T}^{gen} < 0.5$)
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.grid(True)

# Save the scatter plot as a PDF and PNG
plt.savefig(os.path.join(output_dir,f'eta_phi_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'eta_phi_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()


##############################################
## CORRELATION PLOT: PHI vs PHI_DIFF       ##
##############################################

# Create a scatter plot for filtered genpart_pt vs phi_diff
plt.figure(figsize=(12, 8))

# Plot a scatter plot
plt.scatter(particle_phis, phi_diffs, alpha=0.6, color='#4682B4', s=10, edgecolor='k')

# Set plot labels and title
plt.xlabel(r'$\phi gen$')
plt.ylabel(r'$\Delta \phi$ ')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.grid(True)

# Save the scatter plot as a PDF and PNG
plt.savefig(os.path.join(output_dir,f'phi_phi_diff_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'phi_phi_diff_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()


##############################################
## CORRELATION PLOT:PHI RECO/GEN            ##
##############################################

# Create a scatter plot for filtered genpart_pt vs phi_diff
plt.figure(figsize=(12, 8))

# Plot a scatter plot
plt.scatter(particle_phis, reco_phis, alpha=0.6, color='#4682B4', s=10, edgecolor='k')

# Set plot labels and title
plt.xlabel(r'$\phi gen')
plt.ylabel(r'$\phi$ reco ')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.grid(True)

# Save the scatter plot as a PDF and PNG
plt.savefig(os.path.join(output_dir,f'phi_gen_reco_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'phi_gen_reco_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()

##############################################
## CORRELATION PLOT:ETA RECO/GEN            ##
##############################################

# Create a scatter plot for filtered genpart_pt vs phi_diff
plt.figure(figsize=(12, 8))

# Plot a scatter plot
plt.scatter(particle_etas, reco_etas, alpha=0.6, color='#4682B4', s=10, edgecolor='k')

# Set plot labels and title
plt.xlabel(r'$\eta$ gen')
plt.ylabel(r'$\eta$ reco ')
mplhep.cms.label('Private work', data=True, rlabel='Photons PU0')
plt.grid(True)

# Save the scatter plot as a PDF and PNG
plt.savefig(os.path.join(output_dir,f'eta_gen_reco_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.pdf'))
plt.savefig(os.path.join(output_dir,f'eta_gen_reco_correlation_plot_pt_ratio_{args.algo}_{args.subdet}_{args.events}.png'))
plt.close()


# Define the thresholds
eta_diff_threshold = 0.2
phi_diff_threshold = 0.2

# Open two separate files to save the events
'''with open('eta_diff_outliers.txt', 'w') as eta_file, open('phi_diff_outliers.txt', 'w') as phi_file:
    # Write headers to the files
    eta_file.write("best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, reco_eta, reco_phi, event\n")
    phi_file.write("best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, reco_eta, reco_phi, event\n")

    # Iterate through the data and filter events based on the thresholds
    for best_subwindow_pt, genpart_pt, pt_ratio, eta_diff, phi_diff, particle_eta, particle_phi, reco_eta, reco_phi, event in data:
        if eta_diff > eta_diff_threshold or eta_diff < -eta_diff_threshold:
            # Write the event details to the eta_diff file
            eta_file.write(f"{best_subwindow_pt}, {genpart_pt}, {pt_ratio}, {eta_diff}, {phi_diff}, {particle_eta}, {particle_phi}, {reco_eta}, {reco_phi}, {event}\n")
        
        if phi_diff > phi_diff_threshold or phi_diff < -phi_diff_threshold:
            # Write the event details to the phi_diff file
            phi_file.write(f"{best_subwindow_pt}, {genpart_pt}, {pt_ratio}, {eta_diff}, {phi_diff}, {particle_eta}, {particle_phi}, {reco_eta}, {reco_phi}, {event}\n")

print("Events with significant eta_diff or phi_diff have been saved to separate files.")'''