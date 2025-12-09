import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2
import mplhep

# --- CMS style ---
mplhep.style.use("CMS")

# ------------------------------
# 1. Load the input file
# ------------------------------
filename = "merged_Jets_16towersPU200_Delta_Match02_Deltakt_04_200ntuples.txt"
df = pd.read_csv(filename)

df["matched"] = df["matched"].astype(str).str.lower() == "true"
df_matched = df[df["matched"]].copy()


# ------------------------------
# 2. Compute beta(η)
# ------------------------------
def compute_beta(row):
    pt_raw = row["reco_pt"]
    pt_gen = row["gen_pt"]
    rho = row["rho"]
    area = row["jet_area"]

    if rho * area <= 0:
        return np.nan
    return (pt_raw - pt_gen) / (rho * area)

df_matched.loc[:, "beta"] = df_matched.apply(compute_beta, axis=1)
df_matched = df_matched.dropna(subset=["beta"])

eta = df_matched["gen_eta"].values
beta = df_matched["beta"].values

# ------------------------------
# 3. Bin η and compute mean & uncertainty
# ------------------------------
eta_bins = np.linspace(1.6, 2.9, 20)
bin_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])

counts, edges = np.histogram(df_matched["gen_eta"], bins=eta_bins)
for i, c in enumerate(counts):
    print(f"Bin {edges[i]:.2f}-{edges[i+1]:.2f}: {c} jets")

beta_means = []
beta_errors = []

for i in range(len(eta_bins)-1):
    mask = (eta >= eta_bins[i]) & (eta < eta_bins[i+1])
    vals = beta[mask]

    if len(vals) > 3:
        
        mean = np.mean(vals)
        #mean= np.median(vals) 
        std = np.std(vals)
        err = std / np.sqrt(len(vals))

        #print("len(vals)", len(vals), " mean", f"{mean:.3f}", " median", f"{median:.3f}", " err", f"{err:.3f}", " std", f"{std:.3f}" )
        #print("-------------------------------------")
    else:
        mean = np.nan
        err = np.nan

    beta_means.append(mean)
    beta_errors.append(err)

beta_means = np.array(beta_means)
beta_errors = np.array(beta_errors)

mask_valid = ~np.isnan(beta_means)
x = bin_centers[mask_valid]
y = beta_means[mask_valid]
yerr = beta_errors[mask_valid]

# ------------------------------
# 4. Polynomial fit
# ------------------------------
def poly3(x, p0, p1, p2, p3):
    return p0 + p1*x + p2*x**2 + p3*x**3

popt, pcov = curve_fit(poly3, x, y, sigma=yerr, absolute_sigma=True)
p0, p1, p2, p3 = popt
perr = np.sqrt(np.diag(pcov))

# ------------------------------
# 5. Plot β(η) with CMS styling
# ------------------------------
x_fit = np.linspace(1.6, 2.9, 400)
y_fit = poly3(x_fit, *popt)

plt.figure(figsize=(8,6))

plt.errorbar(
    x, y, yerr=yerr,
    fmt='o', markersize=4,
    color='black', ecolor='gray', capsize=3,
    label='β(η) mean per η-bin'
)

plt.plot(x_fit, y_fit, color='red', linewidth=2,
         label='3rd-order polynomial fit')

plt.xlabel("gen jet |η|", fontsize=15)
plt.ylabel(r"$\beta(\eta) = (p_{T,raw} - p_{T,gen}) / (\rho A)$", fontsize=15)
#plt.title("β(η) fit", fontsize=18)
plt.grid(True, alpha=0.3)

# ---- Fit result box ----
chi2_val = np.sum(((y - poly3(x, *popt)) / yerr)**2)
ndf = len(x) - len(popt)
chi2_ndf = chi2_val / ndf

textbox = (
    rf"$\chi^2/\mathrm{{ndf}} = {chi2_ndf:.2f}$" + "\n"
    rf"$p_0 = {p0:.3f} \pm {perr[0]:.3f}$" + "\n"
    rf"$p_1 = {p1:.3f} \pm {perr[1]:.3f}$" + "\n"
    rf"$p_2 = {p2:.3f} \pm {perr[2]:.3f}$" + "\n"
    rf"$p_3 = {p3:.3f} \pm {perr[3]:.3f}$"
)

plt.text(
    0.02, 0.85, textbox,
    transform=plt.gca().transAxes,
    fontsize=9,
    verticalalignment='top',
    bbox=dict(facecolor='white', alpha=0.9)
)

plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=14)
mplhep.cms.text("Preliminary", loc=0, fontsize=16) 
plt.tight_layout()
plt.savefig('Beta_fit_200.pdf')
#plt.show()


# ----------------------------------------------------------
# 6. Apply the η-dependent offset correction
# ----------------------------------------------------------
def beta_eta(eta):
    return poly3(eta, *popt)

eta_vals = df_matched["gen_eta"].values
rho_vals = df_matched["rho"].values
area_vals = df_matched["jet_area"].values
pt_raw_vals = df_matched["reco_pt"].values
pt_gen_vals = df_matched["gen_pt"].values

beta_vals = beta_eta(eta_vals)

C_offset = 1 - (rho_vals * area_vals * beta_vals) / pt_raw_vals
pt_corr = pt_raw_vals * C_offset

df_matched.loc[:, "C_offset"] = C_offset
df_matched.loc[:, "pt_corr"] = pt_corr

print(df_matched[["gen_pt", "reco_pt", "pt_corr", "C_offset"]].head())

# ----------------------------------------------------------
# 7. Jet response before/after
# ----------------------------------------------------------
response_raw = pt_raw_vals / pt_gen_vals
response_corr = pt_corr / pt_gen_vals

# Compute mean and std
mean_raw, std_raw = np.mean(response_raw), np.std(response_raw)
mean_corr, std_corr = np.mean(response_corr), np.std(response_corr)

# Define common bins
bins = np.linspace(0, 5, 51)  # 100 bins from 0 to 5

plt.figure(figsize=(8,6))

plt.hist(response_raw, bins=bins, density=True, alpha=0.55,
         label="Raw response")
plt.hist(response_corr, bins=bins, density=True, alpha=0.55,
         label="Corrected response")

plt.xlabel(r"$p_{T}^{reco} / p_{T}^{gen}$", fontsize=15)
plt.ylabel("Density", fontsize=15)
#plt.title("Response before/after offset correction", fontsize=18)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.xlim(0, 5)
plt.tick_params(axis='both', which='major', labelsize=14)
mplhep.cms.text("Preliminary", loc=0, fontsize=16) 


# ---- Add statistics textbox ----
stats_text = (
    f"Raw: μ = {mean_raw:.3f}, σ = {std_raw:.3f}\n"
    f"Corrected: μ = {mean_corr:.3f}, σ = {std_corr:.3f}"
)
plt.text(
    0.65, 0.75, stats_text,
    transform=plt.gca().transAxes,
    fontsize=10,
    bbox=dict(facecolor='white', alpha=0.9)
)


plt.tight_layout()
plt.savefig('Response_200.pdf')
#plt.show()

# ----------------------------------------------------------
# 8. Resolution vs η
# ----------------------------------------------------------
eta_bins = np.linspace(1.6, 2.9, 15)
eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])

raw_res = []
corr_res = []

for i in range(len(eta_bins)-1):
    mask = (eta_vals >= eta_bins[i]) & (eta_vals < eta_bins[i+1])
    r_raw = response_raw[mask]
    r_corr = response_corr[mask]

    raw_res.append(np.std(r_raw) if len(r_raw) > 5 else np.nan)
    corr_res.append(np.std(r_corr) if len(r_corr) > 5 else np.nan)

plt.figure(figsize=(8,6))
plt.plot(eta_centers, raw_res, "-o", label="Raw resolution")
plt.plot(eta_centers, corr_res, "-o", label="Corrected resolution")

plt.xlabel("gen jet |η|", fontsize=15)
plt.ylabel("RMS", fontsize=15)
#plt.title("Resolution vs η", fontsize=18)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=14)
mplhep.cms.text("Preliminary", loc=0,  fontsize=16) 
plt.tight_layout()
plt.savefig('RMS_200.pdf')
#plt.show()

#--------------------------------------

df_corrected = df_matched.copy()
print(df_corrected.head())

# ----------------------------------------------------------
# 9. Save new file with only matched jets and corrected pt
# ----------------------------------------------------------
# Copy only matched jets
df_corrected = df_matched.copy()

# Add pt_corr column (already computed)
# reco_pt remains untouched
# pt_corr is already in df_matched

# Construct new filename
filename_corrected = filename.replace(".txt", "_OffsetCorrected.txt")

# Save to new txt file
df_corrected.to_csv(filename_corrected, index=False)
print(f"Corrected file with only matched jets saved as: {filename_corrected}")



# ----------------------------------------------------------
# 10. η distribution: all jets vs matched jets
# ----------------------------------------------------------
plt.figure(figsize=(8,6))

# All jets
plt.hist(df["gen_eta"], bins=30, alpha=0.5, label="All jets", color='blue', density=False)

# Matched jets
plt.hist(df_matched["gen_eta"], bins=30, alpha=0.5, label="Matched jets", color='red', density=False)

plt.xlabel("gen jet |η|", fontsize=15)
plt.ylabel("Entries", fontsize=15)
#plt.title("η distribution: All vs Matched jets", fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=14)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.tight_layout()
#plt.savefig('Eta_distribution.pdf')
plt.show()