import pandas as pd
import ROOT
from array import array
import matplotlib.pyplot as plt
import numpy as np
import mplhep

mplhep.style.use("CMS")

# -------------------------
# 1. Load CSV with offset-corrected pT
# -------------------------
df = pd.read_csv("merged_Jets_16towersPU200_Delta_Match02_Deltakt_04_200ntuples_OffsetCorrected.txt")
df = df[df['matched'] == True]  # Only use matched jets
print(f"Total jets loaded: {len(df)}")
print(f"Matched jets used: {len(df)}")

df['response'] = df['pt_corr'] / df['gen_pt']

# --- Plot corrected pT and response ---
plt.figure(figsize=(9,6))
plt.hist(df['pt_corr'], bins=50, alpha=0.6, label='Offset-corrected pT (NO MC correction applied yet)', color='tab:blue')
plt.hist(df['gen_pt'], bins=50, alpha=0.6, label='Gen pT', color='tab:orange')
plt.xlabel('Jet pT [GeV]', fontsize=14)
plt.ylabel('Counts', fontsize=14)
#plt.title('Offset-corrected Jet pT distributions', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.tight_layout()
plt.show()

# -------------------------
# 2. Create TProfile2D (Response Map)
# -------------------------
pt_bins, pt_min, pt_max = 10, 0, 200
eta_bins, eta_min, eta_max = 10, 1.7, 2.8

tprof = ROOT.TProfile2D("tprof", "Jet Response", pt_bins, pt_min, pt_max,
                        eta_bins, eta_min, eta_max, 0, 2)
for index, row in df.iterrows():
    tprof.Fill(row['gen_pt'], row['gen_eta'], row['pt_corr'] / row['gen_pt'])

# --- Plot TProfile2D as 2D heatmap ---
response_array = np.zeros((pt_bins, eta_bins))
pt_centers = [tprof.GetXaxis().GetBinCenter(i) for i in range(1, pt_bins+1)]
eta_centers = [tprof.GetYaxis().GetBinCenter(j) for j in range(1, eta_bins+1)]

for i in range(1, pt_bins+1):
    for j in range(1, eta_bins+1):
        response_array[i-1, j-1] = tprof.GetBinContent(i, j)

plt.figure(figsize=(9,6))
pcm = plt.pcolormesh(eta_centers, pt_centers, response_array, shading='auto', cmap='viridis')
plt.colorbar(pcm, label='Average Response (pt_corr / gen_pt)')
plt.xlabel('Gen jet eta', fontsize=14)
plt.ylabel('Gen jet pT [GeV]', fontsize=14)
#plt.title('2D TProfile2D Response Map (Offset-corrected)', fontsize=16)
plt.tick_params(axis='both', labelsize=12)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.tight_layout()
plt.show()

# -------------------------
# 3. Convert TProfile2D -> TGraph2D
# -------------------------
x_list, y_list, z_list = [], [], []
for i in range(1, tprof.GetNbinsX()+1):
    for j in range(1, tprof.GetNbinsY()+1):
        c = tprof.GetBinContent(i,j)
        if c == 0: continue
        x_list.append(c * tprof.GetXaxis().GetBinCenter(i))
        y_list.append(tprof.GetYaxis().GetBinCenter(j))
        z_list.append(c)

tgraph = ROOT.TGraph2D(len(x_list), array('d', x_list), array('d', y_list), array('d', z_list))

'''plt.figure(figsize=(9,6))
plt.scatter(x_list, y_list, c=z_list, cmap='viridis', s=40)
plt.colorbar(label='Response')
plt.xlabel('Average Offset-corrected pT [GeV]', fontsize=14)
plt.ylabel('Gen jet eta', fontsize=14)
plt.title('TGraph2D Points (Offset-corrected Response)', fontsize=16)
plt.tick_params(axis='both', labelsize=12)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.tight_layout()
plt.show()'''

# -------------------------
# 4. Interpolate and compute MC-corrected pT
# -------------------------
mc_corrected = []
for _, row in df.iterrows():
    response = tgraph.Interpolate(row['pt_corr'], row['gen_eta'])
    if response==0: response=1.0
    mc_corrected.append(row['pt_corr']/response)
df['pT_MC_corrected'] = mc_corrected

# -------------------------
# 5. 2D heatmap axes swapped
# -------------------------
xi = np.linspace(min(y_list), max(y_list), 200)
yi = np.linspace(min(x_list), max(x_list), 200)
X, Y = np.meshgrid(xi, yi)
Z = np.zeros_like(X)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        val = tgraph.Interpolate(Y[i,j], X[i,j])
        Z[i,j] = val if val!=0 else np.nan

plt.figure(figsize=(9,6))
pcm = plt.pcolormesh(X, Y, Z, shading='auto', cmap='jet', vmin=0, vmax=1.25)  # <-- jet colormap
cbar = plt.colorbar(pcm)
cbar.set_label(r'Jet Response ($p_T^{\mathrm{corr}} / p_T^{\mathrm{gen}}$)', fontsize=14)
plt.scatter(y_list, x_list, color='k', s=20, label='')
plt.xlabel(r'gen jet $\eta$', fontsize=14)
plt.ylabel(r'Offset-corrected jet $p_T$ [GeV]', fontsize=14)
#plt.title('MC Jet Response Heatmap with TGraph2D Points (Axes Swapped)', fontsize=16)
plt.legend(fontsize=12)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.tight_layout()
plt.savefig('HeatMap_200.png')
plt.savefig('HeatMap_200.pdf')
plt.show()

# -------------------------
# 6. Plot comparison of responses: raw / offset / MC-corrected
# -------------------------
plt.figure(figsize=(9,6))

# Compute responses
response_raw = df['reco_pt'] / df['gen_pt']          # raw pT
response_offset = df['pt_corr'] / df['gen_pt']       # offset-corrected
response_MC = df['pT_MC_corrected'] / df['gen_pt']   # offset+MC corrected

# Plot histograms
bins = np.linspace(0, 5, 50)
plt.hist(response_raw, bins=bins, alpha=0.5, label='Raw', color='tab:blue', density=True)
plt.hist(response_offset, bins=bins, alpha=0.5, label='Offset-corrected', color='tab:orange', density=True)
plt.hist(response_MC, bins=bins, alpha=0.5, label='Offset+MC-corrected', color='tab:green', density=True)


# Compute mean and std for each
mean_raw, std_raw = np.mean(response_raw), np.std(response_raw)
mean_offset, std_offset = np.mean(response_offset), np.std(response_offset)
mean_MC, std_MC = np.mean(response_MC), np.std(response_MC)

# Text box with statistics
textbox = (
    f"Raw: {mean_raw:.3f} ± {std_raw:.3f}\n"
    f"Offset-corrected:: {mean_offset:.3f} ± {std_offset:.3f}\n"
    f"Offset + MC-corrected: {mean_MC:.3f} ± {std_MC:.3f}"
)
plt.text(0.55, 0.60, textbox, transform=plt.gca().transAxes,
         fontsize=12, bbox=dict(facecolor='white', alpha=0.8))

plt.xlabel(r'$p_T / p_T^{\mathrm{gen}}$', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)
mplhep.cms.text("Preliminary", loc=0, fontsize=16)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('Response_comparison_200.png')
plt.savefig('Response_comparison_200.pdf')
plt.show()

# -------------------------
# 7. Save before/after pT
# -------------------------
before_after_df = df[['pt_corr','pT_MC_corrected']]
before_after_df.to_csv("merged_Jets_16towersPU200_Delta_Match02_Deltakt_04_MC_Corrected_200.txt", index=False)
print("Saved offset-corrected and MC-corrected pT file.")