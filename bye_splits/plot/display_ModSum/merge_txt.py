import glob
import os

# Set your folder path
directory_path = "/home/llr/cms/manoni/CMSSW_12_5_2_patch1/src/Hgcal/bye_splits/bye_splits/plot/display_ModSum/output_Jets_txt_16towersPU200_02antikt_TT1/"
output_file = "merged_Jets_16towersPU200antikt02_TT1.txt"

# Expected header (set manually)
header = "event,gen_eta,gen_phi,gen_pt,reco_eta,reco_phi,reco_pt,eta_diff,phi_diff,pt_ratio,matched\n"

# Find all .txt files
file_list = sorted(glob.glob(os.path.join(directory_path, "*.txt")))
print("Files found:", file_list)

# Open output file for writing
with open(output_file, "w") as outfile:
    outfile.write(header)  # Write header once

    for i, file_path in enumerate(file_list):
        with open(file_path, "r") as infile:
            lines = infile.readlines()

            # Skip header in each file (first line)
            content = lines[1:] if lines[0].strip().startswith("event") else lines

            outfile.writelines(content)

print(f"Merged {len(file_list)} files into {output_file}")

