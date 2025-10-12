## Introduction

This package supports performance studies for the Stage 1 TPG (Trigger Primitive Generation) system of HGCal, aiming to improve object reconstruction by refining the use of Module Sums and Trigger Cells (TCs) / Super Trigger Cells (STCs) in building trigger towers.

Trigger towers are static, geometrically defined objects spanning a fixed 20×72 bin grid in the 𝜂–φ plane, covering the region 𝜂 ∈ [1.5, 3] with a granularity of about 5° in φ and ≈0.08 in 𝜂. This granularity closely matches the existing CMS calorimeter trigger towers.

This project focuses on:

- Performing comprehensive performnace studies (mainly energy and position resolution) across multiple particle types—**photons, pions, jets**—under **zero and 200 pile-up conditions**, along with **minimum bias**, samples for trigger rate analysis.
- Utilizing the full V16 geometry configuration, including both the CEE (silicon layers) and CEH (silicon + scintillator) subdetectors.

By investigating how Module Sums and TCs/STCs can be combined and mapped into trigger towers, this workflow aims to provide a robust foundation for future firmware implementations and performance optimizations in the CMS HGCAL trigger system.

---

## 1  Overview

`mainModuleSums.py` orchestrates the full chain:

1. Ingests ROOT ntuples
2. Assigns Module sums/STCs to fixed η–φ tower grids using several algorithms (baseline, area‑overlap, 4/8/16‑tower splitting logics).
4. Matches reconstructed particles to generator particles, evaluates energy and position resolution in η/φ for performance studies.

The script relies on helper modules (`processingMS`, `resolutionMS`, `geometryMS`, `helperMS`) located in the same package.

---

## 2  Quick Start

```bash
# Example 1 – all photon events, 8‑tower algorithm
python mainModuleSums.py \
    --event -1 \
    --geom V16 \
    --algo 8towers \
    --particle photons \
    --subdet 5 \
    --inputfile root \

# Example 2 – 4 000 pion events, enable STCs, 16‑tower algorithm
python mainModuleSums.py \
    --event -1 --n 4000 \
    --geom V16 \
    --algo 16towers \
    --particle pions \
    --subdet 5 \
    --inputfile root \
    --STCs
```

### 2.1  CLI Arguments

| Flag                 | Default    | Description                                                                         |
| -------------------- | ---------- | ----------------------------------------------------------------------------------- |
| `--subdet`           | `1`        | 1 = CEE (all‑Si), 2 = CEH‑Si, 3 = CEH‑Scint, 4 = CEH (Si + Scint), 5 = CEE + CEH (all layers)   |
| `--event`            | `5492`     | Single event ID, or `-1` to process the full file                                   |
| `--n`                | *None*     | Process *n* randomly selected events                                                |
| `--algo`             | `8towers`  | Tower‑building scheme: `baseline`, `area_overlap`, `4towers`, `8towers`, `16towers` |
| `--particle`         | `photons`  | Choose among `photons`, `pions`, `jets`, `neutrinos`, `jets`                                |
| `--geom`             | `V16`      | CMSSW geometry tag: `V11` (2021) or `V16` (2023)                                    |
| `--inputfile`        | `root`     | Input container: `root` or `hdf5` (hdf5 handling is outdated)                                                   |
| `--STCs / --no-STCs` | *disabled* | Toggle Super Trigger Cells logic                                                |

---

## 3  Prerequisites

* **Python ≥ 3.9**
* Tested on `CMSSW_12_5_2_patch1`release
* To be run on alma9 an el7  container is needed
  * Instructions to run the container:\

    1) Launch the el7 container:\
    `/data_CMS/cms/manoniL1HGCAL/el7_container`\
    This drops you into a Bash shell inside el7.\

    2) Inside the container, run:\
    `source /cvmfs/cms.cern.ch/cmsset_default.sh`\
    `cd /home/llr/cms/manoni/CMSSW_12_5_2_patch1/src`\
    `cmsenv`\

You now have a proper el7 environment with CMSSW ready.

* PyPI packages (listed in `requirements.txt`):

  * pandas, numpy, uproot, awkward, shapely, geojson, matplotlib, tqdm, etc.

---

## 4  Repository Layout

```
HGCalModuleSums/
├── mainModuleSums.py          # Entry point
├── processingMS.py            # I/O + high‑level processing
├── resolutionMS.py            # Gen‑reco matching & resolution computation
├── geometryMS.py              # η–φ grid & GeoJSON helpers
├── helperMS.py                # Misc utilities
├── requirements.txt
└── README.md                  # You are here
```

---
## 5  Inputs

* **Input ROOT files**: 
Input ROOT files can be found in

  * `/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples/` for samples containing only Modules Sums
  * `/data_CMS/cms/manoni/L1HGCAL/final_skimmed_V16ntuples_STCS/` for samples containing Module Sums and STCs.

Samples used at the moment are hardcoded in `mainModulesSums.py`

---

## 6  Outputs

* txt files with the following header:
`event,gen_eta,gen_phi,gen_pt,reco_eta,reco_phi,reco_pt,eta_diff,phi_diff,pt_ratio, matched`

These `.txt` files can then be used as inputs to the plotter scripts to produce the final performance plots.

---
