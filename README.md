# Search for H → φφ → 4γ 

 statistical analysis for Beyond-the-Standard-Model search for
**displaced photons from the Higgs boson**: `gg → H → φφ → 4γ`, using CMS
Run‑2 and Run‑3 data.

The hypothetical decays of the Higgs of light scalars `φ`, each of which decays to two
photons. The `φ` can be long‑lived, the photons are **displaced**. The
search is parametrized in both the scalar mass `m_φ` and lifetime `cτ`.
This repo builds the signal and background models, writes datacards, 
runs the fits, and builds plots for the analysis

---

## Overview
- **Signal:** `H → φφ → 4γ`, scanned over `m_φ ∈ {15, 20, 30, 40, 50, 55} GeV`
  and `cτ ∈ {0, 10, 20, 50, 100, 1000} mm`.
- **Observable:** the reconstructed, vertex‑corrected 4γ mass fit in the window **110–140 GeV** 
- **Categories** (by the transverse decay length of the two `φ` candidates,
  threshold `lxy = 50 cm`):
  - `prompt`: both `φ` prompt
  - `asym`: one prompt, one displaced
  - `displaced`: both displaced
- **Background:** data‑driven, modeled with a Bernstein polynomial fit to the
  mass sidebands; signal modeled with a double Crystal Ball.
- **Eras:** Run‑2 `{2017, 2018}` and Run‑3
  `{2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}`.
---

## Structure
```
.
├── ggHparameters.py        # parameters relevant to the analysis in one place
│                           
├── ggHcuts.py              #all cuts used in analysis in one place
│                           
│
├── run_datacard.py         # builds full datacards and loops over years, masses, lifetimes
├── run_summary.py          # makes summary plots for all years, lifetimes, masses
├── UL_vs_mass.py           
│
├── datacard/               #code related to datacard creation  (fitting machinery, building cards, workspaces, etc)
│   ├── ggHdatacardmaker.py 
│   ├── datacardtools.py    
│   ├── ggHfitter.py       
│   └── ggHdatacardworkspace.py 
│
└── plotting/               # all plotters + plotting infrastructure
    ├── ggHcmsstyle.py      
    ├── CMS_lumi.py         
    ├── tdrstyle.py        
    ├── plottingtools.py   
    ├── plot_summary.py     #blinded sideband plot with S+B model in signal region
    ├── plot_blinded_sidebands.py  #used for SF derivation
    ├── plot_postfit.py     #pre fit and post fit plots from combine
    ├── plot_limits.py      
    ├── plot_optimization.py
    ├── plot_lxy.py         
    └── plot_comparison.py 
```

## Usage

### 1. Build datacards & run the fit — `run_datacard.py`

**To run for a single point: (mass,lifetime,era)**

```bash
python3 run_datacard.py -s /path/to/signal_ggH4g.root -b /path/to/data_EGamma.root -m 30 -ct 100 -y 2018 -c1 prompt -c2 asym -c3 displaced
```

This builds each category card, runs Combine, makes the postfit plots, and combines the cards into
`datacard_ggH_4g_m30_ct100_combined_2018.txt`.

**Process an entire run** (loops all configured mass/lifetime/era points and
combines cards automatically):

```bash
python3 run_datacard.py -process_run2 1   # 2017, 2018
python3 run_datacard.py -process_run3 1   # 2022, 2023, 2024
```

### 2. Summary plots — `run_summary.py`

Loops over all mass / lifetime / year points and writes a signal+background
summary plot for each:

```bash
python3 run_summary.py
```

### 3. Expected-limit scan — `UL_vs_mass.py`

Scans `(mass, lifetime, year)`, building cards, running
`AsymptoticLimits`, and collecting the median expected upper limit on `r`.
Edit the point lists at the bottom of the file, then:

```bash
python3 UL_vs_mass.py
```

### 4. Individual plotters

The plotters under `plotting/` import the package, so run them as **modules**
from this directory:

| Command | Produces |
| --- | --- |
| `python3 -m plotting.plot_summary -m 30 -ct 100 -y 2018` | S+B summary plot |
| `python3 -m plotting.plot_blinded_sidebands -m 30 -ct 100 -y 2018` | blinded sideband comparison |
| `python3 -m plotting.plot_limits` | 1D limit plot (`plot_2D()` for the heatmap) |
| `python3 -m plotting.plot_optimization` | category-boundary optimization scan |
| `python3 -m plotting.plot_lxy` | 2D `L_xy` heatmap |
| `python3 -m plotting.plot_comparison` | ggH vs. VH cross-section comparison |

`plot_postfit.plot(...)` is normally driven by `run_datacard.py` rather than run
standalone.

