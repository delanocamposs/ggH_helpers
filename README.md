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
├── run_datacard.py         # builds datacards for all years, masses, lifetimes
├── run_postfit.py          # builds cards + runs combine + makes postfit plots
├── run_summary.py          # makes summary plots for all years, lifetimes, masses
│
├── datacard/               #code related to datacard creation  (fitting machinery, building cards, workspaces, etc)
│   ├── ggHdatacardmaker.py 
│   ├── datacardtools.py    
│   ├── ggHfitter.py       
│   └── ggHdatacardworkspace.py 
│
└── plotting/               # main plotters + plotting infrastructure
    ├── plottingtools.py    
    ├── plot_summary.py     #blinded sideband plot with S+B model in signal region
    ├── plot_postfit.py     #pre fit and post fit plots from combine
    ├── UL_vs_mass.py       #scan of expected upper limits vs mass/lifetime
    ├── style/              #CMS plotting style helpers
    │   ├── ggHcmsstyle.py  
    │   ├── CMS_lumi.py     
    │   └── tdrstyle.py     
    └── studies/            #one-off scripts from internal studies
        ├── plot_blinded_sidebands.py  #used for SF derivation
        ├── plot_limits.py             #1D/2D limit plots
        ├── plot_optimization.py       #category-boundary optimization scan
        ├── plot_lxy.py                #2D L_xy heatmap
        └── plot_comparison.py         #ggH vs VH cross-section comparison
```

## Usage

The three `run_*.py` scripts share the same interface: a single point via
`-m/-ct/-y` (plus `-c1..-c4` categories where relevant), or a whole run via
`-process_run2 1` / `-process_run3 1`. They are deliberately separated —
`run_datacard.py` only builds cards, and the plotting scripts build whatever
they need on top of that.

### 1. Build datacards — `run_datacard.py`

**Builds datacards only** (no Combine fits, no plots). Run Combine yourself
afterwards.

**Single point** (mass, lifetime, era):

```bash
python3 run_datacard.py -s /path/to/signal_ggH4g.root -b /path/to/data_EGamma.root -m 30 -ct 100 -y 2018 -c1 prompt -c2 asym -c3 displaced
```

This writes each category card plus the combined card
`datacard_ggH_4g_m30_ct100_combined_2018.txt`, then prints a green summary of
what was made.

**Process an entire run** (loops all configured mass/lifetime/era points and
combines cards automatically):

```bash
python3 run_datacard.py -process_run2 1   # 2017, 2018
python3 run_datacard.py -process_run3 1   # 2022, 2023, 2024
```

Then set limits / fit manually, e.g.:

```bash
text2workspace.py datacard_ggH_4g_m30_ct100_combined_2018.txt -o ws.root
combine -M AsymptoticLimits ws.root -m 125
```

### 2. Postfit plots — `run_postfit.py`

**Depends on `run_datacard`.** For each point/category it builds the card, runs
Combine (`MultiDimFit` + `FitDiagnostics`), and draws the pre/post-fit plots.
Same interface as `run_datacard.py`:

```bash
python3 run_postfit.py -s /path/to/signal.root -b /path/to/data.root -m 30 -ct 100 -y 2018 -c1 prompt -c2 asym -c3 displaced
python3 run_postfit.py -process_run2 1
python3 run_postfit.py -process_run3 1
```

### 3. Summary plots — `run_summary.py`

Independent of the datacards — summary plots build their own histograms from the
ntuples. Single point or whole run:

```bash
python3 run_summary.py -m 30 -ct 100 -y 2018
python3 run_summary.py -process_run2 1
python3 run_summary.py -process_run3 1
```

### 4. Expected-limit scan — `UL_vs_mass.py`

Scans `(mass, lifetime, year)`, building cards, running
`AsymptoticLimits`, and collecting the median expected upper limit on `r`.
Edit the point lists at the bottom of the file, then:

```bash
python3 -m plotting.UL_vs_mass
```

### 5. Individual plotters

The plotters under `plotting/` import the package, so run them as **modules**
from this directory:

| Command | Produces |
| --- | --- |
| `python3 -m plotting.plot_summary -m 30 -ct 100 -y 2018` | S+B summary plot |

`plot_postfit.plot(...)` is normally driven by `run_datacard.py` rather than run
standalone.

### 6. Internal-study scripts — `plotting/studies/`

One-off scripts kept for reference; also run as modules from this directory:

| Command | Produces |
| --- | --- |
| `python3 -m plotting.studies.plot_blinded_sidebands -m 30 -ct 100 -y 2018` | blinded sideband comparison (SF derivation) |
| `python3 -m plotting.studies.plot_limits` | 1D limit plot (`plot_2D()` for the heatmap) |
| `python3 -m plotting.studies.plot_optimization` | category-boundary optimization scan |
| `python3 -m plotting.studies.plot_lxy` | 2D `L_xy` heatmap |
| `python3 -m plotting.studies.plot_comparison` | ggH vs. VH cross-section comparison |

