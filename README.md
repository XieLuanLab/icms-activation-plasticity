# ICMS Plasticity — Figure Reproduction Code

Code for reproducing the figures in:

**"Learning induces activation-mechanism-dependent neural plasticity in an intracortical microstimulation task"**
Robin Kim, Roy Lycke, Pavlo Zolotavin, Jon Montes, Chong Xie, Lan Luan.
Manuscript under review (2026).

- **Source code:** https://github.com/XieLuanLab/icms-activation-plasticity
- **Figure data (Zenodo):** https://doi.org/10.5281/zenodo.17727762
- **Processed recordings, NWB (DANDI):** https://dandiarchive.org/dandiset/001868 (in preparation)

## Repository layout

```
icms-activation-plasticity/
├── python/      # Ephys & behavior figure scripts (read data/)
├── matlab/      # Two-photon imaging figure scripts (read data/)
├── processing/  # Analysis pipeline (reference; not needed to make figures)
├── data/        # Figure data — download from Zenodo (not in git)
├── output/      # Generated figure panels (created on run)
└── requirements.txt
```

The repository has two parts.

**Figure generation (reproducible).** Download the figure data from Zenodo
into `data/`, then run the scripts in `python/` (ephys & behavior: Fig 4-6,
S1, S3-S9) and `matlab/` (two-photon imaging: Fig 2-3, S2) to regenerate the
figures. Fig S6 (imaging movement control) is plotted by `python/fig_s6/`.

**Analysis pipeline.** `processing/` contains the upstream code that produced the
processed data from the raw recordings, provided for reference. It is not needed
to regenerate the figures and requires the raw recordings (available on request).

## Setup

### Data
Download the figure data from [Zenodo](https://doi.org/10.5281/zenodo.17727762)
and unpack its contents into the `data/` directory.

### Python (Fig 4-6, S1, S3-S9)
```bash
pip install -r requirements.txt
```
Python 3.12 (the version used for the published figures).

### MATLAB (Fig 2-3, S2)
Requires MATLAB R2023a or later. Run from the repository root after setting
paths via `matlab/config.m` (it points at the imaging data in `data/matlab/`).
Output is written to `output/`.

## Generating figures

Python scripts are run as modules from the repository root:
```bash
python -m python.fig4.modulation_over_time
python -m python.fig5.generate_figure5
python -m python.fig_s9.generate_figure
```

## Notes

- Data tiers: figure data (processed tables read by the scripts) on Zenodo;
  processed neural recordings (NWB) on DANDI; raw recordings (NS5/NEV) are not
  publicly archived due to size and are available on request.
- A few panels need data not on Zenodo (e.g. Fig 1C-E and 1F representative
  traces, Fig 3A ROI overlays); those scripts are included for reference and the
  data is available on request.
- The reference pipelines in `processing/` use lab-server paths and external
  tools and are not required to regenerate the figures.
