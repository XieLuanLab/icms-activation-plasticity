# ICMS Plasticity - Figure Reproduction Code

Code for reproducing the figures in:

**"Learning induces activation-mechanism-dependent neural plasticity in an intracortical microstimulation task"**
Robin Kim, Roy Lycke, Pavlo Zolotavin, Jon Montes, Chong Xie, Lan Luan.
Science Advances (2026), accepted.

- **Source code:** https://github.com/XieLuanLab/icms-activation-plasticity
- **Figure source data (Zenodo):** https://doi.org/10.5281/zenodo.21382755
- **Processed recordings, NWB (DANDI):** https://dandiarchive.org/dandiset/001868

## Repository layout

```
icms-activation-plasticity/
|-- python/      # Ephys and behavior figure scripts (read data/)
|-- matlab/      # Two-photon imaging figure scripts (read data/)
|-- processing/  # Analysis pipeline (reference; not needed to make figures)
|-- data/        # Figure source data - download from Zenodo (not in git)
|-- output/      # Generated figure panels (created on run)
`-- requirements.txt
```

The repository has two parts.

**Figure generation (reproducible).** Download the figure source data from
Zenodo and place the extracted `data/` directory at the repository root. Then run
the scripts in `python/` (ephys and behavior: Fig 4-6, S1, S3-S9) and `matlab/`
(two-photon imaging: Fig 2-3, S2) to regenerate figure panels.

**Analysis pipeline.** `processing/` contains the upstream code that produced the
processed data from the raw recordings, provided for reference. It is not needed
to regenerate the figures and requires the raw recordings (available on request).

## Setup

### Data
Download the figure source data from
[Zenodo](https://doi.org/10.5281/zenodo.21382755). Unpack the archive so that
the repository contains a top-level `data/` directory.

### Python (Fig 4-6, S1, S3-S9)
```bash
pip install -r requirements.txt
```
Python 3.12 was used for the accepted manuscript figures.

### MATLAB (Fig 2-3, S2)
Requires MATLAB R2023a or later. Run from the repository root after setting
paths via `matlab/config.m` (it points at the imaging data in `data/matlab/`).
Output is written to `output/`.

## Generating figures

Python scripts are run as modules from the repository root:
```bash
python -m python.fig4.modulation_over_time
python -m python.fig5.generate_figure5
python -m python.fig_s5.generate_figure
python -m python.fig_s9.generate_figure
```

## Notes

- Data tiers: figure source data (processed tables and arrays read directly by
  these scripts) are on Zenodo; processed neural recordings in NWB format are on
  DANDI; raw recordings (NS5/NEV) are not publicly archived due to size and are
  available on request.
- A few panels need data not on Zenodo (e.g. Fig 1C-E and 1F representative
  traces, Fig 3A ROI overlays); those scripts are included for reference and the
  data is available on request.
- The reference pipelines in `processing/` use lab-server paths and external
  tools and are not required to regenerate the figures.
