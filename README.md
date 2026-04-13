# ICMS Plasticity — Figure Reproduction Code

Code and processed data for reproducing figures in:

**"Learning induces activation-mechanism-dependent neural plasticity in an intracortical microstimulation task"**
Kim et al., *Science Advances* (2026)

## Setup

### Data
Download processed data from Zenodo and place contents in the `data/` directory.

### Python (Figures 4-6, S1, S3-S8)
```bash
pip install -r requirements.txt
```
Python 3.10+ recommended.

### MATLAB (Figures 2-3, S2)
Requires MATLAB R2023a or later.

## Generating figures

Python scripts are run from the repository root:
```bash
python -m python.fig5.generate_figure5
python -m python.fig_s7.generate_figure
```

MATLAB scripts should be run from the repository root after loading data via `matlab/config.m`. Output files are saved to `output/`.

## Notes

- Figure 1C-F data is not included but is available upon request.
- Some MATLAB scripts in `fig3/` reference raw session data on the lab server and are included as reference. The reproducible Fig 3 panels use `config.m` and the saved `.mat` files. Raw data is available upon request.
