"""Fig 6A: Example PL and NPL unit rasters + firing rates (hit vs miss).

4 panels: PL early, NPL early, PL late, NPL late.
Each shows raster (hit=green, miss=orange) and smoothed FR.

Usage:
    python python/fig6/example_raster_and_fr.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

color_map = sns.color_palette("deep", 4)
HIT_C = color_map[2]   # green
MISS_C = color_map[3]   # orange/red
PL_C = color_map[0]     # blue (for title)
NPL_C = color_map[1]    # orange (for title)

PANELS = ['pl_early', 'npl_early', 'pl_late', 'npl_late']


def plot_panel(label):
    d = np.load(DATA_DIR / 'fig6' / f'fig6a_{label}.npz', allow_pickle=True)

    raster = d['raster_array']
    trial_colors = d['trial_colors']
    hit_fr = d['hit_fr']
    miss_fr = d['miss_fr']
    fr_bins = d['fr_bins']
    fr_ylim = float(d['fr_ylim'])
    title = str(d['title'])

    title_color = PL_C if 'pl' in label else NPL_C

    fig, axes = plt.subplots(1, 2, figsize=(3, 2))

    # Raster
    for trial_idx in range(len(raster)):
        spikes = raster[trial_idx]
        color = MISS_C if trial_colors[trial_idx] == 'miss' else HIT_C
        axes[0].vlines(spikes, ymin=trial_idx - 0.5, ymax=trial_idx + 0.5,
                       color=color, linewidth=0.32)

    axes[0].plot([0, 700], [41, 41], color='black', linewidth=1)
    axes[0].set_xlim([-700, 1400])
    axes[0].set_ylim([-1, 43])
    axes[0].set_ylabel('Trial')
    axes[0].set_xlabel('Time (ms)')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Firing rate
    axes[1].plot(fr_bins, hit_fr, c=HIT_C, linewidth=0.32)
    axes[1].plot(fr_bins, miss_fr, c=MISS_C, linewidth=0.32)
    axes[1].plot([0, 700], [fr_ylim, fr_ylim], color='black', linewidth=1)
    axes[1].set_xlim([-700, 1400])
    axes[1].set_ylim([-1, fr_ylim + 2])
    axes[1].set_ylabel('Firing Rate (Hz)')
    axes[1].set_xlabel('Time (ms)')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    # Legend on FR panel (only first PL panel)
    if label == 'pl_early':
        custom_lines = [
            Line2D([0], [0], color=HIT_C, linewidth=0.32),
            Line2D([0], [0], color=MISS_C, linewidth=0.32)
        ]
        axes[1].legend(custom_lines, ['Hit', 'Miss'],
                       fontsize=6, bbox_to_anchor=(0.55, 0.8))

    plt.suptitle(title, color=title_color)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f'{label}_unit.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {label}_unit.svg')


def main():
    for label in PANELS:
        plot_panel(label)


if __name__ == '__main__':
    main()
