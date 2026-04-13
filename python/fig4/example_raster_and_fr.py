"""Fig 4E: Example unit raster and firing rate at 4, 5, 6 uA.

Reproduces the raster + FR panel from a representative unit (ICMS98, unit 15).

Usage:
    python python/fig4/example_raster_and_fr.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE


FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()

_deep3 = sns.color_palette("deep", 3)
STIM_COLORS = {4: _deep3[0], 5: _deep3[1], 6: _deep3[2]}


def main():
    d = np.load(DATA_DIR / 'fig4' /
                'example_raster_and_fr.npz', allow_pickle=True)

    currents = [4, 5, 6]

    n_trials = min(int(d[f'n_trials_{c}uA']) for c in currents)

    fig, axes = plt.subplots(2, 1, figsize=(1.7, 2.6))
    ax_raster = axes[0]
    ax_fr = axes[1]

    # Raster (eventplot)
    all_rasters = []
    all_colors = []
    all_offsets = []
    offset = 0
    linewidth_factor = 40

    for current in currents:
        raster = d[f'raster_{current}uA'][:n_trials]
        line_count = len(raster)
        line_height = np.ceil(line_count / linewidth_factor)
        offsets = np.arange(line_count) + offset

        all_rasters.extend(raster)
        all_colors.extend([STIM_COLORS[current]] * line_count)
        all_offsets.extend(offsets)
        offset += line_count

    ax_raster.eventplot(
        all_rasters,
        orientation="horizontal",
        colors=all_colors,
        linelengths=np.ceil(n_trials / linewidth_factor),
        linewidths=0.4,
        lineoffsets=all_offsets,
    )

    # Stim duration bar
    ax_raster.plot([0, 700], [offset + 2, offset + 2],
                   color='black', linewidth=2)
    ax_raster.set_xticks([0, 700])
    ax_raster.set_xticklabels([])
    ax_raster.set_ylabel('Trial')
    ax_raster.set_xlim(-700, 1400)

    # Firing rate
    for current in currents:
        x = d[f'fr_x_{current}uA']
        y = d[f'fr_y_{current}uA']
        ax_fr.plot(x, y, color=STIM_COLORS[current],
                   linewidth=0.8, label=f'{current} \u00b5A')

    ax_fr.set_xlim(-700, 1400)
    ax_fr.plot([0, 700], [100 + 2, 100 + 2], color='black', linewidth=2)
    ax_fr.set_xticks([0, 700])
    ax_fr.set_xlabel('Time (ms)')
    ax_fr.set_ylabel('Firing rate (Hz)')
    ax_fr.legend(frameon=False, loc='center left', bbox_to_anchor=(0.6, 0.75),
                 handlelength=1.5)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'example_raster_and_fr.svg',
                format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "example_raster_and_fr.svg"}')


if __name__ == '__main__':
    main()
