"""Fig 4F: Longitudinal rasters of all modulated units across sessions.

All modulated units sorted by cortical depth for ICMS92, stim condition
(ch 11, 4 uA), at weeks 0, 1, 2, 3, 4. Spikes shown from -700 to 3400 ms.

Usage:
    python python/fig4/longitudinal_rasters.py
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
from utils.plotting import apply_global_style

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

INTERNEURON_SPACE = 10
MAX_YLIM = 790
XLIM = [-700, 3400]


def main():
    d = np.load(DATA_DIR / 'fig4' / 'longitudinal_rasters_wide.npz', allow_pickle=True)

    session_indices = d['session_indices']
    stim_ch = int(d['stim_channel'])
    stim_cu = int(d['stim_current'])

    # Collect weeks in order
    weeks = []
    for idx in session_indices:
        for key in d.files:
            if key.startswith('wk') and key.endswith('_week'):
                wk = int(d[key])
                if wk not in weeks:
                    weeks.append(wk)
    weeks = sorted(set(weeks))

    n_panels = len(weeks)
    cmap = plt.get_cmap("viridis")
    palette = [cmap(i / 19) for i in range(19)]

    fig, axes = plt.subplots(1, n_panels, figsize=(n_panels * 1.0, 3), sharey=False)
    if n_panels == 1:
        axes = [axes]

    for plot_idx, wk in enumerate(weeks):
        ax = axes[plot_idx]
        n_units = int(d[f'wk{wk}_n_units'])

        offset = 0
        unit_labels = []
        y_positions = []

        for j in range(n_units):
            raster = d[f'wk{wk}_unit{j}_raster']
            num_trials = len(raster)
            trial_offsets = np.arange(offset, offset + num_trials)

            ax.eventplot(
                raster,
                orientation="horizontal",
                colors=[palette[j % len(palette)]],
                linelengths=1.2,
                linewidths=0.5,
                lineoffsets=trial_offsets
            )

            y_positions.append(trial_offsets.mean())
            unit_labels.append(f"{j + 1}")
            offset += num_trials + INTERNEURON_SPACE

            if j < n_units - 1:
                ax.axhline(offset - INTERNEURON_SPACE / 2,
                           color='gray', linewidth=0.5)

        # Stim window shading
        ax.fill_between(
            x=[0, 700], y1=0, y2=offset,
            color=sns.color_palette("deep")[3], alpha=0.15, linewidth=0)

        # Stim duration bar at top
        ax.plot([0, 700], [MAX_YLIM + 10, MAX_YLIM + 10],
                color='black', linewidth=1)

        # Axes
        ax.set_xlim(XLIM)
        ax.set_ylim([-5, MAX_YLIM + 15])
        ax.set_xticks([0, 2000])
        ax.set_xticklabels(['0', '2000'], fontsize=7)
        ax.set_title(f'Week {wk}', fontsize=8)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(unit_labels, rotation=0, fontsize=5)

        if plot_idx == 0:
            ax.set_ylabel('Neuron ID', fontsize=7)
            ax.set_xlabel('Time (ms)', fontsize=7)
        else:
            ax.tick_params(labelleft=True)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        print(f'Week {wk}: {n_units} units')

    plt.tight_layout()
    out = FIG_DIR / 'longitudinal_rasters.svg'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
