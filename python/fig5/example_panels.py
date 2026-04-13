"""Fig 5A-C: Example PL and NPL raw traces, pulse rasters, and spike probability.

Usage:
    python python/fig5/example_panels.py
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

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig5'
FIG_DIR.mkdir(parents=True, exist_ok=True)

STIM_COLORS = sns.color_palette("deep", 3)


def plot_raw_trace(label, title):
    """Plot raw trace with spike markers and stim shading."""
    d = np.load(DATA_DIR / 'fig5' / f'raw_trace_{label}.npz')
    trace = d['raw_trace']
    spikes = d['spike_frames']
    stim_ts = d['stim_ts']
    fs = int(d['fs'])

    time_ms = np.arange(len(trace)) / (fs / 1000)
    stim_ms = stim_ts / (fs / 1000)

    fig, ax = plt.subplots(figsize=(2.2, 1.6))
    ax.plot(time_ms, trace, 'k', linewidth=0.5)

    # Stim markers + blanking
    for i, s in enumerate(stim_ms):
        ax.axvspan(s - 0.5, s + 1.5, color='red', alpha=0.08, zorder=0,
                   label='Blanked' if i == 0 else '')
        ax.axvline(s, ymin=0, ymax=0.85, color=PALETTE[3], linewidth=0.5,
                   label='Stim onset' if i == 0 else '')

    # Spike markers (offset above trace, matching original marker_y_offset=100)
    spike_ms = spikes / (fs / 1000)
    x_start = stim_ms[0] - 0.5
    x_end = stim_ms[-1] + 1.5
    visible = (spike_ms >= x_start) & (spike_ms <= x_end)
    if visible.any():
        spike_y = trace[spikes[visible].astype(int)]
        ax.scatter(spike_ms[visible], spike_y - 100, color=PALETTE[0],
                   marker='^', s=12, zorder=3, label='Detected spike')

    ax.set_xlim(x_start, x_end)
    # Match original ylims from fig_utils.py
    if label == 'pl':
        ax.set_ylim([-4500, -1800])
    elif label == 'npl':
        ax.set_ylim([-1800, -300])
    ax.axis('off')

    # Scale bars
    y_lo = ax.get_ylim()[0] - 50
    ax.plot([x_start, x_start + 5], [y_lo, y_lo], 'k', linewidth=1, clip_on=False)
    ax.text(x_start + 2.5, y_lo - 60, '5 ms', ha='center', fontsize=5)
    ax.plot([x_start, x_start], [y_lo, y_lo + 200], 'k', linewidth=1, clip_on=False)
    ax.text(x_start - 0.5, y_lo + 100, r'200 $\mu$V', ha='right', fontsize=5,
            rotation=90, va='center')

    fig.savefig(FIG_DIR / f'raw_trace_{label}.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved raw_trace_{label}.svg')


def plot_pulse_raster_and_prob(label, title):
    """Plot pulse-aligned raster and spike probability.

    Matches original response_plotting_util.py: plot_pulse_raster + plot_probability.
    """
    d = np.load(DATA_DIR / 'fig5' / f'pulse_raster_{label}.npz', allow_pickle=True)

    import matplotlib as mpl
    mpl.rcParams['lines.linewidth'] = 0.25
    mpl.rcParams['lines.solid_capstyle'] = 'butt'

    fig, (ax_raster, ax_prob) = plt.subplots(1, 2, figsize=(3.75, 1.5))

    # ─── Raster ──────────────────────────────────────────────────
    all_rasters = []
    all_colors = []
    all_offsets = []
    offset = 0

    for ci, current in enumerate([4, 5, 6]):
        rkey = f'raster_{current}uA'
        if rkey not in d:
            continue

        raster = d[rkey]
        n_pulses = int(d[f'n_pulses_{current}uA'])
        color = STIM_COLORS[ci]

        lineoffsets = np.arange(n_pulses) + offset
        all_rasters.extend(raster)
        all_colors.extend([color] * n_pulses)
        all_offsets.extend(lineoffsets)
        offset += n_pulses

    total_pulses = len(all_rasters)
    pulse_linelength = max(1, np.ceil(total_pulses / 40))

    ax_raster.eventplot(
        all_rasters, orientation='horizontal',
        linewidths=0.25, linelengths=pulse_linelength,
        lineoffsets=all_offsets, colors=all_colors)

    ax_raster.set_xlabel('Time (ms)')
    ax_raster.set_ylabel('Pulse')

    # ─── Spike probability with error shading ────────────────────
    pulse_win_ms = float(d['pulse_win_ms'])
    post_blank = float(d['post_stim_blank_ms'])
    pre_blank = float(d['pre_stim_blank_ms'])

    for ci, current in enumerate([4, 5, 6]):
        mean_key = f'prob_mean_{current}uA'
        if mean_key not in d:
            continue

        bin_centers = d[f'prob_bin_centers_{current}uA']
        mean_prob = d[f'prob_mean_{current}uA']
        std_prob = d[f'prob_std_{current}uA']
        n_train = int(d[f'n_train_trials_{current}uA'])
        color = STIM_COLORS[ci]

        alpha_line = 0.2 if n_train < 5 else 1
        alpha_fill = 0.1 if n_train < 5 else 0.3

        ax_prob.plot(bin_centers, mean_prob, color=color,
                     alpha=alpha_line, zorder=2, label=f'{current} \u00b5A')
        ax_prob.fill_between(bin_centers,
                             mean_prob - std_prob, mean_prob + std_prob,
                             color=color, alpha=alpha_fill,
                             linewidth=0, edgecolor='none', zorder=1)

    # Blanking shading
    ax_prob.axvspan(0, post_blank, color='lightgray', alpha=0.5,
                    linewidth=0, edgecolor='none', zorder=0)
    ax_prob.axvspan(pulse_win_ms - pre_blank, pulse_win_ms,
                    color='lightgray', alpha=0.5,
                    linewidth=0, edgecolor='none', zorder=0)

    ax_prob.set_xlabel('Time (ms)')
    ax_prob.set_ylabel('P(spike)')
    ax_prob.legend(fontsize=5)

    plt.tight_layout()
    fig.savefig(FIG_DIR / f'pulse_raster_{label}.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved pulse_raster_{label}.svg')


def main():
    print('=== Raw traces ===')
    plot_raw_trace('pl', 'Pulse-locked')
    plot_raw_trace('npl', 'Non-pulse-locked')

    print('\n=== Pulse rasters + spike probability ===')
    plot_pulse_raster_and_prob('pl', 'Pulse-locked')
    plot_pulse_raster_and_prob('npl', 'Non-pulse-locked')


if __name__ == '__main__':
    main()
