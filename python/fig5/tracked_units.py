"""Fig 5I-M: Tracked PL and NPL unit templates, P(spike), and firing rate.

Horizontal templates (4 channels) across 3 sessions (Weeks 1, 2, 4).
P(spike) and firing rate in 2x3 grid. Matches original tracked_units.py.

Usage:
    python python/fig5/tracked_units.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig5'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def plot_tracked_unit(label):
    d = np.load(DATA_DIR / 'fig5' / f'tracked_unit_{label}.npz', allow_pickle=True)

    n_sessions = int(d['n_sessions'])
    weeks = d['weeks']
    sessions = d['sessions']
    color_idx = int(d['color_idx'])
    color = PALETTE[color_idx]
    pulse_win = float(d['pulse_win_ms'])
    post_blank = float(d['post_stim_blank_ms'])
    pre_blank = float(d['pre_stim_blank_ms'])

    # ─── Horizontal templates (1×3) ─────────────────────────────
    WAVEFORM_SPACING = 120
    fig_t, axes_t = plt.subplots(1, n_sessions, figsize=(4, 1))

    for col in range(n_sessions):
        ax = axes_t[col]
        tmpl = d[f's{col}_template_mean']  # (n_samples, 4)
        n_samples, n_ch = tmpl.shape

        offsets = [(i - 1) * WAVEFORM_SPACING for i in range(n_ch)]

        for i in range(n_ch):
            x = np.arange(n_samples) - 100  # x_offset from original
            alpha = max(0.3, 1.0 - 0.25 * i)
            ax.plot(x + offsets[i], tmpl[:, i], c=color,
                    alpha=alpha, linewidth=0.5)

        ax.set_ylim([-150, offsets[-1] + 50])
        ax.set_title(f'Week {weeks[col]}', fontsize=4, pad=3)
        ax.axis('off')

        # Scale bar on first panel
        if col == 0:
            y_neg = -150
            h_pos_x, h_pos_y = -240, y_neg + 50
            ax.plot([h_pos_x, h_pos_x + 30], [h_pos_y, h_pos_y], 'k', linewidth=0.5)
            ax.plot([h_pos_x, h_pos_x], [h_pos_y, h_pos_y + 100], 'k', linewidth=0.5)
            ax.text(h_pos_x, h_pos_y - 15, '1 ms', fontsize=3)
            ax.text(h_pos_x - 5, h_pos_y, r'100 $\mu$V', fontsize=3, rotation=90)

    fig_t.tight_layout()
    fig_t.savefig(FIG_DIR / f'{label}_template_tracked.svg', format='svg',
                   bbox_inches='tight')
    plt.close(fig_t)
    print(f'Saved {label}_template_tracked.svg')

    # ─── P(spike) and FR (2×3 grid) ─────────────────────────────
    fig, axes = plt.subplots(2, n_sessions, figsize=(4, 2))

    for col in range(n_sessions):
        # P(spike)
        ax_p = axes[0, col]
        key = f's{col}_prob_mean_4uA'
        if key in d:
            bins = d[f's{col}_prob_bins_4uA']
            mean = d[f's{col}_prob_mean_4uA']
            std = d[f's{col}_prob_std_4uA']

            ax_p.plot(bins, mean, color=color, linewidth=1)
            ax_p.fill_between(bins, mean - std, mean + std,
                              color=color, alpha=0.3, linewidth=0, edgecolor='none')

            ax_p.axvspan(0, post_blank, color='lightgray', alpha=0.5,
                         linewidth=0, edgecolor='none', zorder=0)
            ax_p.axvspan(pulse_win - pre_blank, pulse_win,
                         color='lightgray', alpha=0.5,
                         linewidth=0, edgecolor='none', zorder=0)
            ax_p.set_xlim([0, pulse_win])

        ax_p.set_ylim([0, 0.2])
        ax_p.spines['top'].set_visible(False)
        ax_p.spines['right'].set_visible(False)

        # Firing rate
        ax_f = axes[1, col]
        xkey = f's{col}_fr_x_4uA'
        ykey = f's{col}_fr_y_4uA'
        if xkey in d:
            ax_f.plot(d[xkey], d[ykey], color=color, linewidth=1)
            ax_f.axvline(0, color='gray', linestyle='--', linewidth=0.5)
            ax_f.axvline(700, color='gray', linestyle='--', linewidth=0.5)
            ax_f.set_xlim([-700, 1400])
            ax_f.set_xticks([0, 700])

        ax_f.set_ylim([0, 100])
        ax_f.spines['top'].set_visible(False)
        ax_f.spines['right'].set_visible(False)

    axes[0, 0].set_ylabel('P(spike)')
    axes[0, 0].set_xlabel('Time (ms)')
    axes[1, 0].set_ylabel('Firing Rate (Hz)')
    axes[1, 0].set_xlabel('Time (ms)')

    fig.tight_layout()
    fig.savefig(FIG_DIR / f'{label}_metrics_tracked.svg', format='svg',
                 bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {label}_metrics_tracked.svg')


def main():
    print('=== PL tracked unit ===')
    plot_tracked_unit('pl')
    print('\n=== NPL tracked unit ===')
    plot_tracked_unit('npl')


if __name__ == '__main__':
    main()
