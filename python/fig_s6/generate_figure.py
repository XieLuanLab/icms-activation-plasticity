"""Figure S6 -- Control animal ephys (no behavioral training).

Panels A-E: PL mod | NPL mod | PL t2max | NPL t2max | PL P(spike)
All comparisons early vs late, all NS.

Usage (from repo root):
    python -m python.fig_s6.generate_figure
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from datetime import datetime
from scipy.stats import mannwhitneyu

from utils.config import (CONTROL_CSV, OUTPUT_DIR, EARLY_WEEKS, LATE_WEEKS)
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r

apply_global_style()

CURRENTS = [5]
MIN_SPIKES = 15
MIN_TRIALS = 5

FIG_DIR = OUTPUT_DIR / 'fig_s6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

SCATTER_SIZE = 5


def load_control():
    df = pd.read_csv(CONTROL_CSV)
    if 'animal' in df.columns:
        df = df.rename(columns={'animal': 'animal_id'})

    for aid in df['animal_id'].unique():
        mask = df['animal_id'] == aid
        sessions = df.loc[mask, 'session_date'].unique()
        dates = sorted([datetime.strptime(s, '%Y-%m-%d') for s in sessions])
        first_date = dates[0]
        for s in sessions:
            d = datetime.strptime(s, '%Y-%m-%d')
            delta = (d - first_date).days
            df.loc[(df['animal_id'] == aid) & (df['session_date'] == s), 'rel_day'] = delta
            df.loc[(df['animal_id'] == aid) & (df['session_date'] == s), 'rel_week'] = (
                0 if delta == 0 else int(np.ceil(delta / 7)))
    return df


def filter_control(df):
    sub = df[df['stim_current'].isin(CURRENTS)]
    return sub[(sub['modulated'] == True) & (sub['mod_p_val'] < 0.05) &
               (sub['z_score'] > 0) & (sub['rel_week'] < 5) &
               (sub['baseline_too_slow'] == False) & (sub['num_trials'] > MIN_TRIALS) &
               (sub['num_spikes'] > MIN_SPIKES)]


def plot_early_late_panel(ax, df, metric, ylabel, title, color, alternative='two-sided'):
    early = df[df['rel_week'].isin(EARLY_WEEKS)][metric].dropna()
    late = df[df['rel_week'].isin(LATE_WEEKS)][metric].dropna()

    for vals, x in [(early, 0), (late, 0.5)]:
        if len(vals) > 0:
            jitter = np.random.uniform(-0.05, 0.05, size=len(vals))
            ax.scatter(x + jitter, vals, s=SCATTER_SIZE, alpha=0.25, color=color,
                       edgecolors='none', rasterized=False, zorder=1)
            med = vals.median()
            q25, q75 = vals.quantile(0.25), vals.quantile(0.75)
            ax.errorbar(x, med, yerr=[[med - q25], [q75 - med]],
                        fmt='o', color=color, markersize=4, capsize=1.5,
                        elinewidth=0.5, capthick=0.5, zorder=2)

    stat_row = None
    if len(early) > 2 and len(late) > 2:
        stat, p = mannwhitneyu(early, late, alternative=alternative)
        r = rank_biserial_r(stat, len(early), len(late))
        star = sig_text(p)

        y_min, y_max = ax.get_ylim()
        yr = y_max - y_min
        by = y_max + yr * 0.02
        bh = yr * 0.03
        ax.plot([0, 0, 0.5, 0.5], [by, by + bh, by + bh, by], color='k', lw=0.5)
        fs = 6 if star == 'NS' else 8
        ax.text(0.25, by + bh + yr * 0.06, star, ha='center', fontsize=fs, color=color)
        ax.set_ylim(y_min, by + bh + yr * 0.12)

        print(f'  {title} {metric}: e={len(early)} med={early.median():.3f}, '
              f'l={len(late)} med={late.median():.3f}, p={p:.2e} {star}')

        stat_row = {
            'panel': f'{title} {ylabel}', 'comparison': 'early vs late',
            'n1': len(early), 'n2': len(late),
            'median1': f'{early.median():.3f}', 'median2': f'{late.median():.3f}',
            'test': f'Mann-Whitney U ({alternative})', 'U': f'{stat:.0f}',
            'p': f'{p:.2e}', 'r': f'{r:.3f}',
        }

    ax.set_xticks([0, 0.5])
    ax.set_xticklabels(['Early', 'Late'])
    ax.set_xlim(-0.3, 0.8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=8, pad=10)
    return stat_row


def main():
    ctrl = load_control()
    mod = filter_control(ctrl)
    pl = mod[mod['is_pulse_locked'] == True]
    npl = mod[mod['is_pulse_locked'] == False]

    print(f'Control PL: {len(pl)}, NPL: {len(npl)}')

    color = sns.color_palette('deep')[1]
    all_stats = []

    fig, axes = plt.subplots(1, 5, figsize=(12, 2.5))

    for ax, df, metric, ylabel, title in [
        (axes[0], pl, 'z_score', 'Modulation', 'PL'),
        (axes[1], npl, 'z_score', '', 'NPL'),
        (axes[2], pl, 't_to_max_10ms_smoothed', 'Time to max\nFR (ms)', 'PL'),
        (axes[3], npl, 't_to_max_10ms_smoothed', '', 'NPL'),
        (axes[4], pl, 'max_spike_prob', 'P(spike)', 'PL'),
    ]:
        sr = plot_early_late_panel(ax, df, metric, ylabel, title, color)
        if sr: all_stats.append(sr)

    plt.suptitle(f'Control animals - {CURRENTS} \u00b5A', fontsize=10, y=1.05)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / 'control_combined.svg', format='svg', bbox_inches='tight')
    fig.savefig(FIG_DIR / 'control_combined.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\nSaved control_combined.svg')

    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        stats_path = FIG_DIR / 'figure_sControl_stats.csv'
        stats_df.to_csv(stats_path, index=False)
        print(f'Saved {stats_path}')

    print(f'All saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
