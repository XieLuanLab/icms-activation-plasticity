"""Figure S1 — Behavioral characterization.

Panels:
  A: Response time histogram
  B: CDF of response times
  C: Psychometric curves (example session)
  D: Detection thresholds over weeks (individual animals)
  E: Response times over weeks (individual animals)

Usage:
    python python/fig_s1/generate_figure.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import glob
import dill
import math
from datetime import datetime

from utils.config import DATA_DIR, OUTPUT_DIR, MOUSE_NAMES
from utils.plotting import apply_global_style, PALETTE

apply_global_style()

S1_DATA = DATA_DIR / 'fig_s1'
FIG_DIR = OUTPUT_DIR / 'fig_s1'
FIG_DIR.mkdir(parents=True, exist_ok=True)

RT_PKL = S1_DATA / 'rt.pkl'
ANIMAL_MAP = MOUSE_NAMES


def _get_rel_days(sessions):
    fmt = '%d-%b-%Y'
    dates = []
    for s in sessions:
        s = str(s)
        parts = s.replace('\\', '/').split('/')
        for p in reversed(parts):
            try:
                dates.append(datetime.strptime(p, fmt))
                break
            except ValueError:
                continue
    if not dates:
        return [0] * len(sessions)
    first = min(dates)
    return [(d - first).days for d in dates]


def plot_rt_histogram_and_cdf():
    print('Response time histogram and CDF...')
    with open(RT_PKL, 'rb') as f:
        rt_records = dill.load(f)

    df = pd.DataFrame(rt_records)
    df = df[df['week'] < 5]
    rts = df['rt_ms'].values
    sorted_rts = np.sort(rts)
    median_rt = int(np.median(sorted_rts))

    fig, ax = plt.subplots(figsize=(2.5, 2))
    ax.hist(sorted_rts, bins=70, color=PALETTE[0])
    ax.set_xlabel('Response time (ms)')
    ax.set_ylabel('Count')
    ax.set_xlim([0, 700])
    ax.axvline(median_rt, color=PALETTE[3], linestyle='--')
    ax.text(median_rt + 10, ax.get_ylim()[1] * 0.8,
            f'Median: {median_rt} ms', color=PALETTE[3], fontsize=7)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'rt_histogram.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'  Saved rt_histogram.svg (median={median_rt} ms, n={len(rts)})')

    cdf = np.arange(1, len(sorted_rts) + 1) / len(sorted_rts)
    fig, ax = plt.subplots(figsize=(2.5, 2))
    ax.plot(sorted_rts, cdf, color=PALETTE[0], linewidth=1)
    ax.set_xlabel('Response time (ms)')
    ax.set_ylabel('Cumulative probability')
    ax.set_xlim([0, 700])
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'rt_cdf.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('  Saved rt_cdf.svg')


def plot_detection_thresholds():
    print('\nDetection thresholds over weeks...')
    all_rows = []

    csv_files = sorted(glob.glob(str(S1_DATA / '*_blocks.csv')))
    for fp in csv_files:
        df = pd.read_csv(fp)
        animal_id = Path(fp).stem.split('_')[0]
        sessions = df['Session'].unique()
        rel_days = _get_rel_days(sessions)

        for idx, session in enumerate(sessions):
            rel_day = rel_days[idx]
            rel_week = math.ceil(rel_day / 7) if rel_day > 0 else 0
            for ch_label in ['A', 'B', 'C']:
                thr_col = f'Threshold{ch_label}'
                if thr_col in df.columns:
                    thr = df[df['Session'] == session][thr_col].mean()
                    all_rows.append({'animal_id': animal_id, 'week': rel_week,
                                     'channel': ch_label, 'threshold': thr})

    icms83_thr = S1_DATA / 'ICMS83_thresholds.csv'
    if icms83_thr.exists():
        df83 = pd.read_csv(icms83_thr)
        df83 = df83[df83['session'] != '11-Aug-2023']
        ch_map = {3: 'A', 4: 'B', 5: 'C'}
        for _, row in df83.iterrows():
            ch_label = ch_map.get(row['channel'], None)
            if ch_label is None:
                continue
            all_rows.append({'animal_id': 'ICMS83', 'week': int(row['rel_week']),
                             'channel': ch_label, 'threshold': row['threshold']})

    if not all_rows:
        print('  No threshold data found')
        return

    thr_df = pd.DataFrame(all_rows)
    thr_df = thr_df[thr_df['week'] < 5]
    thr_df['mouse'] = thr_df['animal_id'].map(ANIMAL_MAP)

    fig, ax = plt.subplots(figsize=(3, 2))
    ordered_mice = list(ANIMAL_MAP.values())

    for i, mouse in enumerate(ordered_mice):
        mdf = thr_df[thr_df['mouse'] == mouse]
        if len(mdf) == 0:
            continue
        jitter = np.random.uniform(-0.15, 0.15, size=len(mdf))
        ax.scatter(mdf['week'] + jitter, mdf['threshold'],
                   s=8, alpha=0.7, color=PALETTE[i], edgecolors='none')

    for i, mouse in enumerate(ordered_mice):
        mdf = thr_df[thr_df['mouse'] == mouse]
        if len(mdf) == 0:
            continue
        med = mdf.groupby('week')['threshold'].median().reset_index().sort_values('week')
        ax.plot(med['week'], med['threshold'], '-', color=PALETTE[i],
                linewidth=1, label=mouse)

    ax.set_xlabel('Weeks of training')
    ax.set_ylabel('Detection threshold (\u00b5A)')
    ax.set_xticks(sorted(thr_df['week'].unique()))
    ax.set_yticks([2, 4, 6, 8, 10])
    ax.legend(fontsize=5, loc='upper right')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'detection_thresholds.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('  Saved detection_thresholds.svg')


def plot_rt_over_weeks():
    print('\nResponse times over weeks...')
    with open(RT_PKL, 'rb') as f:
        rt_records = dill.load(f)

    df = pd.DataFrame(rt_records)
    df = df[df['week'] < 5]
    df['mouse'] = df['animal_id'].map(ANIMAL_MAP)
    ordered_mice = list(ANIMAL_MAP.values())

    summary = df.groupby(['mouse', 'week'])['rt_ms'].agg(
        median='median',
        q1=lambda x: np.percentile(x, 25),
        q3=lambda x: np.percentile(x, 75)
    ).reset_index()

    fig, ax = plt.subplots(figsize=(3, 2))
    offsets = np.linspace(-0.15, 0.15, len(ordered_mice))

    for i, mouse in enumerate(ordered_mice):
        mdf = summary[summary['mouse'] == mouse]
        ax.errorbar(mdf['week'] + offsets[i], mdf['median'],
                    yerr=[mdf['median'] - mdf['q1'], mdf['q3'] - mdf['median']],
                    fmt='-o', capsize=1, color=PALETTE[i], markersize=4,
                    linewidth=1, label=mouse)

    ax.set_xlabel('Weeks of training')
    ax.set_ylabel('Response time (ms)')
    ax.set_xticks(sorted(df['week'].unique()))
    ax.legend(fontsize=5, loc='upper right')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'rt_over_weeks.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('  Saved rt_over_weeks.svg')


def plot_psych_curves():
    """3 subplots (one per stim channel), each showing 4 blocks overlaid."""
    print('\nPsychometric curves...')
    psych_path = S1_DATA / 'psych_curves.npz'
    if not psych_path.exists():
        print('  psych_curves.npz not found - skipping')
        return

    d = np.load(psych_path)
    channels = list(d['channels'])
    n_blocks = len(d['blocks']) if 'blocks' in d else 4
    ch_labels = ['Stim channel A', 'Stim channel B', 'Stim channel C']
    block_colors = PALETTE[:n_blocks]

    fig, axes = plt.subplots(1, len(channels), figsize=(7, 2), sharey=True)

    for ch_idx, ch_key in enumerate(channels):
        ax = axes[ch_idx]
        ax.set_title(ch_labels[ch_idx], fontsize=7)

        for block in range(n_blocks):
            prefix = f'{ch_key}_block{block}'
            if f'{prefix}_x' not in d:
                continue

            x_data = d[f'{prefix}_x']
            y_data = d[f'{prefix}_y']
            thr = float(d[f'{prefix}_threshold'][0])
            x_fine = d[f'{prefix}_x_fine']
            y_fine = d[f'{prefix}_y_fine']

            y_frac = y_data[:, 0] / y_data[:, 1]
            ax.scatter(x_data, y_frac, color=block_colors[block], s=15,
                       edgecolors='none', alpha=0.7)
            ax.plot(x_fine, y_fine, color=block_colors[block], lw=0.8,
                    label=f'Block {block + 1}: {thr:.1f} \u00b5A')

        ax.set_xlim([0, 10])
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlabel('Stimulation amplitude (\u00b5A)')
        if ch_idx == 0:
            ax.set_ylabel('Response probability')
        ax.legend(fontsize=4, loc='lower right')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'psych_curves.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('  Saved psych_curves.svg')


def main():
    plot_rt_histogram_and_cdf()
    plot_psych_curves()
    plot_detection_thresholds()
    plot_rt_over_weeks()
    print(f'\nAll saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
