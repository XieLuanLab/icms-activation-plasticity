"""Generate comprehensive stats table (CSV) for all Figure 6 comparisons.

Covers:
  - KL divergence (700ms and 120ms windows): PL e/l, NPL e/l, Early PL/NPL, Late PL/NPL
  - Population coupling: PL e/l, NPL e/l, Early PL/NPL, Late PL/NPL
  - Cluster permutation test results
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import json
import pickle
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, entropy
from scipy.stats.mstats import winsorize
from statsmodels.stats.multitest import multipletests

from utils.config import (DF_MERGED_PATH, RAW_DF_PATH, POP_COUPLING_DIR,
                           PSTH_NPZ_PATH, CLUSTERS_SEM_PATH, OUTPUT_DIR,
                           ALL_ANIMALS)
from utils.plotting import rank_biserial_r

FIG_DIR = OUTPUT_DIR / 'fig6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

np.random.seed(42)


def bootstrap_kl(dist1, dist2, bins, n_boot=1000, proportion=0.8):
    dist1 = winsorize(np.asarray(dist1), limits=[0.025, 0.025])
    dist2 = winsorize(np.asarray(dist2), limits=[0.025, 0.025])
    n = int(min(len(dist1), len(dist2)) * proportion)

    kl_vals = []
    for _ in range(n_boot):
        s1 = np.random.choice(dist1, size=n, replace=True)
        s2 = np.random.choice(dist2, size=n, replace=True)
        h1, _ = np.histogram(s1, bins=bins, density=True)
        h2, _ = np.histogram(s2, bins=bins, density=True)
        h1 += 1e-8; h2 += 1e-8
        h1 /= h1.sum(); h2 /= h2.sum()
        kl_vals.append(entropy(h1, h2))
    return np.array(kl_vals)


def collect_kl_stats(df_merged, delta_col, kl_bin_range, window_label):
    """Compute KL divergence stats for one window."""
    pl_df = df_merged[df_merged['is_pulse_locked'] == True]
    npl_df = df_merged[df_merged['is_pulse_locked'] == False]

    bins = np.histogram_bin_edges([], bins=25, range=kl_bin_range)

    pl_early_vals, pl_late_vals = [], []
    npl_early_vals, npl_late_vals = [], []

    for week in range(5):
        pl_hit = pl_df[(pl_df['rel_week'] == week) & (pl_df['trial_outcome'] == 'hit')][delta_col]
        pl_miss = pl_df[(pl_df['rel_week'] == week) & (pl_df['trial_outcome'] == 'miss')][delta_col]
        npl_hit = npl_df[(npl_df['rel_week'] == week) & (npl_df['trial_outcome'] == 'hit')][delta_col]
        npl_miss = npl_df[(npl_df['rel_week'] == week) & (npl_df['trial_outcome'] == 'miss')][delta_col]

        pl_kl = bootstrap_kl(pl_hit, pl_miss, bins)
        npl_kl = bootstrap_kl(npl_hit, npl_miss, bins)

        if week in [0, 1]:
            pl_early_vals.append(pl_kl)
            npl_early_vals.append(npl_kl)
        elif week in [2, 3, 4]:
            pl_late_vals.append(pl_kl)
            npl_late_vals.append(npl_kl)

    pl_e = np.concatenate(pl_early_vals)
    pl_l = np.concatenate(pl_late_vals)
    npl_e = np.concatenate(npl_early_vals)
    npl_l = np.concatenate(npl_late_vals)

    comparisons = [
        ('PL early vs late', pl_e, pl_l),
        ('NPL early vs late', npl_e, npl_l),
        ('Early PL vs NPL', pl_e, npl_e),
        ('Late PL vs NPL', pl_l, npl_l),
    ]

    # Run MWU and collect raw p-values
    raw_results = []
    for name, a, b in comparisons:
        stat, p = mannwhitneyu(a, b, alternative='two-sided')
        raw_results.append((name, a, b, stat, p))

    pvals = [r[4] for r in raw_results]
    _, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')

    rows = []
    for (name, a, b, stat, _), p_corr in zip(raw_results, pvals_corr):
        r = rank_biserial_r(stat, len(a), len(b))
        rows.append({
            'panel': f'KL divergence ({window_label})',
            'comparison': name,
            'n1': len(a),
            'n2': len(b),
            'median1': f'{np.median(a):.4f}',
            'median2': f'{np.median(b):.4f}',
            'test': 'Mann-Whitney U (two-sided, FDR-BH)',
            'U': f'{stat:.1f}',
            'p': f'{p_corr:.2e}',
            'r': f'{r:.3f}',
        })

    return rows


def collect_pop_coupling_stats():
    """Compute pop coupling stats."""
    raw_df = pd.read_pickle(RAW_DF_PATH)

    for ch in raw_df['stim_channel'].unique():
        ch_df = raw_df[raw_df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            raw_df = raw_df[raw_df['stim_channel'] != ch]

    df_mod = raw_df[
        (raw_df['rel_week'] < 5) & (raw_df['baseline_too_slow'] == False) &
        (raw_df['modulated'] == True) & (raw_df['mod_p_val'] < 0.05) &
        (raw_df['num_trials'] > 5) & (raw_df['z_score'] > 0) &
        (raw_df['z_score'] < 100) & (raw_df['num_spikes'] > 50) &
        (raw_df['stim_current'] < 7) & (raw_df['stim_current'] > 3)
    ].rename(columns={'stim_channel': 'stim_ch'})

    def is_early(animal_id, session):
        sel = df_mod[(df_mod['session'] == session) & (df_mod['animal_id'] == animal_id)]
        vals = sel['rel_week'].dropna().unique()
        if len(vals) == 0:
            return 'nan'
        return 'early' if float(vals[0]) < 2 else 'late'

    pl_early = {4: [], 5: [], 6: []}
    pl_late = {4: [], 5: [], 6: []}
    npl_early = {4: [], 5: [], 6: []}
    npl_late = {4: [], 5: [], 6: []}

    for animal_id in ALL_ANIMALS:
        pkl_path = POP_COUPLING_DIR / f'{animal_id}_pop_coupling.pkl'
        if not pkl_path.exists():
            continue
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        for session in data:
            tp = is_early(animal_id, session)
            if tp == 'nan':
                continue
            for (cur, ch), cond in data[session].items():
                if cur not in [4, 5, 6]:
                    continue
                pls = list(cond.get('pl_pc_norm_dict', {}).values())
                npls = list(cond.get('npl_pc_norm_dict', {}).values())
                if tp == 'early':
                    pl_early[cur].extend(pls)
                    npl_early[cur].extend(npls)
                else:
                    pl_late[cur].extend(pls)
                    npl_late[cur].extend(npls)

    def finite(arr):
        a = np.asarray(arr, float)
        return a[np.isfinite(a)]

    def concat_all(d):
        return np.concatenate([finite(d.get(c, [])) for c in [4, 5, 6]])

    groups = {
        'PL early': concat_all(pl_early),
        'PL late': concat_all(pl_late),
        'NPL early': concat_all(npl_early),
        'NPL late': concat_all(npl_late),
    }

    comparisons = [
        ('PL early vs late', groups['PL early'], groups['PL late']),
        ('NPL early vs late', groups['NPL early'], groups['NPL late']),
        ('Early PL vs NPL', groups['PL early'], groups['NPL early']),
        ('Late PL vs NPL', groups['PL late'], groups['NPL late']),
    ]

    raw_results = []
    for name, a, b in comparisons:
        stat, p = mannwhitneyu(a, b, alternative='two-sided')
        raw_results.append((name, a, b, stat, p))

    pvals = [r[4] for r in raw_results]
    _, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')

    rows = []
    for (name, a, b, stat, _), p_corr in zip(raw_results, pvals_corr):
        r = rank_biserial_r(stat, len(a), len(b))
        rows.append({
            'panel': 'Pop coupling',
            'comparison': name,
            'n1': len(a),
            'n2': len(b),
            'median1': f'{np.median(a):.4f}',
            'median2': f'{np.median(b):.4f}',
            'test': 'Mann-Whitney U (two-sided, FDR-BH)',
            'U': f'{stat:.1f}',
            'p': f'{p_corr:.2e}',
            'r': f'{r:.3f}',
        })

    return rows


def collect_cluster_stats():
    """Collect cluster permutation test results from cached JSON."""
    rows = []

    if not CLUSTERS_SEM_PATH.exists():
        print(f'  [SKIP] Cluster cache not found: {CLUSTERS_SEM_PATH}')
        return rows

    with open(CLUSTERS_SEM_PATH) as f:
        cached = json.load(f)

    for key, clusters in cached.items():
        if not clusters:
            rows.append({
                'panel': 'Cluster permutation test',
                'comparison': key,
                'n1': '', 'n2': '',
                'median1': '', 'median2': '',
                'test': 'cluster permutation (N=1000)',
                'U': '',
                'p': 'no clusters',
                'r': '',
            })
        else:
            for i, (start, end, p) in enumerate(clusters):
                rows.append({
                    'panel': 'Cluster permutation test',
                    'comparison': f'{key} cluster {i+1} (bins {start}-{end})',
                    'n1': '', 'n2': '',
                    'median1': '', 'median2': '',
                    'test': 'cluster permutation (N=1000)',
                    'U': '',
                    'p': f'{p:.4f}',
                    'r': '',
                })

    return rows


def pretty_print(df):
    col_widths = {col: max(len(col), df[col].astype(str).str.len().max())
                  for col in df.columns}
    header = ' | '.join(col.ljust(col_widths[col]) for col in df.columns)
    sep = '-+-'.join('-' * col_widths[col] for col in df.columns)
    print(header)
    print(sep)
    for _, row in df.iterrows():
        line = ' | '.join(str(row[col]).ljust(col_widths[col]) for col in df.columns)
        print(line)


def main():
    print('Collecting Figure 6 statistics...\n')
    rows = []

    # KL divergence
    print('KL divergence (700ms)...')
    df_merged = pd.read_pickle(DF_MERGED_PATH)
    rows.extend(collect_kl_stats(df_merged, 'delta_spks', (-30, 50), '700ms'))

    print('KL divergence (120ms)...')
    rows.extend(collect_kl_stats(df_merged, 'delta_spks_120ms', (-10, 15), '120ms'))

    # Pop coupling
    print('Pop coupling...')
    rows.extend(collect_pop_coupling_stats())

    # Cluster test
    print('Cluster test...')
    rows.extend(collect_cluster_stats())

    df = pd.DataFrame(rows)
    save_path = FIG_DIR / 'figure6_stats.csv'
    df.to_csv(save_path, index=False)
    print(f'\nSaved {save_path}\n')
    pretty_print(df)


if __name__ == '__main__':
    main()
