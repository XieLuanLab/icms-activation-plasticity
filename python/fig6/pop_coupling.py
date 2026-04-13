"""Population coupling -- all animals including ICMS83."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from utils.config import RAW_DF_PATH, POP_COUPLING_DIR, OUTPUT_DIR, ALL_ANIMALS
from utils.plotting import apply_global_style, PALETTE, sig_text

FIG_DIR = OUTPUT_DIR / 'fig6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()


def finite(arr):
    a = np.asarray(arr, float)
    return a[np.isfinite(a)]


def main():
    raw_df = pd.read_pickle(RAW_DF_PATH)

    # Filter modulated
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
        if len(vals) == 0: return 'nan'
        return 'early' if float(vals[0]) < 2 else 'late'

    pl_early = {4: [], 5: [], 6: []}
    pl_late = {4: [], 5: [], 6: []}
    npl_early = {4: [], 5: [], 6: []}
    npl_late = {4: [], 5: [], 6: []}

    for animal_id in ALL_ANIMALS:
        pkl_path = POP_COUPLING_DIR / f'{animal_id}_pop_coupling.pkl'
        if not pkl_path.exists():
            print(f'  [SKIP] {pkl_path}')
            continue
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)
        for session in data:
            tp = is_early(animal_id, session)
            if tp == 'nan': continue
            for (cur, ch), cond in data[session].items():
                if cur not in [4, 5, 6]: continue
                pls = list(cond.get('pl_pc_norm_dict', {}).values())
                npls = list(cond.get('npl_pc_norm_dict', {}).values())
                if tp == 'early':
                    pl_early[cur].extend(pls)
                    npl_early[cur].extend(npls)
                else:
                    pl_late[cur].extend(pls)
                    npl_late[cur].extend(npls)

    # Collapsed
    def concat_all(d):
        return np.concatenate([finite(d.get(c, [])) for c in [4, 5, 6]])

    groups = [concat_all(pl_early), concat_all(pl_late),
              concat_all(npl_early), concat_all(npl_late)]
    labels = ['PL\nearly', 'PL\nlate', 'NPL\nearly', 'NPL\nlate']
    colors = [PALETTE[0], PALETTE[0], PALETTE[1], PALETTE[1]]
    alphas = [0.5, 1.0, 0.5, 1.0]

    fig, ax = plt.subplots(figsize=(2.5, 2))
    for x, arr, lab, col, al in zip(range(4), groups, labels, colors, alphas):
        med = np.median(arr)
        iqr = np.percentile(arr, 75) - np.percentile(arr, 25)
        ax.errorbar(x, med, yerr=iqr, fmt='o', color=col, alpha=al,
                    capsize=3, markersize=3, elinewidth=0.5, capthick=0.5)

    tests = []
    for name, x1, x2, a, b in [
        ('PL e/l', 0, 1, groups[0], groups[1]),
        ('NPL e/l', 2, 3, groups[2], groups[3]),
        ('Early PL/NPL', 0, 2, groups[0], groups[2]),
        ('Late PL/NPL', 1, 3, groups[1], groups[3]),
    ]:
        _, p = mannwhitneyu(a, b, alternative='two-sided')
        tests.append({'name': name, 'x1': x1, 'x2': x2, 'p': p})

    pvals = np.array([t['p'] for t in tests])
    _, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')

    y_base = max(np.percentile(arr, 75) for arr in groups) + 1
    for i, (t, pc) in enumerate(zip(tests, pvals_corr)):
        y = y_base + i * 1.5
        stars = sig_text(pc)
        ax.plot([t['x1'], t['x1'], t['x2'], t['x2']], [y, y + 0.3, y + 0.3, y],
                color='k', lw=0.32, clip_on=False)
        offset = 0.1 if stars.lower() in ('ns', 'n.s.') else -0.05
        ax.text((t['x1'] + t['x2']) / 2, y + offset + 0.3, stars,
                ha='center', va='bottom', fontsize=6)
        print(f'  {t["name"]}: p={pc:.2e} {stars}')

    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel('Pop coupling (normalized)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'pop_coupling.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "pop_coupling.svg"}')


if __name__ == '__main__':
    main()
