"""Figure S8 -- Pop coupling panels (700ms, 120ms, control).

Generates standalone pop coupling comparison panels for supplementary figure.

Usage (from repo root):
    python -m python.fig_s8.generate_figure
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from datetime import datetime

from utils.config import (OUTPUT_DIR, POP_COUPLING_DIR, RAW_DF_PATH,
                           CONTROL_DATA_DIR, CONTROL_PC_ANIMALS,
                           EXPERIMENTAL_ANIMALS, ALL_ANIMALS)
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r
from utils.filters import filter_modulated

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig_s8'
FIG_DIR.mkdir(parents=True, exist_ok=True)

ANIMALS_BEHAV = ALL_ANIMALS
ANIMALS_CTRL = CONTROL_PC_ANIMALS


def finite(a):
    a = np.asarray(a, float)
    return a[np.isfinite(a)]


def plot_pop_coupling(groups, labels, colors, alphas, title, save_name, stats_label):
    fig, ax = plt.subplots(figsize=(2.4, 2))

    for x, arr, col, al in zip(range(len(groups)), groups, colors, alphas):
        if len(arr) == 0:
            continue
        med = np.median(arr)
        iqr = np.percentile(arr, 75) - np.percentile(arr, 25)
        ax.errorbar(x, med, yerr=iqr, fmt='o', color=col, alpha=al,
                    capsize=1.5, markersize=3, elinewidth=0.5, capthick=0.5,
                    markerfacecolor=col, markeredgecolor='none')

    tests = []
    all_tests = [
        ('PL e/l', 0, 1, groups[0], groups[1]),
        ('NPL e/l', 2, 3, groups[2], groups[3]),
        ('E PL/NPL', 0, 2, groups[0], groups[2]),
        ('L PL/NPL', 1, 3, groups[1], groups[3]),
    ]
    # Only show brackets for E PL/NPL and L PL/NPL
    bracket_names = {'E PL/NPL', 'L PL/NPL'}
    for name, x1, x2, a, b in all_tests:
        if len(a) > 0 and len(b) > 0:
            stat, p = mannwhitneyu(a, b, alternative='two-sided')
            r = rank_biserial_r(stat, len(a), len(b))
            tests.append({'name': name, 'x1': x1, 'x2': x2, 'p': p, 'r': r,
                          'n1': len(a), 'n2': len(b)})

    if tests:
        pvals = np.array([t['p'] for t in tests])
        _, pc, _, _ = multipletests(pvals, method='fdr_bh')
        # Place brackets above highest error bar cap with padding
        nonempty = [a for a in groups if len(a) > 0]
        max_cap = max(np.median(a) + (np.percentile(a, 75) - np.percentile(a, 25))
                      for a in nonempty)
        y_base = max_cap + 1.0
        bracket_idx = 0
        for i, (t, p) in enumerate(zip(tests, pc)):
            stars = sig_text(p)
            print(f'  {stats_label} {t["name"]}: n=({t["n1"]},{t["n2"]}), '
                  f'p_fdr={p:.2e}, r={t["r"]:.3f} {stars}')
            if t['name'] in bracket_names:
                y = y_base + bracket_idx * 1.2
                bh = 0.5  # vertical part height
                ax.plot([t['x1'], t['x1'], t['x2'], t['x2']],
                        [y, y + bh, y + bh, y], color='k', lw=0.5, clip_on=False)
                fs = 6 if stars == 'NS' else 8
                ax.text((t['x1'] + t['x2']) / 2, y + bh + 0.15, stars,
                        ha='center', fontsize=fs)
                bracket_idx += 1

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_yticks([0, 5, 10])
    ax.set_ylim(-5, 13)
    plt.tight_layout()
    fig.savefig(FIG_DIR / save_name, format='svg', bbox_inches='tight')
    plt.close()
    print(f'  Saved {save_name}')


def get_mod_df():
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df = filter_modulated(raw_df)
    return df.rename(columns={'stim_channel': 'stim_ch'})


def is_early(df_mod, animal_id, session):
    sel = df_mod[(df_mod['session'] == session) & (df_mod['animal_id'] == animal_id)]
    vals = sel['rel_week'].dropna().unique()
    if len(vals) == 0:
        return 'nan'
    return 'early' if float(vals[0]) < 2 else 'late'


def load_behav_pc(pkl_pattern, df_mod):
    pl_e, pl_l, npl_e, npl_l = [], [], [], []
    for aid in ANIMALS_BEHAV:
        pkl = POP_COUPLING_DIR / pkl_pattern.format(aid)
        if not pkl.exists():
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        for session, sd in data.items():
            tp = is_early(df_mod, aid, session)
            if tp == 'nan':
                continue
            for (cur, ch), cond in sd.items():
                if cur not in [4, 5, 6]:
                    continue
                pls = [v for v in cond.get('pl_pc_norm_dict', {}).values()
                       if v is not None and np.isfinite(v)]
                npls = [v for v in cond.get('npl_pc_norm_dict', {}).values()
                        if v is not None and np.isfinite(v)]
                if tp == 'early':
                    pl_e.extend(pls)
                    npl_e.extend(npls)
                else:
                    pl_l.extend(pls)
                    npl_l.extend(npls)
    return [finite(pl_e), finite(pl_l), finite(npl_e), finite(npl_l)]


def load_ctrl_pc():
    from utils.config import CONTROL_SESSIONS
    pl_e, pl_l, npl_e, npl_l = [], [], [], []
    for aid in ANIMALS_CTRL:
        pkl = POP_COUPLING_DIR / f'{aid}_pop_coupling_control.pkl'
        if not pkl.exists():
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)

        # Get first session date from config
        if aid not in CONTROL_SESSIONS or not CONTROL_SESSIONS[aid]:
            continue
        first_date = None
        for s in CONTROL_SESSIONS[aid]:
            for fmt in ('%d-%b-%Y', '%m-%d-%Y', '%Y-%m-%d'):
                try:
                    d = datetime.strptime(s, fmt)
                    if first_date is None or d < first_date:
                        first_date = d
                    break
                except ValueError:
                    continue

        if first_date is None:
            continue

        for session_date, sd in data.items():
            d = None
            for fmt in ('%d-%b-%Y', '%m-%d-%Y', '%Y-%m-%d'):
                try:
                    d = datetime.strptime(session_date, fmt)
                    break
                except ValueError:
                    continue
            if d is None:
                continue
            days = (d - first_date).days
            wk = 0 if days == 0 else int(np.ceil(days / 7))
            if wk >= 5:
                continue
            period = 'early' if wk < 2 else 'late'

            for (cur, ch), cond in sd.items():
                if cur != 5:
                    continue
                pls = [v for v in cond.get('pl_pc_norm_dict', {}).values()
                       if v is not None and np.isfinite(v)]
                npls = [v for v in cond.get('npl_pc_norm_dict', {}).values()
                        if v is not None and np.isfinite(v)]
                if period == 'early':
                    pl_e.extend(pls)
                    npl_e.extend(npls)
                else:
                    pl_l.extend(pls)
                    npl_l.extend(npls)

    return [finite(pl_e), finite(pl_l), finite(npl_e), finite(npl_l)]


def main():
    labels = ['PL\nEarly', 'PL\nLate', 'NPL\nEarly', 'NPL\nLate']
    colors = [PALETTE[0], PALETTE[0], PALETTE[1], PALETTE[1]]
    alphas = [0.5, 1.0, 0.5, 1.0]

    df_mod = get_mod_df()

    print('=== Behavioral 700ms ===')
    groups_700 = load_behav_pc('{}_pop_coupling.pkl', df_mod)
    plot_pop_coupling(groups_700, labels, colors, alphas,
                      'Pop coupling (700ms)', 'pop_coupling_700ms.svg', '700ms')

    print('\n=== Behavioral 120ms ===')
    groups_120 = load_behav_pc('{}_pop_coupling_stpr_120ms.pkl', df_mod)
    plot_pop_coupling(groups_120, labels, colors, alphas,
                      'Pop coupling (120ms)', 'pop_coupling_120ms.svg', '120ms')

    print('\n=== Control ===')
    groups_ctrl = load_ctrl_pc()
    plot_pop_coupling(groups_ctrl, labels, colors, alphas,
                      'Pop coupling (control)', 'pop_coupling_control.svg', 'Control')

    print(f'\nAll saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
