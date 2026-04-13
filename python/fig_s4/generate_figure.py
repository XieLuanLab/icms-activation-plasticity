"""Figure S4 (S-PL) -- PL/NPL detailed supplementary characterization.

Panels:
  A. Filtered traces (PL and NPL)
  B. PL latency vs distance from stim (4/5/6 uA)
  C. PL jitter vs distance from stim (4/5/6 uA)
  D. NPL modulation over weeks
  E. NPL time to max FR over weeks

Usage (from repo root):
    python -m python.fig_s4.generate_figure
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
import statsmodels.api as sm
from matplotlib.collections import PathCollection

from utils.config import (RAW_DF_PATH, POP_COUPLING_DIR, ANIMALS_PC,
                           OUTPUT_DIR, Z_CLIP, TMAX_CLIP, DATA_DIR,
                           FILTERED_TRACES_PATH)
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r
from utils.filters import filter_modulated, filter_npl

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig_s4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

SCATTER_SIZE = 3


# ─── Panels 1-2: Latency and jitter vs distance from stim ───────

def plot_latency_jitter(raw_df):
    print('Latency and jitter vs stim distance...')
    df = filter_modulated(raw_df).copy()
    df['distance_from_stim'] = abs(df['stim_ch_depth_um'] - df['unit_depth_um'])

    X = pd.DataFrame({
        'stim_current': df['stim_current'],
        'distance_from_stim': df['distance_from_stim'],
        'jitter': df['jitter'],
        'latency': df['latency'],
    })
    y = df['is_pulse_locked'].astype(int)
    X = X.dropna()
    y = y.loc[X.index]

    bin_edges = np.arange(60, 800 + 120, 120)
    bin_labels = [f'{s}\u2013{e}' for s, e in zip(bin_edges[:-1], bin_edges[1:])]
    bin_centers = [(s + e) / 2 for s, e in zip(bin_edges[:-1], bin_edges[1:])]
    label_to_center = dict(zip(bin_labels, bin_centers))

    X['stim_distance_bin'] = pd.cut(X['distance_from_stim'], bins=bin_edges,
                                     labels=bin_labels, include_lowest=True)
    X['label'] = y.map({0: 'NPL', 1: 'PL'})

    delta = 0.2 * np.median(np.diff(sorted(label_to_center.values())))
    offsets = {c: off for c, off in zip([4, 5, 6], (-delta, 0.0, delta))}

    stats_rows = []

    for feature in ['latency', 'jitter']:
        fig, ax = plt.subplots(figsize=(3, 2.2))
        for i, cur in enumerate([4, 5, 6]):
            subset = X[(X['stim_current'] == cur) & (X['label'] == 'PL')]
            stats = subset.groupby('stim_distance_bin', observed=False)[
                feature].agg(['mean', 'std'])

            x_vals = np.array([label_to_center[lbl] for lbl in stats.index], dtype=float)
            x_vals += offsets[cur]

            ax.errorbar(x_vals, stats['mean'], yerr=stats['std'],
                        fmt='-o', color=PALETTE[i], ecolor=PALETTE[i],
                        lw=0.5, ms=3, mfc=PALETTE[i], mec='none',
                        elinewidth=0.5, capsize=1.5, capthick=0.5,
                        label=f'{cur} \u00b5A', zorder=2)

        ylabel = 'Jitter ($ms^2$)' if feature == 'jitter' else 'Latency (ms)'
        ax.set_xlabel('Distance from Stim (\u00b5m)')
        ax.set_ylabel(ylabel)
        ax.set_title('Pulse-locked responses', fontsize=8)
        ax.legend(fontsize=6)

        plt.tight_layout()
        fig.savefig(FIG_DIR / f'pl_{feature}.svg', format='svg', bbox_inches='tight')
        plt.close()
        print(f'  Saved pl_{feature}.svg')

        # OLS regression: feature ~ distance_from_stim, per current
        print(f'\n  OLS regression ({feature}):')
        for cur in [4, 5, 6]:
            subset = X[(X['stim_current'] == cur) & (X['label'] == 'PL')].dropna(
                subset=[feature, 'distance_from_stim'])
            if len(subset) < 5:
                continue
            x_ols = sm.add_constant(subset['distance_from_stim'])
            model = sm.OLS(subset[feature], x_ols).fit()
            slope = model.params['distance_from_stim']
            p_val = model.pvalues['distance_from_stim']
            r2 = model.rsquared
            print(f'    {cur} \u00b5A: slope = {slope:.4f}, p = {p_val:.3f}, r\u00b2 = {r2:.2f}, n = {len(subset)}')
            stats_rows.append({
                'panel': f'PL {feature}', 'comparison': f'{cur} \u00b5A OLS',
                'n1': len(subset), 'n2': '',
                'median1': f'{slope:.4f}', 'median2': '',
                'test': 'OLS (feature ~ distance)',
                'U': f'r\u00b2={r2:.2f}', 'p': f'{p_val:.3e}', 'r': f'{slope:.4f}',
            })

        # OLS on binned medians for monotonic segment (latency only)
        if feature == 'latency':
            print(f'\n  OLS on binned medians (120-840 \u00b5m, latency):')
            for cur in [4, 5, 6]:
                subset = X[(X['stim_current'] == cur) & (X['label'] == 'PL')].dropna(
                    subset=[feature, 'distance_from_stim'])
                subset_range = subset[(subset['distance_from_stim'] >= 120) &
                                      (subset['distance_from_stim'] <= 840)]
                if len(subset_range) < 5:
                    continue
                binned = subset_range.groupby('stim_distance_bin', observed=False)[feature].median().dropna()
                bin_x = np.array([label_to_center[lbl] for lbl in binned.index])
                bin_y = binned.values
                if len(bin_x) < 3:
                    continue
                x_ols = sm.add_constant(bin_x)
                model = sm.OLS(bin_y, x_ols).fit()
                slope = model.params[1]
                p_val = model.pvalues[1]
                r2 = model.rsquared
                speed = 1.0 / slope if slope > 0 else float('inf')
                print(f'    {cur} uA: slope = {slope:.3f} ms/um (speed ~ {speed:.2f} m/s), '
                      f'p = {p_val:.3f}, r2 = {r2:.2f}')
                stats_rows.append({
                    'panel': 'PL latency (binned 120-840\u00b5m)',
                    'comparison': f'{cur} \u00b5A OLS binned',
                    'n1': len(bin_x), 'n2': '',
                    'median1': f'{slope:.4f}', 'median2': f'speed={speed:.2f} m/s',
                    'test': 'OLS (binned medians)',
                    'U': f'r\u00b2={r2:.2f}', 'p': f'{p_val:.3e}', 'r': f'{slope:.4f}',
                })

    return stats_rows


# ─── Panels 3-4: NPL modulation and t2max over weeks ────────────

def plot_npl_metrics(raw_df):
    df_npl = filter_npl(raw_df)
    print(f'\nNPL units: {len(df_npl)}')

    stats_rows = []

    for var, ylabel, clip, fname in [
        ('z_score', 'Modulation', Z_CLIP, 'npl_modulation.svg'),
        ('t_to_max_10ms_smoothed', 'Time to max FR (ms)', TMAX_CLIP, 'npl_t2max.svg'),
    ]:
        print(f'\n  NPL {var}:')
        fig, ax = plt.subplots(figsize=(2, 1.7))

        offset = 0.15
        stim_currents = sorted(df_npl['stim_current'].unique())
        offset_map = {s: i * offset - ((len(stim_currents) - 1) / 2) * offset
                      for i, s in enumerate(stim_currents)}
        colors = sns.color_palette("deep", n_colors=len(stim_currents))
        color_map = {s: colors[i] for i, s in enumerate(stim_currents)}

        agg = df_npl.groupby(['rel_week', 'stim_current']).agg(
            median=(var, 'median'),
            q1=(var, lambda x: x.quantile(0.25)),
            q3=(var, lambda x: x.quantile(0.75))
        ).reset_index()
        agg['x'] = agg['rel_week'] + agg['stim_current'].map(offset_map)

        all_y = []
        stars = {}

        for cur in stim_currents:
            ca = agg[agg['stim_current'] == cur]
            cd = df_npl[df_npl['stim_current'] == cur]

            for wk in cd['rel_week'].unique():
                vals = cd[cd['rel_week'] == wk][var].dropna()
                jitter = np.random.uniform(-0.05, 0.05, size=len(vals))
                ax.scatter(wk + offset_map[cur] + jitter, vals,
                           s=1.5, alpha=0.15, color=color_map[cur],
                           edgecolors='none', rasterized=False, zorder=1)

            ax.errorbar(ca['x'], ca['median'],
                         yerr=[ca['median'] - ca['q1'], ca['q3'] - ca['median']],
                         fmt='-o', color=color_map[cur], ms=3, capsize=1.5,
                         elinewidth=0.5, capthick=0.5, lw=0.5, zorder=2)
            all_y.extend(ca['q3'].values)

            early = df_npl[(df_npl['rel_week'].isin([0, 1])) & (df_npl['stim_current'] == cur)][var].dropna()
            late = df_npl[(df_npl['rel_week'].isin([2, 3, 4])) & (df_npl['stim_current'] == cur)][var].dropna()
            alt = 'less' if var == 'z_score' else 'greater'
            stat, p = mannwhitneyu(early, late, alternative=alt)
            r = rank_biserial_r(stat, len(early), len(late))
            print(f'    {cur}\u00b5A: n_e={len(early)}, n_l={len(late)}, p={p:.1e}, r={r:.3f}')
            stars[cur] = sig_text(p)
            stats_rows.append({
                'panel': f'NPL {var}', 'comparison': f'{cur} \u00b5A early vs late',
                'n1': len(early), 'n2': len(late),
                'median1': f'{early.median():.3f}', 'median2': f'{late.median():.3f}',
                'test': 'Mann-Whitney U', 'U': f'{stat:.0f}',
                'p': f'{p:.2e}', 'r': f'{r:.3f}',
            })

        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_xlabel('Weeks of training')
        ax.set_ylabel(ylabel)
        ax.set_title('NPL', fontsize=8)
        ax.legend().set_visible(False)
        for coll in ax.collections:
            if isinstance(coll, PathCollection):
                coll.set_sizes([SCATTER_SIZE])
        if clip:
            ax.set_ylim(top=clip * 1.05)

        # Brackets
        if all_y:
            sv = sorted(all_y)
            yb = sv[-2] * 1.05 if len(sv) > 2 and sv[-1] > 2 * sv[len(sv)//2] else sv[-1] * 1.05
            yr = max(ax.get_ylim()[1] - ax.get_ylim()[0], 1)
            bh = yr * 0.03
            ax.plot([0, 1], [yb, yb], 'k', lw=0.5)
            ax.plot([2, 4], [yb, yb], 'k', lw=0.5)
            ax.plot([0.5, 0.5], [yb, yb + bh], 'k', lw=0.5)
            ax.plot([3, 3], [yb, yb + bh], 'k', lw=0.5)
            ax.plot([0.5, 3], [yb + bh, yb + bh], 'k', lw=0.5)
            base = yb + bh + yr * 0.02
            for i, cur in enumerate([4, 5, 6]):
                if stars.get(cur):
                    fs = 6 if stars[cur] == 'NS' else 8
                    ax.text(1.5, base + i * yr * 0.10, stars[cur],
                            ha='center', fontsize=fs, color=color_map[cur])

        plt.tight_layout()
        fig.savefig(FIG_DIR / fname, format='svg', bbox_inches='tight')
        plt.close()
        print(f'    Saved {fname}')

    return stats_rows


# ─── Panel 5: Pop coupling ──────────────────────────────────────

def plot_pop_coupling(raw_df):
    print('\nPop coupling...')
    df_mod = filter_modulated(raw_df).rename(columns={'stim_channel': 'stim_ch'})

    def is_early(animal_id, session):
        sel = df_mod[(df_mod['session'] == session) & (df_mod['animal_id'] == animal_id)]
        vals = sel['rel_week'].dropna().unique()
        if len(vals) == 0: return 'nan'
        return 'early' if float(vals[0]) < 2 else 'late'

    pl_early, pl_late, npl_early, npl_late = [], [], [], []

    for animal_id in ANIMALS_PC:
        pkl = POP_COUPLING_DIR / f'{animal_id}_pop_coupling.pkl'
        if not pkl.exists():
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        for session, sd in data.items():
            tp = is_early(animal_id, session)
            if tp == 'nan': continue
            for (cur, ch), cond in sd.items():
                if cur not in [4, 5, 6]: continue
                pls = list(cond.get('pl_pc_norm_dict', {}).values())
                npls = list(cond.get('npl_pc_norm_dict', {}).values())
                if tp == 'early':
                    pl_early.extend(pls); npl_early.extend(npls)
                else:
                    pl_late.extend(pls); npl_late.extend(npls)

    def finite(a): return np.asarray(a, float)[np.isfinite(np.asarray(a, float))]
    groups = [finite(pl_early), finite(pl_late), finite(npl_early), finite(npl_late)]
    labels = ['PL\nEarly', 'PL\nLate', 'NPL\nEarly', 'NPL\nLate']
    colors = [PALETTE[0], PALETTE[0], PALETTE[1], PALETTE[1]]
    alphas = [0.5, 1.0, 0.5, 1.0]

    fig, ax = plt.subplots(figsize=(2.5, 2))
    for x, arr, col, al in zip(range(4), groups, colors, alphas):
        if len(arr) == 0: continue
        med = np.median(arr)
        iqr = np.percentile(arr, 75) - np.percentile(arr, 25)
        ax.errorbar(x, med, yerr=iqr, fmt='o', color=col, alpha=al,
                    capsize=1.5, markersize=3, elinewidth=0.5, capthick=0.5)

    tests = []
    stats_rows = []
    for name, x1, x2, a, b in [
        ('PL early vs late', 0, 1, groups[0], groups[1]),
        ('NPL early vs late', 2, 3, groups[2], groups[3]),
        ('Early PL vs NPL', 0, 2, groups[0], groups[2]),
        ('Late PL vs NPL', 1, 3, groups[1], groups[3]),
    ]:
        if len(a) > 0 and len(b) > 0:
            stat, p = mannwhitneyu(a, b, alternative='two-sided')
            r = rank_biserial_r(stat, len(a), len(b))
            tests.append({'name': name, 'x1': x1, 'x2': x2, 'p': p})
            stats_rows.append({
                'panel': 'Pop coupling', 'comparison': name,
                'n1': len(a), 'n2': len(b),
                'median1': f'{np.median(a):.3f}', 'median2': f'{np.median(b):.3f}',
                'test': 'Mann-Whitney U', 'U': f'{stat:.0f}',
                'p': f'{p:.2e}', 'r': f'{r:.3f}',
            })

    if tests:
        pvals = np.array([t['p'] for t in tests])
        _, pc, _, _ = multipletests(pvals, method='fdr_bh')
        y_base = max(np.percentile(a, 75) for a in groups if len(a) > 0) + 0.5
        for i, (t, p) in enumerate(zip(tests, pc)):
            y = y_base + i * 1.2
            star_txt = sig_text(p)
            ax.plot([t['x1'], t['x1'], t['x2'], t['x2']], [y, y+0.2, y+0.2, y],
                    color='k', lw=0.5, clip_on=False)
            fs = 6 if star_txt == 'NS' else 8
            ax.text((t['x1']+t['x2'])/2, y+0.22, star_txt, ha='center', fontsize=fs)
            print(f'  {t["name"]}: p_fdr={p:.2e} {star_txt}')

    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel('Pop coupling (normalized)')
    ax.set_title('Population coupling', fontsize=8)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'pop_coupling.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'  Saved pop_coupling.svg')
    return stats_rows


def plot_filtered_traces():
    """Panel A: Filtered traces for PL and NPL example units."""
    print('Filtered traces...')
    if not FILTERED_TRACES_PATH.exists():
        print('  filtered_traces.npz not found - skipping')
        return
    d = np.load(FILTERED_TRACES_PATH)

    for label, title in [('pl', 'PL'), ('npl', 'NPL')]:
        trace = d[f'{label}_trace']
        spikes = d[f'{label}_spikes']
        stim_ts = d[f'{label}_stim_ts']
        fs = int(d[f'{label}_fs'])

        time_ms = np.arange(len(trace)) / (fs / 1000)
        fig, ax = plt.subplots(figsize=(2.2, 1.6))
        ax.plot(time_ms, trace, linewidth=0.5, color='k')

        stim_ms = stim_ts / (fs / 1000)
        vline_ylims = ax.get_ylim()
        for i, s in enumerate(stim_ms):
            ax.vlines(s, vline_ylims[0], vline_ylims[1], color='C3',
                      linewidth=0.5, label='Stim onset' if i == 0 else '')
            ax.fill_betweenx(vline_ylims, s - 0.5, s + 1.5,
                             color='C3', alpha=0.2, linewidth=0,
                             label='Blanked' if i == 0 else '')

        spike_ms = spikes / (fs / 1000)
        spike_y = trace[spikes.astype(int)]
        marker_offset = 60 if label == 'pl' else 30
        ax.scatter(spike_ms, spike_y - marker_offset, color=PALETTE[0],
                   marker='^', s=12, label='Detected Spike')

        ylims = ax.get_ylim()
        yr = ylims[1] - ylims[0]
        v_len = 200 if label == 'pl' else 100
        sb_y = ylims[0] - yr * 0.05
        ax.plot([-2, 3], [sb_y, sb_y], 'k', linewidth=1)
        ax.text(0.5, sb_y - yr * 0.06, '5 ms', ha='center', fontsize=4)
        ax.plot([-2, -2], [sb_y, sb_y + v_len], 'k', linewidth=1)
        ax.text(-3, sb_y + v_len / 2, f'{v_len} \u00b5V', ha='right', fontsize=4,
                rotation=90, va='center')

        ax.set_ylim([sb_y - yr * 0.1, ylims[1] + yr * 0.1])
        ax.axis('off')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.1),
                  ncol=3, frameon=False, fontsize=4)
        plt.tight_layout(rect=[-0.1, 0, 1, 1])

        fig.savefig(FIG_DIR / f'filtered_trace_{label}.svg',
                    format='svg', bbox_inches='tight')
        plt.close()
        print(f'  Saved filtered_trace_{label}.svg')


def main():
    raw_df = pd.read_pickle(RAW_DF_PATH)
    all_stats = []

    plot_filtered_traces()
    all_stats.extend(plot_latency_jitter(raw_df))
    all_stats.extend(plot_npl_metrics(raw_df))
    # Pop coupling is in S8, not S4

    # Save stats
    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        stats_path = FIG_DIR / 'figure_sPL_stats.csv'
        stats_df.to_csv(stats_path, index=False)
        print(f'\nSaved {stats_path}')

    print(f'\nAll saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
