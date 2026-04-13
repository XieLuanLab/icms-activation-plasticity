"""Figure S5 -- Movement control (symmetric 120ms window analysis).

Layout:
    Row 1 (4 panels): PL mod | NPL mod | PL t2max | NPL t2max

Wheel movement panels saved separately as 1x3.

Usage (from repo root):
    python python/fig_s5/generate_figure.py
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
from matplotlib.collections import PathCollection
from scipy.stats import mannwhitneyu
from scipy.signal import butter, filtfilt
from statsmodels.stats.multitest import multipletests

from utils.config import (RAW_DF_120MS_PATH, POP_COUPLING_DIR, ANIMALS_PC_ALL,
                          OUTPUT_DIR, Z_CLIP, DATA_DIR)
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r
from utils.filters import filter_modulated, filter_pl, filter_npl


apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig_s5'
FIG_DIR.mkdir(parents=True, exist_ok=True)

pl_c = PALETTE[0]
npl_c = PALETTE[1]
SCATTER_SIZE = 3

# Movement control uses the 120ms symmetric window data
TMAX_VAR = 't_to_max_2ms_smoothed'
TMAX_CLIP = 150


def plot_metric(ax, df, metric, ylabel, title, y_clip=None, alternative='two-sided'):
    """Plot a metric over weeks with early-vs-late one-sided test."""
    if len(df) < 5:
        ax.set_title(title, fontsize=8)
        return []

    stats_rows = []
    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {s: i * offset - ((len(stim_currents) - 1) / 2) * offset
                  for i, s in enumerate(stim_currents)}
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {s: colors[i] for i, s in enumerate(stim_currents)}

    agg = df.groupby(['rel_week', 'stim_current']).agg(
        median=(metric, 'median'),
        q1=(metric, lambda x: x.quantile(0.25)),
        q3=(metric, lambda x: x.quantile(0.75))
    ).reset_index()
    agg['x'] = agg['rel_week'] + agg['stim_current'].map(offset_map)

    all_y_vals = []
    stars_arr = {}

    for cur in stim_currents:
        cur_agg = agg[agg['stim_current'] == cur]
        cur_df = df[df['stim_current'] == cur]

        for week in cur_df['rel_week'].unique():
            vals = cur_df[cur_df['rel_week'] == week][metric].dropna()
            jitter = np.random.uniform(-0.05, 0.05, size=len(vals))
            ax.scatter(week + offset_map[cur] + jitter, vals,
                       s=1.5, alpha=0.15, color=color_map[cur],
                       edgecolors='none', rasterized=False, zorder=1)

        ax.errorbar(cur_agg['x'], cur_agg['median'],
                    yerr=[cur_agg['median'] - cur_agg['q1'],
                          cur_agg['q3'] - cur_agg['median']],
                    fmt='-o', label=f"{cur} \u00b5A", color=color_map[cur],
                    ms=3, capsize=1.5, elinewidth=0.5, capthick=0.5, lw=0.5, zorder=2)
        all_y_vals.extend(cur_agg['q3'].values)

        early = df[(df['rel_week'].isin([0, 1])) & (
            df['stim_current'] == cur)][metric].dropna()
        late = df[(df['rel_week'].isin([2, 3, 4])) & (
            df['stim_current'] == cur)][metric].dropna()
        if len(early) > 0 and len(late) > 0:
            stat, p = mannwhitneyu(early, late, alternative=alternative)
            r = rank_biserial_r(stat, len(early), len(late))
            print(
                f'    {cur}\u00b5A: n_e={len(early)}, n_l={len(late)}, U={stat:.0f}, p={p:.1e}, r={r:.3f}')
            stars_arr[cur] = sig_text(p)
            stats_rows.append({
                'panel': title + ' ' + metric, 'comparison': f'{cur} uA early vs late',
                'n1': len(early), 'n2': len(late),
                'median1': f'{early.median():.4f}', 'median2': f'{late.median():.4f}',
                'test': f'Mann-Whitney U ({alternative})',
                'U': f'{stat:.0f}', 'p': f'{p:.2e}', 'r': f'{r:.3f}',
            })

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel('Weeks of training')
    ax.set_title(title, fontsize=8)
    ax.set_ylabel(ylabel)
    ax.legend().set_visible(False)
    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            coll.set_sizes([SCATTER_SIZE])

    if y_clip:
        ax.set_ylim(top=y_clip * 1.05)

    # Brackets
    if all_y_vals and stars_arr:
        sorted_v = sorted(all_y_vals)
        y_bar = sorted_v[-2] * 1.05 if len(sorted_v) > 2 and sorted_v[-1] > 2 * \
            sorted_v[len(sorted_v)//2] else sorted_v[-1] * 1.05
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        bh = y_range * 0.03
        ax.plot([0, 1], [y_bar, y_bar], color='k', lw=0.5)
        ax.plot([2, 4], [y_bar, y_bar], color='k', lw=0.5)
        ax.plot([0.5, 0.5], [y_bar, y_bar + bh], color='k', lw=0.5)
        ax.plot([3, 3], [y_bar, y_bar + bh], color='k', lw=0.5)
        ax.plot([0.5, 3], [y_bar + bh, y_bar + bh], color='k', lw=0.5)

        base_y = y_bar + bh + y_range * 0.02
        for i, cur in enumerate([4, 5, 6]):
            if stars_arr.get(cur):
                star = stars_arr[cur]
                fs = 6 if star == 'NS' else 8
                ax.text(1.5, base_y + i * y_range * 0.08, star,
                        ha='center', fontsize=fs, color=color_map[cur])

    return stats_rows


def plot_pop_coupling(ax, df_mod):
    df_mod = df_mod.rename(columns={'stim_channel': 'stim_ch'})

    def is_early(animal_id, session):
        sel = df_mod[(df_mod['session'] == session) &
                     (df_mod['animal_id'] == animal_id)]
        vals = sel['rel_week'].dropna().unique()
        if len(vals) == 0:
            return 'nan'
        return 'early' if float(vals[0]) < 2 else 'late'

    pl_early, pl_late, npl_early, npl_late = [], [], [], []

    for animal_id in ANIMALS_PC_ALL:
        pkl = POP_COUPLING_DIR / f'{animal_id}_pop_coupling_stpr_120ms.pkl'
        if not pkl.exists():
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        for session, sd in data.items():
            tp = is_early(animal_id, session)
            if tp == 'nan':
                continue
            for (cur, ch), cond in sd.items():
                pls = [v for v in cond.get('pl_pc_norm_dict', {}).values()
                       if v is not None and np.isfinite(v)]
                npls = [v for v in cond.get('npl_pc_norm_dict', {}).values()
                        if v is not None and np.isfinite(v)]
                if tp == 'early':
                    pl_early.extend(pls)
                    npl_early.extend(npls)
                else:
                    pl_late.extend(pls)
                    npl_late.extend(npls)

    def finite(a): return np.asarray(a, float)[
        np.isfinite(np.asarray(a, float))]
    groups = [finite(pl_early), finite(pl_late),
              finite(npl_early), finite(npl_late)]
    labels = ['PL\nEarly', 'PL\nLate', 'NPL\nEarly', 'NPL\nLate']
    colors = [pl_c, pl_c, npl_c, npl_c]
    alphas = [0.5, 1.0, 0.5, 1.0]

    stats_rows = []

    for x, arr, col, al in zip(range(4), groups, colors, alphas):
        if len(arr) == 0:
            continue
        med = np.median(arr)
        iqr = np.percentile(arr, 75) - np.percentile(arr, 25)
        ax.errorbar(x, med, yerr=iqr, fmt='o', color=col, alpha=al,
                    capsize=1.5, markersize=4, elinewidth=0.5, capthick=0.5)

    tests = []
    for name, x1, x2, a, b in [
        ('PL e/l', 0, 1, groups[0], groups[1]),
        ('NPL e/l', 2, 3, groups[2], groups[3]),
        ('E PL/NPL', 0, 2, groups[0], groups[2]),
        ('L PL/NPL', 1, 3, groups[1], groups[3]),
    ]:
        if len(a) > 0 and len(b) > 0:
            stat, p = mannwhitneyu(a, b, alternative='two-sided')
            r = rank_biserial_r(stat, len(a), len(b))
            tests.append({'name': name, 'x1': x1, 'x2': x2, 'p': p})
            stats_rows.append({
                'panel': 'Pop coupling (120ms)', 'comparison': name,
                'n1': len(a), 'n2': len(b),
                'median1': f'{np.median(a):.4f}', 'median2': f'{np.median(b):.4f}',
                'test': 'Mann-Whitney U (two-sided, FDR-BH)',
                'U': f'{stat:.0f}', 'p': '', 'r': f'{r:.3f}',
            })

    if tests:
        pvals = np.array([t['p'] for t in tests])
        _, pc, _, _ = multipletests(pvals, method='fdr_bh')
        for row, p_corr in zip(stats_rows[-len(tests):], pc):
            row['p'] = f'{p_corr:.2e}'

        y_base = max(np.percentile(a, 75) for a in groups if len(a) > 0) + 0.5
        for i, (t, p) in enumerate(zip(tests, pc)):
            y = y_base + i * 0.8
            star_txt = sig_text(p)
            ax.plot([t['x1'], t['x1'], t['x2'], t['x2']], [y, y+0.15, y+0.15, y],
                    color='k', lw=0.5, clip_on=False)
            fs = 6 if star_txt == 'NS' else 8
            ax.text((t['x1']+t['x2'])/2, y+0.17,
                    star_txt, ha='center', fontsize=fs)
            print(f'    {t["name"]}: p={p:.2e} {star_txt}')

    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel('Pop coupling (norm)')
    ax.set_title('')

    return stats_rows


# ─── Wheel movement panels ──────────────────────────────────

FS_WHEEL = 1000.0
ONSET_THRESHOLD_FRAC = 0.5
ONSET_SUSTAINED_MS = 50


def unwrap_by_jump(x, jump_thresh=2000, wrap_val=4096):
    x = np.asarray(x, dtype=float)
    y = x.copy()
    offset = 0
    for k in range(1, len(x)):
        step = x[k] - x[k - 1]
        if step > jump_thresh:
            offset -= wrap_val
        elif step < -jump_thresh:
            offset += wrap_val
        y[k] = x[k] + offset
    return y


def detect_movement_onset(trace, time_ms, threshold,
                          threshold_frac=ONSET_THRESHOLD_FRAC,
                          min_sustained_ms=ONSET_SUSTAINED_MS):
    post_mask = time_ms >= 0
    deviation = np.abs(trace[post_mask])
    post_time = time_ms[post_mask]
    above = deviation > (threshold * threshold_frac)
    min_samples = int(min_sustained_ms * FS_WHEEL / 1000)
    for i in range(len(above) - min_samples):
        if np.all(above[i:i + min_samples]):
            return float(post_time[i])
    return np.nan


def generate_wheel_panels():
    """Generate 1x3 wheel movement figure matching wheel_figure_v3 style."""
    npz_path = DATA_DIR / 'wheel_movement' / 'wheel_traces_raw_icms92.npz'
    if not npz_path.exists():
        print(f'  WARNING: {npz_path} not found, skipping wheel panels')
        return

    _deep = sns.color_palette("deep")
    CLR_HIT = _deep[2]
    CLR_MISS = _deep[3]
    CLR_THRESH = _deep[1]
    CLR_ONSET = _deep[0]
    CLR_RT = _deep[4]

    data = np.load(npz_path)
    time_ms = data['time_ms']
    all_hit = data['hit_traces']
    all_miss = data['miss_traces']
    hit_rt_ms = data['hit_rt_ms']
    hit_current = data['hit_current']
    miss_current = data['miss_current']

    # Downsample 30kHz -> 1kHz
    ds = 30
    hit_ds = all_hit[:, ::ds]
    miss_ds = all_miss[:, ::ds]
    time_ds = time_ms[::ds]

    # Filter to 4-6 uA
    hit_mask = (hit_current >= 4) & (hit_current <= 6)
    miss_mask = (miss_current >= 4) & (miss_current <= 6)
    hit = hit_ds[hit_mask]
    miss = miss_ds[miss_mask]
    hit_rt = hit_rt_ms[hit_mask]

    n_hit = hit.shape[0]
    n_miss = miss.shape[0]
    print(f"  Using {n_hit} hit, {n_miss} miss trials (4-6 uA)")

    # Crop to 0-700ms
    zero_idx = np.argmin(np.abs(time_ds - 0))
    end_idx = np.argmin(np.abs(time_ds - 700))
    hit_post = hit[:, zero_idx:end_idx + 1]
    miss_post = miss[:, zero_idx:end_idx + 1]
    time_post = time_ds[zero_idx:end_idx + 1]

    # Unwrap
    hit_unwrapped = np.zeros_like(hit_post, dtype=float)
    miss_unwrapped = np.zeros_like(miss_post, dtype=float)
    for i in range(n_hit):
        hit_unwrapped[i] = unwrap_by_jump(hit_post[i])
    for i in range(n_miss):
        miss_unwrapped[i] = unwrap_by_jump(miss_post[i])

    # Recenter at t=0
    hit_rel = hit_unwrapped - \
        np.mean(hit_unwrapped[:, 0:3], axis=1, keepdims=True)
    miss_rel = miss_unwrapped - \
        np.mean(miss_unwrapped[:, 0:3], axis=1, keepdims=True)

    # Low-pass filter
    cutoff = 20
    b, a = butter(2, cutoff / (FS_WHEEL / 2), btype='low')
    pad = 100

    hit_filt = np.zeros_like(hit_rel)
    miss_filt = np.zeros_like(miss_rel)
    for i in range(n_hit):
        x_pad = np.pad(hit_rel[i], (pad, pad), mode='edge')
        hit_filt[i] = filtfilt(b, a, x_pad)[pad:-pad]
    for i in range(n_miss):
        x_pad = np.pad(miss_rel[i], (pad, pad), mode='edge')
        miss_filt[i] = filtfilt(b, a, x_pad)[pad:-pad]

    # Flatten hits after RT
    hit_flat = hit_filt.copy()
    for i in range(n_hit):
        rt_idx = np.argmin(np.abs(time_post - hit_rt[i]))
        if 1 <= rt_idx < hit_flat.shape[1]:
            hit_flat[i, rt_idx:] = hit_flat[i, rt_idx - 1]

    # Compute empirical threshold
    rt_displacements = []
    for i in range(n_hit):
        rt_idx = np.argmin(np.abs(time_post - hit_rt[i]))
        if 0 <= rt_idx < hit_filt.shape[1]:
            rt_displacements.append(abs(hit_filt[i, rt_idx]))

    threshold_mv = np.percentile(
        rt_displacements, 25) if rt_displacements else 35.0
    if threshold_mv < 1.0:
        threshold_mv = 35.0
    movement_thresh_mv = threshold_mv * ONSET_THRESHOLD_FRAC

    print(f"  Response threshold (25th pct): {threshold_mv:.1f} mV")
    print(f"  Movement threshold (50%): {movement_thresh_mv:.1f} mV")

    # Detect movement onsets
    hit_onsets = np.array([detect_movement_onset(hit_filt[i], time_post, threshold_mv)
                           for i in range(n_hit)])
    miss_onsets = np.array([detect_movement_onset(miss_filt[i], time_post, threshold_mv)
                            for i in range(n_miss)])

    hit_onsets_valid = hit_onsets[~np.isnan(hit_onsets)]
    miss_onsets_valid = miss_onsets[~np.isnan(miss_onsets)]

    print(f"  Hit onsets: {len(hit_onsets_valid)}/{n_hit}")
    print(f"  Miss onsets: {len(miss_onsets_valid)}/{n_miss}")
    if len(hit_onsets_valid):
        print(f"  Median hit onset: {np.median(hit_onsets_valid):.0f} ms")
    print(f"  Median hit RT: {np.median(hit_rt):.0f} ms")

    # ─── 1x3 figure ──────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(10, 3))

    # === Panel 1: Population overlay ===
    ax = axes[0]
    rng = np.random.default_rng(42)
    plot_n = 200
    hit_plot_idx = rng.choice(n_hit, min(plot_n, n_hit), replace=False)
    miss_plot_idx = rng.choice(n_miss, min(plot_n, n_miss), replace=False)

    for i in miss_plot_idx:
        ax.plot(time_post, miss_filt[i],
                color=CLR_MISS, alpha=0.04, linewidth=0.5)
    for i in hit_plot_idx:
        ax.plot(time_post, hit_flat[i],
                color=CLR_HIT, alpha=0.04, linewidth=0.5)

    ax.plot(time_post, np.mean(miss_filt, axis=0), color=CLR_MISS, linewidth=1.5,
            label=f'Miss mean (n={n_miss})')
    ax.plot(time_post, np.mean(hit_flat, axis=0), color=CLR_HIT, linewidth=1.5,
            label=f'Hit mean (n={n_hit})')

    ax.axhline(threshold_mv, color=CLR_THRESH, linestyle='--', linewidth=0.8, alpha=0.5,
               label=f'Global response threshold ({threshold_mv:.1f} mV)')
    ax.axhline(-threshold_mv, color=CLR_THRESH,
               linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(0, color='k', linewidth=0.3, alpha=0.3)
    ax.set_xlabel('Time from stim onset (ms)')
    ax.set_ylabel('Displacement (mV)')
    ax.set_title('Population overlay')
    ax.legend(fontsize=6, loc='upper left')

    # === Panel 2: Single trial method ===
    ax = axes[1]
    example_idx = 301
    trace = hit_filt[example_idx].copy()
    raw_trace = hit_rel[example_idx]
    onset_ms = hit_onsets[example_idx]
    rt_val = hit_rt[example_idx]

    trace_flat = trace.copy()
    if rt_val <= 700:
        rt_idx = np.argmin(np.abs(time_post - rt_val))
        if 1 <= rt_idx < len(trace_flat):
            trace_flat[rt_idx:] = trace_flat[rt_idx - 1]

    ax.plot(time_post, raw_trace, color=CLR_HIT, alpha=0.2, linewidth=0.8)
    ax.plot(time_post, trace_flat, color=CLR_HIT, linewidth=1.5)

    ax.axhline(movement_thresh_mv, color=CLR_THRESH,
               linestyle='-', linewidth=0.8, alpha=0.7)
    ax.axhline(-movement_thresh_mv, color=CLR_THRESH,
               linestyle='-', linewidth=0.8, alpha=0.7)

    if not np.isnan(onset_ms) and onset_ms <= 700:
        onset_sample = np.argmin(np.abs(time_post - onset_ms))
        ax.axvline(onset_ms, color=CLR_ONSET,
                   linestyle='-', linewidth=1, alpha=0.8)
        ax.scatter([onset_ms], [trace[onset_sample]],
                   color=CLR_ONSET, s=40, zorder=5)

    if rt_val <= 700:
        rt_sample = np.argmin(np.abs(time_post - rt_val))
        ax.axvline(rt_val, color=CLR_RT, linestyle=':', linewidth=1, alpha=0.8)
        ax.scatter([rt_val], [trace[rt_sample]],
                   color=CLR_RT, s=40, zorder=5, marker='D')

    ax.annotate(f'Movement threshold\n(global resp. thresh. x 0.5)',
                xy=(680, movement_thresh_mv), fontsize=6,
                color=CLR_THRESH, ha='right', va='bottom')
    ax.annotate(f'Onset\n{onset_ms:.0f} ms', xy=(onset_ms, trace[onset_sample]),
                xytext=(onset_ms - 60, trace[onset_sample] + 8),
                fontsize=6, color=CLR_ONSET, ha='center',
                arrowprops=dict(arrowstyle='->', color=CLR_ONSET, lw=0.5))
    if rt_val <= 700:
        ax.annotate(f'Response\n{rt_val:.0f} ms', xy=(rt_val, trace[rt_sample]),
                    xytext=(rt_val + 60, trace[rt_sample] - 15),
                    fontsize=6, color=CLR_RT, ha='center',
                    arrowprops=dict(arrowstyle='->', color=CLR_RT, lw=0.5))

    ax.axvline(0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(0, color='k', linewidth=0.3, alpha=0.3)
    ax.set_xlabel('Time from stim onset (ms)')
    ax.set_ylabel('Displacement (mV)')
    ax.set_title('Single trial')

    # === Panel 3: Movement onset distribution (dual y-axis) ===
    ax_left = axes[2]
    ax_right = ax_left.twinx()
    bins = np.arange(0, 750, 25)

    if len(hit_onsets_valid):
        ax_left.hist(hit_onsets_valid, bins=bins, color=CLR_HIT, alpha=0.6,
                     edgecolor='none', label=f'Hit (n={len(hit_onsets_valid)}/{n_hit})')
    ax_left.set_ylabel('Hit count', color=CLR_HIT)
    ax_left.tick_params(axis='y', labelcolor=CLR_HIT)

    if len(miss_onsets_valid):
        ax_right.hist(miss_onsets_valid, bins=bins, color=CLR_MISS, alpha=0.6,
                      edgecolor='none', label=f'Miss (n={len(miss_onsets_valid)}/{n_miss})')
    ax_right.set_ylabel('Miss count', color=CLR_MISS)
    ax_right.tick_params(axis='y', labelcolor=CLR_MISS)
    ax_right.set_ylim(0, ax_left.get_ylim()[1] / 10)

    if len(hit_onsets_valid):
        med_hit_onset = np.median(hit_onsets_valid)
        ax_left.axvline(med_hit_onset, color=CLR_HIT, linestyle='--', linewidth=1.5,
                        label=f'Hit median: {med_hit_onset:.0f} ms')

    if len(miss_onsets_valid):
        med_miss_onset = np.median(miss_onsets_valid)
        ax_left.axvline(med_miss_onset, color=CLR_MISS, linestyle='--', linewidth=1.5,
                        label=f'Miss median: {med_miss_onset:.0f} ms')

    med_rt = np.median(hit_rt)
    ax_left.axvline(med_rt, color=CLR_RT, linestyle=':', linewidth=1.5,
                    label=f'Median RT: {med_rt:.0f} ms')

    ax_left.set_xlabel('Movement onset (ms)')
    ax_left.set_title('Onset latency distribution')

    # Combine legends
    lines1, labels1 = ax_left.get_legend_handles_labels()
    lines2, labels2 = ax_right.get_legend_handles_labels()
    all_lines = lines1 + lines2
    all_labels = labels1 + labels2
    order = []
    for keyword in ['Hit (n=', 'Miss (n=', 'Hit median', 'Miss median', 'Median RT']:
        for i, lab in enumerate(all_labels):
            if keyword in lab and i not in order:
                order.append(i)
                break
    for i in range(len(all_labels)):
        if i not in order:
            order.append(i)
    ax_left.legend([all_lines[i] for i in order], [all_labels[i] for i in order],
                   fontsize=6, loc='upper right')
    ax_right.spines['right'].set_visible(True)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'wheel_movement.svg',
                format='svg', bbox_inches='tight')
    fig.savefig(FIG_DIR / 'wheel_movement.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  Saved wheel_movement.svg')


def main():
    raw_df = pd.read_pickle(RAW_DF_120MS_PATH)
    df_pl = filter_pl(raw_df)
    df_npl = filter_npl(raw_df)
    df_mod = filter_modulated(raw_df)

    print(f'PL: {len(df_pl)}, NPL: {len(df_npl)}, Mod: {len(df_mod)}')

    all_stats = []

    # 4 panels in a single row: PL mod, NPL mod, PL t2max, NPL t2max
    fig, axes = plt.subplots(1, 4, figsize=(13, 2.5))

    print('\n  PL modulation:')
    all_stats.extend(plot_metric(axes[0], df_pl, 'z_score', 'Modulation',
                                 'PL', y_clip=Z_CLIP, alternative='two-sided'))
    print('\n  NPL modulation:')
    all_stats.extend(plot_metric(axes[1], df_npl, 'z_score', '',
                                 'NPL', y_clip=Z_CLIP, alternative='two-sided'))
    print('\n  PL t2max:')
    all_stats.extend(plot_metric(axes[2], df_pl, TMAX_VAR, 'Time to max FR (ms)',
                                 'PL', y_clip=TMAX_CLIP, alternative='two-sided'))
    print('\n  NPL t2max:')
    all_stats.extend(plot_metric(axes[3], df_npl, TMAX_VAR, '',
                                 'NPL', y_clip=175, alternative='two-sided'))

    fig.savefig(FIG_DIR / 'movement_control.svg',
                format='svg', bbox_inches='tight')
    fig.savefig(FIG_DIR / 'movement_control.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\nSaved {FIG_DIR / "movement_control.svg"}')

    # Stats CSV
    stats_df = pd.DataFrame(all_stats)
    stats_path = FIG_DIR / 'movement_control_stats.csv'
    stats_df.to_csv(stats_path, index=False)
    print(f'Saved {stats_path}')

    # --- Wheel movement 1x3 ---
    print('\nWheel movement panels:')
    generate_wheel_panels()


if __name__ == '__main__':
    main()
