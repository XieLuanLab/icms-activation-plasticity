"""Longitudinal analysis and plotting utilities.

Plot metrics over weeks with early vs late comparisons.
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.collections import PathCollection

from .plotting import sig_text, rank_biserial_r


def plot_metric_over_weeks(df, var, ax, alt='two-sided', y_clip=None, iqr_filter=False):
    """Plot metric over weeks with scatter + median +/- IQR, early vs late test."""
    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {s: i * offset - ((len(stim_currents) - 1) / 2) * offset
                  for i, s in enumerate(stim_currents)}
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {s: colors[i] for i, s in enumerate(stim_currents)}

    agg = df.groupby(['rel_week', 'stim_current']).agg(
        median=(var, 'median'),
        q1=(var, lambda x: x.quantile(0.25)),
        q3=(var, lambda x: x.quantile(0.75))
    ).reset_index()
    agg['x'] = agg['rel_week'] + agg['stim_current'].map(offset_map)

    all_y_vals = []
    stars_arr = {}

    for cur in stim_currents:
        cur_agg = agg[agg['stim_current'] == cur]
        cur_df = df[df['stim_current'] == cur]

        for week in cur_df['rel_week'].unique():
            vals = cur_df[cur_df['rel_week'] == week][var].dropna()
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

        test_df = df[df['stim_current'] == cur]
        if iqr_filter:
            vals_all = test_df[var].dropna()
            q1, q3 = vals_all.quantile(0.25), vals_all.quantile(0.75)
            iqr = q3 - q1
            test_df = test_df[(test_df[var] >= q1 - 1.5 * iqr) & (test_df[var] <= q3 + 1.5 * iqr)]

        early = test_df[test_df['rel_week'].isin([0, 1])][var].dropna()
        late = test_df[test_df['rel_week'].isin([2, 3, 4])][var].dropna()
        stat, p = mannwhitneyu(early, late, alternative=alt)
        r = rank_biserial_r(stat, len(early), len(late))
        print(f"  {cur}\u00b5A: n_early={len(early)}, n_late={len(late)}, "
              f"U={stat}, p={p:.1e}, r={r:.3f}")
        stars_arr[cur] = sig_text(p)

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel('Weeks of Training')

    if y_clip is not None:
        ax.set_ylim(top=y_clip * 1.05)

    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            coll.set_sizes([3])

    ax.legend().set_visible(False)
    return ax, stars_arr, all_y_vals


def draw_early_late_brackets(stars_arr, all_y_vals, ax, star_size=8,
                              stars_y_offset=0):
    """Draw early vs late comparison brackets."""
    stim_currents = [4, 5, 6]
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {s: colors[i] for i, s in enumerate(stim_currents)}

    if not all_y_vals:
        return ax

    sorted_vals = sorted(all_y_vals)
    max_val = sorted_vals[-1]
    median_val = sorted_vals[len(sorted_vals) // 2]
    y_bar = sorted_vals[-2] * 1.05 if len(sorted_vals) > 2 and max_val > 2 * median_val else max_val * 1.05

    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else 1

    bracket_h = y_range * 0.03
    ax.plot([0, 1], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([2, 4], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([0.5, 0.5], [y_bar, y_bar + bracket_h], color='black', lw=0.5)
    ax.plot([3, 3], [y_bar, y_bar + bracket_h], color='black', lw=0.5)
    ax.plot([0.5, 3], [y_bar + bracket_h, y_bar + bracket_h], color='black', lw=0.5)

    base_y = y_bar + bracket_h + y_range * 0.02 + stars_y_offset
    line_spacing = y_range * 0.10

    top_y = base_y
    for i, curr in enumerate([4, 5, 6]):
        if stars_arr.get(curr):
            star = stars_arr[curr]
            is_ns = star.upper() in ['NS', 'N.S.']
            fontsize = 6 if is_ns else star_size
            text_y = base_y + i * line_spacing
            ax.text(1.5, text_y, star, ha='center', fontsize=fontsize,
                    color=color_map[curr])
            top_y = max(top_y, text_y)

    needed = top_y + y_range * 0.1
    if needed > ax.get_ylim()[1]:
        ax.set_ylim(ax.get_ylim()[0], needed)

    return ax
