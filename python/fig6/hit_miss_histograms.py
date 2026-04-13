"""Hit vs miss delta_spks histograms and KL divergence (700ms and 120ms)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, entropy
from scipy.stats.mstats import winsorize
from statsmodels.stats.multitest import multipletests

from utils.config import DF_MERGED_PATH, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE, sig_text

FIG_DIR = OUTPUT_DIR / 'fig6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()
np.random.seed(42)


def bootstrap_kl(dist1, dist2, bins, n_boot=1000, proportion=0.4):
    dist1 = winsorize(np.asarray(dist1), limits=[0.01, 0.01])
    dist2 = winsorize(np.asarray(dist2), limits=[0.01, 0.01])
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

    kl_vals = np.array(kl_vals)
    return kl_vals, np.median(kl_vals), np.percentile(kl_vals, 75) - np.percentile(kl_vals, 25)


def _sig_bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], lw=0.5, c='black')
    fontsize = 6 if text in ('n.s.', 'ns', 'NS') else 7
    ax.text((x1 + x2) / 2, y - 0.004, text, ha='center', va='bottom', fontsize=fontsize)


def run_analysis(df_merged, delta_col, kl_bin_range, hist_bin_range, window_label, hist_n_bins=40):
    print(f'\n{"="*60}\n  Window: {window_label}\n{"="*60}')

    pl_df = df_merged[df_merged['is_pulse_locked'] == True]
    npl_df = df_merged[df_merged['is_pulse_locked'] == False]

    pl_early_vals, pl_late_vals = [], []
    npl_early_vals, npl_late_vals = [], []
    pl_early_hits, pl_early_misses = [], []
    pl_late_hits, pl_late_misses = [], []
    npl_early_hits, npl_early_misses = [], []
    npl_late_hits, npl_late_misses = [], []

    bins = np.histogram_bin_edges([], bins=25, range=kl_bin_range)

    for week in range(5):
        pl_hit = pl_df[(pl_df['rel_week'] == week) & (pl_df['trial_outcome'] == 'hit')][delta_col]
        pl_miss = pl_df[(pl_df['rel_week'] == week) & (pl_df['trial_outcome'] == 'miss')][delta_col]
        npl_hit = npl_df[(npl_df['rel_week'] == week) & (npl_df['trial_outcome'] == 'hit')][delta_col]
        npl_miss = npl_df[(npl_df['rel_week'] == week) & (npl_df['trial_outcome'] == 'miss')][delta_col]

        print(f'  Wk{week}: PL hit={len(pl_hit)} miss={len(pl_miss)}, NPL hit={len(npl_hit)} miss={len(npl_miss)}')

        pl_kl, _, _ = bootstrap_kl(pl_hit, pl_miss, bins)
        npl_kl, _, _ = bootstrap_kl(npl_hit, npl_miss, bins)

        if week in [0, 1]:
            pl_early_vals.append(pl_kl); npl_early_vals.append(npl_kl)
            pl_early_hits.extend(pl_hit); pl_early_misses.extend(pl_miss)
            npl_early_hits.extend(npl_hit); npl_early_misses.extend(npl_miss)
        elif week in [2, 3, 4]:
            pl_late_vals.append(pl_kl); npl_late_vals.append(npl_kl)
            pl_late_hits.extend(pl_hit); pl_late_misses.extend(pl_miss)
            npl_late_hits.extend(npl_hit); npl_late_misses.extend(npl_miss)

    # --- Histograms ---
    color_map = sns.color_palette("deep", 4)
    hit_c, miss_c = color_map[2], color_map[3]
    hist_bins = np.histogram_bin_edges(hist_bin_range, bins=hist_n_bins)

    fig, axes = plt.subplots(2, 2, figsize=(4.5, 4))
    def plot_pair(ax, hits, misses, title, title_color):
        ax.hist(winsorize(np.array(hits), limits=[0.01, 0.01]), bins=hist_bins,
                color=hit_c, density=True, edgecolor='none', alpha=0.6)
        ax.hist(winsorize(np.array(misses), limits=[0.01, 0.01]), bins=hist_bins,
                color=miss_c, density=True, edgecolor='none', alpha=0.6)
        ax.set_title(title, color=title_color)

    plot_pair(axes[0, 0], pl_early_hits, pl_early_misses, 'PL early', color_map[0])
    plot_pair(axes[0, 1], pl_late_hits, pl_late_misses, 'PL late', color_map[0])
    plot_pair(axes[1, 0], npl_early_hits, npl_early_misses, 'NPL early', color_map[1])
    plot_pair(axes[1, 1], npl_late_hits, npl_late_misses, 'NPL late', color_map[1])

    y_max_hist = max(ax.get_ylim()[1] for ax in axes.ravel())
    for ax in axes.ravel():
        ax.set_ylim([0, y_max_hist])
        ax.set_xlabel('\u0394 spikes')
        ax.set_ylabel('Density')
    axes[0, 0].legend(['Hit', 'Miss'], fontsize=6)
    plt.tight_layout()
    fig.savefig(FIG_DIR / f'hit_miss_histograms_{window_label}.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved hit_miss_histograms_{window_label}.svg')

    # --- KL divergence ---
    pl_e = np.concatenate(pl_early_vals); pl_l = np.concatenate(pl_late_vals)
    npl_e = np.concatenate(npl_early_vals); npl_l = np.concatenate(npl_late_vals)

    _, pval_pl = mannwhitneyu(pl_e, pl_l, alternative='two-sided')
    _, pval_npl = mannwhitneyu(npl_e, npl_l, alternative='two-sided')
    _, pval_early = mannwhitneyu(pl_e, npl_e, alternative='two-sided')
    _, pval_late = mannwhitneyu(pl_l, npl_l, alternative='two-sided')
    _, pvals_corr, _, _ = multipletests([pval_pl, pval_npl, pval_early, pval_late], method='fdr_bh')

    print(f'  PL early median: {np.median(pl_e):.4f}, late: {np.median(pl_l):.4f}')
    for name, p in zip(['PL e/l', 'NPL e/l', 'Early PL/NPL', 'Late PL/NPL'], pvals_corr):
        print(f'  {name}: p={p:.2e} {sig_text(p)}')

    x = np.array([0, 1, 2, 3])
    medians = [np.median(g) for g in [pl_e, pl_l, npl_e, npl_l]]
    iqrs = [(np.percentile(g, 75) - np.percentile(g, 25)) / 2 for g in [pl_e, pl_l, npl_e, npl_l]]

    fig, ax = plt.subplots(figsize=(2, 1.5))
    for xi, med, iqr, c, alpha in zip(x, medians, iqrs,
            [PALETTE[0], PALETTE[0], PALETTE[1], PALETTE[1]], [0.5, 1, 0.5, 1]):
        ax.errorbar([xi], [med], yerr=[iqr], fmt='o', color=c, alpha=alpha,
                    capsize=3, linewidth=0.32, markeredgecolor='none',
                    elinewidth=0.32, capthick=0.32, markersize=3)

    ax.set_xticks(x)
    ax.set_xticklabels(['PL\nEarly', 'PL\nLate', 'NPL\nEarly', 'NPL\nLate'])
    ax.set_ylabel('KL divergence')
    ax.set_xlim(-1, 4)

    y_max = max(medians) + max(iqrs)
    dr = y_max - min(medians)
    edge_h = dr * 0.08
    tier_sp = dr * 0.35

    tier1 = y_max + dr * 0.15
    _sig_bracket(ax, 0, 1, tier1, edge_h, sig_text(pvals_corr[0]))
    _sig_bracket(ax, 2, 3, tier1, edge_h, sig_text(pvals_corr[1]))
    tier2 = tier1 + tier_sp
    _sig_bracket(ax, 0, 2, tier2, edge_h, sig_text(pvals_corr[2]))
    tier3 = tier2 + tier_sp
    _sig_bracket(ax, 1, 3, tier3, edge_h, sig_text(pvals_corr[3]))
    ax.set_ylim(top=tier3 + dr * 0.15)

    plt.tight_layout()
    fig.savefig(FIG_DIR / f'divergence_{window_label}.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved divergence_{window_label}.svg')


def main():
    df_merged = pd.read_pickle(DF_MERGED_PATH)
    run_analysis(df_merged, 'delta_spks', (-30, 50), (-20, 60), '700ms')


if __name__ == '__main__':
    main()
