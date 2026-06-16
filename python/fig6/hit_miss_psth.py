"""Hit vs miss PSTHs with cluster-based permutation test.

Median +/- IQR traces, gold cluster shading, red onset line.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from utils.config import PSTH_NPZ_PATH, CLUSTERS_SEM_PATH, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE

FIG_DIR = OUTPUT_DIR / 'fig6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()

N_PERM = 1000
CLUSTER_ALPHA = 0.05


def find_clusters(stat_arr, threshold):
    sig_mask = np.abs(stat_arr) > threshold
    clusters = []
    in_cluster = False
    start = 0
    for i in range(len(sig_mask)):
        if sig_mask[i] and not in_cluster:
            start = i; in_cluster = True
        elif not sig_mask[i] and in_cluster:
            clusters.append((start, i - 1, np.sum(stat_arr[start:i])))
            in_cluster = False
    if in_cluster:
        clusters.append((start, len(sig_mask) - 1, np.sum(stat_arr[start:])))
    return clusters


def cluster_permutation_test(hit_arr, miss_arr):
    n_hit, n_miss = len(hit_arr), len(miss_arr)
    combined = np.vstack([hit_arr, miss_arr])
    n_total = n_hit + n_miss
    n_bins = hit_arr.shape[1]

    obs_stats = np.zeros(n_bins)
    for b in range(n_bins):
        _, p = mannwhitneyu(hit_arr[:, b], miss_arr[:, b], alternative='two-sided')
        z = np.sign(np.mean(hit_arr[:, b]) - np.mean(miss_arr[:, b]))
        obs_stats[b] = z * (-np.log10(max(p, 1e-20)))

    threshold = -np.log10(CLUSTER_ALPHA)
    obs_clusters = find_clusters(obs_stats, threshold)
    if not obs_clusters:
        return []

    rng = np.random.default_rng(42)
    null_max = []
    for _ in range(N_PERM):
        idx = rng.permutation(n_total)
        ph, pm = combined[idx[:n_hit]], combined[idx[n_hit:]]
        ps = np.zeros(n_bins)
        for b in range(n_bins):
            _, p = mannwhitneyu(ph[:, b], pm[:, b], alternative='two-sided')
            z = np.sign(np.mean(ph[:, b]) - np.mean(pm[:, b]))
            ps[b] = z * (-np.log10(max(p, 1e-20)))
        pc = find_clusters(ps, threshold)
        null_max.append(max(abs(c[2]) for c in pc) if pc else 0)

    null_max = np.array(null_max)
    return [(s, e, np.mean(null_max >= abs(cs))) for s, e, cs in obs_clusters]


def load_or_run_clusters(hit_arr, miss_arr, key, cache_path):
    cached = {}
    if cache_path.exists():
        with open(cache_path) as f:
            cached = json.load(f)
    if key in cached:
        print(f'  Cached: {key}')
        return [(s, e, p) for s, e, p in cached[key]]

    print(f'  Computing: {key} (hit={len(hit_arr)}, miss={len(miss_arr)})...')
    sig = cluster_permutation_test(hit_arr, miss_arr)
    cached[key] = [[int(s), int(e), float(p)] for s, e, p in sig]
    with open(cache_path, 'w') as f:
        json.dump(cached, f, indent=2)
    return sig


def plot_grid(npz_path=None, cache_path=None):
    if npz_path is None:
        npz_path = PSTH_NPZ_PATH
    if cache_path is None:
        cache_path = CLUSTERS_SEM_PATH

    d = np.load(npz_path)
    bc = d['bin_centers']

    hit_c, miss_c, diff_c = PALETTE[2], PALETTE[3], 'black'

    panels = [
        (0, 0, 'pl_early', 'PL early'),
        (0, 1, 'pl_late', 'PL late'),
        (1, 0, 'npl_early', 'NPL early'),
        (1, 1, 'npl_late', 'NPL late'),
    ]

    # Shared limits
    all_diffs, all_fr_max = [], []
    for _, _, key, _ in panels:
        h, m = d[f'{key}_hit'], d[f'{key}_miss']
        if len(h) > 0:
            all_fr_max.append(np.max(np.percentile(h, 75, axis=0)))
        if len(m) > 0:
            all_fr_max.append(np.max(np.percentile(m, 75, axis=0)))
        if len(h) > 0 and len(m) > 0:
            all_diffs.append(np.median(h, axis=0) - np.median(m, axis=0))

    diff_max = max(np.max(np.abs(a)) for a in all_diffs) * 1.2 if all_diffs else 5
    diff_lim = (-diff_max, diff_max)
    fr_ylim = (0, max(all_fr_max) * 1.15 if all_fr_max else 10)

    fig, axes = plt.subplots(2, 2, figsize=(6, 5))

    for row, col, key, title in panels:
        ax = axes[row, col]
        hit_arr, miss_arr = d[f'{key}_hit'], d[f'{key}_miss']

        if len(hit_arr) == 0 or len(miss_arr) == 0:
            ax.set_title(title, fontsize=7)
            continue

        hit_med = np.median(hit_arr, axis=0)
        miss_med = np.median(miss_arr, axis=0)
        hit_lo, hit_hi = np.percentile(hit_arr, 25, axis=0), np.percentile(hit_arr, 75, axis=0)
        miss_lo, miss_hi = np.percentile(miss_arr, 25, axis=0), np.percentile(miss_arr, 75, axis=0)
        diff = hit_med - miss_med

        # Clusters
        sig = load_or_run_clusters(hit_arr, miss_arr, key, cache_path)
        cluster_onset, cluster_end = None, None
        for start, end, p in sig:
            if p < 0.05:
                cluster_onset, cluster_end = bc[start], bc[end]

        if cluster_onset is not None:
            ax.axvline(cluster_onset, color='red', linestyle='-', linewidth=0.5, alpha=0.7)
            ax.text(cluster_onset + 5, fr_ylim[1] * 0.95, f'{cluster_onset:.0f} ms',
                    fontsize=5, color='red', va='top')

        ax.plot(bc, hit_med, color=hit_c, linewidth=1, label='Hit')
        ax.fill_between(bc, hit_lo, hit_hi, color=hit_c, alpha=0.2)
        ax.plot(bc, miss_med, color=miss_c, linewidth=1, label='Miss')
        ax.fill_between(bc, miss_lo, miss_hi, color=miss_c, alpha=0.2)

        ax_r = ax.twinx()
        ax_r.plot(bc, diff, color=diff_c, linewidth=1, alpha=0.7, label='Hit \u2212 Miss')
        ax_r.axhline(0, color=diff_c, linewidth=0.3, linestyle='--', alpha=0.5)
        ax_r.set_ylim(diff_lim)
        ax_r.spines['top'].set_visible(False)
        ax_r.spines['right'].set_visible(True)
        if col == 1:
            ax_r.set_ylabel('Hit \u2212 Miss (Hz)', color=diff_c)
        else:
            ax_r.set_yticklabels([])
        ax_r.tick_params(axis='y', labelcolor=diff_c)

        ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(122, color='black', linestyle=':', linewidth=0.8, alpha=0.6)
        if row == 0 and col == 0:
            ax.text(125, fr_ylim[1] * 0.5, 'Movement\nonset', fontsize=5,
                    va='center', ha='left', color='black', alpha=0.6)
        ax.set_title(title, fontsize=7)
        ax.set_ylim(fr_ylim)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if col == 0:
            ax.set_ylabel('Firing rate (Hz)')
        if row == 1:
            ax.set_xlabel('Time from stim onset (ms)')

    axes[0, 0].legend(fontsize=7, loc='upper right')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'hit_miss_psth.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'\nSaved {FIG_DIR / "hit_miss_psth.svg"}')


if __name__ == '__main__':
    plot_grid()
