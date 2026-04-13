"""Fig 4D: Stim vs non-stim waveform overlay and centroid distance."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE

FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()


def main():
    d = np.load(DATA_DIR / 'fig4' / 'stim_nonstim_waveforms.npz')

    unit_ids = [2, 17]
    unit_labels = ['A', 'B']
    palette = sns.color_palette('deep', n_colors=2)
    fs = 30000

    # ─── Waveform overlay: mean stim vs nonstim, 3 sparse channels ──
    # Unit 9 from ICMS92 01-Sep-2023, matching manuscript Fig 4D top.
    PRE_PEAK_SPACE = 150

    d_unit = np.load(DATA_DIR / 'fig4' / 'stim_nonstim_unit9.npz')
    stim_wvfs_u9 = d_unit['stim_wvfs']      # (n_spikes, 90, 3)
    nonstim_wvfs_u9 = d_unit['nonstim_wvfs']  # (n_spikes, 90, 3)
    n_ch_sparse = stim_wvfs_u9.shape[2]
    n_samples_u9 = stim_wvfs_u9.shape[1]

    # Sort sparse channels by depth (ch_locs[:, 1])
    sparse_locs = d_unit['ch_locs']
    depth_order = np.argsort(sparse_locs[:, 1])

    fig_wvf, ax = plt.subplots(figsize=(1.25, 1))
    x = np.arange(n_samples_u9)
    offsets = [j * PRE_PEAK_SPACE for j in range(n_ch_sparse)]

    for ch_idx, local_idx in enumerate(depth_order):
        mean_ns = np.mean(nonstim_wvfs_u9[:, :, local_idx], axis=0)
        std_ns = np.std(nonstim_wvfs_u9[:, :, local_idx], axis=0)
        mean_s = np.mean(stim_wvfs_u9[:, :, local_idx], axis=0)
        std_s = np.std(stim_wvfs_u9[:, :, local_idx], axis=0)

        ax.plot(x, mean_ns + offsets[ch_idx], color='k', linewidth=0.75,
                alpha=0.4, label='Non-stim' if ch_idx == 0 else '')
        ax.fill_between(x, mean_ns - std_ns + offsets[ch_idx],
                        mean_ns + std_ns + offsets[ch_idx], color='k', alpha=0.2)

        ax.plot(x, mean_s + offsets[ch_idx], color='r', linewidth=0.75,
                alpha=0.4, label='Stim' if ch_idx == 0 else '')
        ax.fill_between(x, mean_s - std_s + offsets[ch_idx],
                        mean_s + std_s + offsets[ch_idx], color='r', alpha=0.2)

    y_lo = -250
    ax.set_ylim([y_lo, offsets[-1] + 50])
    ax.axis('off')
    ax.legend(loc='lower right', frameon=True, fontsize=3)

    # Scale bars
    ax.plot([0, 30], [y_lo + 10, y_lo + 10], 'k', linewidth=1)
    ax.plot([0, 0], [y_lo + 10, y_lo + 110], 'k', linewidth=1)
    ax.text(0, y_lo - 60, '1 ms', fontsize=5)
    ax.text(-8, y_lo + 10, r'100 $\mu$V', fontsize=5, rotation=90)

    fig_wvf.savefig(FIG_DIR / 'stim_nonstim_waveforms.svg', format='svg',
                     bbox_inches='tight')
    plt.close(fig_wvf)
    print(f'Saved {FIG_DIR / "stim_nonstim_waveforms.svg"}')

    # ─── Centroid distance comparison (10 PCs, all sessions) ────────
    from scipy.stats import mannwhitneyu

    centroid_path = DATA_DIR / 'fig4' / 'centroid_distances.npz'
    if not centroid_path.exists():
        print('\ncentroid_distances.npz not found — skipping centroid plot')
        return

    cd = np.load(centroid_path)
    within_arr = cd['within']
    between_arr = cd['between']

    _, p = mannwhitneyu(within_arr, between_arr, alternative='two-sided')

    def p_stars(p):
        if p < 0.001: return '***'
        if p < 0.01: return '**'
        if p < 0.05: return '*'
        return 'NS'

    print(f'\nCentroid distances (10 PCs, all sessions):')
    print(f'  Within: n={len(within_arr)}, median={np.median(within_arr):.1f}')
    print(f'  Between: n={len(between_arr)}, median={np.median(between_arr):.1f}')
    print(f'  p (two-sided): {p:.2e} {p_stars(p)}')

    fig2, ax2 = plt.subplots(figsize=(1.7, 1.5))

    for vals, x in [(within_arr, 0), (between_arr, 1)]:
        med = np.median(vals)
        q25, q75 = np.percentile(vals, 25), np.percentile(vals, 75)
        ax2.errorbar(x, med, yerr=[[med - q25], [q75 - med]], fmt='o', color='k',
                     capsize=3, markeredgecolor='none', elinewidth=0.5,
                     capthick=0.5, markersize=4)

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(['Within unit\n(stim vs nonstim)',
                         'Between units\n(stim vs stim)'], fontsize=5)
    ax2.set_ylabel('Centroid Distance')
    ax2.set_xlim(-0.5, 1.5)
    ax2.set_ylim(bottom=0)

    y_max = max(np.percentile(within_arr, 75), np.percentile(between_arr, 75)) * 1.1
    bh = y_max * 0.05
    ax2.plot([0, 0, 1, 1], [y_max, y_max + bh, y_max + bh, y_max], color='k', lw=0.5)
    ax2.text(0.5, y_max + bh + 2, p_stars(p), ha='center', fontsize=8)
    ax2.set_ylim(top=y_max + bh + 30)

    fig2.savefig(FIG_DIR / 'centroid_distance.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "centroid_distance.svg"}')


if __name__ == '__main__':
    main()
