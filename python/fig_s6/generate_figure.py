"""Figure S6 -- Imaging movement control (first 5 frames).

Panels:
  A. Example recruited vs not-recruited neuron traces
  B. # neurons recruited in first 5 frames, early vs late phase
  C. Increasing subset: max dF/F in first 5 frames over weeks

Reads data/fig_s6/recruitment_data_first5f.mat (built by matlab/fig_s6/).

Usage (from repo root):
    python -m python.fig_s6.generate_figure
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style

apply_global_style()

WIN_FRAMES = int(sys.argv[1]) if len(sys.argv) > 1 else 5
SUFFIX = f'first{WIN_FRAMES}f'

MAT_PATH = DATA_DIR / 'fig_s6' / f'recruitment_data_{SUFFIX}.mat'
FIG_DIR = OUTPUT_DIR / 'fig_s6'
FIG_DIR.mkdir(parents=True, exist_ok=True)

RED = (0.80, 0.00, 0.00)
BLUE = (0.00, 0.30, 0.85)
GRAY = (0.40, 0.40, 0.40)
LGRAY = (0.90, 0.90, 0.90)


def sig_stars(p):
    if p < 1e-3:
        return '***'
    if p < 1e-2:
        return '**'
    if p < 0.05:
        return '*'
    return 'n.s.'


def main():
    d = loadmat(str(MAT_PATH), squeeze_me=True)
    t_ms = np.asarray(d['t_ms'], dtype=float)
    trace_rec = np.asarray(d['trace_recruited'], dtype=float)
    trace_nrc = np.asarray(d['trace_notrecruited'], dtype=float)
    thr_rec = float(d['thr_recruited'])
    thr_nrc = float(d['thr_notrecruited'])
    first5_ms = float(d['first5_ms'])
    first5Cols = np.atleast_1d(d['first5Cols']).astype(int) - 1  # MATLAB 1-based -> Python 0-based
    earlyVals = np.atleast_1d(d['earlyVals']).astype(float)
    lateVals = np.atleast_1d(d['lateVals']).astype(float)
    p_value = float(d['p_value'])
    amp_weeks = np.atleast_1d(d['amp_weeks']).astype(float)
    amp_mean = np.atleast_1d(d['amp_mean']).astype(float)
    amp_std = np.atleast_1d(d['amp_std']).astype(float)
    amp_p_value = float(d['amp_p_value'])

    fig, axes = plt.subplots(
        1, 3, figsize=(7.2, 1.9),
        gridspec_kw={'width_ratios': [1.5, 0.9, 1.0]},
    )

    # --- Panel A: example neurons (recruited vs not recruited) ---
    ax = axes[0]
    ax.axvspan(0, first5_ms, color=LGRAY, zorder=0, alpha=0.7, lw=0)
    ax.plot(t_ms, trace_rec, color=RED, lw=0.9, zorder=3)
    ax.plot(t_ms, trace_nrc, color=BLUE, lw=0.9, zorder=3)
    ax.scatter(t_ms[first5Cols], trace_rec[first5Cols], s=4, c=[RED], zorder=4, edgecolors='none')
    ax.scatter(t_ms[first5Cols], trace_nrc[first5Cols], s=4, c=[BLUE], zorder=4, edgecolors='none')
    ax.axhline(thr_rec, color=RED, ls='--', lw=0.6, zorder=2)
    ax.axhline(thr_nrc, color=BLUE, ls='--', lw=0.6, zorder=2)
    ax.text(-190, thr_rec - 0.04, r'2$\times$CV', color=RED, fontsize=6, va='top', ha='left')
    ax.text(-190, thr_nrc - 0.04, r'2$\times$CV', color=BLUE, fontsize=6, va='top', ha='left')
    ax.axvline(0, color='k', lw=0.5)
    ax.axvline(first5_ms, color='k', ls=':', lw=0.5)
    ix_r = int(np.nanargmax(trace_rec)); pk_r = float(trace_rec[ix_r])
    ix_n = int(np.nanargmax(trace_nrc)); pk_n = float(trace_nrc[ix_n])
    ax.text(t_ms[ix_r], pk_r + 0.10, 'recruited', color=RED, ha='center', va='bottom',
            fontsize=7, fontweight='bold')
    ax.text(t_ms[ix_n], pk_n + 0.10, 'not recruited', color=BLUE, ha='center', va='bottom',
            fontsize=7, fontweight='bold')
    ax.set_xlim(-200, 800)
    ax.set_ylim(-0.30, 2.80)
    ax.set_xlabel('Time from ICMS onset (ms)')
    ax.set_ylabel(r'$\Delta$F/F')
    ax.set_title('Example neurons')

    # --- Panel B: channel x session count (early vs late) ---
    ax = axes[1]
    rng = np.random.default_rng(0)
    for v, x in zip([earlyVals, lateVals], [1, 2]):
        q1, q2, q3 = np.quantile(v, [0.25, 0.5, 0.75])
        ax.plot([x, x], [q1, q3], color='k', lw=0.75)
        ax.plot([x - 0.08, x + 0.08], [q1, q1], color='k', lw=0.75)
        ax.plot([x - 0.08, x + 0.08], [q3, q3], color='k', lw=0.75)
        ax.plot([x - 0.16, x + 0.16], [q2, q2], color='k', lw=1.3)
        jit = x + 0.06 * rng.standard_normal(len(v))
        ax.scatter(jit, v, s=6, c=[GRAY], alpha=0.7, edgecolors='none', zorder=3)

    allv = np.concatenate([earlyVals, lateVals])
    q1L, q3L = np.quantile(lateVals, [0.25, 0.75])
    fence = q3L + 3 * (q3L - q1L)
    within = allv[allv <= fence]
    visTop = float(within.max()) if within.size else float(allv.max())
    tick = 0.03 * visTop
    yb = visTop * 1.05
    ax.plot([1, 1, 2, 2], [yb - tick, yb, yb, yb - tick], color='k', lw=0.75)
    ax.text(1.5, yb + tick * 1.2, sig_stars(p_value), ha='center', va='bottom', fontsize=10)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Early phase', 'Late phase'])
    ax.set_xlim(0.5, 2.5)
    ax.set_ylim(min(0, float(allv.min())) - 1, yb + 0.18 * visTop)
    ax.set_ylabel(f'# recruited\n(first {WIN_FRAMES} frames)')

    # --- Panel C: increasing subset, weekly trajectory + early/late bracket ---
    ax = axes[2]
    ax.errorbar(amp_weeks, amp_mean, yerr=amp_std, fmt='o-',
                color=RED, markerfacecolor=RED, markeredgecolor='none',
                markersize=2.8, lw=0.9, capsize=4, elinewidth=0.9, capthick=1.0)
    y_top_data = float(np.max(amp_mean + amp_std))
    y_bot_data = float(np.min(amp_mean - amp_std))
    y_range = max(y_top_data - y_bot_data, 1e-6)
    y_bar = y_top_data + 0.10 * y_range
    bracket_h = 0.04 * y_range
    ax.plot([0, 1], [y_bar, y_bar], color='k', lw=0.75)
    ax.plot([2, 4], [y_bar, y_bar], color='k', lw=0.75)
    ax.plot([0.5, 0.5], [y_bar, y_bar + bracket_h], color='k', lw=0.75)
    ax.plot([3, 3], [y_bar, y_bar + bracket_h], color='k', lw=0.75)
    ax.plot([0.5, 3], [y_bar + bracket_h, y_bar + bracket_h], color='k', lw=0.75)
    star_y = y_bar + bracket_h + 0.03 * y_range
    ax.text(1.75, star_y, sig_stars(amp_p_value), ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('Week')
    ax.set_ylabel(r'max $\Delta$F/F' + f'\n(first {WIN_FRAMES} frames)')
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlim(-0.4, 4.4)
    ax.set_ylim(top=star_y + 0.12 * y_range)
    ax.set_title('Increasing subset')

    plt.tight_layout()
    fig.savefig(FIG_DIR / f'recruitment_{SUFFIX}_combined.svg', format='svg', bbox_inches='tight')
    fig.savefig(FIG_DIR / f'recruitment_{SUFFIX}_combined.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
