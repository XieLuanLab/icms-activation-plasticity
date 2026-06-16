"""Shared plotting utilities and global style."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

PALETTE = sns.color_palette("deep")


def apply_global_style():
    """Apply consistent matplotlib styling for all figures."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial'],
        'font.size': 7,
        'axes.titlesize': 7,
        'figure.titlesize': 7,
        'legend.fontsize': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'axes.linewidth': 0.32,
        'lines.linewidth': 0.8,
        'lines.markersize': 2.5,
        'lines.markeredgewidth': 0,
        'errorbar.capsize': 5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'xtick.major.width': 0.32,
        'ytick.major.width': 0.32,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,
        'svg.fonttype': 'none',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.01,
    })


def sig_text(p):
    """Convert p-value to significance text."""
    if p <= 0.001: return '***'
    if p <= 0.01: return '**'
    if p <= 0.05: return '*'
    return 'NS'


def rank_biserial_r(u, n1, n2):
    """Rank-biserial correlation effect size for Mann-Whitney U."""
    return 1 - (2 * u) / (n1 * n2)
