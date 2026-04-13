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


def add_sig_bracket(ax, x1, x2, y, h=0.9, text="ns", lw=0.32):
    """Draw a significance bracket between x1 and x2 at height y."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y],
            color="k", lw=lw, clip_on=False)
    text_offset = 0.2 if text == "ns" else -0.1
    ax.text((x1 + x2) / 2, y + text_offset, text,
            ha="center", va="bottom")


def p_to_stars(p):
    """Convert p-value to significance stars."""
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 5e-2: return "*"
    return "ns"


def sig_text(p):
    """Convert p-value to significance text."""
    if p <= 0.001: return '***'
    if p <= 0.01: return '**'
    if p <= 0.05: return '*'
    return 'NS'


def rank_biserial_r(u, n1, n2):
    """Rank-biserial correlation effect size for Mann-Whitney U."""
    return 1 - (2 * u) / (n1 * n2)


def add_scale_bars_wvf(ax, h_pos, v_pos, h_length, v_length,
                        line_width=1, h_label="", v_label="",
                        v_label_x_offset=0, v_label_y_offset=0,
                        h_label_x_offset=0, h_label_y_offset=0):
    """Add scale bars to a waveform plot."""
    hx, hy = h_pos
    vx, vy = v_pos
    ax.plot([hx, hx + h_length], [hy, hy], 'k', linewidth=line_width)
    ax.plot([vx, vx], [vy, vy + v_length], 'k', linewidth=line_width)
    if h_label:
        ax.text(hx + h_label_x_offset, hy + h_label_y_offset, h_label,
                fontsize=4, ha='center', va='top')
    if v_label:
        ax.text(vx + v_label_x_offset, vy + v_label_y_offset, v_label,
                fontsize=4, ha='right', va='center', rotation=90)
