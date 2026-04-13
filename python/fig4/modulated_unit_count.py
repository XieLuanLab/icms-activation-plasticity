"""Modulated units per stim channel over weeks (ICMS92 at 4 uA)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from utils.filters import filter_modulated
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r
from utils.config import RAW_DF_PATH, OUTPUT_DIR
from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')


FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    apply_global_style()
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df_mod = filter_modulated(raw_df)

    # Filter for 4 uA and ICMS92
    df_sub = df_mod[(df_mod['stim_current'] == 4) &
                    (df_mod['animal_id'] == 'ICMS92')]

    fig, ax = plt.subplots(figsize=(2, 1.6))

    df_unique = df_sub.drop_duplicates(
        subset=['animal_id', 'session', 'unit_id', 'stim_channel'])
    df_counts = (df_unique
                 .groupby(['animal_id', 'rel_week', 'session', 'stim_channel'])
                 .size().reset_index(name='n_modulated'))

    # Scatter
    jitter = np.random.uniform(-0.15, 0.15, len(df_counts))
    ax.scatter(df_counts['rel_week'] + jitter, df_counts['n_modulated'],
               s=15, alpha=0.5, color=PALETTE[0], edgecolors='none')

    # Median line
    med = df_counts.groupby('rel_week')['n_modulated'].median().reset_index()
    ax.plot(med['rel_week'], med['n_modulated'], color=PALETTE[0], linewidth=1)

    # Early vs late
    early = df_counts[df_counts['rel_week'] <= 1]['n_modulated']
    late = df_counts[df_counts['rel_week'] > 1]['n_modulated']
    stat, p = mannwhitneyu(early, late, alternative='less')
    r = rank_biserial_r(stat, len(early), len(late))
    print(f"U={stat}, p={p:.4g}, r={r:.3f}, n_early={len(early)}, n_late={len(late)}")

    # Bracket
    y_bar = df_counts['n_modulated'].max() + 5
    ax.plot([0, 1], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([2, 4], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([0.5, 0.5], [y_bar, y_bar + 2], color='black', lw=0.5)
    ax.plot([3, 3], [y_bar, y_bar + 2], color='black', lw=0.5)
    ax.plot([0.5, 3], [y_bar + 2, y_bar + 2], color='black', lw=0.5)
    ax.text(1.5, y_bar + 2, sig_text(p), ha='center', fontsize=8)

    ax.set_xticks(np.arange(5))
    ax.set_ylim(bottom=0)
    ax.set_title("Modulated units per\nstim channel at 4 \u00b5A (Mouse1)")
    ax.set_xlabel('Weeks of training')
    ax.set_ylabel('Count per session')
    plt.tight_layout()

    fig.savefig(FIG_DIR / 'modulated_unit_count.svg',
                format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "modulated_unit_count.svg"}')


if __name__ == '__main__':
    main()
