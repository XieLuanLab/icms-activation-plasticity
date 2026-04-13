"""Modulation (z-score) and t2max over training weeks — pyramidal and interneuron."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

from utils.config import RAW_DF_PATH, OUTPUT_DIR, Z_CLIP, TMAX_CLIP
from utils.plotting import apply_global_style, sig_text, rank_biserial_r
from utils.filters import filter_modulated
from utils.longitudinal import plot_metric_over_weeks, draw_early_late_brackets

FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)


# TODO: plot_cell_type_metric was in the old figure4/utils.py but is not in
#       utils/longitudinal.py. Inlined here until it is added to the shared utils.
def plot_cell_type_metric(df, var, ax, cell_type='pyramidal', alt='two-sided',
                          star_size=8, stars_y_offset=0, stars_spacing=0.05,
                          y_clip=None, iqr_filter=False):
    """Plot metric over weeks for a specific cell type with early/late brackets."""
    cell_df = df[df['cell_type'] == cell_type]
    ax, stars_arr, all_y_vals = plot_metric_over_weeks(
        cell_df, var=var, alt=alt, ax=ax, y_clip=y_clip, iqr_filter=iqr_filter)
    ax = draw_early_late_brackets(
        stars_arr, all_y_vals, ax, star_size=star_size,
        stars_y_offset=stars_y_offset)
    ax.legend().set_visible(False)
    return ax


def main():
    apply_global_style()
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df_mod = filter_modulated(raw_df, max_z_score=100, min_spikes=50)

    for var, ylabel, fname_prefix, y_clip, alt in [
        ('z_score', 'Modulation', 'mod_z_score', Z_CLIP, 'two-sided'),
        ('t_to_max_10ms_smoothed', 'Time to max FR (ms)', 'mod_t2max', TMAX_CLIP, 'two-sided'),
    ]:
        for cell_type in ['pyramidal', 'interneuron']:
            fig, ax = plt.subplots(figsize=(1.7, 1.2))

            if cell_type == 'interneuron':
                star_size, stars_y_offset, stars_spacing = 6, 0.2, 0.09
            else:
                star_size, stars_y_offset, stars_spacing = 8, 0, 0.05

            if var == 't_to_max_10ms_smoothed' and cell_type == 'interneuron':
                stars_y_offset = 5

            plot_cell_type_metric(
                df_mod, var, ax, cell_type=cell_type, alt=alt,
                star_size=star_size, stars_y_offset=stars_y_offset,
                stars_spacing=stars_spacing, y_clip=y_clip)

            ax.set_xlabel('Weeks of training')
            ax.set_ylabel(ylabel)
            if var == 't_to_max_10ms_smoothed':
                ax.set_ylim(0, 380)
                ax.set_yticks([0, 100, 200, 300])
            plt.tight_layout()

            save_path = FIG_DIR / f'{fname_prefix}_{cell_type}.svg'
            fig.savefig(save_path, format='svg', bbox_inches='tight')
            plt.close()
            print(f'Saved {save_path}')


if __name__ == '__main__':
    main()
