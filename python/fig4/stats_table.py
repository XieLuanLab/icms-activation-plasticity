"""Generate comprehensive stats table (CSV) for all Figure 4 comparisons."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from utils.config import RAW_DF_PATH, OUTPUT_DIR
from utils.plotting import rank_biserial_r
from utils.filters import filter_modulated, filter_pl, filter_npl

FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def collect_stats():
    """Collect all statistical comparisons for Figure 4."""
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df_mod = filter_modulated(raw_df, max_z_score=100, min_spikes=50)

    rows = []

    # --- Cell type counts ---
    df_unique_all = raw_df.drop_duplicates(subset=['animal_id', 'session', 'unit_id'])
    t2p_all = df_unique_all['t2p_ms'].dropna()
    t2p_all = t2p_all[t2p_all > 0]
    n_py_all = (t2p_all >= 0.45).sum()
    n_in_all = (t2p_all < 0.45).sum()

    df_unique_mod = df_mod.drop_duplicates(subset=['animal_id', 'session', 'unit_id'])
    t2p_mod = df_unique_mod['t2p_ms'].dropna()
    t2p_mod = t2p_mod[t2p_mod > 0]
    n_py_mod = (t2p_mod >= 0.45).sum()
    n_in_mod = (t2p_mod < 0.45).sum()

    rows.append({
        'panel': 'Cell type histogram',
        'comparison': 'All units',
        'group1': 'Pyramidal', 'n1': int(n_py_all),
        'group2': 'Interneuron', 'n2': int(n_in_all),
        'test': 'count', 'U': '', 'p': '', 'r': '',
        'n_animals': len(raw_df['animal_id'].unique()),
    })
    rows.append({
        'panel': 'Cell type histogram',
        'comparison': 'Modulated units',
        'group1': 'Pyramidal', 'n1': int(n_py_mod),
        'group2': 'Interneuron', 'n2': int(n_in_mod),
        'test': 'count', 'U': '', 'p': '', 'r': '',
        'n_animals': len(df_mod['animal_id'].unique()),
    })

    # --- Modulation and t2max over weeks (per cell type, per current) ---
    for cell_type in ['pyramidal', 'interneuron']:
        cell_df = df_mod[df_mod['cell_type'] == cell_type]

        for var, var_label, alt in [
            ('z_score', 'Modulation', 'less'),              # early < late
            ('t_to_max_10ms_smoothed', 'Time to max FR', 'greater'),  # early > late
        ]:
            for cur in [4, 5, 6]:
                early = cell_df[(cell_df['rel_week'].isin([0, 1])) &
                                (cell_df['stim_current'] == cur)][var].dropna()
                late = cell_df[(cell_df['rel_week'].isin([2, 3, 4])) &
                               (cell_df['stim_current'] == cur)][var].dropna()

                if len(early) > 0 and len(late) > 0:
                    stat, p = mannwhitneyu(early, late, alternative=alt)
                    r = rank_biserial_r(stat, len(early), len(late))
                else:
                    stat, p, r = np.nan, np.nan, np.nan

                rows.append({
                    'panel': f'{var_label} ({cell_type})',
                    'comparison': f'{cur} \u00b5A early vs late',
                    'group1': f'Early (wk 0-1)',
                    'n1': len(early),
                    'group2': f'Late (wk 2-4)',
                    'n2': len(late),
                    'median1': f'{early.median():.2f}' if len(early) > 0 else '',
                    'median2': f'{late.median():.2f}' if len(late) > 0 else '',
                    'test': f'Mann-Whitney U ({alt})',
                    'U': f'{stat:.1f}' if not np.isnan(stat) else '',
                    'p': f'{p:.2e}' if not np.isnan(p) else '',
                    'r': f'{r:.3f}' if not np.isnan(r) else '',
                    'n_animals': len(cell_df['animal_id'].unique()),
                })

    # --- Modulated unit count (ICMS92, 4 uA) ---
    df_4ua = df_mod[(df_mod['stim_current'] == 4) & (df_mod['animal_id'] == 'ICMS92')]
    df_unique = df_4ua.drop_duplicates(
        subset=['animal_id', 'session', 'unit_id', 'stim_channel'])
    df_counts = (df_unique
                 .groupby(['animal_id', 'rel_week', 'session', 'stim_channel'])
                 .size().reset_index(name='n_modulated'))

    early_counts = df_counts[df_counts['rel_week'] <= 1]['n_modulated']
    late_counts = df_counts[df_counts['rel_week'] > 1]['n_modulated']

    if len(early_counts) > 0 and len(late_counts) > 0:
        stat, p = mannwhitneyu(early_counts, late_counts, alternative='two-sided')
        r = rank_biserial_r(stat, len(early_counts), len(late_counts))
    else:
        stat, p, r = np.nan, np.nan, np.nan

    rows.append({
        'panel': 'Modulated unit count',
        'comparison': 'ICMS92 4 \u00b5A early vs late',
        'group1': 'Early (wk 0-1)',
        'n1': len(early_counts),
        'group2': 'Late (wk 2-4)',
        'n2': len(late_counts),
        'median1': f'{early_counts.median():.1f}' if len(early_counts) > 0 else '',
        'median2': f'{late_counts.median():.1f}' if len(late_counts) > 0 else '',
        'test': 'Mann-Whitney U (two-sided)',
        'U': f'{stat:.1f}' if not np.isnan(stat) else '',
        'p': f'{p:.2e}' if not np.isnan(p) else '',
        'r': f'{r:.3f}' if not np.isnan(r) else '',
        'n_animals': 1,
    })

    # --- Summary counts ---
    n_sessions = raw_df[['animal_id', 'session']].drop_duplicates().shape[0]
    n_units = raw_df[['animal_id', 'session', 'unit_id']].drop_duplicates().shape[0]
    n_mod_units = df_mod[['animal_id', 'session', 'unit_id']].drop_duplicates().shape[0]
    n_animals = len(raw_df['animal_id'].unique())

    rows.append({
        'panel': 'Summary',
        'comparison': 'Dataset overview',
        'group1': f'{n_animals} animals',
        'n1': n_sessions,
        'group2': f'{n_units} total units',
        'n2': n_mod_units,
        'test': 'descriptive',
        'U': '', 'p': '', 'r': '',
        'n_animals': n_animals,
    })

    return pd.DataFrame(rows)


def pretty_print(df):
    """Print stats table as a formatted, readable table."""
    col_widths = {col: max(len(col), df[col].astype(str).str.len().max())
                  for col in df.columns}

    # Header
    header = ' | '.join(col.ljust(col_widths[col]) for col in df.columns)
    sep = '-+-'.join('-' * col_widths[col] for col in df.columns)
    print(header)
    print(sep)

    # Rows
    for _, row in df.iterrows():
        line = ' | '.join(str(row[col]).ljust(col_widths[col]) for col in df.columns)
        print(line)


def main():
    print('Collecting Figure 4 statistics...\n')
    df = collect_stats()

    save_path = FIG_DIR / 'figure4_stats.csv'
    df.to_csv(save_path, index=False)
    print(f'Saved {save_path}\n')
    pretty_print(df)


if __name__ == '__main__':
    main()
