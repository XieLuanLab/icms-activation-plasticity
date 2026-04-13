"""Generate all Figure 5 panels.

Panels:
  1. Pie chart (PL classification pattern across currents)
  2. Pyramidal fraction PL vs NPL
  3. PL modulation over weeks
  4. PL time to max FR over weeks
  5. PL P(spike) over weeks
  6. NPL fraction early vs late

Usage (from repo root):
    python -m python.fig5.generate_figure5
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from utils.config import RAW_DF_PATH, OUTPUT_DIR, Z_CLIP, TMAX_CLIP, EARLY_WEEKS, LATE_WEEKS
from utils.plotting import apply_global_style, PALETTE, sig_text
from utils.filters import filter_modulated, filter_pl
from utils.longitudinal import plot_metric_over_weeks, draw_early_late_brackets

apply_global_style()
FIG_DIR = OUTPUT_DIR / 'fig5'
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ─── Panel 1: Pie chart ──────────────────────────────────────────

def classify_pattern(pl_list):
    diffs = [b - a for a, b in zip(pl_list, pl_list[1:])]
    if all(p == 0 for p in pl_list):
        return 'Not PL'
    elif all(p == 1 for p in pl_list):
        return 'All PL'
    elif all(d >= 0 for d in diffs) and any(d > 0 for d in diffs):
        return 'NPL to PL'
    elif all(d <= 0 for d in diffs) and any(d < 0 for d in diffs):
        return 'PL to NPL'
    else:
        return 'Mixed'


def plot_pie(df_mod):
    grouped = df_mod.groupby(['animal_id', 'session', 'stim_channel', 'unit_id'])
    rows = []
    for (animal_id, session, stim_channel, unit_id), group_df in grouped:
        stim_currents = sorted(group_df['stim_current'].unique())
        if len(stim_currents) < 2:
            continue
        pl_dict = {}
        for current in stim_currents:
            is_locked = group_df[group_df['stim_current'] == current]['is_pulse_locked'].any()
            pl_dict[current] = is_locked
        rows.append({'pl_dict': pl_dict})

    final_df = pd.DataFrame(rows)
    categories = []
    for _, row in final_df.iterrows():
        currents = sorted(row['pl_dict'].keys())
        pl_status = [int(row['pl_dict'][c]) for c in currents]
        categories.append(classify_pattern(pl_status))

    final_df['response_category'] = categories
    counts = final_df['response_category'].value_counts()

    fig, ax = plt.subplots(figsize=(3, 3))
    wedges, texts, autotexts = ax.pie(
        counts, labels=None, colors=PALETTE[:len(counts)],
        autopct='%1.1f%%', startangle=140, pctdistance=0.7)
    for at in autotexts:
        at.set_color('white')
        at.set_fontsize(9)

    ax.legend(labels=counts.index, fontsize=7, loc='center left',
              bbox_to_anchor=(0.45, 0.3), facecolor='white', edgecolor='none',
              framealpha=0.7, frameon=True)
    ax.axis('equal')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'pie_chart.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('Saved pie_chart.svg')
    for cat, n in counts.items():
        print(f'  {cat}: {n} ({n / len(final_df) * 100:.1f}%)')


# ─── Panel 2: Pyramidal fraction PL vs NPL ───────────────────────

def plot_pyr_fraction(df_mod):
    df = df_mod.copy()
    df['is_pulse_locked'] = df['is_pulse_locked'].astype(int)

    group_cols = ['animal_id', 'session', 'unit_id', 'rel_week', 'stim_current']
    unit_current = df.groupby(group_cols, as_index=False).apply(
        lambda g: pd.Series({
            'n_ch': g['stim_channel'].nunique(),
            'pl_ch': g.groupby('stim_channel')['is_pulse_locked'].max().sum(),
        }), include_groups=False).reset_index(drop=True)
    unit_current['is_pl_uc'] = (unit_current['pl_ch'] > unit_current['n_ch'] / 2.0).astype(int)

    keys = ['animal_id', 'session', 'unit_id']
    agg = unit_current.groupby(keys)['is_pl_uc'].agg(
        pl_rows='sum', total_rows='count').reset_index()
    agg['is_pl_unit'] = (agg['pl_rows'] > agg['total_rows'] / 2.0).astype(int)

    ct_map = (df.dropna(subset=['cell_type'])
              .groupby(keys)['cell_type']
              .agg(lambda s: s.mode().iat[0] if not s.mode().empty else s.iloc[0])
              .reset_index())
    unit_level = agg.merge(ct_map, on=keys, how='left')
    unit_level = unit_level[unit_level['cell_type'].isin(['pyramidal', 'interneuron'])]

    ct = pd.crosstab(unit_level['is_pl_unit'], unit_level['cell_type']).reindex(
        index=[0, 1], columns=['interneuron', 'pyramidal'], fill_value=0)

    a, b = int(ct.loc[1, 'interneuron']), int(ct.loc[1, 'pyramidal'])
    c, d = int(ct.loc[0, 'interneuron']), int(ct.loc[0, 'pyramidal'])
    _, p_val = fisher_exact([[a, b], [c, d]], alternative='two-sided')

    pl_total = a + b
    npl_total = c + d
    frac_py_pl = b / pl_total if pl_total else 0
    frac_py_npl = d / npl_total if npl_total else 0

    print(f'Pyr fraction: PL={frac_py_pl:.3f} (n={pl_total}), '
          f'NPL={frac_py_npl:.3f} (n={npl_total}), p={p_val:.2e}')

    fig, ax = plt.subplots(figsize=(2.2, 2.5))
    ax.bar(['NPL', 'PL'], [frac_py_npl, frac_py_pl], color=PALETTE[3])
    ax.set_ylim(0, 1)
    ax.set_ylabel('Pyramidal cell fraction')

    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    stars = sig_text(p_val)
    y_ax = 0.85
    ax.plot([0, 0, 1, 1], [y_ax, y_ax + 0.02, y_ax + 0.02, y_ax],
            lw=0.5, c='k', transform=trans, clip_on=False)
    fontsize = 6 if stars == 'NS' else 8
    ax.text(0.5, y_ax + 0.03, stars, ha='center', va='bottom',
            fontsize=fontsize, transform=trans, clip_on=False)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'pyr_fraction.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('Saved pyr_fraction.svg')


# ─── Panels 3-5: PL features over time ──────────────────────────

def plot_pl_features(raw_df):
    df_pl = filter_pl(raw_df)
    print(f'PL units: {len(df_pl)}')

    features = [
        ('z_score', 'Modulation', Z_CLIP, 0.05, 'pl_modulation.svg', 'two-sided'),
        ('t_to_max_10ms_smoothed', 'Time to max FR (ms)', TMAX_CLIP, 10, 'pl_t2max.svg', 'two-sided'),
        ('max_spike_prob', 'P(spike)', 0.15, 0.002, 'pl_spike_prob.svg', 'two-sided'),
    ]

    for var, ylabel, y_clip, stars_y, fname, alt in features:
        fig, ax = plt.subplots(figsize=(2.3, 1.7))
        print(f'\n  PL {var}:')
        ax, stars_arr, all_y_vals = plot_metric_over_weeks(
            df_pl, var=var, alt=alt, ax=ax, y_clip=y_clip)
        draw_early_late_brackets(stars_arr, all_y_vals, ax, stars_y_offset=stars_y)
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        fig.savefig(FIG_DIR / fname, format='svg', bbox_inches='tight')
        plt.close()
        print(f'  Saved {fname}')


# ─── Panel 6: NPL fraction early vs late ─────────────────────────

def plot_npl_fraction(df_mod):
    df = df_mod.copy()
    df['is_pulse_locked'] = df['is_pulse_locked'].astype(int)

    group_cols = ['animal_id', 'session', 'unit_id', 'rel_week', 'stim_current']
    unit_current = df.groupby(group_cols, as_index=False).apply(
        lambda g: pd.Series({
            'n_ch': g['stim_channel'].nunique(),
            'pl_ch': g.groupby('stim_channel')['is_pulse_locked'].max().sum(),
        }), include_groups=False).reset_index(drop=True)
    unit_current['is_pl_uc'] = (unit_current['pl_ch'] > unit_current['n_ch'] / 2.0).astype(int)

    keys = ['animal_id', 'session', 'unit_id']
    agg = unit_current.groupby(keys)['is_pl_uc'].agg(
        pl_rows='sum', total_rows='count').reset_index()
    agg['is_pl_unit'] = (agg['pl_rows'] > agg['total_rows'] / 2.0).astype(int)

    wk_map = unit_current.groupby(keys)['rel_week'].agg(
        lambda s: s.mode().iat[0] if not s.mode().empty else s.iloc[0]).reset_index()
    unit_level = agg.merge(wk_map, on=keys, how='left')
    conditions = [unit_level['rel_week'].isin(EARLY_WEEKS), unit_level['rel_week'].isin(LATE_WEEKS)]
    unit_level['period'] = np.select(conditions, ['Early', 'Late'], default=None)
    unit_level = unit_level[unit_level['period'].isin(['Early', 'Late'])]

    ct = pd.crosstab(unit_level['period'], unit_level['is_pl_unit']).reindex(
        index=['Early', 'Late'], columns=[0, 1], fill_value=0)
    early_npl, early_pl = int(ct.loc['Early', 0]), int(ct.loc['Early', 1])
    late_npl, late_pl = int(ct.loc['Late', 0]), int(ct.loc['Late', 1])

    frac_npl_early = early_npl / (early_npl + early_pl) if (early_npl + early_pl) else 0
    frac_npl_late = late_npl / (late_npl + late_pl) if (late_npl + late_pl) else 0

    _, p_val = fisher_exact([[early_npl, early_pl], [late_npl, late_pl]], alternative='two-sided')
    stars = sig_text(p_val)

    print(f'NPL fraction: early={frac_npl_early:.3f} (n={early_npl + early_pl}), '
          f'late={frac_npl_late:.3f} (n={late_npl + late_pl}), p={p_val:.2e} {stars}')

    fig, ax = plt.subplots(figsize=(2.2, 2.5))
    ax.bar(['Early', 'Late'], [frac_npl_early, frac_npl_late], color=PALETTE[1])
    ax.set_ylim(0, 0.6)
    ax.set_ylabel('NPL fraction')

    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    y_bracket = 0.88
    ax.plot([0, 0, 1, 1], [y_bracket, y_bracket + 0.03, y_bracket + 0.03, y_bracket],
            lw=0.5, c='k', transform=trans, clip_on=False)
    fontsize = 6 if stars == 'NS' else 8
    ax.text(0.5, y_bracket + 0.04, stars, ha='center', va='bottom',
            fontsize=fontsize, transform=trans, clip_on=False)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'npl_fraction.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('Saved npl_fraction.svg')


# ─── Main ────────────────────────────────────────────────────────

def main():
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df_mod = filter_modulated(raw_df)
    print(f'Total modulated: {len(df_mod)}')

    print('\n=== Pie chart (D) ===')
    plot_pie(df_mod)

    print('\n=== Pyramidal fraction (E) ===')
    plot_pyr_fraction(df_mod)

    print('\n=== PL features over time (F-H) ===')
    plot_pl_features(raw_df)

    # NPL fraction (old Fig 5I, removed from manuscript)
    # plot_npl_fraction(df_mod)

    # Example panels (A-C) and tracked unit panels (I-M)
    import importlib.util
    fig5_dir = Path(__file__).resolve().parent

    print('\n=== Example panels (A-C) ===')
    spec = importlib.util.spec_from_file_location('example_panels', fig5_dir / 'example_panels.py')
    ep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ep)
    ep.main()

    print('\n=== Tracked unit panels (I-M) ===')
    spec = importlib.util.spec_from_file_location('tracked_units', fig5_dir / 'tracked_units.py')
    tu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tu)
    tu.main()


if __name__ == '__main__':
    main()
