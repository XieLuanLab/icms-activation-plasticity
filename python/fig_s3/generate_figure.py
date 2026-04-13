"""Figure S3 (S-Ephys) -- Electrophysiology supplementary.

Panels:
  1. Filtered trace PL example
  2. Filtered trace NPL example
  3. PCA example (PC2 vs PC1, stim vs non-stim)
  4. Modulation all cells over weeks
  5. Time to max FR all cells over weeks
  6. Modulation at detection threshold

Usage (from repo root):
    python -m python.fig_s3.generate_figure
"""
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
from sklearn.decomposition import PCA
from matplotlib.collections import PathCollection

from utils.config import (RAW_DF_PATH, STIM_NONSTIM_PATH, FILTERED_TRACES_PATH,
                           DATA_DIR, OUTPUT_DIR, Z_CLIP, TMAX_CLIP,
                           ICMS83_THRESHOLDS_CSV)
from utils.plotting import apply_global_style, PALETTE, sig_text, rank_biserial_r
from utils.filters import filter_modulated

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig_s3'
FIG_DIR.mkdir(parents=True, exist_ok=True)

SCATTER_SIZE = 3


# ─── Panel 1-2: Filtered traces PL and NPL ──────────────────────

def plot_filtered_traces():
    print('Filtered traces...')
    d = np.load(FILTERED_TRACES_PATH)

    for label, title in [('pl', 'PL'), ('npl', 'NPL')]:
        trace = d[f'{label}_trace']
        spikes = d[f'{label}_spikes']
        stim_ts = d[f'{label}_stim_ts']
        fs = int(d[f'{label}_fs'])

        time_ms = np.arange(len(trace)) / (fs / 1000)

        fig, ax = plt.subplots(figsize=(2.2, 1.6))
        ax.plot(time_ms, trace, linewidth=0.5, color='k')

        # Stim markers -- full-height vlines + pink shading
        stim_ms = stim_ts / (fs / 1000)
        vline_ylims = ax.get_ylim()
        for i, s in enumerate(stim_ms):
            ax.vlines(s, vline_ylims[0], vline_ylims[1], color='C3',
                      linewidth=0.5, label='Stim onset' if i == 0 else '')
            ax.fill_betweenx(vline_ylims, s - 0.5, s + 1.5,
                             color='C3', alpha=0.2, linewidth=0,
                             label='Blanked' if i == 0 else '')

        # Spike markers
        spike_ms = spikes / (fs / 1000)
        spike_y = trace[spikes.astype(int)]
        marker_offset = 60 if label == 'pl' else 30
        ax.scatter(spike_ms, spike_y - marker_offset, color=PALETTE[0],
                   marker='^', s=12, label='Detected Spike')

        # Scale bars
        ylims = ax.get_ylim()
        yr = ylims[1] - ylims[0]
        v_len = 200 if label == 'pl' else 100
        sb_y = ylims[0] - yr * 0.05
        ax.plot([-2, 3], [sb_y, sb_y], 'k', linewidth=1)  # 5 ms
        ax.text(0.5, sb_y - yr * 0.06, '5 ms', ha='center', fontsize=4)
        ax.plot([-2, -2], [sb_y, sb_y + v_len], 'k', linewidth=1)
        ax.text(-3, sb_y + v_len / 2, f'{v_len} \u00b5V', ha='right', fontsize=4,
                rotation=90, va='center')

        ax.set_ylim([sb_y - yr * 0.1, ylims[1] + yr * 0.1])
        ax.axis('off')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.1),
                  ncol=3, frameon=False, fontsize=4)
        plt.tight_layout(rect=[-0.1, 0, 1, 1])

        fig.savefig(FIG_DIR / f'filtered_trace_{label}.svg',
                    format='svg', bbox_inches='tight')
        plt.close()
        print(f'  Saved filtered_trace_{label}.svg')


# ─── Panel 3: PCA example ───────────────────────────────────────

def plot_pca():
    print('\nPCA example...')
    d = np.load(STIM_NONSTIM_PATH)
    unit_ids = [2, 17]
    unit_labels = ['A', 'B']
    palette = sns.color_palette('deep', n_colors=2)

    X_all, labels = [], []
    for i, uid in enumerate(unit_ids):
        stim = d[f'unit{uid}_stim']
        nonstim = d[f'unit{uid}_nonstim']
        X_all.append(nonstim); X_all.append(stim)
        labels += [f'unit{unit_labels[i]}-nonstim'] * len(nonstim)
        labels += [f'unit{unit_labels[i]}-stim'] * len(stim)

    X_all = np.vstack(X_all)
    labels = np.array(labels)
    X_pca = PCA(n_components=2).fit_transform(X_all)

    fig, ax = plt.subplots(figsize=(1.7, 1.5))
    for i, lbl in enumerate(unit_labels):
        color = palette[i]
        for status, shape in zip(['nonstim', 'stim'], ['^', 'o']):
            idx = labels == f'unit{lbl}-{status}'
            ax.scatter(X_pca[idx, 0], X_pca[idx, 1], color=color, marker=shape,
                       alpha=0.5, label=f'unit{lbl}-{status}', s=8, edgecolors='none')

    ax.set_xlabel('PC 1'); ax.set_ylabel('PC 2')
    ax.set_xlim([-900, 900]); ax.set_ylim([-400, 400])
    ax.legend(fontsize=5, loc='lower left', markerscale=1.5)
    fig.savefig(FIG_DIR / 'pca_example.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('  Saved pca_example.svg')


# ─── Panels 4-5: Modulation and t2max all cells ─────────────────

def plot_metric_over_weeks(df, var, ax, y_clip=None, alt='two-sided'):
    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {s: i * offset - ((len(stim_currents) - 1) / 2) * offset
                  for i, s in enumerate(stim_currents)}
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {s: colors[i] for i, s in enumerate(stim_currents)}

    agg = df.groupby(['rel_week', 'stim_current']).agg(
        median=(var, 'median'),
        q1=(var, lambda x: x.quantile(0.25)),
        q3=(var, lambda x: x.quantile(0.75))
    ).reset_index()
    agg['x'] = agg['rel_week'] + agg['stim_current'].map(offset_map)

    stars_arr = {}
    all_y_vals = []

    for cur in stim_currents:
        cur_agg = agg[agg['stim_current'] == cur]
        cur_df = df[df['stim_current'] == cur]

        for week in cur_df['rel_week'].unique():
            vals = cur_df[cur_df['rel_week'] == week][var].dropna()
            jitter = np.random.uniform(-0.05, 0.05, size=len(vals))
            ax.scatter(week + offset_map[cur] + jitter, vals,
                       s=1.5, alpha=0.15, color=color_map[cur],
                       edgecolors='none', rasterized=False, zorder=1)

        ax.errorbar(cur_agg['x'], cur_agg['median'],
                     yerr=[cur_agg['median'] - cur_agg['q1'],
                           cur_agg['q3'] - cur_agg['median']],
                     fmt='-o', label=f"{cur} \u00b5A", color=color_map[cur],
                     ms=3, capsize=1.5, elinewidth=0.5, capthick=0.5, lw=0.5, zorder=2)
        all_y_vals.extend(cur_agg['q3'].values)

        early = df[(df['rel_week'].isin([0, 1])) & (df['stim_current'] == cur)][var].dropna()
        late = df[(df['rel_week'].isin([2, 3, 4])) & (df['stim_current'] == cur)][var].dropna()
        stat, p = mannwhitneyu(early, late, alternative=alt)
        r = rank_biserial_r(stat, len(early), len(late))
        print(f'    {cur}\u00b5A: n_e={len(early)}, n_l={len(late)}, p={p:.1e}, r={r:.3f}')
        stars_arr[cur] = sig_text(p)

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel('Weeks of training')
    if y_clip: ax.set_ylim(top=y_clip * 1.05)
    ax.legend().set_visible(False)
    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            coll.set_sizes([SCATTER_SIZE])

    # Brackets
    if all_y_vals:
        sv = sorted(all_y_vals)
        y_bar = sv[-2] * 1.05 if len(sv) > 2 and sv[-1] > 2 * sv[len(sv)//2] else sv[-1] * 1.05
        yr = max(ax.get_ylim()[1] - ax.get_ylim()[0], 1)
        bh = yr * 0.03
        ax.plot([0, 1], [y_bar, y_bar], 'k', lw=0.5)
        ax.plot([2, 4], [y_bar, y_bar], 'k', lw=0.5)
        ax.plot([0.5, 0.5], [y_bar, y_bar + bh], 'k', lw=0.5)
        ax.plot([3, 3], [y_bar, y_bar + bh], 'k', lw=0.5)
        ax.plot([0.5, 3], [y_bar + bh, y_bar + bh], 'k', lw=0.5)
        base_y = y_bar + bh + yr * 0.02
        for i, cur in enumerate([4, 5, 6]):
            if stars_arr.get(cur):
                star = stars_arr[cur]
                fs = 6 if star == 'NS' else 8
                ax.text(1.5, base_y + i * yr * 0.10, star,
                        ha='center', fontsize=fs, color=color_map[cur])


def plot_modulation_panels(df_mod):
    for var, ylabel, title, clip, fname, alt in [
        ('z_score', 'Modulation', 'All cells', Z_CLIP, 'modulation_all_cells.svg', 'less'),
        ('t_to_max_10ms_smoothed', 'Time to max FR (ms)', 'All cells', TMAX_CLIP, 't2max_all_cells.svg', 'greater'),
    ]:
        print(f'\n  {title} {var}:')
        fig, ax = plt.subplots(figsize=(2, 1.7))
        plot_metric_over_weeks(df_mod, var, ax, y_clip=clip, alt=alt)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=8)
        plt.tight_layout()
        fig.savefig(FIG_DIR / fname, format='svg', bbox_inches='tight')
        plt.close()
        print(f'    Saved {fname}')


# ─── Panel 6: Modulation at threshold ────────────────────────────

def _load_icms83_thresholds():
    """Load ICMS83 thresholds and average per (session, channel)."""
    if not ICMS83_THRESHOLDS_CSV.exists():
        return pd.DataFrame()
    df = pd.read_csv(ICMS83_THRESHOLDS_CSV)
    df = df[df['session'] != '11-Aug-2023']
    avg = df.groupby(['session', 'channel'])['threshold'].mean().reset_index()
    avg['animal_id'] = 'ICMS83'
    avg = avg.rename(columns={'channel': 'stim_channel', 'threshold': 'detection_threshold'})
    return avg[['animal_id', 'session', 'stim_channel', 'detection_threshold']]


def plot_modulation_at_threshold(df_mod, raw_df):
    """Modulation at detection threshold.

    For each (animal, stim_channel, week):
    1. Each session: round detection threshold -> nearest integer current
    2. Get median z-score of all modulated units at that current for that channel
    3. Take median across sessions within that week
    -> One data point per (animal, channel, week)
    """
    print('\n  Modulation at threshold:')
    thresh_map = (raw_df
                  .drop_duplicates(subset=['animal_id', 'session', 'stim_channel'])
                  [['animal_id', 'session', 'stim_channel', 'rel_week', 'detection_threshold']])

    # Fill in ICMS83 thresholds from separate CSV
    icms83_thr = _load_icms83_thresholds()
    if not icms83_thr.empty:
        wk_map = raw_df[raw_df['animal_id'] == 'ICMS83'][['session', 'rel_week']].drop_duplicates()
        icms83_thr = icms83_thr.merge(wk_map, on='session', how='left')
        thresh_map = thresh_map[thresh_map['animal_id'] != 'ICMS83']
        thresh_map = pd.concat([thresh_map, icms83_thr], ignore_index=True)
        print(f'    Added {len(icms83_thr)} ICMS83 threshold entries')

    thresh_map = thresh_map.dropna(subset=['detection_threshold']).copy()
    thresh_map['thr_current'] = thresh_map['detection_threshold'].round().astype(int)
    thresh_map = thresh_map[thresh_map['thr_current'].between(4, 6)]

    df_ = df_mod.merge(thresh_map, on=['animal_id', 'session', 'stim_channel'],
                       how='inner', suffixes=('', '_chan')).copy()
    df_filt = df_[df_['stim_current'] == df_['thr_current']]
    df_filt = df_filt[df_filt['rel_week'].between(0, 4)]

    var = 'z_score'

    session_med = (df_filt.groupby(['animal_id', 'stim_channel', 'session', 'rel_week'])[var]
                   .median().reset_index(name='session_median_z'))

    week_med = (session_med.groupby(['animal_id', 'stim_channel', 'rel_week'])['session_median_z']
                .median().reset_index(name='week_median_z'))

    print(f'    {len(df_filt)} unit-responses -> {len(session_med)} session medians -> {len(week_med)} (animal, channel, week) data points')

    df_phase = week_med.copy()
    df_phase['phase'] = np.where(df_phase['rel_week'] <= 1, 'early', 'late')
    early_vals = df_phase.loc[df_phase['phase'] == 'early', 'week_median_z'].dropna()
    late_vals = df_phase.loc[df_phase['phase'] == 'late', 'week_median_z'].dropna()

    p_val = 1.0
    if len(early_vals) > 0 and len(late_vals) > 0:
        _, p_val = mannwhitneyu(early_vals, late_vals, alternative='less')
    star = sig_text(p_val)
    print(f'    e={len(early_vals)} med={early_vals.median():.2f}, '
          f'l={len(late_vals)} med={late_vals.median():.2f}, p={p_val:.2e} {star}')

    df_var = (week_med.groupby('rel_week', as_index=False)
              .agg(median=('week_median_z', 'median'),
                   q1=('week_median_z', lambda x: x.quantile(0.25)),
                   q3=('week_median_z', lambda x: x.quantile(0.75)))
              .sort_values('rel_week'))

    fig, ax = plt.subplots(figsize=(2, 1.7))

    for wk in week_med['rel_week'].unique():
        vals = week_med[week_med['rel_week'] == wk]['week_median_z'].dropna()
        jitter = np.random.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(wk + jitter, vals, s=SCATTER_SIZE, alpha=0.2, color='k',
                   edgecolors='none', rasterized=False, zorder=1)

    ax.errorbar(df_var['rel_week'], df_var['median'],
                yerr=[df_var['median'] - df_var['q1'], df_var['q3'] - df_var['median']],
                fmt='-o', ms=3, capsize=1.5, elinewidth=0.5, capthick=0.5,
                c='k', lw=0.5, zorder=2)

    ax.set_ylim(0, 15)

    # Bracket
    y_max = ax.get_ylim()[1]
    yr = y_max
    bh = yr * 0.03
    by = yr * 0.85
    ax.plot([0, 1], [by, by], 'k', lw=0.5)
    ax.plot([2, 4], [by, by], 'k', lw=0.5)
    ax.plot([0.5, 0.5], [by, by + bh], 'k', lw=0.5)
    ax.plot([3, 3], [by, by + bh], 'k', lw=0.5)
    ax.plot([0.5, 3], [by + bh, by + bh], 'k', lw=0.5)
    fs = 6 if star == 'NS' else 8
    ax.text(1.75, by + bh + yr * 0.01, star, ha='center', fontsize=fs)

    ax.set_xlabel('Weeks')
    ax.set_ylabel('Modulation')
    ax.set_title('At detection threshold', fontsize=8)
    ax.set_xticks([0, 1, 2, 3, 4])
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'modulation_at_threshold.svg', format='svg', bbox_inches='tight')
    plt.close()
    print('    Saved modulation_at_threshold.svg')


def collect_stats(df_mod, raw_df):
    """Collect all statistical comparisons for Figure S-Ephys."""
    rows = []

    for var, var_label, alt in [
        ('z_score', 'Modulation all cells', 'less'),
        ('t_to_max_10ms_smoothed', 'T2max all cells', 'greater'),
    ]:
        for cur in sorted(df_mod['stim_current'].unique()):
            early = df_mod[(df_mod['rel_week'].isin([0, 1])) &
                           (df_mod['stim_current'] == cur)][var].dropna()
            late = df_mod[(df_mod['rel_week'].isin([2, 3, 4])) &
                          (df_mod['stim_current'] == cur)][var].dropna()

            if len(early) > 0 and len(late) > 0:
                stat, p = mannwhitneyu(early, late, alternative=alt)
                r = rank_biserial_r(stat, len(early), len(late))
            else:
                stat, p, r = np.nan, np.nan, np.nan

            rows.append({
                'panel': var_label,
                'comparison': f'{cur} uA early vs late',
                'n1': len(early),
                'n2': len(late),
                'median1': f'{early.median():.2f}' if len(early) > 0 else '',
                'median2': f'{late.median():.2f}' if len(late) > 0 else '',
                'test': f'Mann-Whitney U ({alt})',
                'U': f'{stat:.1f}' if not np.isnan(stat) else '',
                'p': f'{p:.2e}' if not np.isnan(p) else '',
                'r': f'{r:.3f}' if not np.isnan(r) else '',
            })

    # --- Modulation at threshold ---
    thresh_map = (raw_df
                  .drop_duplicates(subset=['animal_id', 'session', 'stim_channel'])
                  [['animal_id', 'session', 'stim_channel', 'rel_week', 'detection_threshold']])

    icms83_thr = _load_icms83_thresholds()
    if not icms83_thr.empty:
        wk_map = raw_df[raw_df['animal_id'] == 'ICMS83'][['session', 'rel_week']].drop_duplicates()
        icms83_thr = icms83_thr.merge(wk_map, on='session', how='left')
        thresh_map = thresh_map[thresh_map['animal_id'] != 'ICMS83']
        thresh_map = pd.concat([thresh_map, icms83_thr], ignore_index=True)

    thresh_map = thresh_map.dropna(subset=['detection_threshold']).copy()
    thresh_map['thr_current'] = thresh_map['detection_threshold'].round().astype(int)
    thresh_map = thresh_map[thresh_map['thr_current'].between(4, 6)]

    df_ = df_mod.merge(thresh_map, on=['animal_id', 'session', 'stim_channel'],
                       how='inner', suffixes=('', '_chan')).copy()
    df_thr = df_[df_['stim_current'] == df_['thr_current']]
    df_thr = df_thr[df_thr['rel_week'].between(0, 4)]

    session_med = (df_thr.groupby(['animal_id', 'stim_channel', 'session', 'rel_week'])['z_score']
                   .median().reset_index(name='session_median_z'))
    week_med = (session_med.groupby(['animal_id', 'stim_channel', 'rel_week'])['session_median_z']
                .median().reset_index(name='week_median_z'))

    df_phase = week_med.copy()
    df_phase['phase'] = np.where(df_phase['rel_week'] <= 1, 'early', 'late')
    early_vals = df_phase.loc[df_phase['phase'] == 'early', 'week_median_z'].dropna()
    late_vals = df_phase.loc[df_phase['phase'] == 'late', 'week_median_z'].dropna()

    if len(early_vals) > 0 and len(late_vals) > 0:
        stat, p = mannwhitneyu(early_vals, late_vals, alternative='less')
        r = rank_biserial_r(stat, len(early_vals), len(late_vals))
    else:
        stat, p, r = np.nan, np.nan, np.nan

    rows.append({
        'panel': 'Modulation at threshold',
        'comparison': 'early vs late',
        'n1': len(early_vals),
        'n2': len(late_vals),
        'median1': f'{early_vals.median():.2f}' if len(early_vals) > 0 else '',
        'median2': f'{late_vals.median():.2f}' if len(late_vals) > 0 else '',
        'test': 'Mann-Whitney U (two-sided)',
        'U': f'{stat:.1f}' if not np.isnan(stat) else '',
        'p': f'{p:.2e}' if not np.isnan(p) else '',
        'r': f'{r:.3f}' if not np.isnan(r) else '',
    })

    return pd.DataFrame(rows)


def main():
    # Filtered traces moved to fig_s4 (Panel A)
    plot_pca()

    print('\nLoading raw_df...')
    raw_df = pd.read_pickle(RAW_DF_PATH)
    df_mod = filter_modulated(raw_df, max_z_score=50, min_spikes=50)
    print(f'Modulated: {len(df_mod)}')

    plot_modulation_panels(df_mod)
    plot_modulation_at_threshold(df_mod, raw_df)

    # Stats CSV
    print('\nCollecting stats...')
    stats_df = collect_stats(df_mod, raw_df)
    stats_path = FIG_DIR / 'figure_sEphys_stats.csv'
    stats_df.to_csv(stats_path, index=False)
    print(f'Saved {stats_path}')

    print(f'\nAll saved to {FIG_DIR}')


if __name__ == '__main__':
    main()
