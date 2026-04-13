"""Generate Figure S7 -- Electrophysiological Unit Tracking Validation.

Panels:
  A: Example tracked waveforms (overlaid across sessions, rainbow gradient)
  B: Waveform distance -- tracked pairs vs nearby non-matched units
  C: Modulation trends restricted to tracked units only
     - PL z-score, NPL z-score, PL t2max, NPL t2max, PL P(spike)

Usage:
    python -m fig_s7.generate_figure
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import json
import pickle
from datetime import datetime
from collections import Counter
import math

from utils.config import (RAW_DF_PATH, OUTPUT_DIR, TRACK_DIR,
                           TRACKED_TEMPLATES_PATH as CACHE_PATH,
                           EXPERIMENTAL_ANIMALS, ICMS83_SESSIONS,
                           Z_CLIP, TMAX_CLIP, PSPIKE_CLIP,
                           EARLY_WEEKS, LATE_WEEKS,
                           MOUSE_NAMES, ANIMALS_SESSIONS)
from utils.filters import filter_modulated, filter_pl, filter_npl
from utils.plotting import apply_global_style, PALETTE, sig_text
from utils.longitudinal import plot_metric_over_weeks, draw_early_late_brackets

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig_s7'
FIG_DIR.mkdir(parents=True, exist_ok=True)

FS = 30000.0

EXCLUDE_TRACKS = {'ICMS93_6'}  # flagged: amplitude shift + channel drift
MIN_TRACK_SESSIONS = 3


def mouse_name(animal_id):
    return MOUSE_NAMES.get(animal_id, animal_id)


def track_display_name(track_id):
    """Convert e.g. 'ICMS92_7' -> 'M1 #7'."""
    parts = track_id.rsplit('_', 1)
    if len(parts) == 2:
        aid, num = parts
        mn = MOUSE_NAMES.get(aid, aid)
        short = mn.replace('Mouse ', 'M')
        return f'{short} #{num}'
    return track_id


# --- Data loading --------------------------------------------------------

def load_raw_df():
    """Load raw_df_700ms.pkl."""
    return pd.read_pickle(RAW_DF_PATH)


def filter_to_tracked(df, tracked_set):
    """Filter df to only rows matching tracked units."""
    mask = df.apply(lambda r: (r['animal_id'], r['session'], r['unit_id']) in tracked_set, axis=1)
    return df[mask]


def get_tracked_unit_ids():
    """Return set of (animal_id, session, unit_id) for tracked units
    with >= MIN_TRACK_SESSIONS sessions."""
    with open(TRACK_DIR / 'all_tracks.json') as f:
        all_tracks = json.load(f)
    tracked = set()
    for track in all_tracks:
        if track['n_sessions'] < MIN_TRACK_SESSIONS:
            continue
        animal = track['animal_id']
        for session, uid in track['sessions'].items():
            tracked.add((animal, session, uid))
    return tracked


HANDPICKED_EXAMPLES = ['ICMS92_2', 'ICMS92_4', 'ICMS93_11']


def _get_example_tracks(cache, n=3):
    """Return handpicked example tracks, sorted by session count (most first)."""
    result = []
    for tid in HANDPICKED_EXAMPLES:
        if tid in cache:
            result.append((tid, cache[tid]))
    # Sort by number of valid sessions, descending
    result.sort(key=lambda x: -sum(1 for v in x[1]['sessions'].values() if v is not None))
    return result[:n]


# --- Panel A: Example tracked waveforms ----------------------------------

def _rel_week_for_animal(session, animal_id):
    """Compute relative week for a session given animal's first session."""
    first_session = ANIMALS_SESSIONS[animal_id][0]
    d = datetime.strptime(session, '%d-%b-%Y')
    d0 = datetime.strptime(first_session, '%d-%b-%Y')
    days = (d - d0).days
    return 0 if days == 0 else int(math.ceil(days / 7))


def plot_panel_a(fig, gs_a):
    """Overlaid waveforms for 3 example tracked units with rainbow gradient."""
    with open(CACHE_PATH, 'rb') as f:
        cache = pickle.load(f)

    examples = _get_example_tracks(cache, n=3)
    gradient_cmap = plt.cm.rainbow

    gs_inner = gs_a.subgridspec(1, 3, wspace=0.3)

    for ex_idx, (tid, td) in enumerate(examples):
        ax = fig.add_subplot(gs_inner[0, ex_idx])
        animal_id = td['animal_id']
        sessions = td['sessions']
        first_session = ANIMALS_SESSIONS[animal_id][0]
        sessions_chrono = sorted(
            [s for s in sessions.keys() if sessions[s] is not None],
            key=lambda s: datetime.strptime(s, '%d-%b-%Y'))
        n_sess = len(sessions_chrono)

        if n_sess == 0:
            ax.axis('off')
            continue

        # Collect waveforms on peak channel
        waveforms = []
        for session in sessions_chrono:
            info = sessions[session]
            if info is None:
                waveforms.append(None)
                continue
            ch = info['primary_ch']
            t1, t2 = 0, min(181, info['template_mean'].shape[0])
            waveforms.append(info['template_mean'][t1:t2, ch])

        # Plot with rainbow gradient + slight vertical offset per session
        y_offset_step = 5
        total_offset = (n_sess - 1) * y_offset_step
        for i, (session, wvf) in enumerate(zip(sessions_chrono, waveforms)):
            if wvf is None:
                continue
            frac = i / max(n_sess - 1, 1)
            color = gradient_cmap(frac)
            y_off = i * y_offset_step - total_offset / 2
            x = np.arange(len(wvf)) / (FS / 1000)
            ax.plot(x, wvf + y_off, color=color, linewidth=1.2, alpha=0.85)

        ax.set_xlim(0, x[-1])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        def _rel_week(s):
            d = datetime.strptime(s, '%d-%b-%Y')
            d0 = datetime.strptime(first_session, '%d-%b-%Y')
            days = (d - d0).days
            return 0 if days == 0 else int(math.ceil(days / 7))

        unit_label = chr(ord('A') + ex_idx)  # A, B, C
        ax.set_title(f'Unit {unit_label}',
                     fontsize=7, pad=3)

        # Scale bar on each example
        y_min, y_max = ax.get_ylim()
        sb_x, sb_y = 0, y_min + (y_max - y_min) * 0.05
        ax.plot([sb_x, sb_x + 1.0], [sb_y, sb_y], 'k', linewidth=1)
        ax.text(sb_x + 0.5, sb_y - (y_max - y_min) * 0.08, '1 ms',
                ha='center', fontsize=5)
        ax.plot([sb_x, sb_x], [sb_y, sb_y + 100], 'k', linewidth=1)
        ax.text(sb_x - 0.15, sb_y + 50, '100 \u00b5V', ha='right',
                fontsize=5, rotation=90, va='center')

        # Colorbar-style gradient legend on last example
        if ex_idx == len(examples) - 1:
            sm = plt.cm.ScalarMappable(cmap=gradient_cmap,
                                        norm=plt.Normalize(0, 1))
            cbar_ax = fig.add_axes([ax.get_position().x1 + 0.005,
                                     ax.get_position().y0 + 0.02,
                                     0.006, ax.get_position().height * 0.5])
            cb = fig.colorbar(sm, cax=cbar_ax)
            cb.set_ticks([0, 1])
            cb.set_ticklabels(['Early', 'Late'])
            cb.ax.tick_params(labelsize=5, length=1)


# --- Panel B: Track raster -----------------------------------------------

def plot_panel_b(ax):
    """Raster of all tracked units: each row is a track, dots at detected weeks.
    Sorted by track duration (longest on top), then by first week."""
    with open(TRACK_DIR / 'all_tracks.json') as f:
        all_tracks = json.load(f)
    with open(CACHE_PATH, 'rb') as f:
        cache = pickle.load(f)
    outliers = get_outlier_sessions(cache, ratio_threshold=5.0)
    bad_tracks = get_bad_tracks(cache)
    spanning = get_spanning_tracks()

    # Build list of (track_id, weeks_detected, first_week, duration)
    track_info = []
    for track in all_tracks:
        if track['n_sessions'] < MIN_TRACK_SESSIONS:
            continue
        aid = track['animal_id']
        tid = track['track_id']
        if tid in bad_tracks or tid not in spanning:
            continue
        weeks = []
        for session in track['sessions'].keys():
            if (tid, session) in outliers:
                continue
            wk = _rel_week_for_animal(session, aid)
            weeks.append(wk)
        weeks = sorted(set(weeks))
        if not weeks:
            continue
        duration = weeks[-1] - weeks[0]
        track_info.append((tid, weeks, weeks[0], duration))

    # Sort: longest duration first, then earliest start
    track_info.sort(key=lambda x: (-x[3], x[2]))

    # Plot
    for row_idx, (tid, weeks, first_wk, dur) in enumerate(track_info):
        # Horizontal line spanning the track
        ax.plot([min(weeks), max(weeks)], [row_idx, row_idx],
                color='#cccccc', linewidth=0.4, zorder=1)
        # Dots at each detected week
        ax.scatter(weeks, [row_idx] * len(weeks), s=4, color=PALETTE[0],
                   edgecolors='none', zorder=2)

    ax.set_xlabel('Weeks of training')
    ax.set_ylabel('Tracked unit')
    ax.set_xlim(-0.5, 4.5)
    ax.set_xticks(range(5))
    ax.set_ylim(-1, len(track_info))
    ax.set_yticks([0, len(track_info) - 1])
    ax.set_yticklabels([len(track_info), 1])
    ax.invert_yaxis()
    ax.set_title(f'Unit tracking (n = {len(track_info)} tracks)', fontsize=8)

    n_sessions_all = [len(t[1]) for t in track_info]

    print(f'  Total tracks: {len(track_info)}')
    print(f'  Sessions per track: median={np.median(n_sessions_all):.0f}, '
          f'range={min(n_sessions_all)}-{max(n_sessions_all)}')


# --- Panel C: Waveform distance (tracked vs nearby) ----------------------

def plot_panel_c(ax):
    """Waveform distance: tracked pairs vs nearby non-matched units."""
    with open(CACHE_PATH, 'rb') as f:
        cache = pickle.load(f)

    # Exclude outlier sessions and bad tracks (same filter as panel D)
    outliers = get_outlier_sessions(cache, ratio_threshold=5.0)
    bad_tracks = get_bad_tracks(cache)

    def _get_concat_waveform(info):
        tmpl = info['template_mean']
        t1, t2 = 0, min(181, tmpl.shape[0])
        return tmpl[t1:t2, :].flatten()

    # Within-track distances
    within_dists = []
    for tid, td in cache.items():
        if tid in bad_tracks:
            continue
        sessions = td['sessions']
        sessions_chrono = sorted(
            [s for s in sessions.keys() if sessions[s] is not None],
            key=lambda s: datetime.strptime(s, '%d-%b-%Y'))
        # Exclude outlier sessions
        sessions_chrono = [s for s in sessions_chrono if (tid, s) not in outliers]
        if len(sessions_chrono) < 2:
            continue

        wvfs = []
        for session in sessions_chrono:
            info = sessions[session]
            if info is None:
                wvfs.append(None)
                continue
            wvfs.append(_get_concat_waveform(info))

        for i in range(len(wvfs)):
            for j in range(i + 1, len(wvfs)):
                if wvfs[i] is not None and wvfs[j] is not None:
                    min_len = min(len(wvfs[i]), len(wvfs[j]))
                    d = np.linalg.norm(wvfs[i][:min_len] - wvfs[j][:min_len])
                    within_dists.append(d)

    # Between-track distances (nearby units, same session)
    session_units = {}
    for tid, td in cache.items():
        if tid in bad_tracks:
            continue
        aid = td['animal_id']
        for s, info in td['sessions'].items():
            if info is not None and (tid, s) not in outliers:
                session_units.setdefault((aid, s), []).append((tid, info))

    cross_dists = []
    ch_tolerance = 2
    for tid, td in cache.items():
        if tid in bad_tracks:
            continue
        aid = td['animal_id']
        for session in td['sessions']:
            if (tid, session) in outliers:
                continue
            info = td['sessions'][session]
            if info is None:
                continue
            ch = info['primary_ch']
            wvf_tracked = _get_concat_waveform(info)

            for other_tid, other_info in session_units.get((aid, session), []):
                if other_tid == tid:
                    continue
                if abs(other_info['primary_ch'] - ch) > ch_tolerance:
                    continue
                wvf_other = _get_concat_waveform(other_info)
                min_len = min(len(wvf_tracked), len(wvf_other))
                d = np.linalg.norm(wvf_tracked[:min_len] - wvf_other[:min_len])
                cross_dists.append(d)

    within_dists = np.array(within_dists)
    cross_dists = np.array(cross_dists)

    # Normalize relative to median tracked distance
    norm_factor = np.median(within_dists)
    within_norm = within_dists / norm_factor
    cross_norm = cross_dists / norm_factor

    # Histogram
    bins = np.linspace(0, np.percentile(cross_norm, 98), 50)
    ax.hist(within_norm, bins=bins, density=True, histtype='stepfilled',
            alpha=0.45, color=PALETTE[0], edgecolor=PALETTE[0], linewidth=0.8,
            label='Tracked', zorder=2)
    ax.hist(cross_norm, bins=bins, density=True, histtype='stepfilled',
            alpha=0.35, color=PALETTE[1], edgecolor=PALETTE[1], linewidth=0.8,
            label='Nearby units', zorder=1)
    ax.set_xlabel('Normalized waveform distance')
    ax.set_ylabel('Density')
    ax.set_title('Waveform distance', fontsize=8)
    ax.legend(fontsize=5)

    ratio = np.median(cross_norm)
    from scipy.stats import mannwhitneyu
    stat, p_val = mannwhitneyu(within_dists, cross_dists, alternative='less')
    print(f'  Tracked: n={len(within_dists)}, median={np.median(within_dists):.0f}')
    print(f'  Nearby: n={len(cross_dists)}, median={np.median(cross_dists):.0f}')
    print(f'  Ratio (nearby/tracked): {ratio:.1f}x')
    print(f'  Mann-Whitney U p={p_val:.2g}')


# --- Panel D: Modulation trends (tracked only) ---------------------------

def plot_panel_d(fig, gs_d, raw_df, tracked_set):
    """Modulation trends restricted to tracked units.

    1 row x 4 cols: PL z-score | NPL z-score | PL t2max | NPL t2max
    """
    df_pl = filter_pl(raw_df)
    df_npl = filter_npl(raw_df)

    pl_tracked = filter_to_tracked(df_pl, tracked_set)
    npl_tracked = filter_to_tracked(df_npl, tracked_set)
    pl_tracked = pl_tracked[pl_tracked['rel_week'] <= 4]
    npl_tracked = npl_tracked[npl_tracked['rel_week'] <= 4]

    print(f'  PL tracked rows: {len(pl_tracked)}, NPL tracked rows: {len(npl_tracked)}')

    gs_inner = gs_d.subgridspec(1, 4, wspace=0.4)

    panels = [
        (0, 0, pl_tracked, 'z_score', 'Modulation (z)', 'PL', 'two-sided', Z_CLIP),
        (0, 1, npl_tracked, 'z_score', '', 'NPL', 'two-sided', Z_CLIP),
        (0, 2, pl_tracked, 't_to_max_10ms_smoothed', 'Time to max\nfiring rate (ms)', 'PL', 'two-sided', TMAX_CLIP),
        (0, 3, npl_tracked, 't_to_max_10ms_smoothed', '', 'NPL', 'two-sided', TMAX_CLIP),
    ]

    for row, col, df, var, ylabel, title, alt, clip in panels:
        ax = fig.add_subplot(gs_inner[row, col])
        if len(df) < 10:
            ax.set_title(f'{title}\n(insufficient data)', fontsize=6)
            continue

        print(f'\n  {title}:')
        ax_out, stars, y_vals = plot_metric_over_weeks(
            df, var, ax, alt=alt, y_clip=clip)
        draw_early_late_brackets(stars, y_vals, ax)
        ax.set_ylabel(ylabel, fontsize=6)
        ax.set_title(title, fontsize=7, pad=15)


# --- Outlier filtering ---------------------------------------------------

def get_outlier_sessions(cache, ratio_threshold=5.0):
    """Find (track_id, session) pairs with waveform distance > ratio_threshold
    times the track median distance."""
    exclude = set()
    for tid, td in cache.items():
        sessions = td['sessions']
        sc = sorted([s for s in sessions if sessions[s] is not None],
                    key=lambda s: datetime.strptime(s, '%d-%b-%Y'))
        if len(sc) < 3:
            continue
        wvfs = {}
        for s in sc:
            info = sessions[s]
            tmpl = info['template_mean']
            t1, t2 = 0, min(181, tmpl.shape[0])
            wvfs[s] = tmpl[t1:t2, :].flatten()
        min_len = min(len(w) for w in wvfs.values())
        stack = np.array([w[:min_len] for w in wvfs.values()])
        median_wvf = np.median(stack, axis=0)
        dists = {s: np.linalg.norm(wvfs[s][:min_len] - median_wvf) for s in sc}
        med_dist = np.median(list(dists.values()))
        for s in sc:
            ratio = dists[s] / med_dist if med_dist > 0 else 0
            if ratio > ratio_threshold:
                exclude.add((tid, s))
    return exclude


TRACK_DIST_RATIO_THRESHOLD = 3.0


def get_bad_tracks(cache):
    """Return set of track_ids whose median pairwise waveform distance
    exceeds TRACK_DIST_RATIO_THRESHOLD times the global median."""
    def _concat(info):
        tmpl = info['template_mean']
        t1, t2 = 0, min(181, tmpl.shape[0])
        return tmpl[t1:t2, :].flatten()

    track_medians = []
    for tid, td in cache.items():
        sessions = td['sessions']
        sc = sorted([s for s in sessions if sessions[s] is not None],
                    key=lambda s: datetime.strptime(s, '%d-%b-%Y'))
        if len(sc) < 2:
            continue
        wvfs = [_concat(sessions[s]) for s in sc]
        dists = []
        for i in range(len(wvfs)):
            for j in range(i + 1, len(wvfs)):
                ml = min(len(wvfs[i]), len(wvfs[j]))
                dists.append(np.linalg.norm(wvfs[i][:ml] - wvfs[j][:ml]))
        track_medians.append((tid, np.median(dists)))

    global_median = np.median([m for _, m in track_medians])
    bad = set()
    for tid, med in track_medians:
        if med > global_median * TRACK_DIST_RATIO_THRESHOLD:
            bad.add(tid)
    print(f'  Excluding {len(bad)} bad tracks (>{TRACK_DIST_RATIO_THRESHOLD}x global median distance)')
    return bad


def get_spanning_tracks():
    """Return set of track_ids that span both early (weeks 0-1) and late (weeks 2-4)."""
    with open(TRACK_DIR / 'all_tracks.json') as f:
        all_tracks = json.load(f)
    spanning = set()
    for track in all_tracks:
        aid = track['animal_id']
        weeks = set()
        for session in track['sessions'].keys():
            weeks.add(_rel_week_for_animal(session, aid))
        has_early = bool(weeks & {0, 1})
        has_late = bool(weeks & {2, 3, 4})
        if has_early and has_late:
            spanning.add(track['track_id'])
    return spanning


def get_tracked_unit_ids_filtered():
    """Return tracked set with outlier sessions, bad tracks, and non-spanning tracks excluded."""
    with open(CACHE_PATH, 'rb') as f:
        cache = pickle.load(f)
    outliers = get_outlier_sessions(cache, ratio_threshold=5.0)
    bad_tracks = get_bad_tracks(cache)
    spanning = get_spanning_tracks()
    print(f'  Excluding {len(outliers)} outlier (track, session) pairs')
    print(f'  Spanning early+late: {len(spanning)} tracks')

    with open(TRACK_DIR / 'all_tracks.json') as f:
        all_tracks = json.load(f)
    tracked = set()
    for track in all_tracks:
        if track['n_sessions'] < MIN_TRACK_SESSIONS:
            continue
        tid = track['track_id']
        if tid in bad_tracks:
            continue
        if tid not in spanning:
            continue
        animal = track['animal_id']
        for session, uid in track['sessions'].items():
            if (tid, session) in outliers:
                continue
            tracked.add((animal, session, uid))
    return tracked


# --- Main ----------------------------------------------------------------

def main():
    print('Loading data...')
    raw_df = load_raw_df()
    tracked_set = get_tracked_unit_ids_filtered()
    print(f'Tracked (animal, session, unit) entries: {len(tracked_set)}')

    # Create figure: 4 rows
    # Row 1: A -- 3 example tracked waveforms
    # Row 2: B -- track raster  |  C -- waveform distance
    # Row 3: D -- modulation trends (tracked only)
    fig = plt.figure(figsize=(12, 8))
    gs_main = fig.add_gridspec(3, 1, height_ratios=[2.5, 3, 2],
                                hspace=0.5, left=0.06, right=0.94,
                                top=0.97, bottom=0.06)

    # Panel A: example tracked waveforms
    print('\nPanel A: Example tracked waveforms...')
    plot_panel_a(fig, gs_main[0])

    # Panels B + C side by side
    gs_bc = gs_main[1].subgridspec(1, 2, wspace=0.35)

    print('\nPanel B: Track raster...')
    ax_b = fig.add_subplot(gs_bc[0, 0])
    plot_panel_b(ax_b)

    print('\nPanel C: Waveform distance...')
    ax_c = fig.add_subplot(gs_bc[0, 1])
    plot_panel_c(ax_c)

    # Panel D: modulation trends
    print('\nPanel D: Modulation trends (tracked only)...')
    plot_panel_d(fig, gs_main[2], raw_df, tracked_set)

    # Panel labels
    for lbl, gs_pos, x_off in [('A', gs_main[0], 0.02),
                                 ('B', gs_bc[0, 0], 0.02),
                                 ('C', gs_bc[0, 1], 0.50),
                                 ('D', gs_main[2], 0.02)]:
        bbox = gs_pos.get_position(fig)
        fig.text(x_off, bbox.y1 + 0.01, lbl, fontsize=14,
                 fontweight='bold', va='bottom')

    # Save
    out_svg = FIG_DIR / 'figure_s7.svg'
    out_png = FIG_DIR / 'figure_s7.png'
    fig.savefig(out_svg, format='svg', bbox_inches='tight')
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\nSaved {out_svg}')
    print(f'Saved {out_png}')


if __name__ == '__main__':
    main()
