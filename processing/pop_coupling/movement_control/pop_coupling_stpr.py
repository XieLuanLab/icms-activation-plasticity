"""
STPR-based population coupling with configurable window.

Computes spike-triggered population rate (STPR) and normalized PC
using only spikes within a specified window from stim onset.

Saves per-animal pkl files in same format as v4 (stpr_dict, pc_dict, pc_norm_dict).

Usage:
    python -m pop_coupling.movement_control.pop_coupling_stpr --window 120
    python -m pop_coupling.movement_control.pop_coupling_stpr --window 700
"""
import argparse
import ast
import traceback
import pickle
import os
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

import pop_coupling.fig_utils as fig_utils
import batch_process.util.file_util as file_util
from pop_coupling.shared_plotting import PALETTE

SAVE_DIR = Path('pop_coupling')
SAVE_DIR.mkdir(exist_ok=True)

ANIMAL_IDS = ['ICMS92', 'ICMS93', 'ICMS98', 'ICMS100', 'ICMS101', 'ICMS83']
FS = 30000.0
MIN_SPIKES_STPR = 10  # min spikes for reliable STPR (lower for short windows)


def compute_stpr_and_pc(pop_raster, neuron_id, bin_size=0.001, win=0.1,
                        gauss_sigma=0.01, n_shuffles=200, rng=None):
    """Compute STPR and population coupling from concatenated spike times."""
    T = max((spikes.max() if len(spikes) else 0.0) for spikes in pop_raster.values())
    if T <= 0:
        return None, np.nan, np.nan, 0

    edges = np.arange(0.0, T + bin_size, bin_size)
    nbins = edges.size - 1

    pop_counts = np.zeros(nbins, dtype=float)
    for j, spikes in pop_raster.items():
        if j == neuron_id or len(spikes) == 0:
            continue
        pop_counts += np.histogram(spikes, bins=edges)[0]

    pop_rate = pop_counts / bin_size
    sigma_bins = max(1e-9, gauss_sigma / bin_size)
    pop_rate_smooth = gaussian_filter1d(pop_rate, sigma=sigma_bins, mode='nearest')
    pop_mc = pop_rate_smooth - np.mean(pop_rate_smooth)

    spikes_i = np.asarray(pop_raster.get(neuron_id, []), float)
    if spikes_i.size == 0:
        return None, np.nan, np.nan, 0

    win_bins = int(round(win / bin_size))
    spike_bins = np.floor(spikes_i / bin_size).astype(int)

    windows = []
    valid_spike_bins = []
    for sb in spike_bins:
        a, b = sb - win_bins, sb + win_bins
        if a < 0 or b >= nbins:
            continue
        windows.append(pop_mc[a:b + 1])
        valid_spike_bins.append(sb)

    if len(windows) == 0:
        return None, np.nan, np.nan, 0

    stpr = np.vstack(windows).mean(axis=0)
    pc = float(stpr[win_bins])
    n_used = len(windows)

    # Shuffle normalization
    pc_norm = np.nan
    if n_shuffles and n_used > 0:
        if rng is None:
            rng = np.random.default_rng()
        shuf_vals = []
        for _ in range(int(n_shuffles)):
            shift = rng.integers(low=win_bins + 1, high=max(win_bins + 2, nbins - win_bins - 1))
            pop_mc_sh = np.roll(pop_mc, shift)
            vals = [pop_mc_sh[sb - win_bins: sb + win_bins + 1][win_bins]
                    for sb in valid_spike_bins
                    if (sb - win_bins) >= 0 and (sb + win_bins) < nbins]
            if len(vals):
                shuf_vals.append(np.mean(vals))
        if len(shuf_vals):
            den = np.median(np.abs(shuf_vals))
            pc_norm = np.nan if den <= 1e-12 else float(pc / den)

    return stpr, pc, pc_norm, n_used


def concat_segments_rebased(spike_train, segments_s):
    """Concatenate spikes from segments, rebased to continuous time."""
    out = []
    offset = 0.0
    for a, b in segments_s:
        seg = spike_train[(spike_train >= a) & (spike_train < b)]
        if seg.size:
            out.append((seg - a) + offset)
        offset += (b - a)
    return np.concatenate(out) if out else np.array([], dtype=float)


def get_session_path(animal_id, session):
    if animal_id == 'ICMS83':
        return f'E:/ICMS83/{session}'
    root = "C:/data"
    if animal_id in ['ICMS92', 'ICMS93']:
        return os.path.join(root, animal_id, 'Behavior', session)
    return os.path.join(root, animal_id, session)


def get_windowed_segments(trial_df, window_ms, stim_trains=None):
    """Get stim segments restricted to window_ms from first pulse onset.

    Uses stim_timestamps from trial_df if available, otherwise falls back
    to stim_trains list (indexed by trial number).
    """
    window_samples = int(window_ms * FS / 1000)
    currents = trial_df['current'].unique()
    channels = trial_df['channel'].unique()
    currents = currents[(currents >= 4) & (currents < 7)]
    channels = channels[channels > 0]

    has_stim_ts = 'stim_timestamps' in trial_df.columns

    condition_segments = {}
    for current in currents:
        for channel in channels:
            sub_df = trial_df[(trial_df['current'] == current) &
                              (trial_df['channel'] == channel)]
            if sub_df.empty:
                continue
            segments = []
            for _, trial_row in sub_df.iterrows():
                onset = None

                if has_stim_ts:
                    ts = trial_row.get('stim_timestamps')
                    if ts is None:
                        continue
                    if isinstance(ts, str):
                        ts = ast.literal_eval(ts)
                    if isinstance(ts, float) and np.isnan(ts):
                        continue
                    if hasattr(ts, '__len__') and len(ts) > 0:
                        onset = ts[0]
                    elif np.isscalar(ts) and np.isfinite(ts):
                        onset = float(ts)
                elif stim_trains is not None:
                    # Use stim_trains indexed by trial number (1-indexed)
                    trial_num = int(trial_row['trial']) - 1
                    if 0 <= trial_num < len(stim_trains):
                        train = stim_trains[trial_num]
                        if len(train) > 0:
                            onset = train[0]

                if onset is not None:
                    segments.append((onset, onset + window_samples))

            if segments:
                condition_segments[(current, channel)] = segments
    return condition_segments


def lookup_units(df_mod, animal_id, session, ch, cur, is_pl):
    mask = (
        (df_mod['animal_id'] == animal_id) &
        (df_mod['session'] == session) &
        (df_mod['stim_ch'] == ch) &
        (df_mod['stim_current'] == cur) &
        (df_mod['is_pulse_locked'] == is_pl)
    )
    return df_mod.loc[mask, 'unit_id'].tolist()


def process_animal(animal_id, window_ms, df_mod):
    spike_pkl = f'unit_spike_train_data/{animal_id}_spike_data.pkl'
    if not Path(spike_pkl).exists():
        print(f'  [SKIP] {spike_pkl} not found')
        return

    with open(spike_pkl, 'rb') as f:
        spike_data = pickle.load(f)

    session_dict = OrderedDict()

    for session_name, session_data in spike_data.items():
        unit_spike_train_dict = session_data['unit_spike_train_dict']
        session_path = get_session_path(animal_id, session_name)

        # Load trial_df
        sr_pkl = Path(session_path) / 'batch_sort/session_responses_all.pkl'
        if not sr_pkl.exists():
            continue
        with open(sr_pkl, 'rb') as f:
            sr = pickle.load(f)
        trial_df = sr.trial_df

        stim_trains = session_data.get('stim_trains', [])
        condition_segments = get_windowed_segments(trial_df, window_ms, stim_trains=stim_trains)
        condition_results = {}

        for (cur, ch), segments in condition_segments.items():
            segments_s = [(a / FS, b / FS) for a, b in segments]

            pl_ids = lookup_units(df_mod, animal_id, session_name, ch, cur, True)
            npl_ids = lookup_units(df_mod, animal_id, session_name, ch, cur, False)
            mod_ids = sorted(set(pl_ids + npl_ids))

            if len(mod_ids) < 2:
                continue

            # Build pop raster: concatenate segments, rebase to continuous time
            pop_raster = {}
            for uid in mod_ids:
                if uid not in unit_spike_train_dict:
                    continue
                st = np.asarray(unit_spike_train_dict[uid], dtype=float) / FS
                pop_raster[uid] = concat_segments_rebased(st, segments_s)

            if len(pop_raster) < 2:
                continue

            pl_stpr_dict, pl_pc_dict, pl_pc_norm_dict = {}, {}, {}
            npl_stpr_dict, npl_pc_dict, npl_pc_norm_dict = {}, {}, {}

            rng = np.random.default_rng(0)
            for uid in pop_raster:
                stpr, pc, pc_norm, n_used = compute_stpr_and_pc(
                    pop_raster, uid, win=0.05, n_shuffles=200, rng=rng)

                if n_used < MIN_SPIKES_STPR:
                    stpr = None
                    pc = np.nan
                    pc_norm = np.nan

                if uid in pl_ids:
                    pl_stpr_dict[uid] = stpr
                    pl_pc_dict[uid] = pc
                    pl_pc_norm_dict[uid] = pc_norm
                elif uid in npl_ids:
                    npl_stpr_dict[uid] = stpr
                    npl_pc_dict[uid] = pc
                    npl_pc_norm_dict[uid] = pc_norm

            condition_results[(cur, ch)] = {
                "pl_stpr_dict": pl_stpr_dict,
                "pl_pc_dict": pl_pc_dict,
                "pl_pc_norm_dict": pl_pc_norm_dict,
                "npl_stpr_dict": npl_stpr_dict,
                "npl_pc_dict": npl_pc_dict,
                "npl_pc_norm_dict": npl_pc_norm_dict,
            }

        session_dict[session_name] = condition_results
        n_conds = sum(1 for v in condition_results.values()
                      if len(v.get('pl_pc_dict', {})) + len(v.get('npl_pc_dict', {})) > 0)
        print(f'  {session_name}: {n_conds} conditions')

    out_pkl = SAVE_DIR / f'{animal_id}_pop_coupling_stpr_{window_ms}ms.pkl'
    with open(out_pkl, 'wb') as f:
        pickle.dump(session_dict, f)
    print(f'  Saved {out_pkl}')


def aggregate_and_plot(window_ms, df_mod):
    """Load per-animal pkls, aggregate, and generate collapsed figure."""
    from pop_coupling.shared_plotting import apply_global_style
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    apply_global_style()

    pl_c = PALETTE[0]
    npl_c = PALETTE[1]

    def is_early(animal_id, session):
        sel = df_mod[(df_mod['session'] == session) & (df_mod['animal_id'] == animal_id)]
        vals = sel['rel_week'].dropna().unique()
        if len(vals) == 0: return 'nan'
        return 'early' if float(vals[0]) < 2 else 'late'

    def finite(arr):
        a = np.asarray(arr, float)
        return a[np.isfinite(a)]

    pl_early, pl_late = [], []
    npl_early, npl_late = [], []

    for animal_id in ANIMAL_IDS:
        pkl = SAVE_DIR / f'{animal_id}_pop_coupling_stpr_{window_ms}ms.pkl'
        if not pkl.exists():
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        for session, session_data in data.items():
            tp = is_early(animal_id, session)
            if tp == 'nan':
                continue
            for (cur, ch), cond in session_data.items():
                pls = [v for v in cond.get('pl_pc_norm_dict', {}).values()
                       if v is not None and np.isfinite(v)]
                npls = [v for v in cond.get('npl_pc_norm_dict', {}).values()
                        if v is not None and np.isfinite(v)]
                if tp == 'early':
                    pl_early.extend(pls)
                    npl_early.extend(npls)
                else:
                    pl_late.extend(pls)
                    npl_late.extend(npls)

    pl_e = finite(pl_early)
    pl_l = finite(pl_late)
    npl_e = finite(npl_early)
    npl_l = finite(npl_late)

    print(f'\nSTRP Pop coupling {window_ms}ms (collapsed 4-6 uA):')
    print(f'  PL early:  n={len(pl_e)}, median={np.median(pl_e):.4f}')
    print(f'  PL late:   n={len(pl_l)}, median={np.median(pl_l):.4f}')
    print(f'  NPL early: n={len(npl_e)}, median={np.median(npl_e):.4f}')
    print(f'  NPL late:  n={len(npl_l)}, median={np.median(npl_l):.4f}')

    groups = [pl_e, pl_l, npl_e, npl_l]
    labels = ['PL\nEarly', 'PL\nLate', 'NPL\nEarly', 'NPL\nLate']
    colors = [pl_c, pl_c, npl_c, npl_c]
    alphas = [0.5, 1.0, 0.5, 1.0]

    fig, ax = plt.subplots(figsize=(3, 2.5))
    for x, arr, lab, col, al in zip(range(4), groups, labels, colors, alphas):
        if len(arr) == 0:
            continue
        med = np.median(arr)
        iqr = np.percentile(arr, 75) - np.percentile(arr, 25)
        ax.errorbar(x, med, yerr=iqr, fmt='o', color=col, alpha=al,
                    capsize=3, markersize=4, elinewidth=0.5, capthick=0.5)

    tests = []
    if len(pl_e) > 0 and len(pl_l) > 0:
        _, p = mannwhitneyu(pl_e, pl_l, alternative='two-sided')
        tests.append({'x1': 0, 'x2': 1, 'p': p, 'label': 'PL early vs late'})
    if len(npl_e) > 0 and len(npl_l) > 0:
        _, p = mannwhitneyu(npl_e, npl_l, alternative='two-sided')
        tests.append({'x1': 2, 'x2': 3, 'p': p, 'label': 'NPL early vs late'})
    if len(pl_e) > 0 and len(npl_e) > 0:
        _, p = mannwhitneyu(pl_e, npl_e, alternative='two-sided')
        tests.append({'x1': 0, 'x2': 2, 'p': p, 'label': 'Early PL vs NPL'})
    if len(pl_l) > 0 and len(npl_l) > 0:
        _, p = mannwhitneyu(pl_l, npl_l, alternative='two-sided')
        tests.append({'x1': 1, 'x2': 3, 'p': p, 'label': 'Late PL vs NPL'})

    if tests:
        pvals = np.array([t['p'] for t in tests])
        _, pvals_corr, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')

        nonempty = [arr for arr in groups if len(arr) > 0]
        y_base = max(np.percentile(arr, 75) for arr in nonempty) + 0.5
        bracket_spacing = 1.0

        def p_to_stars(p):
            if p < 1e-3: return '***'
            if p < 1e-2: return '**'
            if p < 5e-2: return '*'
            return 'NS'

        for i, (t, pc) in enumerate(zip(tests, pvals_corr)):
            y = y_base + i * bracket_spacing
            stars = p_to_stars(pc)
            h = 0.2
            ax.plot([t['x1'], t['x1'], t['x2'], t['x2']], [y, y + h, y + h, y],
                    color='k', lw=0.5, clip_on=False)
            fontsize = 6 if stars == 'NS' else 8
            ax.text((t['x1'] + t['x2']) / 2, y + h + 0.1, stars,
                    ha='center', va='bottom', fontsize=fontsize)

        ax.set_ylim(bottom=ax.get_ylim()[0],
                    top=y_base + len(tests) * bracket_spacing + 1.5)

        print()
        for t, pc in zip(tests, pvals_corr):
            print(f'  {t["label"]}: p_raw={t["p"]:.2e}, p_fdr={pc:.2e} {p_to_stars(pc)}')

    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel('Pop coupling (normalized)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig_dir = Path(__file__).resolve().parents[2] / 'figures' / 'movement_control'
    fig_dir.mkdir(parents=True, exist_ok=True)
    save_path = fig_dir / f'pop_coupling_stpr_{window_ms}ms.svg'
    fig.savefig(save_path, format='svg', bbox_inches='tight')
    fig.savefig(str(save_path).replace('.svg', '.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\nSaved {save_path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--window', type=int, default=120)
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()

    raw_df = pd.read_pickle('raw_df_700ms_v3.pkl')
    df_mod = fig_utils.get_modulated_only_responses(raw_df)
    df_mod = df_mod.rename(columns={'stim_channel': 'stim_ch'})

    if not args.plot_only:
        for animal_id in ANIMAL_IDS:
            print(f'\n{"="*50}')
            print(f'  {animal_id} — {args.window}ms STPR')
            print(f'{"="*50}')
            try:
                process_animal(animal_id, args.window, df_mod)
            except Exception as e:
                print(f'  [ERROR] {e}')
                traceback.print_exc()

    aggregate_and_plot(args.window, df_mod)


if __name__ == '__main__':
    main()
