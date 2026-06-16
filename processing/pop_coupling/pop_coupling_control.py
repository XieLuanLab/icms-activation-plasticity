"""
Population coupling analysis for control animals.

Adapted from pop_coupling_v2.py for experimental animals.
Computes STPR (spike-triggered population rate) and shuffle-normalized
population coupling (pc_norm) — matching the experimental pipeline output.

Uses session_responses.pkl for PL/NPL classification and
condensed_trials.csv for trial windows (to segment spikes).

Output: pop_coupling/{animal}_pop_coupling_control.pkl
"""

import sys
import os
import ast
import traceback
from pathlib import Path
from collections import OrderedDict

# Project root for imports (needed to unpickle session_responses.pkl)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import dill as pickle

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


# ─── Config ───────────────────────────────────────────────────────────
CONTROL_ROOT = Path("E:/ICMS plasticity control")
ANIMAL_IDS = ["ICMS43", "ICMS45", "ICMS48", "ICMS54", "ICMS56"]
FS = 30000.0
OUT_DIR = Path("pop_coupling")
MIN_SPIKES = 50  # minimum spikes for valid STPR


# ─── STPR + PC computation (from pop_coupling_v2.py) ─────────────────

def compute_stpr_and_pc(pop_raster, neuron_id, bin_size=0.001, win=0.10,
                        gauss_sigma=0.01, n_shuffles=200, rng=None):
    """
    Spike-triggered population rate and population coupling.

    Returns: lags (s), stpr (array), pc (float), n_spikes_used (int), pc_norm (float or None)
        pc_norm is pc normalized by median(|pc_shuffled|) if n_shuffles>0.
    """
    T = max((spikes.max() if len(spikes) else 0.0)
            for spikes in pop_raster.values())
    if T <= 0:
        return None, None, np.nan, 0, None

    edges = np.arange(0.0, T + bin_size, bin_size)
    nbins = edges.size - 1

    # Population rate (exclude target neuron)
    pop_counts = np.zeros(nbins, dtype=float)
    for j, spikes in pop_raster.items():
        if j == neuron_id or len(spikes) == 0:
            continue
        pop_counts += np.histogram(spikes, bins=edges)[0]
    pop_rate = pop_counts / bin_size
    sigma_bins = max(1e-9, gauss_sigma / bin_size)
    pop_rate_smooth = gaussian_filter1d(pop_rate, sigma=sigma_bins, mode='nearest')

    # Mean-center
    pop_mc = pop_rate_smooth - np.mean(pop_rate_smooth)

    # Spike-triggered average around neuron i's spikes
    spikes_i = np.asarray(pop_raster.get(neuron_id, []), float)
    win_bins = int(round(win / bin_size))
    lags = np.arange(-win_bins, win_bins + 1) * bin_size

    if spikes_i.size == 0:
        return lags, np.full(lags.size, np.nan), np.nan, 0, None

    spike_bins = np.floor(spikes_i / bin_size).astype(int)

    windows = []
    for sb in spike_bins:
        a, b = sb - win_bins, sb + win_bins
        if a < 0 or b >= nbins:
            continue
        windows.append(pop_mc[a:b + 1])

    if len(windows) == 0:
        return lags, np.full(lags.size, np.nan), np.nan, 0, None

    stpr = np.vstack(windows).mean(axis=0)
    pc = float(stpr[win_bins])  # value at zero lag
    n_used = len(windows)

    # Shuffle normalization
    pc_norm = None
    if n_shuffles and n_used > 0:
        if rng is None:
            rng = np.random.default_rng()
        shuf_vals = []
        for _ in range(int(n_shuffles)):
            shift = rng.integers(low=win_bins + 1, high=nbins - win_bins - 1)
            pop_mc_sh = np.roll(pop_mc, shift)
            vals = [pop_mc_sh[sb - win_bins: sb + win_bins + 1][win_bins]
                    for sb in spike_bins
                    if (sb - win_bins) >= 0 and (sb + win_bins) < nbins]
            if len(vals):
                shuf_vals.append(np.mean(vals))
        if len(shuf_vals):
            den = np.median(np.abs(shuf_vals))
            pc_norm = np.nan if den <= 1e-12 else float(pc / den)

    return lags, stpr, pc, n_used, pc_norm


# ─── Helper: concatenate spike segments rebased ──────────────────────

def concat_segments_rebased(spike_train, segments):
    """Concatenate spikes from segments, rebasing each to start at 0."""
    out = []
    offset = 0.0
    for a, b in segments:
        seg = spike_train[(spike_train >= a) & (spike_train < b)]
        if seg.size:
            out.append((seg - a) + offset)
        offset += (b - a)
    return np.concatenate(out) if out else np.array([], dtype=float)


# ─── Control-specific helpers ─────────────────────────────────────────

def discover_sessions(animal_dir):
    """Find spike_sorting dirs with session_responses.pkl + condensed_trials.csv."""
    sessions = []
    for spike_dir in sorted(animal_dir.rglob("spike_sorting")):
        sr_path = spike_dir / "session_responses.pkl"
        csv_path = spike_dir / "condensed_trials.csv"
        analyzer_path = spike_dir / "analyzer_final.zarr"
        if sr_path.exists() and csv_path.exists() and analyzer_path.exists():
            sessions.append(spike_dir)
    return sessions


def get_session_date(spike_sorting_dir):
    """Extract session date string from path."""
    exp_dir = spike_sorting_dir.parent
    return exp_dir.parent.name


def parse_stim_timestamps(ts_str):
    """Parse stim_timestamps string from CSV."""
    if isinstance(ts_str, str):
        try:
            return ast.literal_eval(ts_str)
        except (ValueError, SyntaxError):
            return []
    if isinstance(ts_str, (list, np.ndarray)):
        return list(ts_str)
    return []


def get_trial_windows_per_condition(csv_path):
    """
    Build trial windows from condensed_trials.csv.

    Returns: dict[(current, channel)] -> list of (start_s, end_s) tuples
    """
    df = pd.read_csv(csv_path)
    condition_windows = {}

    for (cur, ch), sub_df in df.groupby(['current_uA', 'channel']):
        segments = []
        for _, row in sub_df.iterrows():
            ts = parse_stim_timestamps(row['stim_timestamps'])
            if not ts or len(ts) == 0:
                continue
            segments.append((ts[0], ts[-1]))

        if segments:
            condition_windows[(int(cur), int(ch))] = segments

    return condition_windows


def get_pl_npl_unit_ids(session_responses, channel, current):
    """Get PL and NPL unit IDs for a given (channel, current)."""
    pl_ids = []
    npl_ids = []

    for uid, ur in session_responses.unit_responses.items():
        key = (channel, current)
        if key not in ur.stim_responses:
            continue
        scr = ur.stim_responses[key]
        if scr.pulse_response.is_pulse_locked:
            pl_ids.append(int(uid))
        else:
            npl_ids.append(int(uid))

    return pl_ids, npl_ids


def get_pop_coupling_session(spike_sorting_dir):
    """
    Compute STPR-based population coupling for PL vs NPL units.

    Matches experimental pop_coupling_v2.py output format:
        condition_results[(cur, ch)] = {
            "pl_stpr_dict":     {unit_id: stpr_array},
            "pl_pc_dict":       {unit_id: pc_hz},
            "pl_pc_norm_dict":  {unit_id: pc_normalized},
            "npl_stpr_dict":    {unit_id: stpr_array},
            "npl_pc_dict":      {unit_id: pc_hz},
            "npl_pc_norm_dict": {unit_id: pc_normalized},
        }
    """
    from spikeinterface import full as si

    # Load session responses for PL/NPL labels
    sr_path = spike_sorting_dir / "session_responses.pkl"
    with open(sr_path, 'rb') as f:
        session_responses = pickle.load(f)

    # Load spike trains from analyzer
    analyzer_path = spike_sorting_dir / "analyzer_final.zarr"
    analyzer = si.load_sorting_analyzer(str(analyzer_path))
    unit_spike_train_dict = {
        int(uid): analyzer.sorting.get_unit_spike_train(uid)
        for uid in analyzer.unit_ids
    }

    # Get trial windows per condition from CSV
    csv_path = spike_sorting_dir / "condensed_trials.csv"
    condition_windows = get_trial_windows_per_condition(csv_path)

    condition_results = {}

    for (cur, ch), segments in condition_windows.items():
        if len(segments) < 3:
            continue

        # Get PL/NPL unit IDs for this condition
        pl_unit_ids, npl_unit_ids = get_pl_npl_unit_ids(session_responses, ch, cur)
        mod_ids = sorted(set(pl_unit_ids + npl_unit_ids))

        if len(mod_ids) < 2:
            continue

        # Build pop_raster: concatenate stim segments rebased
        segments_s = [(float(a), float(b)) for a, b in segments]
        pop_raster = {}
        for uid in mod_ids:
            if uid not in unit_spike_train_dict:
                continue
            st = np.asarray(unit_spike_train_dict[uid]).astype(float) / FS
            pop_raster[uid] = concat_segments_rebased(st, segments_s)

        if len(pop_raster) < 2:
            continue

        pl_stpr_dict, pl_pc_dict, pl_pc_norm_dict = {}, {}, {}
        npl_stpr_dict, npl_pc_dict, npl_pc_norm_dict = {}, {}, {}

        rng = np.random.default_rng(0)
        for uid in pop_raster:
            lags, stpr, pc, n_used, pc_norm = compute_stpr_and_pc(
                pop_raster, uid, win=0.1, n_shuffles=200, rng=rng
            )

            if n_used < MIN_SPIKES:
                stpr = None
                pc = np.nan
                pc_norm = np.nan

            if uid in pl_unit_ids:
                pl_stpr_dict[uid] = stpr
                pl_pc_dict[uid] = pc
                pl_pc_norm_dict[uid] = pc_norm
            elif uid in npl_unit_ids:
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

    return condition_results


# ─── Main ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)

    for animal_id in ANIMAL_IDS:
        print(f"\nProcessing {animal_id}...")
        animal_dir = CONTROL_ROOT / animal_id
        sessions = discover_sessions(animal_dir)

        if not sessions:
            print(f"  No sessions found")
            continue

        session_dict = OrderedDict()

        for spike_dir in sessions:
            session_date = get_session_date(spike_dir)
            print(f"  {session_date}...", end=" ", flush=True)

            try:
                condition_results = get_pop_coupling_session(spike_dir)
                session_dict[session_date] = condition_results
                n_conds = len(condition_results)
                n_pl = sum(
                    len(v["pl_pc_dict"]) for v in condition_results.values()
                )
                n_npl = sum(
                    len(v["npl_pc_dict"]) for v in condition_results.values()
                )
                print(f"{n_conds} conditions, {n_pl} PL, {n_npl} NPL responses")
            except Exception as e:
                print(f"[WARN] {e}")
                traceback.print_exc()

        out_pkl = OUT_DIR / f"{animal_id}_pop_coupling_control.pkl"
        with open(out_pkl, "wb") as f:
            pickle.dump(session_dict, f)
        print(f"  -> Saved {out_pkl}")
