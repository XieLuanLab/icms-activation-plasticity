"""
Run UnitMatch on all animals (ICMS83 + experimental).

Generates:
- Per-animal unitmatch results (probability matrix, UIDs, tracks)
- Combined tracking summary across all animals

Usage:
    python -m batch_process.postprocessing.run_unitmatch
    python -m batch_process.postprocessing.run_unitmatch ICMS92
"""
import sys
import ast
import gc
import math
import itertools
from pathlib import Path
from datetime import datetime
from collections import OrderedDict

import numpy as np
import pandas as pd
import dill as pickle
import spikeinterface as si
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

import UnitMatchPy.default_params as default_params
import UnitMatchPy.utils as util
import UnitMatchPy.overlord as ov
import UnitMatchPy.bayes_functions as bf

import batch_process.util.file_util as file_util

ANIMALS = {
    'ICMS83': {
        'root': Path('E:/ICMS83'),
        'sessions': ['20-Jul-2023', '24-Jul-2023', '26-Jul-2023', '28-Jul-2023',
                      '31-Jul-2023', '02-Aug-2023', '04-Aug-2023', '07-Aug-2023',
                      '10-Aug-2023'],
        'analyzer_path': 'batch_sort/stage3/analyzer_final.zarr',
    },
    'ICMS92': {
        'root': Path('C:/data/ICMS92'),
        'sessions': ['30-Aug-2023', '01-Sep-2023', '06-Sep-2023', '08-Sep-2023',
                      '12-Sep-2023', '14-Sep-2023', '19-Sep-2023', '21-Sep-2023',
                      '25-Sep-2023', '27-Sep-2023'],
        'subfolder': 'Behavior',
        'analyzer_path': 'batch_sort/merge/hmerge_analyzer.zarr',
    },
    'ICMS93': {
        'root': Path('C:/data/ICMS93'),
        'sessions': ['30-Aug-2023', '06-Sep-2023', '12-Sep-2023', '14-Sep-2023',
                      '20-Sep-2023', '22-Sep-2023', '26-Sep-2023', '29-Sep-2023',
                      '04-Oct-2023', '06-Oct-2023'],
        'subfolder': 'Behavior',
        'analyzer_path': 'batch_sort/merge/hmerge_analyzer.zarr',
    },
    'ICMS98': {
        'root': Path('C:/data/ICMS98'),
        'sessions': ['20-Oct-2023', '24-Oct-2023', '26-Oct-2023', '31-Oct-2023',
                      '02-Nov-2023', '07-Nov-2023', '14-Nov-2023', '17-Nov-2023',
                      '20-Nov-2023', '22-Nov-2023'],
        'analyzer_path': 'batch_sort/merge/hmerge_analyzer.zarr',
    },
    'ICMS100': {
        'root': Path('C:/data/ICMS100'),
        'sessions': ['26-Oct-2023', '02-Nov-2023', '03-Nov-2023', '09-Nov-2023',
                      '16-Nov-2023', '21-Nov-2023', '05-Dec-2023'],
        'analyzer_path': 'batch_sort/merge/hmerge_analyzer.zarr',
    },
    'ICMS101': {
        'root': Path('C:/data/ICMS101'),
        'sessions': ['27-Oct-2023', '01-Nov-2023', '03-Nov-2023', '16-Nov-2023',
                      '21-Nov-2023', '22-Nov-2023', '29-Nov-2023', '01-Dec-2023',
                      '05-Dec-2023'],
        'analyzer_path': 'batch_sort/merge/hmerge_analyzer.zarr',
    },
}

OUTPUT_ROOT = Path(__file__).resolve().parents[2] / 'unit_tracking'
MIN_CONSECUTIVE_P = 0.5


def get_session_path(animal_id, session_name):
    info = ANIMALS[animal_id]
    if 'subfolder' in info:
        return info['root'] / info['subfolder'] / session_name
    return info['root'] / session_name


def get_analyzer_path(animal_id, session_name):
    return get_session_path(animal_id, session_name) / ANIMALS[animal_id]['analyzer_path']


def get_max_channels(animal_id):
    max_ch = 0
    for session in ANIMALS[animal_id]['sessions']:
        ap = get_analyzer_path(animal_id, session)
        if ap.exists():
            a = si.load_sorting_analyzer(ap)
            max_ch = max(max_ch, a.get_num_channels())
            del a
    return max_ch


def extract_waveforms(animal_id, output_dir, n_channels_target):
    """Extract waveforms for all sessions of one animal."""
    session_dirs = []
    for session_name in ANIMALS[animal_id]['sessions']:
        ap = get_analyzer_path(animal_id, session_name)
        if not ap.exists():
            print(f"  [SKIP] No analyzer for {session_name}")
            continue

        print(f"  Loading {session_name}...")
        analyzer = si.load_sorting_analyzer(ap)

        if not analyzer.has_extension("waveforms"):
            analyzer.compute("random_spikes", method="all")
            analyzer.compute("waveforms")
            analyzer.compute("templates")

        waveforms_ext = analyzer.get_extension("waveforms")
        channel_locations = analyzer.get_channel_locations()
        n_ch = len(channel_locations)

        # Pad channels if needed
        if n_ch < n_channels_target:
            pad_n = n_channels_target - n_ch
            spacing = channel_locations[-1] - channel_locations[-2]
            extra = np.array([channel_locations[-1] + spacing * (i + 1) for i in range(pad_n)])
            channel_locations = np.vstack([channel_locations, extra])

        session_out = output_dir / session_name
        session_out.mkdir(parents=True, exist_ok=True)
        raw_dir = session_out / "RawWaveforms"
        raw_dir.mkdir(exist_ok=True)

        np.save(session_out / "channel_positions.npy", channel_locations)

        good_ids = []
        for uid in analyzer.unit_ids:
            wvf = waveforms_ext.get_waveforms_one_unit(uid)
            if len(wvf) < 2:
                continue
            half = len(wvf) // 2
            avg1 = np.mean(wvf[:half], axis=0)
            avg2 = np.mean(wvf[half:], axis=0)

            if n_ch < n_channels_target:
                pad = ((0, 0), (0, n_channels_target - n_ch))
                avg1 = np.pad(avg1, pad, mode='constant')
                avg2 = np.pad(avg2, pad, mode='constant')

            avg1 -= np.mean(avg1[:15], axis=0, keepdims=True)
            avg2 -= np.mean(avg2[:15], axis=0, keepdims=True)

            np.save(raw_dir / f"Unit{uid}_RawSpikes.npy",
                    np.stack([avg1, avg2], axis=-1))
            good_ids.append(uid)

        pd.DataFrame({'cluster_id': good_ids, 'group': ['good'] * len(good_ids)
                       }).to_csv(session_out / "cluster_group.tsv", sep='\t', index=False)

        print(f"  Saved {len(good_ids)} units ({n_ch} -> {n_channels_target} ch)")
        session_dirs.append(str(session_out))

        del analyzer, waveforms_ext
        gc.collect()

    return session_dirs


def run_unitmatch_animal(animal_id):
    """Run UnitMatch on one animal."""
    output_dir = OUTPUT_ROOT / animal_id
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  {animal_id}")
    print(f"{'='*60}")

    # Extract waveforms
    print("  Getting max channels...")
    n_ch_target = get_max_channels(animal_id)
    print(f"  Target channels: {n_ch_target}")

    session_dirs = extract_waveforms(animal_id, output_dir, n_ch_target)
    if len(session_dirs) < 2:
        print(f"  [SKIP] Need >= 2 sessions")
        return None

    # Run UnitMatch
    param = default_params.get_default_param()
    param['KS_dirs'] = session_dirs

    sample = np.load(list((Path(session_dirs[0]) / 'RawWaveforms').iterdir())[0])
    param['SpikeWidth'] = sample.shape[0]
    param['waveidx'] = np.arange(sample.shape[0])
    param['PeakLoc'] = sample.shape[0] // 2

    wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(session_dirs)
    param = util.get_probe_geometry(channel_pos[0], param)

    waveform, session_id, session_switch, within_session, good_units, param = \
        util.load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only=True)

    flat_unit_ids = np.concatenate([gu.flatten() for gu in good_units]).astype(int)

    clus_info = pd.DataFrame({
        'session_id': session_id,
        'original_ids': flat_unit_ids,
    })

    print(f"  Units: {param['n_units']}, Sessions: {len(session_dirs)}")

    # Extract parameters and scores
    extracted = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    total_score, candidate_pairs, scores_to_include, predictors = \
        ov.extract_metric_scores(extracted, session_switch, within_session, param, niter=2)

    # Bayesian matching
    cond = np.unique(candidate_pairs)
    priors = np.array([0.5, 0.5])
    parameter_kernels = bf.get_parameter_kernels(
        scores_to_include, candidate_pairs, cond, param, add_one=1)
    probability = bf.apply_naive_bayes(
        parameter_kernels, priors, predictors, param, cond)
    prob_matrix = probability[:, 1].reshape(param['n_units'], param['n_units'])

    # Evaluate
    match_threshold = 0.75
    util.evaluate_output(prob_matrix, param, within_session, session_switch,
                         match_threshold=match_threshold)

    # Build tracks
    n_units = len(session_id)
    cross_mask = np.zeros((n_units, n_units), dtype=bool)
    for i in range(n_units):
        for j in range(n_units):
            if session_id[i] != session_id[j]:
                cross_mask[i, j] = True

    matches = (prob_matrix > match_threshold) & cross_mask
    graph = csr_matrix(matches.astype(int))
    _, labels = connected_components(graph, directed=False)

    # Deduplicate and trim
    sessions_list = ANIMALS[animal_id]['sessions']
    tracks = _build_clean_tracks(labels, session_id, flat_unit_ids, prob_matrix,
                                  sessions_list)

    # Save
    np.savez(output_dir / 'unitmatch_results.npz',
             output_prob_matrix=prob_matrix,
             UIDs=labels,
             session_id=session_id,
             flat_unit_ids=flat_unit_ids,
             match_threshold=match_threshold)

    # Save tracks as JSON-friendly dict
    import json
    tracks_serializable = {}
    for tid, tinfo in tracks.items():
        tracks_serializable[str(tid)] = {s: int(u) for s, u in tinfo.items()}
    with open(output_dir / 'tracks.json', 'w') as f:
        json.dump(tracks_serializable, f, indent=2)

    print(f"\n  Tracked >= 2 sessions: {len(tracks)}")
    for tid, smap in sorted(tracks.items(), key=lambda x: -len(x[1])):
        print(f"    Track {tid}: {len(smap)} sessions")

    return {
        'animal_id': animal_id,
        'prob_matrix': prob_matrix,
        'UIDs': labels,
        'session_id': session_id,
        'flat_unit_ids': flat_unit_ids,
        'tracks': tracks,
    }


def _build_clean_tracks(UIDs, session_id, flat_unit_ids, prob_matrix, sessions_list):
    """Build deduplicated and trimmed tracks."""
    raw_tracks = {}
    for uid in np.unique(UIDs):
        mask = UIDs == uid
        indices = np.where(mask)[0]
        sessions = session_id[mask]
        units = flat_unit_ids[mask]
        unique_sessions = np.unique(sessions)
        if len(unique_sessions) < 2:
            continue

        clean_map = {}
        clean_indices = {}
        for s in unique_sessions:
            s_mask = sessions == s
            s_indices = indices[s_mask]
            s_units = units[s_mask]
            if len(s_indices) == 1:
                clean_map[s] = int(s_units[0])
                clean_indices[s] = int(s_indices[0])
            else:
                best_score = -1
                for idx, u in zip(s_indices, s_units):
                    other_indices = indices[sessions != s]
                    avg_p = np.mean(prob_matrix[idx, other_indices]) if len(other_indices) > 0 else 0
                    if avg_p > best_score:
                        best_score = avg_p
                        clean_map[s] = int(u)
                        clean_indices[s] = int(idx)
        raw_tracks[uid] = {'map': clean_map, 'indices': clean_indices}

    # Trim weak links
    trimmed = {}
    for uid, info in raw_tracks.items():
        smap = info['map']
        sidx = info['indices']
        sorted_s = sorted(smap.keys())
        chains = []
        current_chain = [sorted_s[0]]
        for i in range(len(sorted_s) - 1):
            p = prob_matrix[sidx[sorted_s[i]], sidx[sorted_s[i + 1]]]
            if p >= MIN_CONSECUTIVE_P:
                current_chain.append(sorted_s[i + 1])
            else:
                chains.append(current_chain)
                current_chain = [sorted_s[i + 1]]
        chains.append(current_chain)
        best = max(chains, key=len)
        if len(best) >= 2:
            trimmed[uid] = {sessions_list[s]: smap[s] for s in best}

    return trimmed


def _get_rel_week(session_name, first_session_name):
    """Compute relative week using ceil(days / 7), week 0 for first session."""
    fmt = '%d-%b-%Y'
    d = datetime.strptime(session_name, fmt)
    d0 = datetime.strptime(first_session_name, fmt)
    days = (d - d0).days
    return math.ceil(days / 7) if days > 0 else 0


def filter_tracks_by_week(tracks, sessions_list, max_week=4):
    """Filter track session maps to only include sessions within [0, max_week].

    Drops tracks that end up with fewer than 2 sessions after filtering.
    """
    first_session = sessions_list[0]
    filtered = {}
    for tid, session_map in tracks.items():
        trimmed = {s: u for s, u in session_map.items()
                   if _get_rel_week(s, first_session) <= max_week}
        if len(trimmed) >= 2:
            filtered[tid] = trimmed
    return filtered


def combine_all_tracks(max_week=4):
    """Load all per-animal tracks, filter to max_week, and combine into one summary."""
    import json
    all_tracks = []

    for animal_id in ANIMALS:
        tracks_path = OUTPUT_ROOT / animal_id / 'tracks.json'
        if not tracks_path.exists():
            continue
        with open(tracks_path) as f:
            tracks = json.load(f)

        if max_week is not None:
            sessions_list = ANIMALS[animal_id]['sessions']
            tracks = filter_tracks_by_week(tracks, sessions_list, max_week=max_week)

        for track_id, session_map in tracks.items():
            all_tracks.append({
                'animal_id': animal_id,
                'track_id': f"{animal_id}_{track_id}",
                'n_sessions': len(session_map),
                'sessions': session_map,
            })

    # Summary
    df = pd.DataFrame([{
        'animal_id': t['animal_id'],
        'track_id': t['track_id'],
        'n_sessions': t['n_sessions'],
    } for t in all_tracks])

    if len(df) > 0:
        df.to_csv(OUTPUT_ROOT / 'all_tracks_summary.csv', index=False)

        print(f"\n{'='*60}")
        print("COMBINED TRACKING SUMMARY")
        print(f"{'='*60}")
        print(f"Total tracked units: {len(df)}")
        print(f"Tracked >= 3 sessions: {(df['n_sessions'] >= 3).sum()}")
        print(f"Mean track length: {df['n_sessions'].mean():.1f}")
        print(f"Max track length: {df['n_sessions'].max()}")
        print(f"\nPer animal:")
        for aid in sorted(df['animal_id'].unique()):
            a = df[df['animal_id'] == aid]
            print(f"  {aid}: {len(a)} tracks, lengths: {sorted(a['n_sessions'].values, reverse=True)}")

    # Save full tracks dict
    with open(OUTPUT_ROOT / 'all_tracks.json', 'w') as f:
        json.dump(all_tracks, f, indent=2)

    return all_tracks


if __name__ == '__main__':
    animals_to_run = list(ANIMALS.keys())
    if len(sys.argv) > 1:
        animals_to_run = [a for a in sys.argv[1:] if a in ANIMALS]

    for animal_id in animals_to_run:
        try:
            run_unitmatch_animal(animal_id)
        except Exception as e:
            print(f"  [ERROR] {animal_id}: {e}")
            import traceback
            traceback.print_exc()

    combine_all_tracks()
