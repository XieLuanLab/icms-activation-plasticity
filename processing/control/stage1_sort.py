"""
Control Spike Sorting Pipeline - Stage 1 (Batch)
=================================================
Runs preprocessing, sorting, auto-curation, and generates PPTs for ALL sessions.
Skips sessions that are already sorted. Everything saved to disk.

After this completes, review the PPTs and run stage2_finalize.py.

Usage (from terminal for parallel processing):
    cd c:\\Users\\robin\\Desktop\\Python\\research\\icms-plasticity
    c:\\Users\\robin\\Desktop\\Python\\research\\si_v103_env\\Scripts\\python.exe pipeline\\stage1_sort.py
"""

import sys
import os
import shutil
import traceback
import time
from pathlib import Path

# Project root for imports — must be before project-relative imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import matplotlib
matplotlib.use("Agg")

import json  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import spikeinterface as si  # noqa: E402
import spikeinterface.sorters as ss  # noqa: E402
import spikeinterface.preprocessing as sp  # noqa: E402
from numcodecs import Blosc  # noqa: E402

from batch_process.util.plotting import plot_units_in_batches  # noqa: E402
from preprocessing.curate import curate_units  # noqa: E402
from preprocessing.icms_pipeline import preprocess_icms, stim_timestamps_seconds_to_samples  # noqa: E402
from util.load_control import (  # noqa: E402
    discover_sessions,
    extract_trial_parameters_from_rep,
)




def process_session(exp_dir, reps, job_kwargs, max_reps=None, cleanup_zarr=False, probe_json_path=None, force_resort=False, survey_current_uA=5.0):
    """
    Process a single session: preprocess → sort → curate → PPTs.

    Parameters
    ----------
    exp_dir : str
        Path to the session's ExperimentData folder.
    reps : list
        List of rep dicts from discover_sessions.
    job_kwargs : dict
        SpikeInterface job kwargs (n_jobs, chunk_duration, etc).
    max_reps : int or None
        Limit number of reps to load. None = all reps.
    cleanup_zarr : bool
        If True, delete preprocessed.zarr after sorting to save disk.
    probe_json_path : str or None
        Path to a custom probe layout JSON. If None, default is used.
    force_resort : bool
        If True, delete old results and re-sort even if already done.
    """
    save_dir = Path(exp_dir) / "spike_sorting"

    # --- Skip if already sorted (unless force) ---
    if not force_resort and (save_dir / "analyzer_raw.zarr").exists() and (save_dir / "curation_labels.npz").exists():
        print(f"  [SKIP] Already sorted: {save_dir}")
        return True

    # Clean old results if forcing re-sort
    if force_resort and save_dir.exists():
        for sub in ("mountainsort5_output", "analyzer_raw.zarr", "curation_labels.npz",
                     "curated_units.pptx", "rejected_units.pptx", "preprocessed.zarr",
                     "condensed_trials.csv"):
            p = save_dir / sub
            if p.exists():
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
        print(f"  [FORCE] Cleaned old results")

    save_dir.mkdir(exist_ok=True)

    if max_reps is not None:
        reps = reps[:max_reps]

    preproc_zarr = save_dir / "preprocessed.zarr"

    # --- RESUME: skip expensive load/condense/preprocess if zarr exists ---
    if preproc_zarr.exists():
        print("  [RESUME] preprocessed.zarr found — skipping to sorting step.")
    else:
        # 1. Extract trial parameters & Load recordings (Concatenated)
        import probeinterface as pi
        if probe_json_path is None:
            probe_json_path = os.path.join(PROJECT_ROOT, 'util', 'control_probe.json')
        probegroup = pi.read_probeinterface(probe_json_path)

        all_trials_dfs = []
        recordings = []
        cumulative_samples = 0
        fs = None

        MAX_NS5_SIZE_GB = 5.0  # Skip abnormally large NS5 files (likely corrupted)
        for rep in reps:
            # Check file size before loading
            ns5_size_gb = os.path.getsize(rep['ns5_path']) / (1024**3)
            if ns5_size_gb > MAX_NS5_SIZE_GB:
                print(f"  [!] SKIP Rep {rep['rep_num']}: NS5 is {ns5_size_gb:.1f} GB (max {MAX_NS5_SIZE_GB} GB)")
                continue
            print(f"  [+] Loading Rep {rep['rep_num']} (.ns5 and .nev)")
            try:
                # Load segment
                rec_rep = si.extractors.read_blackrock(rep['ns5_path'])
                rec_rep = rec_rep.set_probegroup(probegroup)
                if fs is None:
                    fs = rec_rep.get_sampling_frequency()
                
                # Extract trials
                df = extract_trial_parameters_from_rep(rep, survey_current_uA=survey_current_uA)
                
                # Shift timestamps by cumulative offset
                if cumulative_samples > 0:
                    offset_s = cumulative_samples / fs
                    df['stim_start_s'] = df['stim_start_s'].apply(lambda x: x + offset_s if x is not None else None)
                    df['stim_timestamps'] = df['stim_timestamps'].apply(lambda ts: [t + offset_s for t in ts])
                
                all_trials_dfs.append(df)
                recordings.append(rec_rep)
                cumulative_samples += rec_rep.get_num_samples(0)
            except Exception as e:
                print(f"      [!] Failed to load Rep {rep['rep_num']}: {e}")

        if not recordings:
            print("  [!] No recordings loaded, skipping.")
            return False

        session_df = pd.concat(all_trials_dfs, ignore_index=True)
        
        # First, concatenate all original recordings into one segment
        multi_rec = si.concatenate_recordings(recordings)
        
        # --- Option: Condensed Sorting ---
        # Instead of sorting 2 hours of mostly empty air, we extract a window 
        # around each trial and stitch them into a shorter, dense recording.
        PRE_WINDOW_S = 5.0  # seconds before stim start
        POST_WINDOW_S = 5.0 # seconds after stim start
        
        print(f"\n  [!] CONDENSING: Stitching {len(session_df)} trials into a dense recording...")
        
        trial_recordings = []
        condensed_trials = []
        current_condensed_start = 0.0
        
        # Sort session_df to be absolutely sure order is preserved
        session_df = session_df.sort_values('stim_start_s').reset_index(drop=True)

        for _, row in session_df.iterrows():
            t_start = row['stim_start_s'] - PRE_WINDOW_S
            t_end = row['stim_start_s'] + POST_WINDOW_S
            
            # Ensure we stay within recording bounds
            t_start = max(0.0, t_start)
            t_end = min(multi_rec.get_total_duration(), t_end)
            
            # Slice the recording
            slice_rec = multi_rec.frame_slice(start_frame=int(t_start*fs), end_frame=int(t_end*fs))
            trial_recordings.append(slice_rec)
            
            # Calculate new shifted timestamps for this trial in the condensed timeline
            condensed_offset = current_condensed_start - t_start
            
            new_row = row.copy()
            new_row['stim_start_s'] = row['stim_start_s'] + condensed_offset
            new_row['stim_timestamps'] = [t + condensed_offset for t in row['stim_timestamps']]
            condensed_trials.append(new_row)
            
            # Update condensed clock
            current_condensed_start += (t_end - t_start)

        # Re-stitch into a single dense recording
        condensed_rec = si.concatenate_recordings(trial_recordings)
        condensed_df = pd.DataFrame(condensed_trials)
        
        print(f"  Condensed Duration: {condensed_rec.get_total_duration():.2f}s (Saved {(multi_rec.get_total_duration() - condensed_rec.get_total_duration())/60:.1f} minutes)")
        
        # Overwrite originals so the rest of the pipeline uses the condensed data
        multi_rec = condensed_rec
        session_df = condensed_df
        
        # Save parameters for later verification/plotting
        session_df.to_csv(save_dir / "condensed_trials.csv", index=False)
        
        # 3. Preprocess
        print("\n  [1/6] Preprocessed Condensed Signal...")
        stim_samples = stim_timestamps_seconds_to_samples(session_df, fs=fs)
        # Controls were stimulated at 50 Hz (vs 100 Hz for the stim cohort).
        _, rec_wvf = preprocess_icms(multi_rec, stim_samples, stim_freq_hz=50)

        # 4. Materialize to zarr
        print("  [2/6] Saving preprocessed recording to zarr...")

        compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
        
        # Save directly to zarr using standard SI chunking!
        rec_disk = rec_wvf.save(folder=preproc_zarr, format="zarr", overwrite=True,
                                compressor=compressor, **job_kwargs)

    # 5. MountainSort5
    print("  [3/6] Running MountainSort5...")

    # --- Windows fix: release memory-mapped file handles before sorting ---
    # The zarr save above may leave mmap handles open via rec_wvf / multi_rec.
    # On Windows this causes PermissionError when the sorter tries to rmtree
    # its output folder.  Deleting refs + gc.collect() releases the handles.
    import gc
    # Clean up intermediate variables (may not exist on resume path)
    try:
        del rec_wvf, multi_rec, condensed_rec, trial_recordings, recordings
    except NameError:
        pass
    gc.collect()

    # Read the materialized recording from disk (no dangling mmap refs)
    rec_disk = si.read_zarr(preproc_zarr)

    ms5_output_folder = save_dir / "mountainsort5_output"
    ms5_params = {
        "scheme": "2",
        "detect_threshold": 5.5,
        "npca_per_channel": 3,
        "npca_per_subdivision": 10,
        "snippet_mask_radius": 0,
        "scheme2_detect_channel_radius": 200,
        "scheme2_training_duration_sec": 300,
        "filter": False,
        "whiten": True,
        "skip_alignment": True,
    }
    t_sort = time.time()

    # Suppress isosplit5 spam via Python-level stream filter.
    # (Do NOT redirect OS-level fds — that hides real errors.)
    import warnings
    import io
    warnings.filterwarnings("ignore", module=".*isosplit.*")
    warnings.filterwarnings("ignore", module=".*mountainsort.*")

    class _SpamFilter(io.TextIOBase):
        _SPAM = ("size did not change after splitting",
                 "new parcels having no points",
                 "warning: size did not change",
                 "warning: new parcels")
        def __init__(self, stream):
            self._stream = stream
        def write(self, msg):
            if msg and not any(s in msg.lower() for s in self._SPAM):
                return self._stream.write(msg)
            return 0
        def flush(self):
            self._stream.flush()
        def fileno(self):
            return self._stream.fileno()
        def __getattr__(self, name):
            return getattr(self._stream, name)

    _orig_stdout, _orig_stderr = sys.stdout, sys.stderr
    sys.stdout = _SpamFilter(_orig_stdout)
    sys.stderr = _SpamFilter(_orig_stderr)

    sorting = ss.run_sorter(
        sorter_name="mountainsort5",
        recording=sp.scale(rec_disk, dtype="float"),
        folder=ms5_output_folder,
        remove_existing_folder=True,
        verbose=False,
        **ms5_params,
    )

    # Restore Python-level streams
    sys.stdout = _orig_stdout
    sys.stderr = _orig_stderr
    warnings.resetwarnings()



    print(f"  Found {len(sorting.get_unit_ids())} raw units (sort took {time.time() - t_sort:.1f}s)")

    # 6. Create raw analyzer
    print("  [4/6] Creating SortingAnalyzer...")
    analyzer = si.create_sorting_analyzer(
        sorting=sorting,
        recording=rec_disk,
        format="zarr",
        folder=save_dir / "analyzer_raw.zarr",
        sparse=False,
        overwrite=True,
    )
    analyzer.compute(
        ["random_spikes", "waveforms", "templates", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )

    # 7. Auto-curate
    print("  [5/6] Auto-curating...")
    good_ids, bad_ids, labels = curate_units(analyzer, min_spikes=100)

    np.savez(
        save_dir / "curation_labels.npz",
        good_ids=good_ids,
        bad_ids=bad_ids,
        label_keys=list(labels.keys()),
        label_values=list(labels.values()),
    )

    # Auto-generate curation_config.json for user editing between stages
    config_path = save_dir / "curation_config.json"
    if not config_path.exists():
        config = {
            "rescue_ids": [],
            "remove_ids": [],
            "merge_pairs": [],
        }
        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

    # 8. Generate PPTs
    print("  [6/6] Generating PPTs...")
    sparse_raw = analyzer.copy()
    sparse_raw.sparsity = si.compute_sparsity(
        analyzer, method="radius", radius_um=60
    )
    sparse_raw.compute(
        ["random_spikes", "waveforms", "templates", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )

    # Good units
    plot_units_in_batches(
        sparse_raw, save_dir=save_dir, ppt_name="curated_units",
        unit_ids=good_ids
    )

    # Rejected units (color-coded)
    unit_colors = {}
    for uid in bad_ids:
        reason = labels.get(uid, "unknown")
        if reason == "artifact":
            unit_colors[uid] = "red"
        elif reason == "noise":
            unit_colors[uid] = "orange"
        elif reason == "low_count":
            unit_colors[uid] = "gray"
        else:
            unit_colors[uid] = "purple"

    plot_units_in_batches(
        sparse_raw, save_dir=save_dir, ppt_name="rejected_units",
        unit_ids=bad_ids, unit_colors=unit_colors
    )

    # 9. Keep preprocessed.zarr for stage2 (PPTs + sortingview need it)
    # Delete manually after stage2 is done, or set cleanup_zarr=True
    if cleanup_zarr and preproc_zarr.exists():
        shutil.rmtree(preproc_zarr)
        print(f"  Cleaned up: {preproc_zarr}")

    print(f"\n  [OK] Complete! Good: {len(good_ids)} | Rejected: {len(bad_ids)}")
    print(f"    PPTs in: {save_dir}")
    return True


def run_stage1(base_folder="E:/ICMS plasticity control/ICMS56", max_reps=None, probe_json_path=None, allow_missing_mat=True, force_resort=False, survey_current_uA=5.0):
    """
    Batch process all control animal sessions.

    Parameters
    ----------
    base_folder : str
        Root data directory.
    max_reps : int or None
        Max reps per session. None = use all reps.
    probe_json_path : str or None
        Path to custom probe layout JSON.
    allow_missing_mat : bool
        If True, synthesizes trials natively from .nev when .mat is missing.
    """
    # n_jobs=2: Limit parallelism to reduce memory pressure (OOM crash prevention)
    job_kwargs = dict(n_jobs=2, chunk_duration="1s", progress_bar=True)

    sessions = discover_sessions(base_folder, allow_missing_mat=allow_missing_mat)
    print(f"Found {len(sessions)} session folders\n")

    results = {}
    t_start = time.time()

    for i, (exp_dir, reps) in enumerate(sessions.items(), 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(sessions)}] {exp_dir}")
        print(f"  Reps: {len(reps)}")
        print(f"{'='*60}")

        try:
            success = process_session(exp_dir, reps, job_kwargs, max_reps=max_reps, probe_json_path=probe_json_path, force_resort=force_resort, survey_current_uA=survey_current_uA)
            results[exp_dir] = "OK" if success else "FAILED"
        except Exception as e:
            print(f"\n  [ERROR] {e}")
            traceback.print_exc()
            results[exp_dir] = f"ERROR: {e}"

    # Final summary
    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"BATCH COMPLETE — {len(sessions)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for exp_dir, status in results.items():
        icon = "[OK]" if status == "OK" else "[-]" if "SKIP" in status else "[X]"
        print(f"  {icon} {Path(exp_dir).parent.name}/{Path(exp_dir).name}: {status}")


if __name__ == "__main__":
    run_stage1(base_folder="E:/ICMS plasticity control/ICMS56", force_resort=False)
