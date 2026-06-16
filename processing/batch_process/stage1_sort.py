"""
Stage 1: preprocess + spike sort (MountainSort5) + auto-curate, per session.

Loads via the session-loader dispatch (standard DataLoader for Mouse 1-5, the
Mouse 6 loader for ICMS83), runs preprocess_icms, sorts, and auto-curates.
Output: <session>/batch_sort/stage1/

Usage:
    python -m batch_process.stage1_sort <session_folder> [<session_folder> ...]
"""

import sys
import os
import gc
import shutil
import time
import json
import traceback
from pathlib import Path

import numpy as np
import spikeinterface as si
import spikeinterface.sorters as ss
import spikeinterface.preprocessing as sp
from numcodecs import Blosc

from batch_process.util.dataloaders import load_session
from batch_process.util.plotting import plot_units_in_batches
from preprocessing.icms_pipeline import preprocess_icms
from preprocessing.curate import curate_units




def process_session(data_folder, job_kwargs, force_resort=False):
    """Run stage 1 on a single session."""
    data_folder = Path(data_folder)
    session_name = data_folder.name

    save_dir = data_folder / "batch_sort" / "stage1"

    # Skip if already done
    if not force_resort and (save_dir / "stage1_analyzer_curated_sparse.zarr").exists():
        print(f"  [SKIP] Already sorted: {save_dir}")
        return True

    # Clean old results if forcing re-sort
    if force_resort and save_dir.exists():
        for sub in ("mountainsort5_output", "stage1_analyzer_raw.zarr",
                     "stage1_analyzer_raw_sparse.zarr", "stage1_analyzer_curated",
                     "stage1_analyzer_curated.zarr", "stage1_analyzer_curated_sparse.zarr",
                     "curation_labels.npz",
                     "curated_units.pptx", "rejected_units.pptx", "preprocessed.zarr"):
            p = save_dir / sub
            if p.exists():
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
        print(f"  [FORCE] Cleaned old results")

    save_dir.mkdir(parents=True, exist_ok=True)

    # ─── Load data ────────────────────────────────────────
    print("  Loading session data...")
    dataloader = load_session(str(data_folder), get_recording=True)
    rec = dataloader.recording

    # Save trial_df
    dataloader.trial_df.to_csv(
        data_folder / "batch_sort" / "trial_df.csv", index=False)

    n_trials = len(dataloader.trial_df)
    channels = sorted(dataloader.trial_df['channel'].unique())
    print(f"  {n_trials} trials, channels={channels}, "
          f"duration={rec.get_num_samples() / dataloader.fs:.0f}s")

    # Stim timestamps are already in samples from the adapter
    all_stim_samples = [[int(ts) for ts in trial_ts]
                        for trial_ts in dataloader.all_stim_timestamps
                        if len(trial_ts) > 0]

    if not all_stim_samples:
        print("  [SKIP] No stim timestamps")
        return False

    # ─── Preprocess ───────────────────────────────────────
    preproc_zarr = save_dir / "preprocessed.zarr"

    if preproc_zarr.exists():
        print("  [RESUME] preprocessed.zarr found — skipping to sorting.")
    else:
        print("  Preprocessing...")
        _, rec_wvf = preprocess_icms(rec, all_stim_samples)

        print("  Saving preprocessed to zarr...")
        compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
        rec_wvf.save(folder=preproc_zarr, format="zarr",
                     overwrite=True, compressor=compressor, **job_kwargs)

        # Release memory
        del rec_wvf, rec
        gc.collect()

    # ─── Sort ─────────────────────────────────────────────
    print("  Running MountainSort5...")
    rec_disk = si.read_zarr(preproc_zarr)

    ms5_output = save_dir / "mountainsort5_output"
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
    }

    t_sort = time.time()
    sorting = ss.run_sorter(
        sorter_name="mountainsort5",
        recording=sp.scale(rec_disk, dtype="float"),
        folder=ms5_output,
        remove_existing_folder=True,
        verbose=True,
        **ms5_params,
    )
    print(f"  Found {len(sorting.get_unit_ids())} raw units "
          f"({time.time() - t_sort:.0f}s)")

    # ─── Raw analyzer ─────────────────────────────────────
    print("  Creating raw analyzer...")
    analyzer = si.create_sorting_analyzer(
        sorting=sorting,
        recording=rec_disk,
        format="zarr",
        folder=save_dir / "stage1_analyzer_raw.zarr",
        sparse=False,
        overwrite=True,
    )
    analyzer.compute(
        ["random_spikes", "waveforms", "templates", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )

    # ─── Auto-curate ──────────────────────────────────────
    print("  Auto-curating...")
    good_ids, bad_ids, labels = curate_units(analyzer, min_spikes=100)
    print(f"  {len(good_ids)} good, {len(bad_ids)} bad")

    np.savez(
        save_dir / "curation_labels.npz",
        good_ids=good_ids,
        bad_ids=bad_ids,
        label_keys=list(labels.keys()),
        label_values=list(labels.values()),
    )

    config_path = save_dir / "curation_config.json"
    if not config_path.exists():
        with open(config_path, "w") as f:
            json.dump({"rescue_ids": [], "remove_ids": []}, f, indent=4)

    # ─── PPTs ─────────────────────────────────────────────
    print("  Generating PPTs...")
    sparse_raw_path = save_dir / "stage1_analyzer_raw_sparse.zarr"
    if sparse_raw_path.exists():
        shutil.rmtree(sparse_raw_path)
    sparse_raw = analyzer.save_as(
        folder=sparse_raw_path, format="zarr")
    sparse_raw.sparsity = si.compute_sparsity(
        analyzer, method="radius", radius_um=60)
    sparse_raw.compute(
        ["random_spikes", "waveforms", "templates", "template_similarity", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )

    if len(good_ids) > 0:
        plot_units_in_batches(
            sparse_raw, save_dir=save_dir, ppt_name="curated_units",
            unit_ids=good_ids)

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

    if len(bad_ids) > 0:
        plot_units_in_batches(
            sparse_raw, save_dir=save_dir, ppt_name="rejected_units",
            unit_ids=bad_ids, unit_colors=unit_colors)

    # ─── Save curated analyzer ────────────────────────────
    if len(good_ids) == 0:
        print("  [WARN] No good units — skipping curated analyzer.")
        # Still create the marker file so we don't re-sort
        sparse_curated_path = save_dir / "stage1_analyzer_curated_sparse.zarr"
        sparse_curated_path.mkdir(parents=True, exist_ok=True)
        (sparse_curated_path / "_NO_GOOD_UNITS").touch()
    else:
        curated_path = save_dir / "stage1_analyzer_curated"
        if curated_path.exists():
            shutil.rmtree(curated_path)
        curated_analyzer = analyzer.select_units(
            unit_ids=good_ids, format="zarr",
            folder=curated_path)
        curated_analyzer.compute(
            ["random_spikes", "waveforms", "templates",
             "template_similarity", "correlograms",
             "spike_amplitudes", "unit_locations"],
            extension_params={
                "unit_locations": {"method": "center_of_mass"},
                "correlograms": {"window_ms": 100},
            },
            **job_kwargs,
        )

        # Sparse curated — save directly instead of .copy() to avoid
        # provenance resolution issues with relative sorter paths
        sparse_curated_path = save_dir / "stage1_analyzer_curated_sparse.zarr"
        if sparse_curated_path.exists():
            shutil.rmtree(sparse_curated_path)
        sparse_curated = curated_analyzer.save_as(
            folder=sparse_curated_path, format="zarr")
        sparse_curated.sparsity = si.compute_sparsity(
            curated_analyzer, method="radius", radius_um=60)
        sparse_curated.compute(
            ["random_spikes", "waveforms", "templates", "correlograms"],
            extension_params={"correlograms": {"window_ms": 100}},
            **job_kwargs,
        )

    # ─── Cleanup ──────────────────────────────────────────
    # Keep raw analyzer (needed for stage 2 rescue/reject)
    sparse_raw_path = save_dir / "stage1_analyzer_raw_sparse.zarr"
    if sparse_raw_path.exists():
        shutil.rmtree(sparse_raw_path)
    if ms5_output.exists():
        shutil.rmtree(ms5_output)

    print(f"  [OK] Stage 1 complete!")
    return True


def main():
    job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

    force = "--force" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--force"]

    sessions = args
    if not sessions:
        print("Usage: python -m batch_process.stage1_sort <session_folder> [<session_folder> ...]")
        return

    t_start = time.time()
    results = {}

    for i, session in enumerate(sessions, 1):
        session_name = Path(session).name
        print(f"\n{'='*60}")
        print(f"[{i}/{len(sessions)}] {session_name}")
        print(f"{'='*60}")

        try:
            t0 = time.time()
            success = process_session(session, job_kwargs, force_resort=force)
            elapsed = time.time() - t0
            results[session_name] = f"OK ({elapsed/60:.1f}min)"
        except Exception as e:
            print(f"\n  [ERROR] {e}")
            traceback.print_exc()
            results[session_name] = f"ERROR: {e}"

    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"BATCH COMPLETE — {len(sessions)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for name, status in results.items():
        icon = "[OK]" if "OK" in status else "[X]"
        print(f"  {icon} {name}: {status}")


if __name__ == "__main__":
    main()
