"""
Control Spike Sorting Pipeline - Stage 3 (Merge + Finalize)
============================================================
Run AFTER reviewing merge_candidates.pptx from stage2_curate.py
and editing merge_pairs in curation_config.json.

This stage:
1. Loads analyzer_curated.zarr (stage 2 output)
2. Reads merge_pairs from curation_config.json
3. Applies merges (if any)
4. Computes all final extensions
5. Saves analyzer_final.zarr
6. Generates final_units.pptx

If merge_pairs is empty, it just copies curated -> final with full extensions.

Usage (from terminal):
    cd C:\\Users\\robin\\Desktop\\Python\\research\\icms-plasticity
    C:\\Users\\robin\\Desktop\\Python\\research\\si_v103_env\\Scripts\\python.exe pipeline\\stage3_merge.py
"""

import sys
import os
import json
import shutil
import traceback
import time
import types
from pathlib import Path

# Project root for imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402
import spikeinterface as si  # noqa: E402

from batch_process.util.plotting import plot_units_in_batches  # noqa: E402
from util.load_control import discover_sessions  # noqa: E402


def process_session_stage3(sorting_dir, job_kwargs, force=False):
    """
    Merge + finalize a single session.

    Parameters
    ----------
    sorting_dir : Path
        Path to the session's spike_sorting folder.
    job_kwargs : dict
        SpikeInterface job kwargs.
    force : bool
        Re-run even if analyzer_final.zarr already exists.

    Returns
    -------
    bool
        True if successful.
    """
    sorting_dir = Path(sorting_dir)
    curated_path = sorting_dir / "analyzer_curated.zarr"
    final_path = sorting_dir / "analyzer_final.zarr"
    preproc_zarr = sorting_dir / "preprocessed.zarr"
    config_path = sorting_dir / "curation_config.json"

    # --- Checks ---
    if not curated_path.exists():
        print(f"  [SKIP] No analyzer_curated.zarr — run stage2 first")
        return False

    if not force and final_path.exists():
        print(f"  [SKIP] Already finalized: {final_path}")
        return True

    # --- Load config ---
    merge_pairs = []
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
        merge_pairs = config.get("merge_pairs", [])
    else:
        print("  [!] No curation_config.json — proceeding without merges")

    # --- Load curated analyzer ---
    print("  Loading analyzer_curated.zarr...")
    analyzer = si.load_sorting_analyzer(curated_path)

    has_recording = preproc_zarr.exists()
    if has_recording:
        rec = si.load(preproc_zarr)
        analyzer.set_temporary_recording(rec)
        print(f"  Recording loaded: {preproc_zarr}")
    else:
        print("  [!] No preprocessed.zarr — recordingless mode")

    # Workaround: SI provenance bug
    analyzer.get_sorting_provenance = types.MethodType(
        lambda self: None, analyzer)

    # --- Apply merges ---
    if merge_pairs:
        print(f"\n  Applying {len(merge_pairs)} merge(s): {merge_pairs}")

        typed_merge_groups = [
            [analyzer.unit_ids.dtype.type(uid) for uid in group]
            for group in merge_pairs
        ]

        merged_path = sorting_dir / "_analyzer_merged_tmp.zarr"
        if merged_path.exists():
            shutil.rmtree(merged_path)

        analyzer = analyzer.merge_units(
            merge_unit_groups=typed_merge_groups,
            merging_mode="soft",
            new_id_strategy="take_first",
            format="zarr",
            folder=merged_path,
        )

        if has_recording:
            analyzer.set_temporary_recording(rec)

        print(f"  Merged! Unit count: {len(analyzer.unit_ids)}")
    else:
        print("  No merge pairs — finalizing as-is")

    # --- Save as analyzer_final.zarr ---
    print("  Saving analyzer_final.zarr...")
    if final_path.exists():
        shutil.rmtree(final_path)

    # Workaround for merged analyzer provenance
    analyzer.get_sorting_provenance = types.MethodType(
        lambda self: None, analyzer)

    final = analyzer.save_as(format="zarr", folder=final_path)

    if has_recording:
        final.set_temporary_recording(rec)

    # --- Compute extensions ---
    print("  Computing extensions...")
    if has_recording:
        final.compute(
            ["random_spikes", "waveforms", "templates", "correlograms",
             "template_similarity", "unit_locations", "spike_amplitudes"],
            extension_params={
                "random_spikes": {"max_spikes_per_unit": 2000},
                "correlograms": {"window_ms": 100, "bin_ms": 0.5},
                "unit_locations": {"method": "center_of_mass"},
            },
            **job_kwargs,
        )
    else:
        final.compute(
            ["correlograms", "template_similarity"],
            extension_params={"correlograms": {"window_ms": 100, "bin_ms": 0.5}},
            **job_kwargs,
        )

    # --- Generate final PPT ---
    if has_recording:
        print("  Generating final_units.pptx...")
        sparse = final.copy()
        sparse.sparsity = si.compute_sparsity(
            final, method="radius", radius_um=60)
        sparse.compute(
            ["random_spikes", "waveforms", "templates", "correlograms"],
            extension_params={"correlograms": {"window_ms": 100}},
            **job_kwargs,
        )
        plot_units_in_batches(
            sparse, save_dir=sorting_dir, ppt_name="final_units")
    else:
        print("  [!] Skipping PPT (no recording)")

    # --- Cleanup temp merged zarr ---
    merged_tmp = sorting_dir / "_analyzer_merged_tmp.zarr"
    if merged_tmp.exists():
        shutil.rmtree(merged_tmp)

    print(f"\n  Stage 3 complete! {len(final.unit_ids)} final units")
    print(f"  -> {final_path}")
    return True


def run_stage3(base_folder="E:/ICMS plasticity control/ICMS48",
               session_path=None, force=False, allow_missing_mat=True):
    """
    Batch merge + finalize for control sessions.

    Parameters
    ----------
    base_folder : str
        Root data directory for the animal.
    session_path : str or None
        If provided, only process this session's spike_sorting folder.
    force : bool
        Re-run even if analyzer_final.zarr already exists.
    """
    job_kwargs = dict(n_jobs=0.5, chunk_duration="1s", progress_bar=True)

    # Find all spike_sorting dirs
    sessions = discover_sessions(base_folder, allow_missing_mat=allow_missing_mat)

    # Filter to single session if requested
    if session_path is not None:
        session_path = Path(session_path).resolve()
        matched = {k: v for k, v in sessions.items()
                   if Path(k).resolve() == session_path}
        if not matched:
            print(f"Session not found: {session_path}")
            print("Available sessions:")
            for k in sessions:
                print(f"  {k}")
            return
        sessions = matched

    print(f"Processing {len(sessions)} session(s)\n")

    results = {}
    t_start = time.time()

    for i, exp_dir in enumerate(sessions, 1):
        sorting_dir = Path(exp_dir) / "spike_sorting"
        print(f"\n{'='*60}")
        print(f"[{i}/{len(sessions)}] {sorting_dir}")
        print(f"{'='*60}")

        try:
            ok = process_session_stage3(sorting_dir, job_kwargs, force=force)
            results[str(sorting_dir)] = "OK" if ok else "SKIPPED"
        except Exception as e:
            print(f"\n  [ERROR] {e}")
            traceback.print_exc()
            results[str(sorting_dir)] = f"ERROR: {e}"

    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"STAGE 3 COMPLETE -- {len(sessions)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for path, status in results.items():
        icon = "[OK]" if status == "OK" else "[-]" if "SKIP" in status else "[X]"
        print(f"  {icon} {Path(path).parent.name}: {status}")


if __name__ == "__main__":
    control_animals = [
        "ICMS43",
        "ICMS45",
        "ICMS48",
        "ICMS54",
        "ICMS56",
    ]
    for animal in control_animals:
        base = f"E:/ICMS plasticity control/{animal}"
        print(f"\n{'='*60}")
        print(f"  STAGE 3 — {animal}")
        print(f"{'='*60}")
        run_stage3(base_folder=base, force=True)
