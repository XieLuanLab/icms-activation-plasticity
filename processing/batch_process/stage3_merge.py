"""
Stage 3: apply curated merges + finalize, per session.

Reads batch_sort/stage2/stage2_analyzer.zarr + curation_config.json, applies the
selected merges, and writes batch_sort/stage3/analyzer_final.zarr.

Usage:
    python -m batch_process.stage3_merge <session_folder> [<session_folder> ...] [--force]
"""

import sys
import os
import json
import shutil
import traceback
import time
import types
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import spikeinterface as si

from batch_process.util.plotting import plot_units_in_batches




def process_session(data_folder, job_kwargs, force=False):
    """Run stage 3 merge + finalize on a single session."""
    data_folder = Path(data_folder)
    session_name = data_folder.name

    stage2_dir = data_folder / "batch_sort" / "stage2"
    stage3_dir = data_folder / "batch_sort" / "stage3"
    curated_path = stage2_dir / "stage2_analyzer.zarr"
    config_path = stage2_dir / "curation_config.json"
    preproc_zarr = data_folder / "batch_sort" / "stage1" / "preprocessed.zarr"
    final_path = stage3_dir / "analyzer_final.zarr"

    # --- Checks ---
    if not curated_path.exists():
        print(f"  [SKIP] No stage2_analyzer.zarr — run stage2 first")
        return False

    if not force and final_path.exists():
        print(f"  [SKIP] Already finalized: {final_path}")
        return True

    stage3_dir.mkdir(parents=True, exist_ok=True)

    # --- Load config ---
    merge_pairs = []
    remove_ids = []
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
        merge_pairs = config.get("merge_pairs", [])
        remove_ids = config.get("remove_ids", [])
    else:
        print("  [!] No curation_config.json — proceeding without merges")

    # --- Load curated analyzer ---
    print("  Loading stage2_analyzer.zarr...")
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

    # --- Remove units ---
    if remove_ids:
        print(f"  Removing {len(remove_ids)} unit(s): {remove_ids}")
        keep_ids = [uid for uid in analyzer.unit_ids if uid not in remove_ids]
        removed_path = stage3_dir / "_analyzer_removed_tmp.zarr"
        if removed_path.exists():
            shutil.rmtree(removed_path)
        analyzer = analyzer.select_units(
            unit_ids=keep_ids, format="zarr", folder=removed_path)
        if has_recording:
            analyzer.set_temporary_recording(rec)
        analyzer.get_sorting_provenance = types.MethodType(
            lambda self: None, analyzer)
        print(f"  After removal: {len(analyzer.unit_ids)} units")

    # --- Apply merges ---
    if merge_pairs:
        print(f"\n  Applying {len(merge_pairs)} merge(s): {merge_pairs}")

        typed_merge_groups = [
            [analyzer.unit_ids.dtype.type(uid) for uid in group]
            for group in merge_pairs
        ]

        merged_path = stage3_dir / "_analyzer_merged_tmp.zarr"
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

        analyzer.get_sorting_provenance = types.MethodType(
            lambda self: None, analyzer)

        print(f"  Merged! Unit count: {len(analyzer.unit_ids)}")
    else:
        print("  No merge pairs — finalizing as-is")

    # --- Save as analyzer_final.zarr ---
    print("  Saving analyzer_final.zarr...")
    if final_path.exists():
        shutil.rmtree(final_path)

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
            sparse, save_dir=stage3_dir, ppt_name="final_units")
    else:
        print("  [!] Skipping PPT (no recording)")

    # --- Cleanup temp zarrs ---
    for tmp_name in ("_analyzer_merged_tmp.zarr", "_analyzer_removed_tmp.zarr"):
        tmp_path = stage3_dir / tmp_name
        if tmp_path.exists():
            shutil.rmtree(tmp_path)

    print(f"\n  Stage 3 complete! {len(final.unit_ids)} final units")
    print(f"  -> {final_path}")
    return True


def main():
    job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

    force = "--force" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--force"]

    sessions = args
    if not sessions:
        print("Usage: python -m batch_process.stage3_merge <session_folder> [<session_folder> ...]")
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
            success = process_session(session, job_kwargs, force=force)
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
