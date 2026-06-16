"""
Experimental Stage 2 — PPT-based waveform curation (replaces sortingview).

Uses the same curation pipeline as the control animals:
  1. Load curation_config.json for rescue_ids / remove_ids
  2. Select accepted units from stage1_analyzer_raw.zarr
  3. Run waveform curation:
     a. Remove bad waveforms (blanking artifacts, amplitude outliers)
     b. Classify waveforms as spike-like vs non-spike-like (raw waveforms)
     c. Subcluster spike-like waveforms (UMAP + Isosplit6)
     d. Evaluate subclusters (artifact, SNR, spike count)
  4. Create stage2_analyzer.zarr with accepted subclusters only
  5. Generate PPT reports and curation diagnostics

Works for any experimental animal (ICMS83, ICMS92, ICMS93, ICMS98, etc.).
Requires stage1 to have been run first.

Usage:
    cd C:\\Users\\robin\\Desktop\\Python\\research\\icms-plasticity
    python -m batch_process.stage2_curate E:/ICMS83/20-Jul-2023/batch_sort/stage1
    python -m batch_process.stage2_curate C:/data/ICMS92/behavior/21-Sep-2023/batch_sort/stage1
"""

import sys
import os
import json
import shutil
import traceback
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import spikeinterface as si
import spikeinterface.preprocessing as sp
import matplotlib.pyplot as plt

from batch_process.util.curate_util import (
    remove_bad_waveforms_A,
    remove_bad_waveforms_B,
)
from batch_process.visualization.spike_classifier import (
    get_MAD_threshold_primary_ch,
    extract_waveform_features,
    classify_waveform,
)
from batch_process.util.subcluster_util import (
    cluster_waveforms,
    accept_subcluster,
)
from batch_process.util.plotting import plot_units_in_batches
from batch_process.util.template_util import get_unit_primary_ch_wvfs


# ─── Helpers ──────────────────────────────────────────────────────

def load_curation_config(stage1_dir):
    """Load curation_config.json if it exists."""
    config_path = Path(stage1_dir) / "curation_config.json"
    if config_path.exists():
        with open(config_path) as f:
            return json.load(f)
    print("  [!] No curation_config.json found, using auto-labels as-is")
    return {"rescue_ids": [], "remove_ids": [], "merge_pairs": []}


def get_accepted_ids(stage1_dir, config):
    """Apply rescue/remove from config to auto-curation labels."""
    labels_file = Path(stage1_dir) / "curation_labels.npz"
    if not labels_file.exists():
        raise FileNotFoundError(f"No curation_labels.npz in {stage1_dir}")

    data = np.load(labels_file, allow_pickle=True)
    good_ids = list(data["good_ids"])
    bad_ids = list(data["bad_ids"])

    for uid in config.get("rescue_ids", []):
        if uid in bad_ids:
            bad_ids.remove(uid)
            good_ids.append(uid)
            print(f"    Rescued unit {uid}")
        elif uid not in good_ids:
            print(f"    [!] Unit {uid} not found, skipping rescue")

    for uid in config.get("remove_ids", []):
        if uid in good_ids:
            good_ids.remove(uid)
            bad_ids.append(uid)
            print(f"    Removed unit {uid}")
        elif uid not in bad_ids:
            print(f"    [!] Unit {uid} not found, skipping removal")

    return sorted(good_ids)


def load_raw_recording(stage1_dir):
    """
    Load the raw (unfiltered) recording for spike classification.

    Tries to reconstruct from original NS5 files based on the animal type.
    Falls back to preprocessed.zarr if raw reconstruction fails.
    """
    stage1_dir = Path(stage1_dir)
    batch_sort_dir = stage1_dir.parent  # batch_sort/
    session_dir = batch_sort_dir.parent  # e.g. 20-Jul-2023/

    # Try ICMS83 loader
    if "ICMS83" in str(session_dir):
        try:
            from batch_process.util.dataloader_mouse6 import load_mouse6
            print("  Loading raw recording via Mouse 6 loader...")
            dl = load_mouse6(str(session_dir), get_recording=True)
            return dl.recording
        except Exception as e:
            print(f"  [!] ICMS83 raw load failed: {e}")

    # Try standard experimental loader
    try:
        from util.load_data import load_data
        print("  Loading raw recording via standard loader...")
        dl = load_data(
            str(session_dir), make_folder=False,
            save_folder_name="batch_sort", first_N_files=4,
            server_mount_drive="S:")
        return dl.recording
    except Exception as e:
        print(f"  [!] Standard raw load failed: {e}")

    # Fallback: use preprocessed recording (less ideal for classification
    # but still functional)
    preproc_zarr = stage1_dir / "preprocessed.zarr"
    if preproc_zarr.exists():
        print("  [FALLBACK] Using preprocessed recording for classification")
        return si.read_zarr(preproc_zarr)

    return None


# ─── Waveform curation (mirrors control pipeline) ────────────────

def curate_unit_waveforms(analyzer, raw_analyzer, unit_id):
    """
    Run waveform curation pipeline for a single unit.

    Returns
    -------
    waveform_labels : np.ndarray of int, shape (n_spikes,)
        -1 = bad waveform, 0 = non-spike-like, 1+ = subcluster ID
    subcluster_verdicts : list of (int, str)
        (subcluster_id, verdict) pairs
    """
    # Step 1: Remove bad waveforms (blanking artifacts, amplitude outliers)
    good_A, bad_A = remove_bad_waveforms_A(analyzer, unit_id)
    good_B, bad_B = remove_bad_waveforms_B(analyzer, unit_id)
    total_num = len(good_A) + len(bad_A)

    bad_set = set(bad_A) | set(bad_B)
    good_indices = sorted(set(range(total_num)) - bad_set)

    waveform_labels = np.full(total_num, -1, dtype=int)
    print(f"    Step 1: {len(bad_set)}/{total_num} bad waveforms removed")

    if len(good_indices) < 20:
        print(f"    Too few good waveforms ({len(good_indices)}), skipping")
        return waveform_labels, []

    # Step 2: Classify spike-like vs non-spike-like using raw waveforms
    # Clip good_indices to raw waveform count (raw recording may be shorter)
    raw_wvf_ext = raw_analyzer.get_extension("waveforms")
    raw_n_wvfs = raw_wvf_ext.get_waveforms_one_unit(unit_id=unit_id).shape[0]
    if max(good_indices) >= raw_n_wvfs:
        n_clipped = sum(1 for i in good_indices if i >= raw_n_wvfs)
        good_indices = [i for i in good_indices if i < raw_n_wvfs]
        print(f"    [!] Clipped {n_clipped} indices beyond raw waveform count ({raw_n_wvfs})")
        if len(good_indices) < 20:
            print(f"    Too few good waveforms after clipping ({len(good_indices)}), skipping")
            return waveform_labels, []

    threshold = get_MAD_threshold_primary_ch(analyzer, unit_id)
    centered_raw_wvfs, _ = extract_waveform_features(
        raw_analyzer, analyzer, unit_id, good_indices)

    spike_like = []
    non_spike_like = []
    for idx, wvf in enumerate(centered_raw_wvfs):
        is_spike, _ = classify_waveform(wvf, threshold)
        (spike_like if is_spike else non_spike_like).append(idx)

    good_indices = np.array(good_indices)
    spike_like = np.array(spike_like)
    non_spike_like = np.array(non_spike_like)

    if len(non_spike_like) > 0:
        waveform_labels[good_indices[non_spike_like]] = 0
    print(f"    Step 2: {len(spike_like)} spike-like, "
          f"{len(non_spike_like)} non-spike-like")

    if len(spike_like) < 20:
        print(f"    Too few spike-like waveforms, skipping")
        return waveform_labels, []

    # Step 3: Subcluster (UMAP + Isosplit6) on samples 20-42 of primary ch
    spike_wvfs = centered_raw_wvfs[spike_like]
    samples_to_use = np.arange(20, 42)
    subcluster_labels, _, _ = cluster_waveforms(
        spike_wvfs[:, samples_to_use])

    final_spike_indices = good_indices[spike_like]
    waveform_labels[final_spike_indices] = subcluster_labels

    unique_sc = np.unique(subcluster_labels)
    print(f"    Step 3: {len(unique_sc)} subcluster(s)")

    # Step 4: Evaluate each subcluster
    verdicts = []
    for sc_id in unique_sc:
        sc_mask = subcluster_labels == sc_id
        sc_indices = final_spike_indices[sc_mask]
        _, _, verdict = accept_subcluster(analyzer, unit_id, sc_indices)
        verdicts.append((int(sc_id), verdict))
        print(f"      SC {sc_id}: {np.sum(sc_mask)} spikes -> {verdict}")

    return waveform_labels, verdicts


def make_curation_diagnostic_figure(analyzer, unit_id, labels, verdicts):
    """
    Create a before/after figure for a single unit's waveform curation.

    3 panels on the primary channel:
      - Left:   All waveforms (before curation)
      - Middle: Accepted waveforms (kept subclusters)
      - Right:  Rejected waveforms (bad + non-spike + rejected subclusters)
    """
    primary_wvfs = get_unit_primary_ch_wvfs(analyzer, unit_id)
    n_total = primary_wvfs.shape[0]

    accepted_sc_ids = {sc_id for sc_id, v in verdicts if v == "accept"}
    accepted_mask = np.isin(labels, list(accepted_sc_ids))
    rejected_mask = ~accepted_mask

    acc_wvfs = primary_wvfs[accepted_mask]
    rej_wvfs = primary_wvfs[rejected_mask]

    max_traces = 500
    rng = np.random.default_rng(42)

    def subsample(wvfs):
        if len(wvfs) > max_traces:
            idx = rng.choice(len(wvfs), max_traces, replace=False)
            return wvfs[idx]
        return wvfs

    all_sub = subsample(primary_wvfs)
    acc_sub = subsample(acc_wvfs)
    rej_sub = subsample(rej_wvfs)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)

    ylim_min = np.percentile(primary_wvfs, 0.5)
    ylim_max = np.percentile(primary_wvfs, 99.5)

    titles = [
        f"All ({n_total})",
        f"Accepted ({len(acc_wvfs)})",
        f"Rejected ({len(rej_wvfs)})",
    ]
    colors = ["gray", "#2196F3", "#F44336"]
    datasets = [all_sub, acc_sub, rej_sub]

    for ax, title, color, data in zip(axes, titles, colors, datasets):
        if len(data) > 0:
            ax.plot(data.T, color=color, alpha=0.05, linewidth=0.5)
            ax.plot(np.mean(data, axis=0), color="black", linewidth=1.5)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(ylim_min, ylim_max)
        ax.set_xlabel("Samples")
    axes[0].set_ylabel("Amplitude")

    verdict_str = ", ".join(
        f"SC{sc_id}={v}" for sc_id, v in verdicts) if verdicts else "no subclusters"
    fig.suptitle(f"Unit {unit_id} — {verdict_str}", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


# ─── Session processing ──────────────────────────────────────────

def process_session(stage1_dir, job_kwargs, force=False):
    """
    Run stage 2 waveform curation on a single experimental session.

    Parameters
    ----------
    stage1_dir : str or Path
        Path to the stage1 output folder, e.g.
        E:/ICMS83/20-Jul-2023/batch_sort/stage1
    """
    stage1_dir = Path(stage1_dir)

    # Check for existing output
    stage2_path = stage1_dir.parent / "stage2"
    if not force and (stage2_path / "stage2_analyzer.zarr").exists():
        print(f"  [SKIP] Already curated: {stage2_path}")
        return True

    # ─── 1. Find the raw analyzer (has ALL units) ─────────────
    raw_analyzer_path = stage1_dir / "stage1_analyzer_raw.zarr"
    curated_path = stage1_dir / "stage1_analyzer_curated_sparse.zarr"
    preproc_zarr = stage1_dir / "preprocessed.zarr"

    use_raw_analyzer = raw_analyzer_path.exists()

    if use_raw_analyzer:
        print("  Loading stage1_analyzer_raw.zarr...")
        analyzer_raw = si.load_sorting_analyzer(raw_analyzer_path)
    elif curated_path.exists():
        print("  [!] No raw analyzer — using curated (rescue not available)")
        analyzer_raw = si.load_sorting_analyzer(curated_path)
    else:
        print("  [ERROR] No analyzer found in stage1 dir")
        return False

    # Attach preprocessed recording
    if preproc_zarr.exists():
        rec = si.read_zarr(preproc_zarr)
        analyzer_raw.set_temporary_recording(rec)
    else:
        print("  [ERROR] No preprocessed.zarr found")
        return False

    # Workaround: SI provenance bug with deleted sorter output
    import types
    analyzer_raw.get_sorting_provenance = types.MethodType(
        lambda self: None, analyzer_raw)

    # ─── 2. Load config + get accepted IDs ────────────────────
    config = load_curation_config(stage1_dir)

    if use_raw_analyzer:
        accepted_ids = get_accepted_ids(stage1_dir, config)
    else:
        # Fallback: all units in the curated analyzer are accepted
        accepted_ids = list(analyzer_raw.unit_ids)
        if config.get("remove_ids"):
            accepted_ids = [uid for uid in accepted_ids
                           if uid not in config["remove_ids"]]

    print(f"  Accepted units: {len(accepted_ids)}")
    if not accepted_ids:
        print("  [!] No accepted units")
        return False

    # ─── 3. Select accepted units + compute all waveforms ─────
    accepted_typed = np.array(accepted_ids, dtype=analyzer_raw.unit_ids.dtype)
    temp_path = stage1_dir / "_temp_accepted.zarr"
    if temp_path.exists():
        shutil.rmtree(temp_path)

    print("  Selecting accepted units + computing waveforms...")
    acc_analyzer = analyzer_raw.select_units(
        unit_ids=accepted_typed, format="zarr", folder=temp_path)
    acc_analyzer.compute(
        ["random_spikes", "waveforms", "templates", "spike_amplitudes"],
        extension_params={"random_spikes": {"max_spikes_per_unit": 500_000}},
        **job_kwargs,
    )

    # ─── 4. Load raw (unfiltered) recording for classification ─
    print("  Loading raw recording for waveform classification...")
    raw_rec = load_raw_recording(stage1_dir)
    if raw_rec is None:
        print("  [ERROR] Could not load any recording for raw waveforms")
        if temp_path.exists():
            shutil.rmtree(temp_path)
        return False

    # Verify sample alignment
    n_preproc = rec.get_num_samples(0)
    n_raw = raw_rec.get_num_samples(0)
    if n_preproc != n_raw:
        print(f"  [WARNING] Sample mismatch: preprocessed={n_preproc}, "
              f"raw={n_raw} (diff={abs(n_preproc - n_raw)})")

    # ─── 5. Create raw analyzer (in-memory, temporary) ────────
    print("  Computing raw waveforms...")
    raw_analyzer = si.create_sorting_analyzer(
        sorting=acc_analyzer.sorting, recording=raw_rec, sparse=False)
    raw_analyzer.compute(
        ["random_spikes", "waveforms", "templates"],
        extension_params={"random_spikes": {"max_spikes_per_unit": 500_000}},
        **job_kwargs,
    )

    # ─── 6. Run waveform curation for each unit ───────────────
    print("\n  === Waveform Curation ===")
    fs = rec.get_sampling_frequency()
    kept_spikes = {}
    diagnostic_figs = []
    n_kept = 0
    n_dropped = 0
    n_split = 0

    for unit_id in acc_analyzer.unit_ids:
        print(f"\n  Unit {unit_id}:")
        spike_train = acc_analyzer.sorting.get_unit_spike_train(unit_id)
        labels, verdicts = curate_unit_waveforms(
            acc_analyzer, raw_analyzer, unit_id)

        # Diagnostic figure
        try:
            fig = make_curation_diagnostic_figure(
                acc_analyzer, unit_id, labels, verdicts)
            diagnostic_figs.append((unit_id, fig))
        except Exception as e:
            print(f"    [!] Diagnostic figure failed: {e}")

        accepted = [(sc_id, v) for sc_id, v in verdicts if v == "accept"]
        if not accepted:
            print(f"    -> DROPPED")
            n_dropped += 1
            continue

        if len(accepted) == 1:
            sc_id = accepted[0][0]
            mask = labels == sc_id
            kept_spikes[int(unit_id)] = spike_train[mask]
            print(f"    -> KEPT {np.sum(mask)}/{len(spike_train)} spikes")
            n_kept += 1
        else:
            for sc_id, _ in accepted:
                mask = labels == sc_id
                n = np.sum(mask)
                if n < 50:
                    print(f"    -> SC {sc_id}: too few spikes ({n}), skipping")
                    continue
                new_id = int(unit_id) * 1000 + sc_id
                kept_spikes[new_id] = spike_train[mask]
                print(f"    -> Split: unit {new_id} ({n} spikes)")
            n_split += 1

    if not kept_spikes:
        print("\n  [!] No units survived curation!")
        if temp_path.exists():
            shutil.rmtree(temp_path)
        return False

    print(f"\n  Summary: {n_kept} kept, {n_split} split, {n_dropped} dropped "
          f"-> {len(kept_spikes)} curated units")

    # ─── 7. Create stage2_analyzer.zarr ───────────────────────
    stage2_path.mkdir(parents=True, exist_ok=True)

    print("\n  Creating stage2_analyzer.zarr...")
    curated_sorting = si.NumpySorting.from_unit_dict(
        [kept_spikes], sampling_frequency=fs)

    analyzer_path = stage2_path / "stage2_analyzer.zarr"
    if analyzer_path.exists():
        shutil.rmtree(analyzer_path)

    curated = si.create_sorting_analyzer(
        sorting=curated_sorting, recording=rec,
        format="zarr", folder=analyzer_path, sparse=False, overwrite=True)

    curated.compute(
        ["random_spikes", "waveforms", "templates", "correlograms",
         "template_similarity", "unit_locations", "spike_amplitudes"],
        extension_params={
            "random_spikes": {"max_spikes_per_unit": 2000},
            "correlograms": {"window_ms": 100, "bin_ms": 0.5},
            "unit_locations": {"method": "center_of_mass"},
        },
        **job_kwargs,
    )

    # ─── 8. Generate curated PPT ──────────────────────────────
    print("  Generating curated PPT...")
    sparse_path = stage2_path / "stage2_analyzer_sparse.zarr"
    if sparse_path.exists():
        shutil.rmtree(sparse_path)
    sparse = curated.save_as(folder=sparse_path, format="zarr")
    sparse.sparsity = si.compute_sparsity(
        curated, method="radius", radius_um=60)
    sparse.compute(
        ["random_spikes", "waveforms", "templates", "template_similarity", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )
    plot_units_in_batches(
        sparse, save_dir=stage2_path, ppt_name="curated_units_stage2")

    # ─── 9. Save diagnostic PPT ──────────────────────────────
    if diagnostic_figs:
        from batch_process.util.ppt_image_inserter import PPTImageInserter
        diag_ppt = PPTImageInserter(
            grid_dims=(1, 1), spacing=(0.05, 0.05), title_font_size=14)
        diag_img = stage2_path / "_diag_tmp.png"
        for uid, fig in diagnostic_figs:
            fig.savefig(str(diag_img), dpi=150, bbox_inches="tight")
            diag_ppt.add_image(str(diag_img))
            plt.close(fig)
        diag_ppt.save(stage2_path / "curation_diagnostic.pptx")
        if diag_img.exists():
            diag_img.unlink()
        print(f"  Diagnostic PPT: {stage2_path / 'curation_diagnostic.pptx'}")

    # ─── 10. Merge candidate analysis ─────────────────────────
    try:
        from curation.merge_utils import (
            compute_merge_candidates,
            print_merge_report,
            generate_merge_comparison_figure,
        )
        session_name = f"{stage1_dir.parent.parent.name}/{stage1_dir.parent.name}"
        candidates = compute_merge_candidates(curated, list(curated.unit_ids))
        print_merge_report(candidates, session_name=session_name)

        if candidates:
            from batch_process.util.ppt_image_inserter import PPTImageInserter
            ppt = PPTImageInserter(
                grid_dims=(1, 1), spacing=(0.05, 0.05), title_font_size=14)
            merge_img = stage2_path / "_merge_tmp.png"
            for c in candidates:
                fig = generate_merge_comparison_figure(
                    curated, c["unit_a"], c["unit_b"])
                fig.savefig(str(merge_img), dpi=150)
                ppt.add_image(str(merge_img))
                plt.close(fig)
            ppt.save(stage2_path / "merge_candidates.pptx")
            if merge_img.exists():
                merge_img.unlink()
            print(f"  Merge PPT: {stage2_path / 'merge_candidates.pptx'}")
    except Exception as e:
        print(f"  [!] Merge analysis failed: {e}")

    # ─── 11. Generate stage2 curation_config.json ──────────────
    stage2_config_path = stage2_path / "curation_config.json"
    if not stage2_config_path.exists():
        with open(stage2_config_path, "w") as f:
            json.dump({"remove_ids": [], "merge_pairs": []}, f, indent=4)
        print(f"  Created {stage2_config_path} — edit before running stage 3")

    # ─── 12. Cleanup ──────────────────────────────────────────
    if temp_path.exists():
        shutil.rmtree(temp_path)

    print(f"\n  Stage 2 complete! {len(kept_spikes)} curated units")
    print(f"  -> {analyzer_path}")
    return True




# ─── CLI ──────────────────────────────────────────────────────────

def main():
    if len(sys.argv) > 1:
        job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)
        force = "--force" in sys.argv
        process_session(sys.argv[1], job_kwargs, force=force)
    else:
        print("Usage: python -m batch_process.stage2_curate <stage1_dir> [--force]")


if __name__ == "__main__":
    main()
