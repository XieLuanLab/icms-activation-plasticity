"""
Control Spike Sorting Pipeline - Stage 2 (Waveform Curation)
=============================================================
Run AFTER reviewing PPTs from stage1_sort.py and editing curation_config.json.

This stage:
1. Reads curation_config.json for rescue_ids / remove_ids
2. Selects accepted units from analyzer_raw.zarr
3. Reconstructs raw (unfiltered) condensed recording from NS5 files
4. Runs waveform curation (same methods as experimental batch_process/stage3_v3.py):
   a. Remove bad waveforms (blanking artifacts, amplitude outliers)
   b. Classify waveforms as spike-like vs non-spike-like (raw waveform features)
   c. Subcluster spike-like waveforms (UMAP + Isosplit6)
   d. Evaluate subclusters (artifact check, SNR, spike count)
5. Creates analyzer_curated.zarr with accepted subclusters only
6. Generates PPT and merge candidate report

Usage (from terminal):
    cd C:\\Users\\robin\\Desktop\\Python\\research\\icms-plasticity
    C:\\Users\\robin\\Desktop\\Python\\research\\si_v103_env\\Scripts\\python.exe pipeline\\stage2_curate.py
"""

import sys
import os
import json
import shutil
import traceback
import time
from pathlib import Path

# Project root for imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import matplotlib
matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import spikeinterface as si  # noqa: E402

from batch_process.util.curate_util import (  # noqa: E402
    remove_bad_waveforms_A,
    remove_bad_waveforms_B,
)
from batch_process.visualization.spike_classifier import (  # noqa: E402
    get_MAD_threshold_primary_ch,
    extract_waveform_features,
    classify_waveform,
)
from batch_process.util.subcluster_util import (  # noqa: E402
    cluster_waveforms,
    accept_subcluster,
)
from batch_process.util.plotting import plot_units_in_batches  # noqa: E402
from batch_process.util.template_util import get_unit_primary_ch_wvfs  # noqa: E402
from util.load_control import (  # noqa: E402
    discover_sessions,
    extract_trial_parameters_from_rep,
)

# Must match stage 1 condensation parameters
PRE_WINDOW_S = 5.0
POST_WINDOW_S = 5.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_curation_config(sorting_dir):
    """Load curation_config.json if it exists, else return defaults."""
    config_path = Path(sorting_dir) / "curation_config.json"
    if config_path.exists():
        with open(config_path) as f:
            return json.load(f)
    print("  [!] No curation_config.json found, using auto-labels as-is")
    return {"rescue_ids": [], "remove_ids": [], "merge_pairs": []}


def get_accepted_ids(sorting_dir, config):
    """Apply rescue/remove from config to auto-curation labels."""
    labels_file = Path(sorting_dir) / "curation_labels.npz"
    if not labels_file.exists():
        raise FileNotFoundError(f"No curation_labels.npz in {sorting_dir}")

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


def reconstruct_raw_condensed_recording(exp_dir, reps, probe_json_path=None,
                                        survey_current_uA=5.0, max_reps=None):
    """
    Re-load NS5 files and reconstruct the condensed recording WITHOUT any
    preprocessing (no filtering, no artifact removal). Uses the same
    condensation windows (PRE/POST_WINDOW_S) as stage 1 so spike indices
    align with the preprocessed recording.
    """
    import probeinterface as pi
    if probe_json_path is None:
        probe_json_path = os.path.join(PROJECT_ROOT, 'util', 'control_probe.json')
    probegroup = pi.read_probeinterface(probe_json_path)

    if max_reps is not None:
        reps = reps[:max_reps]

    recordings = []
    all_trials_dfs = []
    cumulative_samples = 0
    fs = None

    for rep in reps:
        print(f"    Loading Rep {rep['rep_num']} (raw NS5)...")
        try:
            rec_rep = si.extractors.read_blackrock(rep['ns5_path'])
            rec_rep = rec_rep.set_probegroup(probegroup)
            if fs is None:
                fs = rec_rep.get_sampling_frequency()

            df = extract_trial_parameters_from_rep(
                rep, survey_current_uA=survey_current_uA)

            if cumulative_samples > 0:
                offset_s = cumulative_samples / fs
                df['stim_start_s'] = df['stim_start_s'].apply(
                    lambda x: x + offset_s if x is not None else None)
                df['stim_timestamps'] = df['stim_timestamps'].apply(
                    lambda ts: [t + offset_s for t in ts])

            all_trials_dfs.append(df)
            recordings.append(rec_rep)
            cumulative_samples += rec_rep.get_num_samples(0)
        except Exception as e:
            print(f"    [!] Failed to load Rep {rep['rep_num']}: {e}")

    if not recordings:
        return None

    session_df = pd.concat(all_trials_dfs, ignore_index=True)
    multi_rec = si.concatenate_recordings(recordings)

    # Apply same condensation as stage 1
    session_df = session_df.sort_values('stim_start_s').reset_index(drop=True)
    trial_recordings = []
    for _, row in session_df.iterrows():
        t_start = max(0.0, row['stim_start_s'] - PRE_WINDOW_S)
        t_end = min(multi_rec.get_total_duration(),
                    row['stim_start_s'] + POST_WINDOW_S)
        slice_rec = multi_rec.frame_slice(
            start_frame=int(t_start * fs), end_frame=int(t_end * fs))
        trial_recordings.append(slice_rec)

    condensed_raw = si.concatenate_recordings(trial_recordings)
    print(f"    Raw condensed: {condensed_raw.get_total_duration():.1f}s, "
          f"{condensed_raw.get_num_samples(0)} samples")
    return condensed_raw


# ---------------------------------------------------------------------------
# Waveform curation (mirrors batch_process/stage3_v3.py)
# ---------------------------------------------------------------------------

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
    # Step 1: Remove bad waveforms (blanking artifacts, wrong primary ch,
    #         amplitude outliers)
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

    Parameters
    ----------
    analyzer : SortingAnalyzer
        The accepted-units analyzer (has all waveforms for this unit).
    unit_id : int
        Unit ID.
    labels : np.ndarray
        Per-spike labels: -1=bad, 0=non-spike-like, 1+=subcluster ID.
    verdicts : list of (int, str)
        (subcluster_id, verdict) pairs from accept_subcluster.

    Returns
    -------
    fig : matplotlib.Figure
    """
    import matplotlib.pyplot as plt

    primary_wvfs = get_unit_primary_ch_wvfs(analyzer, unit_id)
    n_total = primary_wvfs.shape[0]

    # Build accepted mask from labels + verdicts
    accepted_sc_ids = {sc_id for sc_id, v in verdicts if v == "accept"}
    accepted_mask = np.isin(labels, list(accepted_sc_ids))
    rejected_mask = ~accepted_mask

    acc_wvfs = primary_wvfs[accepted_mask]
    rej_wvfs = primary_wvfs[rejected_mask]

    # Subsample for plotting (max 500 traces per panel)
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

    # Consistent y-limits
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

    # Verdict summary in subtitle
    verdict_str = ", ".join(
        f"SC{sc_id}={v}" for sc_id, v in verdicts) if verdicts else "no subclusters"
    fig.suptitle(f"Unit {unit_id} — {verdict_str}", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Session processing
# ---------------------------------------------------------------------------

def process_session_stage2(exp_dir, reps, job_kwargs, survey_current_uA=5.0,
                           force=False, probe_json_path=None, max_reps=None):
    """Process a single session through waveform curation."""
    sorting_dir = Path(exp_dir) / "spike_sorting"

    if not (sorting_dir / "analyzer_raw.zarr").exists():
        print(f"  [SKIP] No analyzer_raw.zarr")
        return False
    if not force and (sorting_dir / "analyzer_curated.zarr").exists():
        print(f"  [SKIP] Already curated")
        return True

    preproc_zarr = sorting_dir / "preprocessed.zarr"
    if not preproc_zarr.exists():
        print(f"  [!] preprocessed.zarr not found — cannot curate waveforms")
        return False

    # 1. Load config + accepted IDs
    config = load_curation_config(sorting_dir)
    accepted_ids = get_accepted_ids(sorting_dir, config)
    print(f"  Accepted units: {len(accepted_ids)}")
    if not accepted_ids:
        print("  [!] No accepted units")
        return False

    # 2. Load analyzer + recording
    print("  Loading analyzer...")
    analyzer = si.load_sorting_analyzer(sorting_dir / "analyzer_raw.zarr")
    rec = si.load(preproc_zarr)
    analyzer.set_temporary_recording(rec)

    # Workaround: SI provenance bug with deleted sorter output
    import types
    analyzer.get_sorting_provenance = types.MethodType(
        lambda self: None, analyzer)

    # 3. Select accepted units -> temp analyzer with ALL waveforms
    accepted_typed = np.array(accepted_ids, dtype=analyzer.unit_ids.dtype)
    temp_path = sorting_dir / "_temp_accepted.zarr"
    if temp_path.exists():
        shutil.rmtree(temp_path)

    print("  Selecting accepted units + computing all waveforms...")
    acc_analyzer = analyzer.select_units(
        unit_ids=accepted_typed, format="zarr", folder=temp_path)
    acc_analyzer.compute(
        ["random_spikes", "waveforms", "templates", "spike_amplitudes"],
        extension_params={"random_spikes": {"max_spikes_per_unit": 500_000}},
        **job_kwargs,
    )

    # 4. Reconstruct raw condensed recording (no filtering)
    print("  Reconstructing raw condensed recording...")
    raw_rec = reconstruct_raw_condensed_recording(
        exp_dir, reps, probe_json_path=probe_json_path,
        survey_current_uA=survey_current_uA, max_reps=max_reps)
    if raw_rec is None:
        print("  [!] Failed to reconstruct raw recording")
        return False

    # Verify sample alignment
    n_preproc = rec.get_num_samples(0)
    n_raw = raw_rec.get_num_samples(0)
    if n_preproc != n_raw:
        print(f"  [WARNING] Sample mismatch: preprocessed={n_preproc}, "
              f"raw={n_raw} (diff={abs(n_preproc - n_raw)})")

    # 5. Create raw analyzer (in-memory, temporary)
    print("  Computing raw waveforms (this may take a while)...")
    raw_analyzer = si.create_sorting_analyzer(
        sorting=acc_analyzer.sorting, recording=raw_rec, sparse=False)
    raw_analyzer.compute(
        ["random_spikes", "waveforms", "templates"],
        extension_params={"random_spikes": {"max_spikes_per_unit": 500_000}},
        **job_kwargs,
    )

    # 6. Run waveform curation for each unit
    print("\n  === Waveform Curation ===")
    import matplotlib.pyplot as plt
    fs = rec.get_sampling_frequency()
    kept_spikes = {}  # new_unit_id -> spike_sample_indices
    diagnostic_figs = []  # (unit_id, fig) for before/after PPT
    n_kept = 0
    n_dropped = 0
    n_split = 0

    for unit_id in acc_analyzer.unit_ids:
        print(f"\n  Unit {unit_id}:")
        spike_train = acc_analyzer.sorting.get_unit_spike_train(unit_id)
        labels, verdicts = curate_unit_waveforms(
            acc_analyzer, raw_analyzer, unit_id)

        # Generate before/after diagnostic figure
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
        return False

    print(f"\n  Summary: {n_kept} kept, {n_split} split, {n_dropped} dropped "
          f"-> {len(kept_spikes)} curated units")

    # 7. Create analyzer_curated.zarr
    print("\n  Creating analyzer_curated.zarr...")
    curated_sorting = si.NumpySorting.from_unit_dict(
        [kept_spikes], sampling_frequency=fs)

    curated_path = sorting_dir / "analyzer_curated.zarr"
    if curated_path.exists():
        shutil.rmtree(curated_path)

    curated = si.create_sorting_analyzer(
        sorting=curated_sorting, recording=rec,
        format="zarr", folder=curated_path, sparse=False, overwrite=True)

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

    # 8. Generate PPT
    print("  Generating curated PPT...")
    sparse = curated.copy()
    sparse.sparsity = si.compute_sparsity(
        curated, method="radius", radius_um=60)
    sparse.compute(
        ["random_spikes", "waveforms", "templates", "correlograms"],
        extension_params={"correlograms": {"window_ms": 100}},
        **job_kwargs,
    )
    plot_units_in_batches(
        sparse, save_dir=sorting_dir, ppt_name="curated_units_stage2")

    # 9. Save diagnostic PPT (accepted vs rejected waveforms)
    if diagnostic_figs:
        from batch_process.util.ppt_image_inserter import PPTImageInserter
        diag_ppt = PPTImageInserter(
            grid_dims=(1, 1), spacing=(0.05, 0.05), title_font_size=14)
        diag_img = sorting_dir / "_diag_tmp.png"
        for uid, fig in diagnostic_figs:
            fig.savefig(str(diag_img), dpi=150, bbox_inches="tight")
            diag_ppt.add_image(str(diag_img))
            plt.close(fig)
        diag_ppt.save(sorting_dir / "curation_diagnostic.pptx")
        if diag_img.exists():
            diag_img.unlink()
        print(f"  Diagnostic PPT: {sorting_dir / 'curation_diagnostic.pptx'}")

    # 10. Merge candidate analysis
    try:
        from curation.merge_utils import (
            compute_merge_candidates,
            print_merge_report,
            generate_merge_comparison_figure,
        )
        session_name = (f"{sorting_dir.parent.parent.name}/"
                        f"{sorting_dir.parent.name}")
        candidates = compute_merge_candidates(
            curated, list(curated.unit_ids))
        print_merge_report(candidates, session_name=session_name)

        if candidates:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from batch_process.util.ppt_image_inserter import PPTImageInserter

            ppt = PPTImageInserter(
                grid_dims=(1, 1), spacing=(0.05, 0.05), title_font_size=14)
            for c in candidates:
                fig = generate_merge_comparison_figure(
                    curated, c["unit_a"], c["unit_b"])
                img_path = sorting_dir / "_merge_tmp.png"
                fig.savefig(str(img_path), dpi=150)
                ppt.add_image(str(img_path))
                plt.close(fig)
            ppt.save(sorting_dir / "merge_candidates.pptx")
            print(f"  Merge PPT: {sorting_dir / 'merge_candidates.pptx'}")
    except Exception as e:
        print(f"  [!] Merge analysis failed: {e}")

    # 11. Cleanup temp files
    if temp_path.exists():
        shutil.rmtree(temp_path)
    merge_tmp = sorting_dir / "_merge_tmp.png"
    if merge_tmp.exists():
        merge_tmp.unlink()

    print(f"\n  Stage 2 complete! {len(kept_spikes)} curated units")
    print(f"  -> {curated_path}")
    return True


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run_stage2(base_folder="E:/ICMS plasticity control/ICMS48",
               session_path=None,
               survey_current_uA=5.0, force=False, probe_json_path=None,
               allow_missing_mat=True, max_reps=None):
    """
    Waveform curation for control sessions.

    Parameters
    ----------
    base_folder : str
        Root data directory for the animal.
    session_path : str or None
        Path to a specific session's ExperimentData (or SurveyData) folder.
        If provided, only that session is processed.
        E.g. "E:/ICMS plasticity control/ICMS45/2-26-2022/26-Feb-2022/ExperimentData"
    survey_current_uA : float
        Current to assign to survey sessions lacking a .mat file.
    force : bool
        Re-run even if analyzer_curated.zarr already exists.
    """
    job_kwargs = dict(n_jobs=0.5, chunk_duration="1s", progress_bar=True)

    sessions = discover_sessions(
        base_folder, allow_missing_mat=allow_missing_mat)

    # Filter to single session if requested
    if session_path is not None:
        session_path = Path(session_path).resolve()
        matched = {k: v for k, v in sessions.items()
                   if Path(k).resolve() == session_path}
        if not matched:
            print(f"Session not found: {session_path}")
            print(f"Available sessions:")
            for k in sessions:
                print(f"  {k}")
            return
        sessions = matched

    print(f"Processing {len(sessions)} session(s)\n")

    results = {}
    t_start = time.time()

    for i, (exp_dir, reps) in enumerate(sessions.items(), 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(sessions)}] {exp_dir}")
        print(f"  Reps: {len(reps)}")
        print(f"{'='*60}")

        try:
            ok = process_session_stage2(
                exp_dir, reps, job_kwargs,
                survey_current_uA=survey_current_uA,
                force=force, probe_json_path=probe_json_path,
                max_reps=max_reps)
            results[exp_dir] = "OK" if ok else "FAILED"
        except Exception as e:
            print(f"\n  [ERROR] {e}")
            traceback.print_exc()
            results[exp_dir] = f"ERROR: {e}"

    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"STAGE 2 COMPLETE -- {len(sessions)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for exp_dir, status in results.items():
        icon = "[OK]" if status == "OK" else "[X]"
        print(f"  {icon} {Path(exp_dir).parent.name}/"
              f"{Path(exp_dir).name}: {status}")


def generate_diagnostics_only(base_folder="E:/ICMS plasticity control/ICMS48",
                              session_path=None, allow_missing_mat=True):
    """
    Generate curation_diagnostic.pptx for sessions that already have
    analyzer_raw.zarr and analyzer_curated.zarr. Does NOT re-run curation.

    For each original unit, compares spike trains between raw and curated
    to identify accepted vs rejected spikes, then plots waveforms.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from batch_process.util.ppt_image_inserter import PPTImageInserter

    sessions = discover_sessions(base_folder, allow_missing_mat=allow_missing_mat)

    if session_path is not None:
        session_path = Path(session_path).resolve()
        sessions = {k: v for k, v in sessions.items()
                    if Path(k).resolve() == session_path}

    for exp_dir in sessions:
        sorting_dir = Path(exp_dir) / "spike_sorting"
        raw_zarr = sorting_dir / "analyzer_raw.zarr"
        curated_zarr = sorting_dir / "analyzer_curated.zarr"
        preproc_zarr = sorting_dir / "preprocessed.zarr"

        if not raw_zarr.exists() or not curated_zarr.exists():
            continue

        print(f"\n{'='*60}")
        print(f"Diagnostics: {Path(exp_dir).parent.name}/{Path(exp_dir).name}")

        # Load both analyzers
        raw_analyzer = si.load_sorting_analyzer(raw_zarr)
        curated_analyzer = si.load_sorting_analyzer(curated_zarr)

        if preproc_zarr.exists():
            rec = si.load(preproc_zarr)
            raw_analyzer.set_temporary_recording(rec)
            curated_analyzer.set_temporary_recording(rec)

        # Get auto-curation labels to know which units were accepted into stage2
        curation_npz = sorting_dir / "curation_labels.npz"
        config = load_curation_config(sorting_dir)

        if curation_npz.exists():
            accepted_ids = get_accepted_ids(sorting_dir, config)
        else:
            accepted_ids = list(raw_analyzer.unit_ids)

        curated_unit_ids = set(curated_analyzer.unit_ids)

        # Build curated spike sets for fast lookup
        curated_spikes = {}
        for cuid in curated_analyzer.unit_ids:
            train = curated_analyzer.sorting.get_unit_spike_train(cuid)
            curated_spikes[int(cuid)] = set(train.tolist())

        diag_ppt = PPTImageInserter(
            grid_dims=(1, 1), spacing=(0.05, 0.05), title_font_size=14)
        diag_img = sorting_dir / "_diag_tmp.png"
        n_figs = 0

        for unit_id in accepted_ids:
            # Get all spikes for this unit from raw
            try:
                raw_train = raw_analyzer.sorting.get_unit_spike_train(unit_id)
            except Exception:
                continue

            raw_wvfs = get_unit_primary_ch_wvfs(raw_analyzer, unit_id)

            # Waveforms are stored for a random subset of spikes.
            # Get which spike indices (into the spike train) have waveforms.
            rs_ext = raw_analyzer.get_extension("random_spikes")
            wvf_indices = rs_ext.get_selected_indices_in_spike_train(
                unit_id=unit_id, segment_index=0)
            sampled_spike_samples = raw_train[wvf_indices]

            # Figure out which spikes ended up in curated
            kept_spike_set = set()
            if int(unit_id) in curated_spikes:
                kept_spike_set |= curated_spikes[int(unit_id)]
            for sc in range(1, 10):
                split_id = int(unit_id) * 1000 + sc
                if split_id in curated_spikes:
                    kept_spike_set |= curated_spikes[split_id]

            # Build masks over the sampled waveforms only
            accepted_mask = np.array(
                [s in kept_spike_set for s in sampled_spike_samples])
            rejected_mask = ~accepted_mask

            # Count totals from full spike train for display
            n_total = len(raw_train)
            n_accepted_total = sum(1 for s in raw_train if s in kept_spike_set)

            acc_wvfs = raw_wvfs[accepted_mask]
            rej_wvfs = raw_wvfs[rejected_mask]

            # Determine verdict string
            if n_accepted_total == 0:
                verdict_str = "DROPPED"
            elif int(unit_id) in curated_spikes:
                verdict_str = f"KEPT {n_accepted_total}/{n_total}"
            else:
                split_ids = [int(unit_id) * 1000 + sc for sc in range(1, 10)
                             if int(unit_id) * 1000 + sc in curated_spikes]
                verdict_str = f"SPLIT -> {split_ids}"

            # Plot
            max_traces = 500
            rng = np.random.default_rng(42)

            def subsample(wvfs):
                if len(wvfs) > max_traces:
                    return wvfs[rng.choice(len(wvfs), max_traces, replace=False)]
                return wvfs

            fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
            ylim_min = np.percentile(raw_wvfs, 0.5)
            ylim_max = np.percentile(raw_wvfs, 99.5)

            panels = [
                (subsample(raw_wvfs), f"All ({n_total})", "gray"),
                (subsample(acc_wvfs) if len(acc_wvfs) > 0 else np.empty((0, raw_wvfs.shape[1])),
                 f"Accepted ({n_accepted_total})", "#2196F3"),
                (subsample(rej_wvfs) if len(rej_wvfs) > 0 else np.empty((0, raw_wvfs.shape[1])),
                 f"Rejected ({n_total - n_accepted_total})", "#F44336"),
            ]

            for ax, (data, title, color) in zip(axes, panels):
                if len(data) > 0:
                    ax.plot(data.T, color=color, alpha=0.05, linewidth=0.5)
                    ax.plot(np.mean(data, axis=0), color="black", linewidth=1.5)
                ax.set_title(title, fontsize=11)
                ax.set_ylim(ylim_min, ylim_max)
                ax.set_xlabel("Samples")
            axes[0].set_ylabel("Amplitude")
            fig.suptitle(f"Unit {unit_id} — {verdict_str}",
                         fontsize=12, fontweight="bold")
            fig.tight_layout()

            fig.savefig(str(diag_img), dpi=150, bbox_inches="tight")
            diag_ppt.add_image(str(diag_img))
            plt.close(fig)
            n_figs += 1

        if n_figs > 0:
            diag_ppt.save(sorting_dir / "curation_diagnostic.pptx")
            print(f"  Saved curation_diagnostic.pptx ({n_figs} units)")
        if diag_img.exists():
            diag_img.unlink()


if __name__ == "__main__":
    control_animals = [
        "ICMS43",
        "ICMS45",
        "ICMS48",
        "ICMS54",
        "ICMS56",
        "ICMS57",
    ]
    for animal in control_animals:
        base = f"E:/ICMS plasticity control/{animal}"
        print(f"\n{'='*60}")
        print(f"  STAGE 2 — {animal}")
        print(f"{'='*60}")
        run_stage2(base_folder=base, force=True)
