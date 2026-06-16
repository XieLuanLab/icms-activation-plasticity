"""
Stage 4: build session responses + extract modulation / pulse-locking metrics.

Generates session_responses (pkl), stim_condition_results (CSV), and per-unit
figures from batch_sort/stage3/analyzer_final.zarr + trial_df.csv.

Usage:
    python -m batch_process.stage4_postprocess <session_folder> [<session_folder> ...] [--force]
"""

import sys
import os
import ast
import time
import traceback
import itertools
import shutil
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dill as pickle
import spikeinterface as si

from batch_process.postprocessing.responses_v2.session_responses_v2 import SessionResponses
from batch_process.postprocessing.responses_v2.unit_response_v2 import UnitResponse
from batch_process.postprocessing.responses_v2.stim_condition_response_v2 import (
    StimConditionResponse, PulseResponse, TrainResponse,
)
import batch_process.postprocessing.stim_response_util as stim_response_util
import batch_process.postprocessing.responses_v2.response_plotting_util as rpu
import batch_process.util.file_util as file_util



StimData = stim_response_util.StimData


def load_trial_df(data_folder):
    """Load trial_df.csv and parse list columns from string representation."""
    csv_path = Path(data_folder) / "batch_sort" / "trial_df.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"No trial_df.csv at {csv_path}")

    df = pd.read_csv(csv_path)

    # Parse stringified lists back to actual lists
    for col in ("stim_timestamps", "frame_timestamps"):
        if col in df.columns:
            df[col] = df[col].apply(ast.literal_eval)

    return df


def build_stim_data(trial_df, trial_type="all"):
    """Build a StimData namedtuple from trial_df for a given trial type.

    trial_type: 'all', 'hit', or 'miss'
    """
    if trial_type == "hit":
        mask = (trial_df["current"] > 0) & (trial_df["response"] == 1)
    elif trial_type == "miss":
        mask = (trial_df["current"] > 0) & (trial_df["response"] == 0)
    else:  # all
        mask = pd.Series(True, index=trial_df.index)

    sub = trial_df[mask]

    timestamps = sub["stim_timestamps"].tolist()
    currents = sub["current"].tolist()
    depths = sub["channel"].tolist()

    return StimData(timestamps, currents, depths)


def store_result(stim_response, results, pre_stim_threshold=0.5):
    """Extract metrics from a single stim condition response."""
    pr = stim_response.pulse_response
    tr = stim_response.train_response

    session_path = stim_response.session_responses.session_path
    animal_id = file_util.get_animal_id(session_path)
    session = file_util.get_date_str(session_path)

    result = {
        "animal_id": animal_id,
        "session": session,
        "rel_day": None,
        "rel_week": None,
        "unit_id": stim_response.unit_id,
        "stim_channel": stim_response.stim_channel,
        "stim_current": stim_response.stim_current,
        "unit_location": stim_response.unit_response.unit_location,
        "cell_type": stim_response.unit_response.cell_type,
        "include": True,
    }

    if pr and tr:
        result.update({
            "is_pulse_locked": pr.is_pulse_locked,
            "max_spike_prob": np.max(pr.stim_metrics["mean_prob_spike"]),
            "pli": pr.pli,
            "latency": pr.stim_metrics["mean_latency"],
            "latency_std": pr.stim_metrics["std_latency"],
            "jitter": pr.stim_metrics["std_latency"]**2,

            "modulated": tr.paired_p_val < 0.05,
            "mod_p_val": tr.paired_p_val,
            "z_score": tr.z_score,

            "train_mean_fr": tr.stim_mean_fr,
            "pre_stim_mean_fr": tr.pre_stim_mean_fr,
            "pre_stim_std_fr": tr.pre_stim_sigma_fr,

            "max_fr_10ms_smoothed": np.max(tr.ten_ms_smoothed_data['firing_rate']) if tr.ten_ms_smoothed_data['firing_rate'].size > 0 else np.nan,
            "t_to_crossing_10ms_smoothed": tr.ten_ms_smoothed_data.get('time_to_threshold_crossing', np.nan),
            "t_to_max_10ms_smoothed": tr.ten_ms_smoothed_data.get('time_to_max_fr', np.nan),

            "max_fr_2ms_smoothed": np.max(tr.two_ms_smoothed_data['firing_rate']) if tr.two_ms_smoothed_data['firing_rate'].size > 0 else np.nan,
            "t_to_crossing_2ms_smoothed": tr.two_ms_smoothed_data.get('time_to_threshold_crossing', np.nan),
            "t_to_max_2ms_smoothed": tr.two_ms_smoothed_data.get('time_to_max_fr', np.nan),

            "t_to_crossing_2ms_binned": tr.two_ms_binned_data.get('time_to_threshold_crossing', np.nan),
            "t_to_max_2ms_binned": tr.two_ms_binned_data.get('time_to_max_ms', np.nan),

            "train_binned_fr_2ms": tr.two_ms_binned_data.get('binned_firing_rate', np.nan),
            "train_smoothed_fr_2ms": tr.two_ms_smoothed_data.get('firing_rate', np.nan),
            "train_smoothed_fr_10ms": tr.ten_ms_smoothed_data['firing_rate'],

            "num_spikes": len(pr.rel_spike_timestamps),
            "num_trials": len(tr.raster_array),
            "baseline_too_slow": tr.pre_stim_mean_fr < pre_stim_threshold,
        })

    results.append(result)


def extract_metrics_to_csv(session_responses, stim_trial_type):
    """Extract modulation metrics from session responses and save to CSV."""
    results = []

    for unit_id in session_responses.unit_ids:
        unit_response = session_responses.get_unit_response(unit_id)
        stim_channels = sorted(
            {ch for ch, _ in unit_response.stim_responses.keys()})

        for stim_channel in stim_channels:
            if stim_channel == 0:
                continue
            currents = sorted(
                [curr for ch, curr in unit_response.stim_responses.keys()
                 if ch == stim_channel])

            for stim_current in currents:
                if stim_current == 0:
                    continue
                stim_response = unit_response.get_stim_response(
                    stim_channel=stim_channel, stim_current=stim_current)
                store_result(stim_response, results)

    save_folder = Path(session_responses.session_path) / "batch_sort"
    rpu.save_results_to_csv(results, save_folder, stim_trial_type)


def process_session(data_folder, job_kwargs, force=False,
                    stim_trial_types=("all",)):
    """Run stage 4 postprocessing on a single session."""
    data_folder = Path(data_folder)
    save_folder = data_folder / "batch_sort"
    stage3_dir = save_folder / "stage3"
    analyzer_path = stage3_dir / "analyzer_final.zarr"

    # Check prerequisites
    if not analyzer_path.exists():
        print(f"  [SKIP] No analyzer_final.zarr — run stage3 first")
        return False

    # Check if already done (CSV exists for each trial type)
    if not force:
        all_done = all(
            (save_folder / f"stim_condition_results_{tt}.csv").exists()
            for tt in stim_trial_types
        )
        if all_done:
            print(f"  [SKIP] Already processed (CSVs exist)")
            return True

    # Load trial data
    print("  Loading trial_df.csv...")
    trial_df = load_trial_df(data_folder)

    # Load timing params
    timing_params = file_util.load_timing_params(
        path="batch_process/postprocessing/timing_params.json")

    # Load analyzer
    print("  Loading analyzer_final.zarr...")
    analyzer = si.load_sorting_analyzer(analyzer_path)

    # Compute extensions if needed
    if not analyzer.has_extension("unit_locations"):
        print("  Computing extensions...")
        analyzer.compute(
            ["random_spikes", "waveforms", "templates",
             "correlograms", "unit_locations"],
            extension_params={
                "random_spikes": {"method": "all"},
                "correlograms": {"window_ms": 100, "bin_ms": 0.5},
            },
            **job_kwargs,
        )

    unit_ids = analyzer.unit_ids

    # Build stim data for each trial type
    for stim_key in stim_trial_types:
        print(f"\n  --- Processing trial type: {stim_key.upper()} ---")

        stim_data = build_stim_data(trial_df, trial_type=stim_key)
        if len(stim_data.timestamps) == 0:
            print(f"  [SKIP] No {stim_key} trials")
            continue

        unique_depths = list(set(stim_data.depths))
        unique_currents = list(set(stim_data.currents))
        unique_stim_params = list(
            itertools.product(unique_depths, unique_currents))

        # Subset of trial_df for SessionResponses
        trial_df_subset = trial_df[[
            'trial', 'response', 'response_time',
            'current', 'channel']].copy()

        # Initialize session responses
        session_responses = SessionResponses(
            str(data_folder),
            str(analyzer_path),
            trial_df_subset,
            stim_data.timestamps,
            timing_params,
        )

        # Process each unit
        for unit_id in unit_ids:
            spike_timestamps = analyzer.sorting.get_unit_spike_train(
                unit_id=unit_id)
            unit_response = UnitResponse(
                unit_id, spike_timestamps, session_responses)

            print(f"  [{stim_key.upper()}] Processing unit {unit_id}")

            for stim_params in unique_stim_params:
                depth, current = stim_params
                stim_timestamps = stim_response_util.get_stim_ts_indices(
                    stim_data.timestamps, stim_data.currents,
                    stim_data.depths, current, depth)

                if len(stim_timestamps) > 0:
                    stim_condition_response = StimConditionResponse(
                        stim_timestamps, depth, current, unit_response)
                    pulse_response = PulseResponse(
                        timing_params, stim_condition_response)
                    train_response = TrainResponse(
                        timing_params, stim_condition_response)

                    if current > 0:
                        stim_condition_response.check_spike_counts_match()

        # Save session responses pkl
        pkl_path = save_folder / f"session_responses_{stim_key}.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(session_responses, f)
        print(f"  Saved {pkl_path}")

        # Extract metrics to CSV
        print(f"  Extracting metrics...")
        session_responses._sorting_analyzer = analyzer
        extract_metrics_to_csv(session_responses, stim_key)

    # Generate figures (using 'all' trial type)
    all_pkl = save_folder / "session_responses_all.pkl"
    if all_pkl.exists():
        print("\n  Generating stim response figures...")
        with open(all_pkl, "rb") as f:
            session_responses = pickle.load(f)
        session_responses._sorting_analyzer = analyzer

        figures_dir = save_folder / "figures" / "stim_condition_figures"
        if figures_dir.exists():
            shutil.rmtree(figures_dir)
        figures_dir.mkdir(parents=True, exist_ok=True)

        for unit_id in session_responses.unit_ids:
            ur = session_responses.get_unit_response(unit_id)
            stim_conditions = ur.show_stim_conditions()
            responses_by_channel = defaultdict(list)

            for stim_condition in stim_conditions:
                stim_ch = stim_condition[0]
                current = stim_condition[1]
                sr = ur.get_stim_response(stim_ch, current)
                responses_by_channel[stim_ch].append((current, sr))

            for stim_ch, responses in responses_by_channel.items():
                fig, axs = rpu.initialize_plotting()
                rpu.plot_train_rasters(responses, axs)
                rpu.plot_train_firing_rate(responses, axs)
                rpu.plot_pulse_raster(responses, axs)
                rpu.plot_probability(responses, axs)
                rpu.plot_z_scores(responses, axs)
                rpu.plot_template(session_responses, unit_id, axs)

                plt.suptitle(f"Unit {unit_id} - Stim Channel D{stim_ch}")
                plt.tight_layout()
                rpu.save_figure(fig, figures_dir, unit_id, stim_ch)

        # Ensure analysis dir exists for PPT save
        analysis_dir = data_folder.parent / "analysis"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        (save_folder / "figures").mkdir(parents=True, exist_ok=True)

        rpu.make_ppt(str(data_folder), figures_dir=figures_dir)

    print(f"\n  Stage 4 complete!")
    return True


def main():
    job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

    force = "--force" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--force"]

    sessions = args
    if not sessions:
        print("Usage: python -m batch_process.stage4_postprocess <session_folder> [<session_folder> ...]")
        return

    # Only process sessions that have analyzer_final.zarr
    sessions = [s for s in sessions
                if (Path(s) / "batch_sort" / "stage3" / "analyzer_final.zarr").exists()]

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
