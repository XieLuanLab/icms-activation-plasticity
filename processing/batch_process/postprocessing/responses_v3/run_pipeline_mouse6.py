"""
Run the v3 response pipeline on ICMS83 sessions.

ICMS83 uses a different data layout than experimental animals:
- Data on E:/ICMS83/ instead of C:/data/
- trial_df from CSV instead of DataLoader
- Analyzer at batch_sort/stage3/analyzer_final.zarr instead of merge/hmerge_analyzer.zarr

Usage:
    python -m batch_process.postprocessing.responses_v3.run_pipeline_mouse6
    python -m batch_process.postprocessing.responses_v3.run_pipeline_mouse6 E:/ICMS83/20-Jul-2023
"""
import sys
import ast
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import dill as pickle
import spikeinterface as si

from batch_process.postprocessing.responses_v3.session_responses import SessionResponses
from batch_process.postprocessing.responses_v3.unit_response import UnitResponse
from batch_process.postprocessing.responses_v3.stim_condition_response import (
    StimConditionResponse, PulseResponse, TrainResponse)
from batch_process.postprocessing.responses_v3.window_config import WindowConfig
from batch_process.postprocessing.responses_v3.process_stim_responses import process_stim_responses

import batch_process.postprocessing.stim_response_util as stim_response_util
import batch_process.util.file_util as file_util

StimData = stim_response_util.StimData

ALL_SESSIONS = [
    "20-Jul-2023", "24-Jul-2023", "26-Jul-2023", "28-Jul-2023",
    "31-Jul-2023", "02-Aug-2023", "04-Aug-2023", "07-Aug-2023",
    "10-Aug-2023",
]

ICMS83_ROOT = Path("E:/ICMS83")


def load_trial_df(data_folder):
    """Load trial_df.csv and parse list columns."""
    csv_path = Path(data_folder) / "batch_sort" / "trial_df.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"No trial_df.csv at {csv_path}")

    df = pd.read_csv(csv_path)
    for col in ("stim_timestamps", "frame_timestamps"):
        if col in df.columns:
            df[col] = df[col].apply(ast.literal_eval)
    return df


def build_stim_data(trial_df, trial_type="all"):
    """Build StimData from trial_df for a given trial type."""
    if trial_type == "hit":
        mask = (trial_df["current"] > 0) & (trial_df["response"] == 1)
    elif trial_type == "miss":
        mask = (trial_df["current"] > 0) & (trial_df["response"] == 0)
    else:
        mask = pd.Series(True, index=trial_df.index)

    sub = trial_df[mask]
    return StimData(
        sub["stim_timestamps"].tolist(),
        sub["current"].tolist(),
        sub["channel"].tolist())


def run_mouse6_session(data_folder, window_config=None,
                       stim_trial_types=None, force=False):
    """Run v3 pipeline on a single ICMS83 session."""
    if window_config is None:
        window_config = WindowConfig.standard_700ms()
    if stim_trial_types is None:
        stim_trial_types = ['all']

    data_folder = Path(data_folder)
    save_folder = data_folder / "batch_sort"
    analyzer_path = save_folder / "stage3" / "analyzer_final.zarr"

    if not analyzer_path.exists():
        print(f"  [SKIP] No analyzer_final.zarr at {analyzer_path}")
        return

    # Check if already done
    if not force:
        all_done = all(
            (save_folder / f"stim_condition_results_{tt}_{window_config.label}.csv").exists()
            for tt in stim_trial_types)
        if all_done:
            print(f"  [SKIP] Already processed for {window_config.label}")
            return

    job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)
    timing_params = file_util.load_timing_params(
        path=str(Path(__file__).resolve().parents[2] / "postprocessing" / "timing_params.json"))

    print(f"  Loading trial_df.csv...")
    trial_df = load_trial_df(data_folder)

    print(f"  Loading analyzer...")
    analyzer = si.load_sorting_analyzer(analyzer_path)

    if not analyzer.has_extension("unit_locations"):
        print("  Computing extensions...")
        analyzer.compute(
            ["random_spikes", "waveforms", "templates",
             "correlograms", "unit_locations"],
            extension_params={
                "random_spikes": {"method": "all"},
                "correlograms": {"window_ms": 100, "bin_ms": 0.5}},
            **job_kwargs)

    unit_ids = analyzer.unit_ids

    for stim_key in stim_trial_types:
        print(f"\n  --- {stim_key.upper()} [{window_config.label}] ---")

        stim_data = build_stim_data(trial_df, trial_type=stim_key)
        if len(stim_data.timestamps) == 0:
            print(f"  [SKIP] No {stim_key} trials")
            continue

        unique_depths = list(set(stim_data.depths))
        unique_currents = list(set(stim_data.currents))
        unique_stim_params = list(
            itertools.product(unique_depths, unique_currents))

        trial_df_subset = trial_df[[
            'trial', 'response', 'response_time',
            'current', 'channel']].copy()

        session_responses = SessionResponses(
            str(data_folder), str(analyzer_path), trial_df_subset,
            stim_data.timestamps, timing_params,
            window_config=window_config)

        for unit_id in unit_ids:
            spike_timestamps = analyzer.sorting.get_unit_spike_train(
                unit_id=unit_id)
            unit_response = UnitResponse(
                unit_id, spike_timestamps, session_responses)

            print(f"  [{stim_key.upper()}] Processing unit {unit_id}")

            for depth, current in unique_stim_params:
                stim_timestamps = stim_response_util.get_stim_ts_indices(
                    stim_data.timestamps, stim_data.currents,
                    stim_data.depths, current, depth)

                if len(stim_timestamps) > 0:
                    scr = StimConditionResponse(
                        stim_timestamps, depth, current, unit_response)
                    PulseResponse(
                        timing_params=timing_params,
                        stim_condition_response=scr)
                    TrainResponse(
                        timing_params=timing_params,
                        stim_condition_response=scr,
                        window_config=window_config)

                    if current > 0 and not window_config.symmetric:
                        scr.check_spike_counts_match()

        # Save session responses
        pkl_path = save_folder / f"session_responses_{stim_key}_{window_config.label}.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(session_responses, f)
        print(f"  Saved {pkl_path}")

        # Extract metrics to CSV
        session_responses._sorting_analyzer = analyzer
        process_stim_responses(session_responses, stim_key,
                               window_label=window_config.label)


def run_all_mouse6(window_configs=None, stim_trial_types=None, force=False):
    """Run v3 pipeline on all ICMS83 sessions."""
    if window_configs is None:
        window_configs = [WindowConfig.standard_700ms()]
    if stim_trial_types is None:
        stim_trial_types = ['all']

    for session_name in ALL_SESSIONS:
        data_folder = ICMS83_ROOT / session_name
        if not data_folder.exists():
            print(f"[SKIP] {data_folder} does not exist")
            continue

        for window_config in window_configs:
            print(f"\n{'='*60}")
            print(f"ICMS83 {session_name} [{window_config.label}]")
            print(f"{'='*60}")

            run_mouse6_session(
                data_folder, window_config=window_config,
                stim_trial_types=stim_trial_types, force=force)

    print("\nAll ICMS83 sessions complete.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and not sys.argv[1].startswith("--"):
        run_mouse6_session(
            sys.argv[1],
            window_config=WindowConfig.standard_700ms(),
            stim_trial_types=['all'],
            force="--force" in sys.argv)
    else:
        run_all_mouse6(
            stim_trial_types=['all'],
            force="--force" in sys.argv)
