"""
Control-group session response creation using the shared response classes.

Loads control data from condensed_trials.csv + analyzer_final.zarr,
then constructs the same SessionResponses → UnitResponse → StimConditionResponse
hierarchy used by the experimental pipeline.

Differences from experimental:
- No behavioral component (no hit/miss separation) — all trials used
- 50 Hz stim (20ms pulse window) vs 100 Hz (10ms)
- Data loaded from condensed_trials.csv instead of DataLoader + .mat
- Stim conditions keyed by (channel, current) instead of (depth, current)
"""

import ast
import json
import itertools
import numpy as np
import pandas as pd
import dill as pickle
from pathlib import Path
from collections import namedtuple

from spikeinterface import full as si

# Reuse the experimental v2 response class hierarchy
from batch_process.postprocessing.responses_v2.session_responses_v2 import SessionResponses
from batch_process.postprocessing.responses_v2.unit_response_v2 import UnitResponse
from batch_process.postprocessing.responses_v2.stim_condition_response_v2 import (
    StimConditionResponse, PulseResponse, TrainResponse,
)
import batch_process.postprocessing.stim_response_util as stim_response_util

# StimData matches the experimental namedtuple interface
StimData = namedtuple("StimData", ["timestamps", "currents", "depths"])

TIMING_PARAMS_PATH = Path(__file__).parent / "control_timing_params.json"


def load_timing_params(path=None):
    """Load control timing parameters from JSON."""
    path = path or TIMING_PARAMS_PATH
    with open(path) as f:
        return json.load(f)


def load_session_data(session_dir):
    """Load analyzer and trial CSV from a control session directory.

    Returns:
        (df, analyzer) where df has columns: channel, current_uA, stim_timestamps
    """
    spike_dir = Path(session_dir) / "spike_sorting"

    csv_path = spike_dir / "condensed_trials.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"No condensed_trials.csv in {spike_dir}")

    df = pd.read_csv(csv_path)
    df['stim_timestamps'] = df['stim_timestamps'].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Prefer finalized analyzer
    analyzer_path = spike_dir / "analyzer_final.zarr"
    if not analyzer_path.exists():
        analyzer_path = spike_dir / "analyzer_raw.zarr"
    if not analyzer_path.exists():
        raise FileNotFoundError(f"No analyzer found in {spike_dir}")

    analyzer = si.load_sorting_analyzer(analyzer_path)
    return df, analyzer


def build_stim_data_from_csv(df, fs=30000):
    """Convert control CSV into a StimData namedtuple matching experimental format.

    Each row in the CSV is one trial with (channel, current_uA, stim_timestamps).
    We map channel → "depth" to match the experimental interface where stim conditions
    are keyed by (depth, current).

    Returns:
        StimData with per-trial timestamps (in samples), currents, and depths (=channels).
    """
    all_timestamps = []
    all_currents = []
    all_depths = []  # For control, "depth" = stim channel ID

    for _, row in df.iterrows():
        ch = int(row['channel'])
        cur = int(row['current_uA'])
        ts = (np.array(row['stim_timestamps']) * fs).astype(np.int64)

        all_timestamps.append(ts)
        all_currents.append(cur)
        all_depths.append(ch)

    return StimData(all_timestamps, all_currents, all_depths)


def ensure_extensions(analyzer, job_kwargs=None):
    """Compute required extensions if missing."""
    job_kwargs = job_kwargs or dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

    required = ["random_spikes", "waveforms", "templates", "correlograms", "unit_locations"]
    extension_params = {
        "random_spikes": {"method": "all"},
        "correlograms": {"window_ms": 100, "bin_ms": 0.5},
    }

    missing = [ext for ext in required if not analyzer.has_extension(ext)]
    if missing:
        print(f"  Computing missing extensions: {missing}")
        analyzer.compute(missing, extension_params=extension_params, **job_kwargs)


def create_session_responses(session_dir, timing_params=None):
    """Create a SessionResponses object for a control session.

    This is the main entry point. It:
    1. Loads the analyzer and trial data
    2. Builds StimData from the CSV
    3. Constructs SessionResponses with all trials (no hit/miss)
    4. Processes each unit × stim condition

    Returns:
        SessionResponses object (same class as experimental pipeline)
    """
    timing_params = timing_params or load_timing_params()
    df, analyzer = load_session_data(session_dir)

    ensure_extensions(analyzer)

    stim_data = build_stim_data_from_csv(df, fs=timing_params['fs'])

    # Build unique stim conditions: (channel, current) pairs
    unique_depths = list(set(stim_data.depths))
    unique_currents = list(set(stim_data.currents))
    unique_stim_params = list(itertools.product(unique_depths, unique_currents))

    # Create session responses — no hit/miss for control
    spike_dir = Path(session_dir) / "spike_sorting"
    analyzer_path = spike_dir / "analyzer_final.zarr"
    if not analyzer_path.exists():
        analyzer_path = spike_dir / "analyzer_raw.zarr"

    session_responses = SessionResponses(
        session_path=str(session_dir),
        sorting_analyzer_path=str(analyzer_path),
        trial_df=df,
        all_stim_timestamps=stim_data.timestamps,
        timing_params=timing_params,
    )
    # Pre-load the analyzer so it's available for processing
    session_responses._sorting_analyzer = analyzer

    # Process each unit
    for unit_id in analyzer.unit_ids:
        print(f"  Processing unit {unit_id}...")
        spike_timestamps = analyzer.sorting.get_unit_spike_train(unit_id)
        unit_response = UnitResponse(unit_id, spike_timestamps, session_responses)

        for stim_params in unique_stim_params:
            depth, current = stim_params
            stim_timestamps = stim_response_util.get_stim_ts_indices(
                stim_data.timestamps, stim_data.currents, stim_data.depths,
                current, depth
            )

            if len(stim_timestamps) > 0:
                scr = StimConditionResponse(
                    stim_timestamps, depth, current, unit_response)
                PulseResponse(timing_params, scr)
                TrainResponse(timing_params, scr)

    return session_responses


def save_session_responses(session_responses, save_dir, suffix=''):
    """Save SessionResponses to pickle."""
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    tag = f"_{suffix}" if suffix else ""
    save_path = save_dir / f"session_responses{tag}.pkl"
    with open(save_path, "wb") as f:
        pickle.dump(session_responses, f)
    print(f"  Saved to {save_path}")
    return save_path


def main(session_dirs=None, base_folder=None):
    """Batch process control sessions.

    Args:
        session_dirs: List of session directory paths. If None, discovers from base_folder.
        base_folder: Root folder to discover sessions from (used if session_dirs is None).
    """
    if session_dirs is None:
        if base_folder is None:
            raise ValueError("Provide either session_dirs or base_folder")
        from util.load_control import discover_sessions
        session_dirs = discover_sessions(base_folder)

    timing_params = load_timing_params()

    for session_dir in session_dirs:
        print(f"\n{'='*60}")
        print(f"Processing: {session_dir}")

        spike_dir = Path(session_dir) / "spike_sorting"
        if not spike_dir.exists():
            print("  -> No spike_sorting directory. Skipping.")
            continue

        try:
            sr = create_session_responses(session_dir, timing_params)
            save_session_responses(sr, spike_dir)
        except FileNotFoundError as e:
            print(f"  -> ERROR: {e}. Skipping.")
        except Exception as e:
            print(f"  -> ERROR: {e}")
            raise


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main(session_dirs=[sys.argv[1]])
    else:
        print("Usage: python create_session_responses.py <session_dir>")
        print("   or: import and call main(base_folder=...)")
