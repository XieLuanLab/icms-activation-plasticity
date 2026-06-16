"""
Orchestration: create SessionResponses objects and save to disk.

Changes from v2:
- Accepts WindowConfig parameter
- Uses v3 response classes with keyword arguments
- Saves with window label in filename
"""
from spikeinterface import full as si
import itertools
import numpy as np
import dill as pickle
import pandas as pd
from pathlib import Path

from batch_process.postprocessing.responses_v3.session_responses import SessionResponses
from batch_process.postprocessing.responses_v3.unit_response import UnitResponse
from batch_process.postprocessing.responses_v3.stim_condition_response import (
    StimConditionResponse, PulseResponse, TrainResponse)
from batch_process.postprocessing.responses_v3.window_config import WindowConfig

import batch_process.postprocessing.stim_response_util as stim_response_util
import batch_process.util.file_util as file_util
from util.load_data import load_data


def process_unit(unit_id, unique_stim_params, stim_data, timing_params,
                 unit_response, window_config):
    """Process a single unit's response to the stimulus."""
    for stim_params in unique_stim_params:
        depth, current = stim_params
        stim_timestamps = stim_response_util.get_stim_ts_indices(
            stim_data.timestamps, stim_data.currents, stim_data.depths,
            current, depth)

        if len(stim_timestamps) > 0:
            stim_condition_response = StimConditionResponse(
                stim_timestamps, depth, current, unit_response)

            pulse_response = PulseResponse(
                timing_params=timing_params,
                stim_condition_response=stim_condition_response)
            train_response = TrainResponse(
                timing_params=timing_params,
                stim_condition_response=stim_condition_response,
                window_config=window_config)

            if current > 0 and not window_config.symmetric:
                stim_condition_response.check_spike_counts_match()


def get_session_responses(data_folder, window_config=None,
                          stim_trial_types=None, server_mount_drive="S:"):
    """
    Build and save SessionResponses for a single data folder.

    Args:
        data_folder: Path to the session data folder.
        window_config: WindowConfig specifying analysis windows.
        stim_trial_types: List of trial types to process (e.g., ['all', 'hit', 'miss']).
        server_mount_drive: Drive letter for server mount.

    Returns:
        dict: Mapping of stim_trial_type -> SessionResponses for the last
              processed folder. Also saves pkl files to disk.
    """
    if window_config is None:
        window_config = WindowConfig.standard_700ms()
    if stim_trial_types is None:
        stim_trial_types = ['all']

    job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

    timing_params = file_util.load_timing_params(
        path=str(Path(__file__).resolve().parents[2] / "postprocessing" / "timing_params.json"))
    save_folder = Path(data_folder) / "batch_sort"

    dataloader = load_data(
        data_folder, make_folder=False, save_folder_name="batch_sort",
        first_N_files=4, server_mount_drive=server_mount_drive)
    trial_df = dataloader.trial_df[[
        'trial', 'response', 'response_time', 'trial_start', 'trial_end',
        'current', 'channel']]

    analyzer = si.load_sorting_analyzer(
        Path(save_folder) / "merge/hmerge_analyzer.zarr")
    sorting_analyzer_path = Path(data_folder) / \
        "batch_sort/merge/hmerge_analyzer.zarr"

    stim_data_all, hit_stim_data, miss_stim_data = stim_response_util.get_stim_data(
        data_folder)
    if hit_stim_data is None:
        print(f"Skipping session {data_folder}...")
        return {}

    unit_ids = analyzer.unit_ids

    if not analyzer.has_extension("unit_locations"):
        extensions_to_compute = [
            "random_spikes", "waveforms", "templates",
            "correlograms", "unit_locations"]
        extension_params = {
            "random_spikes": {"method": "all"},
            "correlograms": {"window_ms": 100, "bin_ms": 0.5}}
        analyzer.compute(extensions_to_compute,
                         extension_params=extension_params, **job_kwargs)

    stim_data_map = {
        "all": stim_data_all,
        "hit": hit_stim_data,
        "miss": miss_stim_data
    }

    results = {}
    for stim_key in stim_trial_types:
        stim_data = stim_data_map[stim_key]
        unique_depths = list(set(stim_data.depths))
        unique_currents = list(set(stim_data.currents))
        unique_stim_params = list(
            itertools.product(unique_depths, unique_currents))

        session_responses = SessionResponses(
            data_folder, sorting_analyzer_path, trial_df,
            stim_data.timestamps, timing_params, window_config=window_config)

        for unit_id in unit_ids:
            spike_timestamps = analyzer.sorting.get_unit_spike_train(
                unit_id=unit_id)
            unit_response = UnitResponse(
                unit_id, spike_timestamps, session_responses)
            print(f"[{stim_key.upper()}] Processing unit id: {unit_id}")
            process_unit(unit_id, unique_stim_params, stim_data,
                         timing_params, unit_response, window_config)

        # Save with window label
        save_path = save_folder / f'session_responses_{stim_key}_{window_config.label}.pkl'
        with open(save_path, "wb") as file:
            pickle.dump(session_responses, file)
        print(f"Session responses [{stim_key}] saved to {save_path}")
        results[stim_key] = session_responses

    return results
