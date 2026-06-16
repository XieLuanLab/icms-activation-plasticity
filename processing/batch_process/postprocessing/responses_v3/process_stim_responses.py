"""
Post-processing: extract features from SessionResponses into CSV/pkl.

Changes from v2:
- Output filenames include window label
- Reuses response_plotting_util from v2 for save logic
"""
import numpy as np
import pandas as pd
from pathlib import Path
import dill as pickle

import batch_process.util.file_util as file_util


def store_result(stim_response, results, include=True, pre_stim_threshold=0.5):
    """Extract features from a single stim condition response."""
    pr, tr = stim_response.pulse_response, stim_response.train_response

    animal_id = file_util.get_animal_id(
        stim_response.session_responses.session_path)
    session = file_util.get_date_str(
        stim_response.session_responses.session_path)

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
        "include": include,
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


def process_stim_responses(session_responses, stim_trial_type, window_label=None):
    """
    Extract features from session_responses and save to CSV/pkl.

    Args:
        session_responses: SessionResponses object.
        stim_trial_type: Trial type string (e.g., 'all', 'hit', 'miss').
        window_label: Optional window label for filename. If None, uses
            session_responses.window_config.label.
    """
    if window_label is None:
        window_label = getattr(
            getattr(session_responses, '_window_config', None),
            'label', '700ms')

    results = []

    for unit_id in session_responses.unit_ids:
        print(f"Processing unit id {unit_id}...")
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

    # Save
    save_folder = Path(session_responses.session_path) / "batch_sort"
    csv_path = save_folder / f"stim_condition_results_{stim_trial_type}_{window_label}.csv"
    pkl_path = save_folder / f"stim_condition_results_{stim_trial_type}_{window_label}.pkl"

    df = pd.DataFrame(results)
    all_columns = set(col for result in results for col in result.keys())
    for col in all_columns:
        if col not in df.columns:
            df[col] = np.nan

    numeric_columns = df.select_dtypes(include=['number']).columns.tolist()
    df[numeric_columns] = df[numeric_columns].apply(
        pd.to_numeric, errors='coerce')

    df.to_csv(csv_path, index=False)
    df.to_pickle(pkl_path)
    print(f"Results saved to {csv_path}")

    return df
