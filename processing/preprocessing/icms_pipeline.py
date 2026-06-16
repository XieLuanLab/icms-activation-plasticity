"""
ICMS preprocessing pipeline using stock SpikeInterface + custom artifact removal.

This adapts the proven 'new_pipeline3' from the old SI fork to work with any
SI version by importing the custom preprocessors from this package.

The pipeline assumes stim timestamps are provided in SAMPLES (int).
"""

import numpy as np
import spikeinterface.preprocessing as sp
from pathlib import Path

# Import our custom preprocessors (standalone, no SI fork needed)
from preprocessing.custom_preprocessors import mean_artifact_subtract, trend_subtract


def preprocess_icms(rec, all_stim_timestamps_samples, save_folder=None,
                    stim_freq_hz=100):
    """
    Full ICMS artifact removal + spike extraction preprocessing pipeline.

    Parameters
    ----------
    rec : BaseRecording
        Raw recording (e.g. from se.read_blackrock).
    all_stim_timestamps_samples : list of list of int
        Stim pulse timestamps in SAMPLES, grouped by trial.
        E.g. [[pulse1, pulse2, ...], [pulse1, pulse2, ...], ...]
    save_folder : str or Path, optional
        If provided, saves the pre-whitened recording JSON for waveform extraction.
    stim_freq_hz : float
        Stim pulse frequency. Determines post-stim window for median
        subtraction and trend removal (must fit within inter-pulse gap).
        100 Hz (ICMS83) → 10 ms window, 50 Hz → 20 ms window.

    Returns
    -------
    rec_for_sorting : BaseRecording
        Whitened recording ready for spike sorting.
    rec_for_waveforms : BaseRecording
        Pre-whitened recording for waveform extraction (preserves amplitude).
    """

    # Flatten all pulse timestamps across trials
    flattened_stim = [ts for trial in all_stim_timestamps_samples for ts in trial]

    # Post-stim window scales with inter-pulse interval
    inter_pulse_ms = 1000.0 / stim_freq_hz  # 10 ms @ 100 Hz, 20 ms @ 50 Hz

    # Artifact removal timing parameters (in ms)
    art_remove_pre = 0.5      # ms before stim pulse to blank
    art_remove_post1 = 1.4    # ms after pulse for initial blank & cubic interp
    art_remove_post2 = 1.5    # ms after pulse for final blank (slightly wider)
    trend_remove_start = 1.4  # ms after pulse to start trend fit
    trend_remove_end = inter_pulse_ms  # fit within inter-pulse gap

    # ---- Step 1: Zero-blank the raw stim artifact ----
    # Zeros out the sharp transient so downstream steps don't get confused
    rec1 = sp.remove_artifacts(
        rec, flattened_stim,
        ms_before=art_remove_pre, ms_after=art_remove_post1,
        mode='zeros'
    )

    # ---- Step 2: Median artifact template subtraction ----
    # Computes the median waveform across all pulses in each trial,
    # then subtracts it. Removes the stereotyped artifact component.
    rec2 = mean_artifact_subtract(
        rec1,
        list_triggers=all_stim_timestamps_samples,
        post_stim_window_ms=inter_pulse_ms,
        mode='median'
    )

    # ---- Step 3: Polynomial trend subtraction ----
    # Fits a 3rd-order polynomial to the post-pulse residual (1.4-10 ms)
    # and subtracts. Removes slow exponential/polynomial decay.
    rec3 = trend_subtract(
        rec2,
        list_triggers=all_stim_timestamps_samples,
        post_pulse_start_ms=trend_remove_start,
        post_pulse_end_ms=trend_remove_end,
        mode='poly',
        poly_order=3
    )

    # ---- Step 4: Common median reference ----
    # Removes correlated noise across all channels
    rec4 = sp.common_reference(
        rec3, operator="median", reference="global"
    )

    # ---- Step 5: Cubic interpolation over artifact zone ----
    # Smoothly interpolates through the blanked zone
    rec5 = sp.remove_artifacts(
        rec4, flattened_stim,
        ms_before=art_remove_pre, ms_after=art_remove_post1,
        mode='cubic'
    )

    # ---- Step 6: Bandpass filter for spikes ----
    rec6 = sp.bandpass_filter(rec5, freq_min=300, freq_max=5000)

    # ---- Step 7: Final zero-blank ----
    # Clean final blank (slightly wider) to remove any filter ringing
    rec7 = sp.remove_artifacts(
        rec6, flattened_stim,
        ms_before=art_remove_pre, ms_after=art_remove_post2,
        mode='zeros'
    )

    # No amplitude clipping — matches experimental pipeline.

    # This is the recording to use for waveform extraction (preserves amplitude)
    rec_for_waveforms = rec7

    # ---- Step 8: Whiten for spike sorting ----
    rec_for_sorting = sp.whiten(rec7, dtype='float32')

    # Optionally save the pre-whitened recording for later waveform extraction
    if save_folder is not None:
        save_path = Path(save_folder) / "recording.json"
        rec_for_waveforms.dump(save_path)
        print(f"  Saved recording to: {save_path}")

    return rec_for_sorting, rec_for_waveforms


def stim_timestamps_seconds_to_samples(session_df, fs=30000.0):
    """
    Convert the stim_timestamps column from seconds to samples (int),
    grouped by trial.

    Parameters
    ----------
    session_df : pd.DataFrame
        Must have a 'stim_timestamps' column with lists of float (seconds).
    fs : float
        Sampling frequency in Hz.

    Returns
    -------
    all_stim_timestamps_samples : list of list of int
        Timestamps in samples, one list per trial.
    """
    all_stim = []
    for _, row in session_df.iterrows():
        trial_ts = [int(round(t * fs)) for t in row['stim_timestamps']]
        if len(trial_ts) > 0:
            all_stim.append(trial_ts)
    return all_stim
