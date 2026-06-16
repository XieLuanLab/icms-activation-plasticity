"""
Custom SpikeInterface-compatible preprocessing modules for ICMS artifact removal.

These are standalone versions of the custom preprocessors originally added to
a local SpikeInterface fork. They work with any SI version (0.100+) as external
classes by importing the base classes directly.

Contains:
  - MeanArtifactSubtractedRecording / mean_artifact_subtract
      Computes the mean/median waveform across all stim pulses in a trial,
      then subtracts it from each pulse window.

  - TrendSubtractedRecording / trend_subtract
      Fits a polynomial or spline to the post-stimulus window for each pulse
      and subtracts the slow drift/trend.
"""

import numpy as np
from spikeinterface.preprocessing.basepreprocessor import BasePreprocessor, BasePreprocessorSegment
from spikeinterface.core.core_tools import define_function_from_class


# =============================================================================
# Mean Artifact Subtract
# =============================================================================

class MeanArtifactSubtractedRecording(BasePreprocessor):
    """
    Subtract the mean or median of post-stimulus pulse windows in a trial.

    For each trial (list of pulse timestamps), compute the average waveform
    across all pulses within the post-stimulus window, then subtract that
    average from every individual pulse window. This removes the stereotyped
    artifact component while preserving unique neural activity.

    Parameters
    ----------
    recording : RecordingExtractor
        The recording extractor to be processed.
    list_triggers : list of list of int
        List of lists containing pulse times IN SAMPLES for each trial.
    post_stim_window_ms : float
        Post-stimulus window in milliseconds for averaging.
    mode : str, optional
        'mean' or 'median', by default 'mean'.
    """

    name = "mean_artifact_subtract"
    installed = True

    def __init__(self, recording, list_triggers, post_stim_window_ms, mode="mean"):
        self._list_triggers = list_triggers
        self._post_stim_window_ms = post_stim_window_ms
        self._mode = mode
        self._sf = recording.get_sampling_frequency()

        BasePreprocessor.__init__(self, recording)
        for parent_segment in recording._recording_segments:
            self.add_recording_segment(
                MeanArtifactSubtractedRecordingSegment(
                    parent_segment, list_triggers, self._post_stim_window_ms, mode
                )
            )

        self._kwargs = {
            "recording": recording,
            "list_triggers": list_triggers,
            "post_stim_window_ms": post_stim_window_ms,
            "mode": mode,
        }


class MeanArtifactSubtractedRecordingSegment(BasePreprocessorSegment):
    def __init__(self, parent_recording_segment, list_triggers, post_stim_window_ms, mode):
        BasePreprocessorSegment.__init__(self, parent_recording_segment)
        self._list_triggers = list_triggers
        self._post_stim_window_samples = int(post_stim_window_ms * 30)  # 30 kHz
        self._mode = mode

    def get_traces(self, start_frame, end_frame, channel_indices):
        traces = self.parent_recording_segment.get_traces(start_frame, end_frame, channel_indices)
        traces = traces.copy()  # avoid modifying upstream data
        
        if not hasattr(self, '_trial_templates'):
            self._trial_templates = {}

        for trial_idx, trial_triggers in enumerate(self._list_triggers):
            # Check if any triggers fall in this chunk before doing expensive work
            triggers_in_chunk = [ts for ts in trial_triggers if start_frame <= ts < end_frame]
            if not triggers_in_chunk:
                continue

            # If the template for this trial hasn't been computed yet, compute it!
            if trial_idx not in self._trial_templates:
                post_stim_traces = []
                for ts in trial_triggers:
                    # Fetch unconditionally from parent across chunk boundaries!
                    window = self.parent_recording_segment.get_traces(
                        ts, ts + self._post_stim_window_samples, channel_indices
                    )
                    if window.shape[0] == self._post_stim_window_samples:
                        post_stim_traces.append(window)
                        
                if len(post_stim_traces) > 0:
                    post_stim_traces = np.array(post_stim_traces)
                    if self._mode == "mean":
                        self._trial_templates[trial_idx] = np.mean(post_stim_traces, axis=0)
                    else:  # median
                        self._trial_templates[trial_idx] = np.median(post_stim_traces, axis=0)
                else:
                    self._trial_templates[trial_idx] = None

            avg_trace = self._trial_templates.get(trial_idx, None)
            if avg_trace is None:
                continue

            # Subtract the average from each pulse window in THIS chunk
            for ts in triggers_in_chunk:
                ts_rel = ts - start_frame
                end_ts_rel = ts_rel + self._post_stim_window_samples
                
                # Handle edge cases where window overlaps chunk boundary
                if end_ts_rel > traces.shape[0]:
                    valid_len = traces.shape[0] - ts_rel
                    traces[ts_rel:, :] -= avg_trace[:valid_len].astype(traces.dtype)
                else:
                    traces[ts_rel:end_ts_rel, :] -= avg_trace.astype(traces.dtype)

        return traces


# API function
mean_artifact_subtract = define_function_from_class(
    source_class=MeanArtifactSubtractedRecording, name="mean_artifact_subtract"
)


# =============================================================================
# Trend Subtract
# =============================================================================

class TrendSubtractedRecording(BasePreprocessor):
    """
    Subtract a polynomial or spline fit of post-stimulus windows.

    For each pulse, fits a polynomial (or spline) to the slow trend in the
    post-stimulus period and subtracts it. This removes the slow exponential
    or polynomial decay artifact that persists after blanking.

    Parameters
    ----------
    recording : RecordingExtractor
        The recording extractor to be processed.
    list_triggers : list of list of int
        List of lists containing pulse times IN SAMPLES for each trial.
    post_pulse_start_ms : float
        Time in ms after each pulse to start the fit window.
    post_pulse_end_ms : float
        Time in ms after each pulse to end the fit window.
    mode : str, optional
        'poly' or 'spline', by default 'poly'.
    poly_order : int, optional
        Polynomial order, by default 3.
    spline_order : int, optional
        Spline order, by default 3.
    dspline : int, optional
        Samples between spline knots, by default 1000.
    """

    name = "trend_subtract"
    installed = True

    def __init__(
        self,
        recording,
        list_triggers,
        post_pulse_start_ms,
        post_pulse_end_ms=10.0,
        mode="poly",
        poly_order=3,
        spline_order=3,
        dspline=1000,
    ):
        assert mode in ["poly", "spline"], "Mode must be either 'poly' or 'spline'"

        self._list_triggers = list_triggers
        self._sf = recording.get_sampling_frequency()

        BasePreprocessor.__init__(self, recording)
        for parent_segment in recording._recording_segments:
            recording_segment = TrendSubtractedRecordingSegment(
                parent_segment,
                list_triggers,
                post_pulse_start_ms,
                post_pulse_end_ms,
                mode,
                poly_order,
                spline_order,
                dspline,
            )
            self.add_recording_segment(recording_segment)

        self._kwargs = {
            "recording": recording,
            "list_triggers": list_triggers,
            "post_pulse_start_ms": post_pulse_start_ms,
            "post_pulse_end_ms": post_pulse_end_ms,
            "mode": mode,
            "poly_order": poly_order,
            "spline_order": spline_order,
            "dspline": dspline,
        }


class TrendSubtractedRecordingSegment(BasePreprocessorSegment):
    def __init__(
        self,
        parent_recording_segment,
        list_triggers,
        post_pulse_start_ms,
        post_pulse_end_ms,
        mode,
        poly_order,
        spline_order,
        dspline,
    ):
        BasePreprocessorSegment.__init__(self, parent_recording_segment)
        self._list_triggers = [[np.int64(ts) for ts in trial] for trial in list_triggers]
        self._sf = parent_recording_segment.sampling_frequency
        self._post_pulse_start_samples = int(self._sf * (post_pulse_start_ms / 1000.0))
        self._post_pulse_end_samples = int(self._sf * (post_pulse_end_ms / 1000.0))
        self._mode = mode
        self._poly_order = poly_order
        self._spline_order = spline_order
        self._dspline = dspline

    def get_traces(self, start_frame, end_frame, channel_indices):
        import warnings

        traces = self.parent_recording_segment.get_traces(start_frame, end_frame, channel_indices)
        traces = traces.copy()

        num_channels = traces.shape[1]

        for trial_triggers in self._list_triggers:
            for ts in trial_triggers:
                start_pulse = ts + self._post_pulse_start_samples
                end_pulse = ts + self._post_pulse_end_samples
                
                # Check if this pulse's trend window overlaps with the current chunk
                overlap_start = max(start_pulse, start_frame)
                overlap_end = min(end_pulse, end_frame)
                
                if overlap_start < overlap_end:  # Overlaps this chunk
                    # Fetch FULL window from parent so polyfit is never truncated
                    full_window = self.parent_recording_segment.get_traces(
                        start_pulse, end_pulse, channel_indices
                    )
                    
                    if full_window.shape[0] < (end_pulse - start_pulse):
                        continue # Near end of recording, too unstable to fit
                        
                    x_full = np.arange(full_window.shape[0])
                    for channel_idx in range(num_channels):
                        y_full = full_window[:, channel_idx]
                        
                        if self._mode == "poly":
                            try:
                                with warnings.catch_warnings():
                                    warnings.simplefilter("ignore")
                                    fit_full = np.polyval(np.polyfit(x_full, y_full, deg=self._poly_order), x_full)
                            except Exception:
                                continue
                        elif self._mode == "spline":
                            from scipy.interpolate import LSQUnivariateSpline
                            splknots = np.arange(self._dspline / 2.0, len(x_full) - self._dspline / 2.0 + 2, self._dspline)
                            if len(splknots) == 0:
                                continue
                            spl = LSQUnivariateSpline(x=x_full, y=y_full, t=splknots, k=self._spline_order)
                            fit_full = spl(x_full)

                        # Extract ONLY the portion of the fit that lands in this chunk
                        rel_start = overlap_start - start_pulse
                        rel_end = overlap_end - start_pulse
                        fit_chunk = fit_full[rel_start:rel_end]
                        
                        # Indices within the current traces array to subtract from
                        t1 = overlap_start - start_frame
                        t2 = overlap_end - start_frame
                        
                        traces[t1:t2, channel_idx] -= fit_chunk.astype(traces.dtype)

        return traces


# API function
trend_subtract = define_function_from_class(
    source_class=TrendSubtractedRecording, name="trend_subtract"
)
