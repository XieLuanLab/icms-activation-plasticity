# responses/stim_condition_response.py
import numpy as np
import spikeinterface.full as si
import batch_process.postprocessing.stim_response_util as stim_response_util
import batch_process.postprocessing.responses_v2.pulse_locked_response_metrics as plrm
import batch_process.util.template_util as template_util


class StimConditionResponse:
    '''
    Stim condition response consists of two responses: train response and pulse response 
    Each is subclass of base response
    '''

    def __init__(self, stim_timestamps, stim_channel, stim_current, unit_response):
        self._unit_id = unit_response.unit_id
        self._stim_channel = stim_channel
        self._stim_current = stim_current
        # stim times for stim condition
        self._stim_timestamps = np.array(stim_timestamps)
        self._all_unit_spike_timestamps = unit_response.spike_timestamps
        self._unit_response = unit_response  # Reference to parent UnitResponse
        # Access to SessionResponses
        self._session_responses = unit_response.session_responses

        self._pulse_response = None
        self._train_response = None

        self._unit_response._add_stim_response(
            stim_channel, stim_current, self)

        self._stim_spikes = None
        self._stim_spike_indices = None
        self._stim_waveforms = None

    @property
    def stim_spikes(self):
        if self._stim_spikes is None:
            window_samples = self._session_responses.timing_params['pulse_post_samples']
            stim_spikes = []
            stim_spike_indices = []

            # Check which spikes occur within the stimulus window
            for idx, spike_ts in enumerate(self.all_unit_spike_timestamps):
                is_stim_related = np.any(
                    (spike_ts > self.stim_timestamps) & (spike_ts <= self.stim_timestamps + window_samples))

                if is_stim_related:
                    stim_spikes.append(spike_ts)
                    stim_spike_indices.append(idx)

            self._stim_spikes = np.array(stim_spikes)
            self._stim_spike_indices = np.array(stim_spike_indices)

        return self._stim_spikes

    @property
    def stim_spike_indices(self):
        if self._stim_spike_indices is None:
            self.stim_spikes  # Ensure stim spikes are calculated
        return self._stim_spike_indices

    @property
    def stim_waveforms(self):
        if self._stim_waveforms is None:
            primary_ch_index = template_util.get_dense_primary_ch_index(
                self._session_responses.sorting_analyzer, self.unit_response.unit_id)
            wvf_ext = self._session_responses.sorting_analyzer.get_extension(
                "waveforms")

            # Get waveforms for stim spikes
            if self.stim_spike_indices.size == 0:
                return []
            self._stim_waveforms = wvf_ext.get_waveforms_one_unit(
                self.unit_response.unit_id)[self.stim_spike_indices, :, primary_ch_index]

        return self._stim_waveforms

    def _add_pulse_response(self, pulse_response):
        self._pulse_response = pulse_response

    def _add_train_response(self, train_response):
        self._train_response = train_response

    def transform_tr_features(self):
        tr = self.train_response
        stim_spike_indices = tr._spike_indices_dict['stim_spike_indices']
        valid_arrays = [arr for arr in stim_spike_indices if len(arr) > 0]
        flattened = np.concatenate(
            valid_arrays).tolist() if valid_arrays else []

        raster = tr.raster_array
        stim_raster = [trial[(trial >= 0) & (trial < 710)] for trial in raster]

        return flattened, stim_raster

    def check_spike_counts_match(self):
        pr = self.pulse_response
        tr = self.train_response
        stim_spike_indices = tr._spike_indices_dict['stim_spike_indices']

        valid_arrays = [arr for arr in stim_spike_indices if len(arr) > 0]
        flattened = np.concatenate(
            valid_arrays).tolist() if valid_arrays else []

        pr_ts_count = len(pr.rel_spike_timestamps)
        tr_ts_count = sum(len(a)
                          for a in tr.spike_times_dict['train_spike_times'])

        if pr_ts_count != tr_ts_count:
            flattened, stim_raster = self.transform_tr_features()
            print(f"Mismatch in spike counts:")
            print(f"  Pulse response count: {pr_ts_count}")
            print(f"  Train response count: {tr_ts_count}")
            print(f"  Pulse response indices: {pr._original_indices}")
            print(f"  Flattened train indices: {flattened}")
            print(f"  Raster stim: {stim_raster}")

        assert pr_ts_count == tr_ts_count, "Mismatch between pulse and train spike counts"

    def __repr__(self):
        return f"{self.__class__.__name__}(stim_condition=D{self._stim_channel}, {self._stim_current} uA)"

    # Properties for accessing private attributes
    @property
    def unit_id(self):
        return self._unit_id

    @property
    def stim_channel(self):
        return self._stim_channel

    @property
    def stim_current(self):
        return self._stim_current

    @property
    def stim_timestamps(self):
        return self._stim_timestamps

    @property
    def all_unit_spike_timestamps(self):
        return self._all_unit_spike_timestamps

    @property
    def unit_response(self):
        return self._unit_response

    @property
    def session_responses(self):
        return self._session_responses

    @property
    def pulse_response(self):
        return self._pulse_response

    @property
    def train_response(self):
        return self._train_response


class BaseResponse:
    def __init__(self, timing_params, stim_condition_response):
        self._all_unit_spike_timestamps = stim_condition_response.all_unit_spike_timestamps
        self._stim_timestamps = stim_condition_response.stim_timestamps
        self._timing_params = timing_params
        self._stim_condition_response = stim_condition_response

        # Register this response with the parent stim_condition_response
        self._register_with_parent()

    def _register_with_parent(self):
        raise NotImplementedError(
            "This method should be implemented by subclasses")

    @property
    def stim_timestamps(self):
        return self._stim_timestamps

    @property
    def all_unit_spike_timestamps(self):
        return self._all_unit_spike_timestamps

    @property
    def stim_channel(self):
        return self._stim_condition_response.stim_channel

    @property
    def stim_current(self):
        return self._stim_condition_response.stim_current

    @property
    def timing_params(self):
        return self._timing_params

    @property
    def stim_condition_response(self):
        return self._stim_condition_response


class PulseResponse(BaseResponse):
    def __init__(self, stim_condition_response, timing_params):
        super().__init__(stim_condition_response, timing_params)
        self.stim_metrics = {}
        self.baseline_metrics = {}
        self._compute_pulse_response()
        self._compute_pulse_locked_response_metrics()

    @property
    def non_stim_waveforms(self):
        all_spike_ts = self._all_unit_spike_timestamps
        stim_waveform_indices = self._original_indices
        pass

    @property
    def stim_waveforms(self):
        pass

    def _adjust_firing_rate_for_blank_regions(self):
        # Extend the firing rate and fr_times arrays to include t=0 and t=pulse_win_ms with fr=0
        pulse_win_ms = self._timing_params["pulse_win_ms"]
        post_stim_blank_ms = self._timing_params["post_stim_blank_ms"]
        pre_stim_blank_ms = self._timing_params["pre_stim_blank_ms"]
        bin_size = self._fr_times[1] - \
            self._fr_times[0]  # Assuming uniform binning

        # Create a full time axis from 0 to pulse_win_ms with consistent binning
        fr_times_full = np.arange(0, pulse_win_ms, bin_size) + bin_size / 2
        # Avoid floating-point issues
        fr_times_full = np.round(fr_times_full, decimals=8)

        # Initialize the full firing rate array with zeros
        firing_rate_full = np.zeros_like(fr_times_full)
        binned_firing_rate_full = np.zeros_like(fr_times_full)

        # Create a mask for the non-blanked regions
        non_blank_mask = (fr_times_full >= post_stim_blank_ms) & (
            fr_times_full <= (pulse_win_ms - pre_stim_blank_ms)
        )

        # Verify that the non-blanked times align with the original firing rate times
        fr_times_non_blank = fr_times_full[non_blank_mask]

        # Due to potential floating-point precision issues, use np.isclose for comparison
        if len(fr_times_non_blank) != len(self._fr_times) or not np.allclose(fr_times_non_blank, self._fr_times):
            raise ValueError(
                "Time bins do not align between fr_times_non_blank and self._fr_times")

        # Assign the computed firing rates to the appropriate positions in the full array
        firing_rate_full[non_blank_mask] = self._firing_rate
        binned_firing_rate_full[non_blank_mask] = self._binned_firing_rate

        # increase size of firing rate array to include blanked portions of window
        self._fr_times = fr_times_full
        self._firing_rate = firing_rate_full
        self._binned_firing_rate = binned_firing_rate_full

    def _compute_pulse_response(self):

        # TODO: check if pulse number same as train
        self._rel_spike_timestamps, self._pulse_indices, self._original_indices = (
            stim_response_util.get_relative_spike_data(
                self._all_unit_spike_timestamps, self._stim_timestamps, self._timing_params
            )
        )
        self._raster_array = stim_response_util.get_pulse_interval_raster_array(
            self._rel_spike_timestamps, self._pulse_indices, self._stim_timestamps
        )

        bin_min = self._timing_params['post_stim_blank_ms']
        bin_max = self._timing_params['pulse_win_ms'] - \
            self._timing_params['pre_stim_blank_ms']
        self._firing_rate, self._binned_firing_rate, self._fr_times = stim_response_util.compute_smoothed_firing_rate(
            self._raster_array, bin_min, bin_max, bin_size=0.1, sigma=0.25
        )

        if len(self._fr_times) == 0:
            return

        self._adjust_firing_rate_for_blank_regions()

    def _compute_pulse_locked_response_metrics(self):
        # Compute stimulation metrics
        # rel_stim_spike_times, stim_trial_indices, _ = stim_response_util.get_relative_spike_data(
        #     self.all_unit_spike_timestamps, self.stim_timestamps, self.timing_params
        # )
        # stim_raster_array = stim_response_util.get_pulse_interval_raster_array(
        #     rel_stim_spike_times, stim_trial_indices)

        stim_raster_array = self._raster_array

        # Compute baseline metrics
        baseline_interval_start_times = plrm.get_baseline_interval_start_times(
            self._stim_timestamps, self._timing_params)
        rel_baseline_spike_times, baseline_pulse_indices, _ = stim_response_util.get_relative_spike_data(
            self._all_unit_spike_timestamps, baseline_interval_start_times, self._timing_params
        )
        baseline_raster_array = stim_response_util.get_pulse_interval_raster_array(
            rel_baseline_spike_times, baseline_pulse_indices, self._stim_timestamps
        )

        # Get bootstrapped metrics for baseline
        (
            self.baseline_metrics["bin_centers"],
            self.baseline_metrics["mean_prob_spike"],
            self.baseline_metrics["std_prob_spike"],
            self.baseline_metrics["mean_latency"],
            self.baseline_metrics["std_latency"],
        ) = plrm.get_bootstrapped_metrics(baseline_raster_array, self._timing_params, subsample_fraction=0.2, num_iterations=500)

        # Get bootstrapped metrics for stimulation
        (
            self.stim_metrics["bin_centers"],
            self.stim_metrics["mean_prob_spike"],
            self.stim_metrics["std_prob_spike"],
            self.stim_metrics["mean_latency"],
            self.stim_metrics["std_latency"],
        ) = plrm.get_bootstrapped_metrics(stim_raster_array, self._timing_params, subsample_fraction=0.2, num_iterations=500)

        # Calculate the Pulse-Locked Index (PLI)
        self._pli = plrm.calculate_pulse_locked_index(
            self.stim_metrics["bin_centers"], self.stim_metrics["mean_prob_spike"], self._timing_params
        )

        # Calculate the null distribution and determine if pulse-locked
        null_dist = plrm.get_null_distribution(
            stim_raster_array, self._timing_params, plot_flag=False)
        self._is_pulse_locked = plrm.is_pulse_locked(self._pli, null_dist)

    def _register_with_parent(self):
        self._stim_condition_response._add_pulse_response(self)

    @property
    def rel_spike_timestamps(self):
        return self._rel_spike_timestamps

    @property
    def trial_indices(self):
        return self._trial_indices

    @property
    def raster_array(self):
        return self._raster_array

    @property
    def firing_rate(self):
        return self._firing_rate

    @property
    def binned_firing_rate(self):
        return self._binned_firing_rate

    @property
    def fr_times(self):
        return self._fr_times

    @property
    def pli(self):
        return self._pli

    @property
    def is_pulse_locked(self):
        return self._is_pulse_locked


# %%


class TrainResponse(BaseResponse):
    def __init__(self, stim_condition_response, timing_params):
        super().__init__(stim_condition_response, timing_params)
        self._compute_train_response()
        self._register_with_parent()

    def _compute_train_response(self):
        # Group stimulus pulses into trains
        self._stim_trains = stim_response_util.group_stim_pulses_into_trains(
            self._stim_timestamps, self._timing_params)

        # Analyze the response to the stimulus train
        # (
        #     self._spike_times_dict,
        #     self._spike_indices_dict,
        #     self._t_test_dict,
        #     self._z_score_dict,
        # ) = stim_response_util.analyze_stimulus_train_response(
        #     self._stim_trains, self._all_unit_spike_timestamps, self._timing_params
        # )

        # Todo: added 6/25
        # (
        #     self._spike_times_dict,
        #     _,
        #     _,
        #     _,
        # ) = stim_response_util.analyze_stimulus_train_response(
        #     self._stim_trains, self._all_unit_spike_timestamps, self._timing_params, larger_window_flag=True
        # )
        # Todo: Added 7/23
        (
            self._spike_times_dict,
            self._spike_indices_dict,
            self._t_test_dict,
            self._z_score_dict,
        ) = stim_response_util.analyze_stimulus_train_response(
            self._stim_trains, self._all_unit_spike_timestamps, self._timing_params, different_window_flag=True
        )

        # Generate a raster array for visualization or further analysis
        self._raster_array = stim_response_util.get_train_win_raster_arr(
            self._spike_times_dict,
            self._stim_trains,
            self._timing_params,
        )

        # Compute the smoothed firing rate across the train window
        bin_min = -self._timing_params['stim_duration_ms']  # -700ms
        bin_max = self._timing_params['window_size_ms']  # 1400ms
        pulse_win_ms = self._timing_params['pulse_win_ms']

        # bin and don't smooth, use 2 ms bin size
        bin_size = 2  # ms
        bin_centers, hist, time_to_max_ms = stim_response_util.get_time_to_max_spikes(
            self._raster_array, bin_min=0, bin_max=700 + pulse_win_ms, bin_size=bin_size)

        # time_to_threshold_crossing = stim_response_util.get_time_to_fr_threshold_crossing(
        #     bin_centers, firing_rate, self.pre_stim_mean_fr, self.pre_stim_sigma_fr)

        # Normalize histogram to get average spike count per trial per bin
        num_trials = len(self._raster_array)
        hist = hist / num_trials
        # Convert bin size to seconds
        binned_firing_rate = hist / (bin_size / 1000)
        time_to_threshold_crossing = stim_response_util.get_time_to_fr_threshold_crossing(
            bin_centers, binned_firing_rate, self.pre_stim_mean_fr, self.pre_stim_sigma_fr)

        self._two_ms_binned_data = {
            'bin_centers': bin_centers,
            'binned_firing_rate': binned_firing_rate,
            'hist': hist,
            'time_to_max_ms': time_to_max_ms,
            'time_to_threshold_crossing': time_to_threshold_crossing
        }

        # bin and smooth, use 2 ms bin size with 5 ms std kernel
        firing_rate, hist, bin_centers = stim_response_util.compute_smoothed_firing_rate(
            self._raster_array, bin_min=0, bin_max=700 + pulse_win_ms, bin_size=2, sigma=5
        )
        time_to_max_fr = stim_response_util.get_time_max_fr(
            bin_centers, firing_rate)
        time_to_threshold_crossing = stim_response_util.get_time_to_fr_threshold_crossing(
            bin_centers, firing_rate, self.pre_stim_mean_fr, self.pre_stim_sigma_fr)

        self._two_ms_smoothed_data = {
            'bin_centers': bin_centers,
            'firing_rate': firing_rate,
            'hist': hist,
            'time_to_max_fr': time_to_max_fr,
            'time_to_threshold_crossing': time_to_threshold_crossing
        }

        # bin and smooth, use 10 ms bin size with 20 ms std kernel
        firing_rate, hist, bin_centers = stim_response_util.compute_smoothed_firing_rate(
            self._raster_array, bin_min=bin_min, bin_max=bin_max, bin_size=10, sigma=20
        )
        time_to_max_fr = stim_response_util.get_time_max_fr(
            bin_centers, firing_rate, bin_min=0, bin_max=700)
        time_to_threshold_crossing_10ms = stim_response_util.get_time_to_fr_threshold_crossing(
            bin_centers, firing_rate, self.pre_stim_mean_fr, self.pre_stim_sigma_fr, bin_min=0, bin_max=700)

        self._ten_ms_smoothed_data = {
            'bin_centers': bin_centers,
            'firing_rate': firing_rate,
            'hist': hist,
            'time_to_max_fr': time_to_max_fr,
            'time_to_threshold_crossing': time_to_threshold_crossing_10ms
        }

    def _register_with_parent(self):
        self._stim_condition_response._add_train_response(self)

    # Properties to access the computed data
    @property
    def stim_trains(self):
        return self._stim_trains

    @property
    def spike_times_dict(self):
        return self._spike_times_dict

    @property
    def raster_array(self):
        return self._raster_array

    @property
    def two_ms_binned_data(self):
        return self._two_ms_binned_data

    @property
    def two_ms_smoothed_data(self):
        return self._two_ms_smoothed_data

    @property
    def ten_ms_smoothed_data(self):
        return self._ten_ms_smoothed_data

    @property
    def t_test_dict(self):
        return self._t_test_dict

    @property
    def paired_t_val(self):
        return self._t_test_dict['paired_t_val']

    @property
    def paired_p_val(self):
        return self._t_test_dict['paired_p_val']

    @property
    def ones_samp_t_val(self):
        return self._t_test_dict['ones_samp_t_val']

    @property
    def one_samp_p_val(self):
        return self._t_test_dict['one_samp_p_val']

    @property
    def z_score_dict(self):
        return self._z_score_dict

    @property
    def z_score(self):
        return self._z_score_dict['z_score']

    @property
    def stim_mean_fr(self):
        return self._z_score_dict['stim_mean_fr']

    @property
    def pre_stim_mean_fr(self):
        return self._z_score_dict['pre_stim_mean_fr']

    @property
    def pre_stim_sigma_fr(self):
        return self._z_score_dict['pre_stim_sigma_fr']
