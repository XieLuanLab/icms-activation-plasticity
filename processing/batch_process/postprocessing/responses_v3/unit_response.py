"""UnitResponse — stores per-unit spike data and stim condition responses."""

from .stim_condition_response import StimConditionResponse
import numpy as np
import batch_process.util.template_util as template_util


class UnitResponse:
    def __init__(self, unit_id, spike_timestamps, session_responses):
        self._session_responses = session_responses
        self._unit_id = unit_id
        self._cell_type = None
        self._unit_location = None
        self._stim_conditions = []
        self._stim_responses = {}
        self._timing_params = self._session_responses.timing_params
        self._spike_timestamps = spike_timestamps

        self._primary_channel_template = None
        self._non_stim_spikes = None
        self._non_stim_spike_indices = None
        self._non_stim_waveforms = None

        self._stim_spikes = None
        self._stim_spike_indices = None
        self._stim_waveforms = None

        self._session_responses._add_unit_response(unit_id, self)

    def is_sorting_analyzer_loaded(self):
        return self._session_responses.is_sorting_analyzer_loaded()

    @property
    def cell_type(self):
        if self._session_responses.cell_type_df is None:
            return
        self._cell_type = self._session_responses.cell_type_df.loc[self._unit_id, "cell_type"]
        return self._cell_type

    @property
    def unit_location(self):
        if not self.is_sorting_analyzer_loaded():
            return
        if self._unit_location is None:
            ul_ext = self._session_responses.sorting_analyzer.get_extension(
                "unit_locations")
            if ul_ext is None:
                self._session_responses.sorting_analyzer.compute(
                    'unit_locations')
                ul_ext = self._session_responses.sorting_analyzer.get_extension(
                    "unit_locations")
            unit_index = np.where(
                self._session_responses.sorting_analyzer.unit_ids == self._unit_id)
            self._unit_location = ul_ext.get_data()[unit_index][0]
        return self._unit_location

    @property
    def primary_channel_template(self):
        if not self.is_sorting_analyzer_loaded():
            return
        if self._primary_channel_template is None:
            self._primary_channel_template = template_util.get_unit_primary_ch_template(
                self._session_responses.sorting_analyzer, self._unit_id)
        return self._primary_channel_template

    def _get_nonstim_spike_data(self):
        window_samples = self._timing_params['pulse_post_samples']
        all_stim_timestamps = self._session_responses.all_stim_timestamps

        if len(all_stim_timestamps) == 0:
            self._non_stim_spikes = np.array(self._spike_timestamps)
            self._non_stim_indices = np.arange(len(self._spike_timestamps))
            return

        concat_stim_timestamps = np.concatenate(all_stim_timestamps)
        non_stim_spikes = []
        non_stim_spike_indices = []
        for idx, spike_ts in enumerate(self._spike_timestamps):
            is_stim_related = np.any((spike_ts > concat_stim_timestamps) & (
                spike_ts <= concat_stim_timestamps + window_samples))
            if not is_stim_related:
                non_stim_spikes.append(spike_ts)
                non_stim_spike_indices.append(idx)

        self._non_stim_spikes = np.array(non_stim_spikes)
        self._non_stim_spike_indices = np.array(non_stim_spike_indices)

    @property
    def non_stim_spikes(self):
        if self._non_stim_spikes is None:
            self._get_nonstim_spike_data()
        return self._non_stim_spikes

    @property
    def non_stim_spike_indices(self):
        if self._non_stim_spike_indices is None:
            self._get_nonstim_spike_data()
        return self._non_stim_spike_indices

    @property
    def non_stim_waveforms(self):
        if not self.is_sorting_analyzer_loaded():
            return
        if self._non_stim_waveforms is None:
            wvf_ext = self._session_responses._sorting_analyzer.get_extension(
                "waveforms")
            if self._session_responses._sorting_analyzer.get_extension('random_spikes').params['max_spikes_per_unit'] == 500:
                non_stim_spike_indices = self.non_stim_spike_indices[
                    self.non_stim_spike_indices < 500]
            else:
                non_stim_spike_indices = self.non_stim_spike_indices
            primary_ch_index = template_util.get_dense_primary_ch_index(
                self._session_responses._sorting_analyzer, self._unit_id)
            self._non_stim_waveforms = wvf_ext.get_waveforms_one_unit(
                self._unit_id)[non_stim_spike_indices, :, primary_ch_index]
        return self._non_stim_waveforms

    @property
    def stim_spikes(self):
        if self._non_stim_spikes is None:
            self._get_nonstim_spike_data()
        stim_spike_indices = np.setdiff1d(
            np.arange(len(self._spike_timestamps)), self._non_stim_spike_indices)
        self._stim_spikes = self._spike_timestamps[stim_spike_indices]
        return self._stim_spikes

    @property
    def stim_spike_indices(self):
        if self._non_stim_spike_indices is None:
            self._get_nonstim_spike_data()
        self._stim_spike_indices = np.setdiff1d(
            np.arange(len(self._spike_timestamps)), self._non_stim_spike_indices)
        return self._stim_spike_indices

    @property
    def stim_waveforms(self):
        if not self.is_sorting_analyzer_loaded():
            return
        if self._stim_waveforms is None:
            primary_ch_index = template_util.get_dense_primary_ch_index(
                self._session_responses._sorting_analyzer, self._unit_id)
            wvf_ext = self._session_responses.sorting_analyzer.get_extension(
                "waveforms")
            self._stim_waveforms = wvf_ext.get_waveforms_one_unit(
                self._unit_id)[self.stim_spike_indices, :, primary_ch_index]
        return self._stim_waveforms

    def _add_stim_response(self, stim_channel, stim_current, stim_response):
        stim_condition = (stim_channel, stim_current)
        self._stim_conditions.append(stim_condition)
        self._stim_responses[stim_condition] = stim_response

    def get_stim_response(self, stim_channel, stim_current):
        stim_condition = (stim_channel, stim_current)
        return self._stim_responses.get(stim_condition, None)

    def show_stim_conditions(self):
        for stim_channel, stim_current in self._stim_conditions:
            print(f"Stim channel depth: {stim_channel}, current: {stim_current} uA")
        return self._stim_conditions

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(unit_id={self.unit_id}, "
            f"stim_responses={len(self.stim_responses)} conditions)")

    @property
    def session_responses(self):
        return self._session_responses

    @property
    def unit_id(self):
        return self._unit_id

    @property
    def stim_responses(self):
        return self._stim_responses

    @property
    def timing_params(self):
        return self._timing_params

    @property
    def spike_timestamps(self):
        return self._spike_timestamps
