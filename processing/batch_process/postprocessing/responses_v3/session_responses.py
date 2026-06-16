"""SessionResponses — top-level container for a session's response data."""

from .unit_response import UnitResponse
from .window_config import WindowConfig
from spikeinterface import full as si
from pathlib import Path


class SessionResponses:
    def __init__(self, session_path, sorting_analyzer_path, trial_df,
                 all_stim_timestamps, timing_params, window_config=None):
        self._session_path = session_path
        self._sorting_analyzer_path = sorting_analyzer_path
        self._sorting_analyzer = None
        self._trial_df = trial_df
        self._unit_responses = {}
        self._all_stim_timestamps = all_stim_timestamps
        self._timing_params = timing_params
        self._cell_type_df = None
        self._window_config = window_config or WindowConfig.standard_700ms()

    def load_sorting_analyzer(self):
        if self._sorting_analyzer is None:
            print("Loading sorting analyzer from path...")
            self._sorting_analyzer = si.load_sorting_analyzer(
                self._sorting_analyzer_path)
        return self._sorting_analyzer

    def is_sorting_analyzer_loaded(self):
        if self._sorting_analyzer is None:
            print("Sorting analyzer not loaded! Call `load_sorting_analyzer()` first.")
            return False
        return True

    def _add_unit_response(self, unit_id, unit_response):
        self._unit_responses[unit_id] = unit_response

    def get_unit_response(self, unit_id):
        if unit_id not in self._unit_responses:
            raise KeyError(f"UnitResponse for unit_id {unit_id} not found!")
        return self._unit_responses[unit_id]

    def __repr__(self):
        return (f"{self.__class__.__name__}(session_path={self._session_path!r}, "
                f"units={len(self._unit_responses)}, window={self._window_config.label})")

    # --- Backward compat for unpickling old files without window_config ---
    def __getattr__(self, name):
        if name == '_window_config':
            self._window_config = WindowConfig.standard_700ms()
            return self._window_config
        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")

    @property
    def cell_type_df(self):
        if not self.is_sorting_analyzer_loaded():
            return
        if self._cell_type_df is None:
            import batch_process.util.classify_cell_type as cell_classifier
            print("Computing cell type classifications...")
            self._cell_type_df = cell_classifier.classify_units_into_cell_types(
                self.sorting_analyzer)["df"]
        return self._cell_type_df

    @property
    def session_path(self):
        return self._session_path

    @property
    def trial_df(self):
        return self._trial_df

    @property
    def unit_responses(self):
        return self._unit_responses

    @property
    def unit_ids(self):
        return list(self._unit_responses.keys())

    @property
    def sorting_analyzer(self):
        if self._sorting_analyzer is None:
            print("Sorting analyzer not loaded.")
            return
        return self._sorting_analyzer

    @property
    def sorting_analyzer_path(self):
        return self._sorting_analyzer_path

    @property
    def all_stim_timestamps(self):
        return self._all_stim_timestamps

    @property
    def timing_params(self):
        return self._timing_params

    @property
    def window_config(self):
        return self._window_config
