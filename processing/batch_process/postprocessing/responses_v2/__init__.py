"""Legacy (v2) response model.

Superseded by responses_v3 for the main feature table (raw_df). Retained because
pop_coupling and the per-unit QC figures in stage4_postprocess still build on it.
"""
from .session_responses_v2 import SessionResponses
from .unit_response_v2 import UnitResponse
from .stim_condition_response_v2 import StimConditionResponse, PulseResponse, TrainResponse

__all__ = [
    "SessionResponses",
    "UnitResponse",
    "StimConditionResponse",
    "PulseResponse",
    "TrainResponse",
]
