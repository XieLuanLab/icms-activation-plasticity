"""
Single entry point to run the full response analysis pipeline. This is the
current pipeline; its output (stim_condition_results) builds the figure raw_df.

Usage:
    from batch_process.postprocessing.responses_v3.run_pipeline import run_pipeline
    from batch_process.postprocessing.responses_v3.window_config import WindowConfig

    # Run with default 700ms window
    run_pipeline("C:/data/ICMS92/Behavior/30-Aug-2023")

    # Run with multiple windows
    run_pipeline(
        "C:/data/ICMS92/Behavior/30-Aug-2023",
        window_configs=[
            WindowConfig.standard_700ms(),
            WindowConfig.short_150ms(),
            WindowConfig.extended_2700ms(),
        ]
    )
"""
import dill as pickle
from pathlib import Path

from batch_process.postprocessing.responses_v3.window_config import WindowConfig
from batch_process.postprocessing.responses_v3.get_session_responses import get_session_responses
from batch_process.postprocessing.responses_v3.process_stim_responses import process_stim_responses


def run_pipeline(data_folder, window_configs=None, stim_trial_types=None,
                 server_mount_drive="S:"):
    """
    Run the full response analysis pipeline for a session.

    Args:
        data_folder: Path to the session data folder.
        window_configs: List of WindowConfig objects. Defaults to [standard_700ms].
        stim_trial_types: List of trial types (e.g., ['all', 'hit', 'miss']).
        server_mount_drive: Drive letter for server mount.
    """
    if window_configs is None:
        window_configs = [WindowConfig.standard_700ms()]
    if stim_trial_types is None:
        stim_trial_types = ['all']

    for window_config in window_configs:
        print(f"\n{'='*60}")
        print(f"Processing with window: {window_config.label}")
        print(f"  pre_stim_ms={window_config.pre_stim_ms}, "
              f"post_stim_ms={window_config.post_stim_ms}")
        print(f"{'='*60}\n")

        # Step 1: Build session responses (heavy — spike extraction)
        session_map = get_session_responses(
            data_folder,
            window_config=window_config,
            stim_trial_types=stim_trial_types,
            server_mount_drive=server_mount_drive)

        # Step 2: Extract features to CSV/pkl (lightweight)
        for stim_key, session_responses in session_map.items():
            session_responses.load_sorting_analyzer()
            process_stim_responses(
                session_responses, stim_key,
                window_label=window_config.label)

    print("\nPipeline complete.")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        data_folder = sys.argv[1]
    else:
        data_folder = None

    if data_folder:
        run_pipeline(data_folder)
    else:
        print("Usage: python -m batch_process.postprocessing.responses_v3.run_pipeline <data_folder>")
