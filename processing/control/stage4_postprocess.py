"""
Control Spike Sorting Pipeline - Stage 4 (Postprocessing)
=========================================================
Run AFTER stage3_merge.py has produced analyzer_final.zarr for all sessions.

This stage:
1. Creates SessionResponses for each session (responses to stim conditions)
2. Extracts modulation metrics → CSV
3. Generates per-unit stim response figures → PNG + PPT
4. (Optional) Runs longitudinal analysis across sessions

Outputs per session (in spike_sorting/):
  - session_responses.pkl
  - stim_condition_results_control.csv
  - figures/stim_condition_figures/*.png
  - figures/stim_responses_{session}.pptx

Usage (from terminal):
    cd C:\\Users\\robin\\Desktop\\Python\\research\\icms-plasticity
    C:\\Users\\robin\\Desktop\\Python\\research\\si_v103_env\\Scripts\\python.exe pipeline\\stage4_postprocess.py
"""

import sys
import os
import io
import traceback
import time
from pathlib import Path

# Project root for imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

# Fix Windows encoding for unicode characters
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use("Agg")

from util.load_control import discover_sessions  # noqa: E402
from control.postprocessing.process_control_responses import run as run_postprocessing  # noqa: E402


def process_session_stage4(exp_dir, force=False):
    """
    Run postprocessing for a single session.

    Parameters
    ----------
    exp_dir : str
        Path to the session's ExperimentData folder.
    force : bool
        Re-run even if outputs already exist.

    Returns
    -------
    bool
        True if successful.
    """
    spike_dir = Path(exp_dir) / "spike_sorting"
    analyzer_path = spike_dir / "analyzer_final.zarr"
    csv_path = spike_dir / "condensed_trials.csv"

    # --- Checks ---
    if not analyzer_path.exists():
        print(f"  [SKIP] No analyzer_final.zarr -- run stage3 first")
        return False

    if not csv_path.exists():
        print(f"  [SKIP] No condensed_trials.csv")
        return False

    # Skip if already processed (check for CSV output)
    results_csv = spike_dir / "stim_condition_results_control.csv"
    if not force and results_csv.exists():
        print(f"  [SKIP] Already processed: {results_csv}")
        return True

    # Run the full postprocessing pipeline
    run_postprocessing(exp_dir, load_pickle=False)
    return True


def run_stage4(base_folder, force=False, allow_missing_mat=True):
    """
    Batch postprocessing for all control sessions in a base folder.

    Parameters
    ----------
    base_folder : str
        Root data directory for the animal.
    force : bool
        Re-run even if already processed.
    """
    sessions = discover_sessions(base_folder, allow_missing_mat=allow_missing_mat)
    print(f"Found {len(sessions)} session(s)\n")

    results = {}
    t_start = time.time()

    for i, exp_dir in enumerate(sessions, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(sessions)}] {exp_dir}")
        print(f"{'='*60}")

        try:
            ok = process_session_stage4(exp_dir, force=force)
            results[exp_dir] = "OK" if ok else "SKIPPED"
        except Exception as e:
            print(f"\n  [ERROR] {e}")
            traceback.print_exc()
            results[exp_dir] = f"ERROR: {e}"

    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"STAGE 4 COMPLETE -- {len(sessions)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for exp_dir, status in results.items():
        icon = "[OK]" if status == "OK" else "[-]" if "SKIP" in status else "[X]"
        print(f"  {icon} {Path(exp_dir).parent.name}/{Path(exp_dir).name}: {status}")


if __name__ == "__main__":
    for animal in ["ICMS43", "ICMS45", "ICMS48", "ICMS54", "ICMS56"]:
        print(f"\n\n{'#'*60}")
        print(f"  STAGE 4: {animal}")
        print(f"{'#'*60}")
        run_stage4(base_folder=f"E:/ICMS plasticity control/{animal}", force=True)
