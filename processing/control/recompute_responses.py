"""
Batch recompute session_responses.pkl for all control sessions
that have analyzer_final.zarr after stage 3 merge.
"""
import sys
import os
import time
import traceback
from pathlib import Path

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import matplotlib
matplotlib.use("Agg")

from control.postprocessing.process_control_responses import run

CONTROL_ROOT = Path("E:/ICMS plasticity control")
ANIMAL_IDS = ["ICMS43", "ICMS45", "ICMS48", "ICMS54", "ICMS56"]


def discover_experiment_dirs(animal_dir):
    """Find ExperimentData dirs that have analyzer_final.zarr."""
    sessions = []
    for spike_dir in sorted(animal_dir.rglob("spike_sorting")):
        if (spike_dir / "analyzer_final.zarr").exists():
            # session_dir = ExperimentData (parent of spike_sorting)
            sessions.append(spike_dir.parent)
    return sessions


if __name__ == "__main__":
    t0 = time.time()
    results = {}

    for animal_id in ANIMAL_IDS:
        animal_dir = CONTROL_ROOT / animal_id
        sessions = discover_experiment_dirs(animal_dir)
        print(f"\n{'='*60}")
        print(f"  {animal_id}: {len(sessions)} session(s)")
        print(f"{'='*60}")

        for i, session_dir in enumerate(sessions, 1):
            label = f"{animal_id}/{session_dir.parent.name}"
            print(f"\n  [{i}/{len(sessions)}] {label}")
            try:
                run(str(session_dir), load_pickle=False)
                results[label] = "OK"
            except Exception as e:
                print(f"  [ERROR] {e}")
                traceback.print_exc()
                results[label] = f"ERROR: {e}"

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"POSTPROCESSING COMPLETE — {len(results)} sessions in {elapsed/60:.1f} min")
    print(f"{'='*60}")
    for path, status in results.items():
        icon = "[OK]" if status == "OK" else "[X]"
        print(f"  {icon} {path}: {status}")
