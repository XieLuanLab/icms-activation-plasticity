"""
Run v3 pipeline on all experimental animals and ICMS83 with specified window configs.

Usage:
    # 700ms on all animals
    python -m batch_process.postprocessing.responses_v3.run_all_animals --window 700ms

    # 150ms on all animals
    python -m batch_process.postprocessing.responses_v3.run_all_animals --window 150ms

    # Both windows on all animals
    python -m batch_process.postprocessing.responses_v3.run_all_animals --window 700ms 150ms

    # Force re-run
    python -m batch_process.postprocessing.responses_v3.run_all_animals --window 700ms --force
"""
import sys
import pickle
from pathlib import Path

from batch_process.postprocessing.responses_v3.window_config import WindowConfig
from batch_process.postprocessing.responses_v3.run_pipeline import run_pipeline
from batch_process.postprocessing.responses_v3.run_pipeline_mouse6 import run_all_mouse6
import batch_process.util.file_util as file_util

PROJECT_ROOT = Path(__file__).resolve().parents[3]

WINDOW_MAP = {
    '700ms': WindowConfig.standard_700ms(),
    'sym120ms': WindowConfig.symmetric_120ms(),
    'sym100ms': WindowConfig.symmetric_100ms(),
    'sym150ms': WindowConfig.symmetric_150ms(),
    'sym700ms': WindowConfig.symmetric_700ms(),
    '2700ms': WindowConfig.extended_2700ms(),
}


def get_experimental_folders():
    """Load all experimental animal data folders from pickle."""
    pkl_path = PROJECT_ROOT / "all_animal_folder_list.pkl"
    with open(pkl_path, "rb") as f:
        all_animal_folders_list = pickle.load(f)
    all_folders = [item for sublist in all_animal_folders_list for item in sublist]
    return file_util.sort_data_folders(all_folders)


def run_experimental(window_configs, stim_trial_types=None, force=False):
    """Run v3 on all experimental animals (92, 93, 98, 100, 101)."""
    if stim_trial_types is None:
        stim_trial_types = ['all']

    folders = get_experimental_folders()
    total = len(folders)

    for i, data_folder in enumerate(folders):
        for wc in window_configs:
            # Check if already done
            save_folder = Path(data_folder) / "batch_sort"
            all_done = all(
                (save_folder / f"stim_condition_results_{tt}_{wc.label}.csv").exists()
                for tt in stim_trial_types)

            if all_done and not force:
                print(f"[{i+1}/{total}] {data_folder} [{wc.label}] — SKIP (exists)")
                continue

            print(f"\n[{i+1}/{total}] {data_folder} [{wc.label}]")
            try:
                run_pipeline(
                    data_folder,
                    window_configs=[wc],
                    stim_trial_types=stim_trial_types)
            except Exception as e:
                print(f"  [ERROR] {e}")
                continue


def run_all(window_configs, stim_trial_types=None, force=False):
    """Run v3 on all animals: ICMS83 + experimental."""
    if stim_trial_types is None:
        stim_trial_types = ['all']

    print("=" * 60)
    print("ICMS83")
    print("=" * 60)
    run_all_mouse6(
        window_configs=window_configs,
        stim_trial_types=stim_trial_types,
        force=force)

    print("\n" + "=" * 60)
    print("EXPERIMENTAL ANIMALS (92, 93, 98, 100, 101)")
    print("=" * 60)
    run_experimental(
        window_configs=window_configs,
        stim_trial_types=stim_trial_types,
        force=force)

    print("\nAll animals complete.")


if __name__ == "__main__":
    force = "--force" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]

    # Parse --window flag
    window_labels = []
    if "--window" in sys.argv:
        idx = sys.argv.index("--window")
        for j in range(idx + 1, len(sys.argv)):
            if sys.argv[j].startswith("--"):
                break
            window_labels.append(sys.argv[j])

    if not window_labels:
        window_labels = ['700ms']

    window_configs = []
    for label in window_labels:
        if label in WINDOW_MAP:
            window_configs.append(WINDOW_MAP[label])
        else:
            print(f"Unknown window: {label}. Options: {list(WINDOW_MAP.keys())}")
            sys.exit(1)

    print(f"Windows: {[wc.label for wc in window_configs]}")
    print(f"Force: {force}")

    run_all(window_configs=window_configs, stim_trial_types=['all'], force=force)
