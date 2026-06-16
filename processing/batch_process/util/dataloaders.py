"""Session loader dispatch for the stim cohort.

Mouse 1-5 (ICMS92/93/98/100/101) use the standard DataLoader. Mouse 6 (ICMS83)
was acquired on a different rig with a different on-disk layout (separate
key-file format, E:/ drive) so it needs its own loader; both return the same
DataLoader interface (recording, trial_df, all_stim_timestamps) and everything
downstream is identical.
"""
import os
from pathlib import Path

from util.load_data import load_data
from batch_process.util.dataloader_mouse6 import load_mouse6


def load_session(data_folder, get_recording=True):
    if "ICMS83" in str(data_folder):
        dl = load_mouse6(str(data_folder), get_recording=get_recording)
        dl.save_folder = Path(data_folder) / "batch_sort"
        os.makedirs(dl.save_folder, exist_ok=True)
        return dl
    return load_data(
        str(data_folder), make_folder=True, save_folder_name="batch_sort",
        first_N_files=4, server_mount_drive="S:", get_recording=get_recording,
    )
