"""
ICMS83 DataLoader adapter.

Produces the same trial_df and recording as the standard DataLoader
but handles ICMS83's unique data layout:

  - Key files at S:/ICMS83/Keys/{date}/Ch{intan}_D{depth}_H{H}M{M}.mat
    with dataCellArr format (N x 1 cell, each entry 8 x 1):
      [0] stimLevel (µA)     [1] response (0/1)
      [2] response_time (s)  [3] unknown flag
      [4] ImageStartStop     [5] FrameTimestamps
      [6] StimStartStop      [7] StimTimestamps (70 pulses @ 100 Hz)

  - NS5/NEV files in E:/ICMS83/{date}/experiment/
    Named D{depth}_H{H}M{M}_ephys{XXXX}.ns5 (one stim channel per file)

  - behavior.mat named D{depth}_H{H}M{M}_behavior.mat (top-level key "s")

  - Depth-to-channel mapping (consistent across all sessions):
      D3: depthIndex=3, ripple_ch=28, intan_ch=29
      D4: depthIndex=4, ripple_ch=11, intan_ch=10
      D5: depthIndex=5, ripple_ch=24, intan_ch=27

Usage:
    from batch_process.util.dataloader_mouse6 import load_mouse6

    dl = load_mouse6("E:/ICMS83/20-Jul-2023")
    dl.trial_df   # same schema as DataLoader.trial_df
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import scipy.io as sio
from pathlib import Path
from spikeinterface import full as si, extractors as se
from probeinterface.io import read_probeinterface
from batch_process.util.impedance_analyzer import ImpedanceAnalyzer


# ─── Constants ────────────────────────────────────────────────
FS = 30000.0
SERVER_KEY_ROOT = Path("S:/ICMS83/Keys")

# Filename → depth index mapping
DEPTH_MAP = {3: 3, 4: 4, 5: 5}  # D{X} in filename → depthIndex

# depth → (ripple_ch, intan_ch)  — consistent across all ICMS83 sessions
CHANNEL_MAP = {
    3: {"ripple": 28, "intan": 29},
    4: {"ripple": 11, "intan": 10},
    5: {"ripple": 24, "intan": 27},
}


class ICMS83DataLoader:
    """DataLoader adapter for ICMS83 data."""

    def __init__(self, data_folder, server_key_root=None):
        self.data_folder = Path(data_folder)
        self.experiment_folder = self.data_folder / "experiment"
        self.animalID = "ICMS83"
        self.fs = FS

        # Determine date from folder name
        self.date = self.data_folder.name

        # Server key folder
        key_root = Path(server_key_root) if server_key_root else SERVER_KEY_ROOT
        self.key_folder = key_root / self.date

        if not self.experiment_folder.exists():
            raise FileNotFoundError(
                f"No experiment/ folder in {self.data_folder}")

    # ─── File discovery ───────────────────────────────────────

    def get_ns5_and_mat_files(self):
        """Find behavior.mat files and their matching NS5 files."""
        all_files = os.listdir(self.experiment_folder)

        # Find behavior files: D{X}_H{H}M{M}_behavior.mat
        self.mat_files = sorted(
            f for f in all_files
            if re.match(r"D\d+_H\d+M\d+_behavior\.mat", f)
        )

        # Extract base names (D{X}_H{H}M{M})
        self.base_names = []
        for mf in self.mat_files:
            m = re.match(r"(D\d+_H\d+M\d+)_behavior\.mat", mf)
            self.base_names.append(m.group(1))

        # Find matching NS5 files
        self.ns5_files = []
        for base in self.base_names:
            matches = [f for f in all_files
                       if f.startswith(base + "_") and f.endswith(".ns5")]
            if matches:
                self.ns5_files.append(matches[0])
            else:
                self.ns5_files.append(None)

    def _find_key_file(self, base_name):
        """Find the server key file matching a behavior base name.

        base_name: e.g. "D3_H11M7"
        Returns: full path to Ch{X}_D{Y}_H{H}M{M}.mat, or None
        """
        if not self.key_folder.exists():
            return None

        key_files = list(self.key_folder.glob("*.mat"))
        for kf in key_files:
            m = re.match(r"Ch\d+_(D\d+_H\d+M\d+)\.mat", kf.name)
            if m and m.group(1) == base_name:
                return kf
        return None

    # ─── Key file parsing ─────────────────────────────────────

    def get_key_files(self):
        """Parse dataCellArr from server key files.

        Populates the same attributes as DataLoader.get_key_files():
            all_currents, all_depths, all_stim_timestamps,
            all_stim_start_stop, all_image_start_stop,
            all_frame_timestamps, all_behave_responses,
            all_response_times
        """
        (
            self.all_currents,
            self.all_depths,
            self.all_stim_timestamps,
            self.all_stim_start_stop,
            self.all_image_start_stop,
            self.all_behave_responses,
            self.all_response_times,
            self.all_frame_timestamps,
        ) = ([], [], [], [], [], [], [], [])

        self.blocks_used = []  # indices of blocks that have key files

        for i, base_name in enumerate(self.base_names):
            key_path = self._find_key_file(base_name)
            if key_path is None:
                print(f"  [SKIP] No key file for {base_name}")
                continue

            self.blocks_used.append(i)
            self._extract_dataCellArr(key_path, i)

    def _extract_dataCellArr(self, key_path, block_index):
        """Parse one dataCellArr key file and append to all_* lists."""
        mat = sio.loadmat(str(key_path))
        dca = mat["dataCellArr"]
        n_trials = dca.shape[0]

        # Get depth from filename
        m = re.match(r"D(\d+)_", self.base_names[block_index])
        depth_index = int(m.group(1))

        # Time offset for this block (seconds → samples)
        offset_samples = self.start_times[block_index] * self.fs

        for t in range(n_trials):
            trial = dca[t, 0]

            current = float(trial[0, 0].flat[0])
            response = int(trial[1, 0].flat[0])
            response_time = float(trial[2, 0].flat[0])

            # ImageStartStop (field 4) — seconds → samples
            img_ss = trial[4, 0].flatten()
            img_ss_samples = (img_ss * self.fs + offset_samples).tolist()

            # FrameTimestamps (field 5) — seconds → samples
            frame_ts = trial[5, 0].flatten()
            frame_ts_samples = (frame_ts * self.fs + offset_samples).tolist()

            # StimStartStop (field 6) — seconds → samples
            stim_ss = trial[6, 0].flatten()
            stim_ss_samples = (stim_ss * self.fs + offset_samples).tolist()

            # StimTimestamps (field 7) — seconds → samples
            stim_ts = trial[7, 0].flatten()
            stim_ts_samples = (stim_ts * self.fs + offset_samples).tolist()

            self.all_currents.append(current)
            self.all_depths.append(depth_index)
            self.all_stim_timestamps.append(stim_ts_samples)
            self.all_stim_start_stop.append(stim_ss_samples)
            self.all_image_start_stop.append(img_ss_samples)
            self.all_frame_timestamps.append(frame_ts_samples)
            self.all_behave_responses.append(response)
            self.all_response_times.append(response_time)

    # ─── Recording handling ───────────────────────────────────

    def get_start_times(self):
        """Compute cumulative start times for each NS5 block."""
        self.start_times = []
        cumulative = 0.0
        for ns5_file in self.ns5_files:
            self.start_times.append(cumulative)
            if ns5_file is not None:
                fpath = str(self.experiment_folder / ns5_file)
                rec = se.read_blackrock(file_path=fpath, stream_id="5")
                cumulative += rec.get_num_frames() / self.fs

    def append_recordings(self):
        """Load NS5 recordings into list and compute start times."""
        self.recordings_list = []
        for ns5_file in self.ns5_files:
            if ns5_file is None:
                continue
            fpath = str(self.experiment_folder / ns5_file)
            rec = se.read_blackrock(file_path=fpath, stream_id="5")
            self.recordings_list.append(rec)

        if not self.recordings_list:
            raise RuntimeError("No NS5 recordings found")

        # Compute start times (cumulative durations)
        durations = [r.get_num_frames() / self.fs for r in self.recordings_list]
        self.start_times = [0.0] + list(np.cumsum(durations[:-1]))

    def get_concatenated_recording(self):
        """Concatenate recordings into single segment and slice to 32 channels."""
        self.recording = si.concatenate_recordings(self.recordings_list)
        self.channel_ids = self.recording.get_channel_ids()
        if hasattr(self.recording, 'select_channels'):
            self.recording = self.recording.select_channels(self.channel_ids[:32])
        else:
            self.recording = self.recording.channel_slice(self.channel_ids[:32])

    def attach_probe(self):
        probe_path = os.path.join(
            os.path.dirname(__file__), "net32Ch.json")
        pi = read_probeinterface(probe_path)
        probe = pi.probes[0]
        self.recording.set_probe(probe, in_place=True)

    def remove_high_impedance_ch(self):
        ia = ImpedanceAnalyzer()
        ia.get_intan_impedances(
            animal_id=self.animalID, server_mount_drive="S:",
            session_date=self.date)
        _, good_channels = ia.get_good_impedances(threshold=2e6)
        good_ripple_channels = ia.intan_to_ripple(good_channels)
        good_indices = sorted(good_ripple_channels - 1)
        channel_ids = self.recording.channel_ids
        if hasattr(self.recording, 'select_channels'):
            self.recording = self.recording.select_channels(
                channel_ids=channel_ids[good_indices])
        else:
            self.recording = self.recording.channel_slice(
                channel_ids=channel_ids[good_indices])

    # ─── Trial dataframe ─────────────────────────────────────

    def get_trial_dataframe(self):
        """Build trial_df matching DataLoader schema."""
        num_trials = len(self.all_currents)
        self.trial_df = pd.DataFrame({
            "trial": range(1, num_trials + 1),
            "response": self.all_behave_responses,
            "response_time": self.all_response_times,
            "current": self.all_currents,
            "channel": self.all_depths,
            "stim_timestamps": self.all_stim_timestamps,
            "frame_timestamps": self.all_frame_timestamps,
            "image_start": [ss[0] if len(ss) >= 2 else np.nan
                            for ss in self.all_image_start_stop],
            "image_stop": [ss[1] if len(ss) >= 2 else np.nan
                           for ss in self.all_image_start_stop],
        })
        return self.trial_df

    # ─── Convenience methods (match DataLoader API) ──────────

    def get_hit_trials_stim_timestamps(self):
        self.hit_trials_stim_timestamps = []
        for _, trial in self.trial_df.iterrows():
            if trial["response"] == 1 and trial["current"] > 0:
                self.hit_trials_stim_timestamps.append(
                    trial["stim_timestamps"])
        self.hit_trials_df = self.trial_df[
            (self.trial_df["current"] > 0) &
            (self.trial_df["response"] == 1)
        ].reset_index(drop=True)

    def get_miss_trials_stim_timestamps(self):
        self.miss_trials_stim_timestamps = []
        for _, trial in self.trial_df.iterrows():
            if trial["response"] == 0 and trial["current"] > 0:
                self.miss_trials_stim_timestamps.append(
                    trial["stim_timestamps"])
        self.miss_trials_df = self.trial_df[
            (self.trial_df["current"] > 0) &
            (self.trial_df["response"] == 0)
        ].reset_index(drop=True)

    def get_save_folder(self):
        """Set save folder to batch_sort/ inside the session."""
        self.save_folder = self.data_folder / "batch_sort"
        print(f"Using save folder: {self.save_folder}")


# ─── Top-level loader function ────────────────────────────────

def load_mouse6(
    data_folder,
    server_key_root=None,
    get_recording=True,
    make_folder=False,
    save_folder_name=None,
):
    """
    Load ICMS83 session — drop-in replacement for load_data().

    Parameters
    ----------
    data_folder : str
        Path to ICMS83 session, e.g. "E:/ICMS83/20-Jul-2023"
    server_key_root : str, optional
        Override server key root (default: S:/xl_stimulation/ICMS83/Keys)
    get_recording : bool
        Whether to load and concatenate NS5 recordings
    """
    dl = ICMS83DataLoader(data_folder, server_key_root=server_key_root)
    dl.get_ns5_and_mat_files()

    if get_recording:
        dl.append_recordings()
        dl.get_concatenated_recording()
        dl.attach_probe()
        try:
            dl.remove_high_impedance_ch()
        except Exception as e:
            print(f"  [WARN] Could not filter impedance: {e}")
    else:
        # Still need start_times for timestamp rebasing
        dl.get_start_times()

    dl.get_key_files()
    dl.get_trial_dataframe()

    if save_folder_name:
        dl.save_folder = Path(data_folder) / save_folder_name
    elif make_folder:
        dl.save_folder = Path(data_folder) / "batch_sort"
        os.makedirs(dl.save_folder, exist_ok=True)

    return dl


# ─── CLI test ─────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    session = sys.argv[1] if len(sys.argv) > 1 else "E:/ICMS83/20-Jul-2023"
    print(f"Loading ICMS83 session: {session}")

    dl = load_mouse6(session, get_recording=False)
    print(f"\nTrial DataFrame: {len(dl.trial_df)} trials")
    print(dl.trial_df.head(10))
    print(f"\nUnique currents: {sorted(dl.trial_df['current'].unique())}")
    print(f"Unique channels: {sorted(dl.trial_df['channel'].unique())}")
    print(f"\nTrials per channel:")
    print(dl.trial_df.groupby('channel').size())
    print(f"\nHit rate by current:")
    for cur in sorted(dl.trial_df['current'].unique()):
        sub = dl.trial_df[dl.trial_df['current'] == cur]
        hr = sub['response'].mean()
        print(f"  {cur} uA: {hr:.1%} ({sub['response'].sum()}/{len(sub)})")
