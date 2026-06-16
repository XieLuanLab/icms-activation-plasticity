import os
from pathlib import Path
import glob
import re
import sys

import numpy as np
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.preprocessing as sp


def discover_sessions(base_dir, allow_missing_mat=False):
    """
    Crawls the directory looking for sessions containing RepX.mat, repX.nev, repX.ns5 files.
    Groups them by their parent session folder and sorts them chronologically by repetition.
    """
    base_path = Path(base_dir)

    # Recursively hunt down all .ns5 files
    ns5_files = list(base_path.rglob("*.ns5"))

    sessions = {}
    for ns5 in ns5_files:
        exp_dir = ns5.parent
        basename = ns5.stem
        basename_lower = basename.lower()

        # Parse the 'rep' or 'channel' index to use as an identifier
        if basename_lower.startswith('rep'):
            match = re.search(r'rep(\d+)', basename_lower)
            if not match: continue
            rep_idx = int(match.group(1))
            mat_exact = ns5.with_suffix('.mat')
            mat_cap = exp_dir / f"Rep{rep_idx}.mat"
            
        elif basename_lower.startswith('ch'):
            match = re.search(r'ch(\d+)', basename_lower)
            if not match: continue
            rep_idx = int(match.group(1)) + 1000  # Offset to prevent collision
            mat_exact = ns5.with_suffix('.mat')
            mat_cap = exp_dir / f"Ch{match.group(1)}.mat"
            
        else:
            continue

        nev_file = ns5.with_suffix('.nev')

        has_mat = mat_exact.exists()
        if not has_mat and mat_cap.exists():
            has_mat = True
            mat_file = mat_cap
        else:
            mat_file = mat_exact

        # Verify that the necessary trio of files actually exist (or allow missing mat)
        if nev_file.exists() and (has_mat or allow_missing_mat):
            if exp_dir not in sessions:
                sessions[exp_dir] = []

            sessions[exp_dir].append({
                "rep_num": rep_idx,
                "ns5_path": str(ns5),
                "nev_path": str(nev_file),
                "mat_path": str(mat_file) if has_mat else None
            })

    # Sort the repetitions within each session so that Rep1 comes before Rep2
    for exp_dir in sessions:
        sessions[exp_dir].sort(key=lambda x: x["rep_num"])

    return sessions


def load_session_reps_as_multisegment(session_reps, probe_json_path=None):
    """
    Takes a list of valid Rep dictionary objects for a single session,
    loads them individually using SpikeInterface's Neo/Blackrock extractor,
    attaches the probe layout, and appends them into a single multi-segment recording.
    """
    # Default to icms48_probe.json in the same directory as this script
    if probe_json_path is None:
        probe_json_path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), 'control_probe.json')
    import probeinterface as pi

    # Load the probe group layout from JSON once
    probegroup = pi.read_probeinterface(probe_json_path)

    recordings = []

    for rep in session_reps:
        print(f"  [+] Loading Rep {rep['rep_num']} (.ns5 and .nev)")

        try:
            # Note: you may need to specify stream_id optionally depending on how Blackrock encoded your analog signals!
            # Typically stream_id='6' or '5' for raw 30kHz un-filtered spiking broadband
            rec = se.read_blackrock(rep['ns5_path'])

            # Attach the probe! Doing this at the segment level ensures correct channel wiring early
            rec = rec.set_probegroup(probegroup)

            recordings.append(rec)
        except Exception as e:
            print(f"      [!] Failed to load Rep {rep['rep_num']}: {e}")

    if not recordings:
        return None

    # Combines multiple file repetitions into a single Multi-Segment Recording natively
    multiseg_recording = si.append_recordings(recordings)
    return multiseg_recording


def plot_segment_traces(recording, session_df, segment_index=None, channel=None, current_uA=None, trial_index=None):
    """
    Plots a short trace window centered on a stim trial, with vertical lines at each pulse.

    Parameters
    ----------
    channel : int, optional
        Filter to a specific electrode channel.
    current_uA : float, optional
        Filter to a specific current amplitude.
    trial_index : int, optional
        Directly pick a row index from the (filtered) dataframe.
        Defaults to 0 (first matching trial).
    """
    import matplotlib.pyplot as plt
    import spikeinterface.widgets as sw

    df = session_df.copy()

    # Apply optional filters
    if channel is not None:
        df = df[df['channel'] == channel]
    if current_uA is not None:
        df = df[df['current_uA'] == current_uA]

    # Only keep trials that actually have pulses
    df = df[df['stim_timestamps'].apply(len) > 0].reset_index(drop=True)

    if df.empty:
        print(
            f"No trials with pulses found for channel={channel}, current={current_uA}.")
        print(f"Available conditions:")
        summary = session_df[session_df['stim_timestamps'].apply(
            len) > 0].groupby(['channel', 'current_uA']).size()
        print(summary.to_string())
        return

    idx = trial_index if trial_index is not None else 0
    if idx >= len(df):
        print(
            f"trial_index={idx} out of range, only {len(df)} matching trials.")
        return

    trial = df.iloc[idx]
    chan = trial['channel']
    cur = trial['current_uA']
    pulses = trial['stim_timestamps']

    # Auto-select the correct segment from the rep number
    # Each rep is a separate segment: rep 1 = segment 0, rep 2 = segment 1, etc.
    seg_idx = int(trial['rep_num']) - \
        1 if segment_index is None else segment_index

    # Plot 0.3 seconds total to clearly show a healthy number of pulses
    start_time = pulses[0] - 0.05
    end_time = pulses[0] + 0.25

    # Check if probe is attached for depth ordering
    has_probe = recording.get_property('contact_vector') is not None
    if has_probe:
        print(f"  Probe attached — using SI-native depth ordering (shallow -> deep).")
    else:
        print(f"  [!] No probe attached — channels will be in default order.")

    print(f"Plotting: Ch {chan}, {cur} uA, rep {trial['rep_num']} (segment {seg_idx}), trial #{idx} "
          f"(Time: {start_time:.3f}s -> {end_time:.3f}s, {len(pulses)} pulses)")

    fig, ax = plt.subplots(figsize=(12, 8))
    sw.plot_traces(
        recording,
        segment_index=seg_idx,
        time_range=(start_time, end_time),
        order_channel_by_depth=has_probe,
        return_in_uV=True,
        backend="matplotlib",
        ax=ax
    )

    # Overlay the individual extracted pulses as proof!
    for p in pulses:
        ax.axvline(x=p, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)

    plt.title(
        f"Rep {trial['rep_num']} (seg {seg_idx}): {cur} uA on Ch {chan}\n(Blue lines = {len(pulses)} extracted pulses)")
    plt.show()


def extract_stim_packets_from_nev(nev_path):
    """
    Reads stim pulse packets directly from a Blackrock/Ripple NEV binary file.
    No cache files, no dependencies beyond numpy + struct.

    Returns
    -------
    timestamps_s : np.ndarray (float64)
        Pulse times in seconds (30 kHz clock / 30000).
    packet_ids : np.ndarray (uint16)
        Raw PacketIDs.  Stim channels are electrode + 5120.
    """
    import struct as _struct
    import os

    # Check if the file is empty or too small to even contain a basic header
    if os.path.getsize(nev_path) < 20:
        return np.array([], dtype=np.float64), np.array([], dtype=np.uint16)

    with open(nev_path, 'rb') as f:
        # NEV basic header layout (after 8-byte FileTypeID + 4-byte FileSpec/Flags):
        #   Offset 12: Bytes in Headers   (uint32) - total size of basic + extended headers
        #   Offset 16: Bytes in Data Packet (uint32) - size of each data packet
        f.seek(12)
        # offset 12: total header bytes
        header_bytes = _struct.unpack('<I', f.read(4))[0]
        # offset 16: bytes per data packet)
        bytes_in_data = _struct.unpack('<I', f.read(4))[0]

        f.seek(header_bytes)
        raw = f.read()

    num_packets = len(raw) // bytes_in_data
    if num_packets == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.uint16)

    # Structured dtype: only first 6 bytes matter (timestamp u4 + packet_id u2)
    dt = np.dtype([
        ('timestamp', '<u4'),
        ('packet_id', '<u2'),
        ('_pad', f'V{bytes_in_data - 6}')
    ])
    packets = np.frombuffer(raw[:num_packets * bytes_in_data], dtype=dt)

    timestamps_s = packets['timestamp'].astype(np.float64) / 30000.0
    packet_ids = packets['packet_id']
    return timestamps_s, packet_ids


def extract_trial_parameters_from_rep(rep_dict, survey_current_uA=5.0):
    """
    Parses the .nev events and .mat replicant array to build a synchronized trial dataframe
    for a single repetition segment.

    Parameters
    ----------
    rep_dict : dict
        Repetition dictionary from discover_sessions.
    survey_current_uA : float
        Current to assign to survey sessions that lack a .mat file (default 5.0).

    Pipeline
    --------
    1. Load the replicant matrix from the .mat file  -> (current, channel) per trial.
    2. Use Neo to pull the TTL envelope rising-edges  -> trial start times.
    3. Directly parse the NEV binary (no cache files)  -> per-pulse stim timestamps.
    4. For each trial, filter stim packets by electrode PacketID (ch + 5120)
       within a time window around the trial start to get the individual 50 Hz pulses.
    """
    import pandas as pd
    
    if rep_dict.get('mat_path') is not None:
        import scipy.io as sio
        import neo

        mat = sio.loadmat(rep_dict['mat_path'], simplify_cells=True)
        replicant = mat.get('replicant')
        
        # replicant is an (N, 2) array: [Channel, Current]
        channels = [int(row[0]) for row in replicant]
        currents = [float(row[1]) for row in replicant]

        # ----- 1. Trial start times from Neo (TTL envelope on analog_input_channel_1) -----
        reader = neo.io.BlackrockIO(filename=rep_dict['nev_path'])
        blk = reader.read_block()
        seg = blk.segments[0]

        ev_chan = next((e for e in seg.events if e.name == 'analog_input_channel_1'), None)
        stim_start_times = []
        if ev_chan is not None:
            for time, label in zip(ev_chan.times, ev_chan.labels):
                if label == '32767':  # Rising edge of trial envelope
                    stim_start_times.append(float(time))

        if len(stim_start_times) != len(currents):
            print(f"      [!] Warning in Rep {rep_dict['rep_num']}: "
                  f"{len(stim_start_times)} TTLs vs {len(currents)} MAT trials.")

        # ----- 2. Per-pulse stim timestamps from direct NEV binary parse -----
        all_timestamps_s, all_packet_ids = extract_stim_packets_from_nev(rep_dict['nev_path'])

        pulse_times_list = []
        for current, chan, start_t in zip(currents, channels, stim_start_times):
            if current == 0 or start_t is None:
                pulse_times_list.append([])
                continue

            # Blackrock/Ripple stim channels use PacketID = electrode + 5120
            elec_packet_id = chan + 5120

            # Time window: 50 ms before trial start to 1.5 s after
            mask = (
                (all_packet_ids == elec_packet_id)
                & (all_timestamps_s >= start_t - 0.05)
                & (all_timestamps_s <= start_t + 1.5)
            )
            pulse_times_list.append(all_timestamps_s[mask].tolist())

        # ----- 3. Assemble DataFrame -----
        n_trials = len(currents)
        padded_starts = (
            stim_start_times[:n_trials]
            if len(stim_start_times) >= n_trials
            else stim_start_times + [None] * (n_trials - len(stim_start_times))
        )

        df = pd.DataFrame({
            'rep_num':         rep_dict['rep_num'],
            'trial':           range(1, n_trials + 1),
            'current_uA':      currents,
            'channel':         channels,
            'stim_start_s':    padded_starts,
            'stim_timestamps': pulse_times_list,
        })

        # Drop the 0 uA catch trials
        df = df[df['current_uA'] > 0].reset_index(drop=True)
        return df
        
    else:
        # ----- FALLBACK: Synthesis from just .nev gaps -----
        print(f"      [!] No .mat found for Rep {rep_dict['rep_num']}. Synthesizing trials from .nev gaps...")
        all_timestamps_s, all_packet_ids = extract_stim_packets_from_nev(rep_dict['nev_path'])

        if len(all_timestamps_s) == 0:
            return pd.DataFrame(columns=['rep_num', 'trial', 'current_uA', 'channel', 'stim_start_s', 'stim_timestamps'])

        diffs = np.diff(all_timestamps_s)
        # Assuming >0.5s gaps represent separate trials
        trial_gaps_idx = np.where(diffs > 0.5)[0]
        trial_starts_idx = np.concatenate(([0], trial_gaps_idx + 1))

        channels = []
        currents = []
        stim_start_times = []
        pulse_times_list = []

        from collections import Counter
        
        for i, start_idx in enumerate(trial_starts_idx):
            end_idx = trial_starts_idx[i+1] if i + 1 < len(trial_starts_idx) else len(all_timestamps_s)

            trial_ts = all_timestamps_s[start_idx:end_idx]
            trial_pids = all_packet_ids[start_idx:end_idx]

            stim_pids = trial_pids[trial_pids >= 5120]
            if len(stim_pids) > 0:
                chan = Counter(stim_pids).most_common(1)[0][0] - 5120
            else:
                chan = 1

            channels.append(int(chan))
            currents.append(survey_current_uA)
            stim_start_times.append(trial_ts[0])
            pulse_times_list.append(trial_ts.tolist())

        df = pd.DataFrame({
            'rep_num':         rep_dict['rep_num'],
            'trial':           range(1, len(currents) + 1),
            'current_uA':      currents,
            'channel':         channels,
            'stim_start_s':    stim_start_times,
            'stim_timestamps': pulse_times_list,
        })

        return df


if __name__ == "__main__":
    import pandas as pd

    # Test execution
    base_folder = "E:/ICMS plasticity control/ICMS48"
    sessions = discover_sessions(base_folder)

    print(f"Found {len(sessions)} valid session folders!")

    for exp_dir, reps in sessions.items():
        print(f"\nProcessing Folder: {exp_dir}")
        print(f"  -> Found {len(reps)} complete Reps (NS5, NEV, MAT)")

        # Only load first rep for faster iteration
        reps = reps[:1]

        # 1. Parse trial parameters from .nev and .mat into a DataFrame
        all_trials_dfs = []
        for rep in reps:
            df = extract_trial_parameters_from_rep(rep)
            all_trials_dfs.append(df)

        session_df = pd.concat(all_trials_dfs, ignore_index=True)
        print("\n  -> Trial Extractor Output (First 5 Trials):")
        print(session_df.head(5).to_string())

        # 2. Build the recording object for SpikeInterface (single rep)
        multi_rec = load_session_reps_as_multisegment(reps)

        if multi_rec:
            print(f"\n  -> Successfully created recording:")
            print(f"     Channels: {multi_rec.get_num_channels()}")
            print(f"     Segments: {multi_rec.get_num_segments()}")
            print(
                f"     Sampling Rate: {multi_rec.get_sampling_frequency()} Hz")

            # Plot raw traces
            print("\n  -> Plotting RAW traces:")
            plot_segment_traces(multi_rec, session_df)

            # Add project root to path for preprocessing imports
            _project_root = os.path.dirname(
                os.path.dirname(os.path.abspath(__file__)))
            if _project_root not in sys.path:
                sys.path.insert(0, _project_root)

            from preprocessing.icms_pipeline import preprocess_icms, stim_timestamps_seconds_to_samples
            from preprocessing.curate import curate_units

            print("\n  -> Running preprocessing pipeline...")
            stim_samples = stim_timestamps_seconds_to_samples(
                session_df, fs=multi_rec.get_sampling_frequency())
            # Controls were stimulated at 50 Hz (vs 100 Hz for the stim cohort),
            # so the artifact-subtraction window must use the 50 Hz inter-pulse gap.
            rec_sorted, rec_wvf = preprocess_icms(multi_rec, stim_samples, stim_freq_hz=50)

            print(f"     Preprocessed recording: {rec_sorted}")
            print("\n  -> Plotting PREPROCESSED traces:")
            plot_segment_traces(rec_wvf, session_df)

            print("\n  -> Materializing preprocessed traces...")
            save_dir = Path(exp_dir) / "spike_sorting"
            save_dir.mkdir(exist_ok=True)
            preproc_zarr = save_dir / "preprocessed.zarr"

            job_kwargs = dict(n_jobs=-1, chunk_duration="1s",
                              progress_bar=True)

            from numcodecs import Blosc
            compressor = Blosc(cname='zstd', clevel=5, shuffle=Blosc.BITSHUFFLE)
            rec_disk = rec_wvf.save(folder=preproc_zarr, format="zarr", overwrite=True,
                                    compressor=compressor, **job_kwargs)
            print(f"     Saved to: {preproc_zarr}")

            # 5. Spike sort with MountainSort5
            print("\n  -> Running MountainSort5...")
            ms5_params = {
                "scheme": "2",
                "detect_threshold": 4,
                "npca_per_channel": 3,
                "npca_per_subdivision": 10,
                "snippet_mask_radius": 0,
                "scheme2_detect_channel_radius": 200,
                "scheme2_training_duration_sec": 300,
                "filter": False,   # already filtered
                "whiten": True,    # MS5 handles whitening internally
            }

            sorting = ss.run_sorter(
                sorter_name="mountainsort5",
                recording=sp.scale(rec_disk, dtype="float"),
                remove_existing_folder=True,
                verbose=True,
                **ms5_params,
            )
            print(f"     Found {len(sorting.get_unit_ids())} units")

            # 6. Create SortingAnalyzer
            print("\n  -> Creating SortingAnalyzer...")
            save_dir = Path(exp_dir) / "spike_sorting"
            save_dir.mkdir(exist_ok=True)

            job_kwargs = dict(n_jobs=1, chunk_duration="1s", progress_bar=True)

            analyzer = si.create_sorting_analyzer(
                sorting=sorting,
                recording=rec_disk,
                format="zarr",
                folder=save_dir / "analyzer_raw.zarr",
                sparse=False,
                overwrite=True,
            )

            extensions_to_compute = [
                "random_spikes",
                "waveforms",
                "templates",
                "correlograms",
            ]
            extension_params = {"correlograms": {"window_ms": 100}}
            analyzer.compute(
                extensions_to_compute,
                extension_params=extension_params,
                **job_kwargs,
            )

            print(f"\n  -> Sorting complete! Results saved to: {save_dir}")
            print(f"     Units: {len(sorting.get_unit_ids())}")
            print(f"     Analyzer: {analyzer.folder}")

            # 7. Curate units (remove stim artifacts and noise)
            print("\n  -> Curating units...")
            good_ids, bad_ids, labels = curate_units(analyzer, min_spikes=100)

            # Save curated analyzer
            curated_analyzer = analyzer.select_units(
                unit_ids=good_ids,
                format="zarr",
                folder=save_dir / "analyzer_curated.zarr",
            )
            curated_analyzer.compute(
                ["random_spikes", "waveforms", "templates", "correlograms"],
                extension_params={"correlograms": {"window_ms": 100}},
                **job_kwargs,
            )
            print(f"     Curated analyzer saved: {curated_analyzer.folder}")

            # 8. Plot units to PPT (sparse analyzer required)
            print("\n  -> Creating sparse analyzer and plotting to PPT...")
            sparse_curated = curated_analyzer.copy()
            sparse_curated.sparsity = si.compute_sparsity(
                curated_analyzer, method="radius", radius_um=60
            )
            sparse_curated.compute(
                ["random_spikes", "waveforms", "templates", "correlograms"],
                extension_params={"correlograms": {"window_ms": 100}},
                **job_kwargs,
            )

            from batch_process.util.plotting import plot_units_in_batches
            plot_units_in_batches(
                sparse_curated, save_dir=save_dir, ppt_name="curated_units"
            )
            print(f"     PPT saved to: {save_dir / 'curated_units.pptx'}")

            # 9. View curated results in browser with sortingview
            print("\n  -> Launching sortingview (curated) in browser...")
            import spikeinterface.widgets as sw
            os.environ['KACHERY_ZONE'] = 'default'
            
            sw.plot_sorting_summary(curated_analyzer, backend="sortingview")

        # Break after testing the first session folder
        break
