from batch_process.util.ppt_image_inserter import PPTImageInserter
import numpy as np
from spikeinterface import full as si
from pathlib import Path
from spikeinterface.postprocessing import compute_correlograms

# from batch_process.util.plotting import *
from batch_process.util.file_util import *
from batch_process.util.plotting import plot_units_in_batches
import batch_process.util.template_util as template_util
# Lazy import: load_data pulls in psignifit which isn't in all envs.
# Only needed by load_raw_analyzer(), imported there instead.
import umap
from isosplit6 import isosplit6

import matplotlib.pyplot as plt
# Lazy import: WaveformClassifier needs xgboost, only used by remove_bad_waveforms_C

# %%
# analyzer.sparsity = None
# analyzer.compute("templates")
# a = analyzer.get_extension("templates")
# plt.plot(a.get_unit_template(unit_id=2))


# TEST
# unit_id = 1
# analyzer.sparsity = None
# analyzer.compute("random_spikes")
# analyzer.compute("waveforms")
# analyzer.compute("templates")

# dense_analyzer = analyzer.copy()

# sparse_analyzer = analyzer.copy()
# sparse_analyzer.sparsity = si.compute_sparsity(
#     sparse_analyzer, method="radius", radius_um=60)
# sparse_analyzer.compute("random_spikes")
# sparse_analyzer.compute("waveforms")
# sparse_analyzer.compute("templates")
# sparse_templates = sparse_analyzer.get_extension("templates")

# plt.plot(sparse_templates.get_unit_template(unit_id=unit_id))

# plt.figure()
# plt.title("Dense template")
# dense_templates = dense_analyzer.get_extension("templates")

# plt.plot(dense_templates.get_unit_template(unit_id=unit_id))

# # Test extremum ch id
# template_extremum_ch = si.get_template_extremum_channel(
#     sparse_analyzer, outputs="id")  # 2
# unit_id_to_ch_ids = sparse_analyzer.sparsity.unit_id_to_channel_ids
# primary_ch = template_extremum_ch[unit_id]
# ch_ids = unit_id_to_ch_ids[unit_id]
# all_ch_ids = analyzer.channel_ids
# np.where(ch_ids == str(primary_ch))[0][0]
# np.where(all_ch_ids == str(primary_ch))[0][0]

# plt.plot(templates.get_unit_template(unit_id=unit_id))


def sparse_template_max_amplitude(template):
    max_amplitudes = np.min(template, axis=0)
    ts = np.argmin(template, axis=0)
    return max_amplitudes, ts


# TODO
def get_ts_offsets(sparse_analyzer):
    # analyzer sparsity needs to have radius of 60 um
    unit_ids = sparse_analyzer.unit_ids
    templates_ext = sparse_analyzer.get_extension("templates")

    ts_offsets = {}
    for unit_id in unit_ids:
        template = templates_ext.get_unit_template(unit_id)
        _, ts = sparse_template_max_amplitude(template)
        primary_ts = ts[0]  # Timestamp of the primary channel
        offsets = ts - primary_ts  # Calculate offsets relative to the primary channel
        ts_offsets[unit_id] = np.mean(np.abs(offsets[1:]))
    return ts_offsets


def classify_waveform(new_wvf):
    point1 = new_wvf[6]
    point2 = new_wvf[74]
    return -60 <= point1 <= 100 and -21 <= point2 <= 35


# TODO
def sort_units_by_primary_channel(analyzer, unit_ids_subset):
    ch_loc = analyzer.get_channel_locations()[:, 1]
    ch_ids = analyzer.channel_ids
    ch_idx = analyzer.channel_ids_to_indices(ch_ids)
    ch_idx_to_loc_dict = {idx: loc for idx, loc in zip(ch_idx, ch_loc)}
    unit_id_to_ch_idx_dict = si.get_template_extremum_channel(
        analyzer, outputs="index")
    # Get the extremum channel location for each unit in the subset
    unit_id_to_extremum_loc_dict = {
        unit_id: ch_idx_to_loc_dict[unit_id_to_ch_idx_dict[unit_id]] for unit_id in unit_ids_subset
    }

    # Sort the unit IDs based on the extremum channel location
    sorted_unit_ids = sorted(
        unit_id_to_extremum_loc_dict.keys(),
        key=lambda unit_id: unit_id_to_extremum_loc_dict[unit_id],
    )[::-1]

    return sorted_unit_ids


def mean_norm_peak(dense_template, sparse_ch_indices, min_idx=30, range_offset=10):
    max_peaks = np.max(
        np.abs(dense_template[min_idx -
               range_offset: min_idx + range_offset + 1, :]),
        axis=0,
    )
    all_ch_indices = np.arange(dense_template.shape[1])
    # Exclude sparse template channels
    non_sparse_ch_indices = np.setdiff1d(all_ch_indices, sparse_ch_indices)
    max_peaks_non_sparse = max_peaks[non_sparse_ch_indices]
    template_peak = np.min(dense_template[min_idx, :])
    normalized_peaks_non_sparse = max_peaks_non_sparse / np.abs(template_peak)
    top_peaks = np.sort(normalized_peaks_non_sparse)[-5:]
    mnr = np.mean(top_peaks)

    return mnr


def exclude_artifact_unit(unit_id, analyzer):
    # analyzer must be dense
    num_channels = analyzer.get_num_channels()
    sparsity = si.compute_sparsity(analyzer, method="radius", radius_um=120)
    unit_id_to_ch_ind = sparsity.unit_id_to_channel_indices
    sparse_ch_indices = unit_id_to_ch_ind[unit_id]
    templates = analyzer.get_extension("templates")
    dense_template = templates.get_unit_template(unit_id)
    assert dense_template.shape[1] == num_channels, f"Dense template ch num is not {num_channels}!"
    mnr = mean_norm_peak(dense_template, sparse_ch_indices)
    return mnr > 0.3


def curate_units(analyzer, sparse_analyzer, min_spikes=100, job_kwargs=None):
    """
    Curate units after sorting to remove stim/imaging artifacts.
    Look at:
        1) acg 30 ms component
        2) shape of waveform
        3) mean of normalized peaks far from primary channel
        4) spike count

    Returns good and bad unit ids
    """
    if job_kwargs is None:
        job_kwargs = {}
    unit_ids = analyzer.unit_ids
    # ccgs, bins = analyzer.compute("correlograms", window_ms=100).get_data()

    ccgs, bins = analyzer.get_extension("correlograms").get_data()

    # ccgs, bins = compute_correlograms(we, window_ms=100)
    bins = bins[50:]  # get only positive values since symmetric
    ccgs = ccgs[:, :, 50:]
    acg_range = np.where((bins > 31) & (bins < 35))[0]
    ts_offset_threshold = 10  # samples, max mean sample offset from primary ch peak
    ts_offsets = get_ts_offsets(sparse_analyzer)

    good_ids1 = []
    artifacts = []

    for i, unit_id in enumerate(unit_ids):
        acg = ccgs[i, i, :]
        sort_idx = np.argsort(acg)[::-1]
        # Exclude artifact unit check
        if exclude_artifact_unit(unit_id, analyzer):
            artifacts.append(unit_id)
            continue

        if np.any(np.isin(sort_idx[:5], acg_range)) & (ts_offsets[unit_id] > ts_offset_threshold):
            artifacts.append(unit_id)
        else:
            good_ids1.append(unit_id)

    good_ids2 = []
    noise = []

    # classify by primary waveform
    unit_primary_ch_wvfs_dict = template_util.get_unit_primary_ch_templates(
        analyzer)
    for unit_id in good_ids1:
        template = unit_primary_ch_wvfs_dict[unit_id]
        if classify_waveform(template):
            good_ids2.append(unit_id)
        else:
            noise.append(unit_id)

    good_ids = sort_units_by_primary_channel(analyzer, good_ids2)
    bad_ids = sort_units_by_primary_channel(analyzer, artifacts + noise)

    return good_ids, bad_ids


def curate_and_plot(analyzer, parent_folder, subfolder_path, save_path=None, curate_only=False):
    sorting = analyzer.sorting

    good_ids, bad_ids = curate_units(analyzer)
    good_sorting = sorting.select_units(
        unit_ids=good_ids, renamed_unit_ids=good_ids)
    bad_sorting = sorting.select_units(
        unit_ids=bad_ids, renamed_unit_ids=bad_ids)

    if save_path:
        folder = Path(parent_folder) / save_path
        if folder.exists() and folder.is_dir():
            shutil.rmtree(folder)
        good_we = analyzer.select_units(
            unit_ids=good_ids, format="zarr", folder=folder)
    else:
        good_analyzer = analyzer.select_units(unit_ids=good_ids)
    bad_analyzer = analyzer.select_units(unit_ids=bad_ids)

    if not curate_only:
        # Check if subfolder_path is provided and non-empty
        if not subfolder_path:
            raise ValueError("subfolder_path cannot be empty.")

        results_path = Path(parent_folder) / subfolder_path
        if results_path.exists() and results_path.is_dir():
            if results_path == Path(parent_folder):
                raise ValueError(
                    "Cannot delete the entire parent folder. Please provide a valid subfolder path.")
            # Ensure only subfolder contents are deleted
            shutil.rmtree(results_path)
        results_path.mkdir(parents=True, exist_ok=True)

        plot_units_in_batches(
            good_analyzer, save_dir=results_path, ppt_name="curated")

    return bad_analyzer, bad_ids, good_analyzer


def remove_bad_waveforms_A(analyzer, unit_id, plot_flag=False):
    """
    Removes bad waveforms based on multiple criteria.

    Parameters:
    - analyzer: SpikeInterface sorting analyzer.
    - unit_id: ID of the unit to process.
    - plot_flag (bool): If True, validates and plots results.

    Returns:
    - good_indices (np.ndarray): Indices of good waveforms.
    - bad_indices (np.ndarray): Indices of bad waveforms.
    """
    # Load waveforms and channel info
    wvf_ext = analyzer.get_extension("waveforms")
    wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)
    num_wvfs = wvfs.shape[0]
    nbefore = wvf_ext.nbefore

    primary_wvfs = template_util.get_unit_primary_ch_wvfs(analyzer, unit_id)
    primary_ch_idx_dense = si.get_template_extremum_channel(
        analyzer, outputs="index")[unit_id]

    # Check if the max amplitude channel matches the template primary channel
    max_indices = np.argmin(wvfs[:, nbefore, :], axis=1)
    has_incorrect_primary_ch = max_indices != primary_ch_idx_dense

    # Criteria 2: Three consecutive zeros in specific ranges
    def detect_consecutive_zeros(waveforms, start, end):
        is_zero = waveforms[:, start:end] == 0
        return np.array([
            np.any(np.convolve(row.astype(int), [1, 1, 1], mode="valid") == 3)
            for row in is_zero
        ])

    three_consecutive_zeros_24_30 = detect_consecutive_zeros(
        primary_wvfs, 24, 31)
    three_consecutive_zeros_32_38 = detect_consecutive_zeros(
        primary_wvfs, 32, 38)

    # Criteria 3: High amplitude in pre-peak range
    high_amplitude_mask = np.any(primary_wvfs[:, 29:32] > -40, axis=1)

    # Combine masks
    combined_mask = (
        # has_incorrect_primary_ch |
        three_consecutive_zeros_24_30
        | three_consecutive_zeros_32_38
        | high_amplitude_mask
    )

    good_indices = np.where(~combined_mask)[0]
    bad_indices = np.where(combined_mask)[0]

    # Optional validation and plotting
    if plot_flag:
        validate_remove_bad_waveforms(
            analyzer, unit_id, primary_ch_idx_dense, good_indices, bad_indices)

    return good_indices, bad_indices


def remove_bad_waveforms_B(analyzer, unit_id, hard_threshold=-50, num_stds=10):
    """
    Reject outliers using the median as the central tendency and 3 SD for range,
    with a hard upper threshold. Optionally plot the histogram with bounds.
    """
    amp_ext = analyzer.get_extension(
        "spike_amplitudes").get_data(outputs="by_unit")[0]
    amplitudes = amp_ext[unit_id]
    # Calculate median and standard deviation
    median = np.median(amplitudes)
    std = np.std(amplitudes)

    # Define bounds: median ± num_stds SD
    lower_bound = median - num_stds * std
    upper_bound = median + num_stds * std

    # Combine hybrid filtering with hard threshold
    good_indices = np.where(
        (amplitudes >= lower_bound) &
        (amplitudes <= upper_bound) &
        (amplitudes <= hard_threshold)
    )[0]

    bad_indices = np.setdiff1d(np.arange(len(amplitudes)), good_indices)

    return good_indices, bad_indices


def get_unit_id_waveform_dict(raw_analyzer, analyzer):
    unit_id_waveform_dict = {}
    for unit_id in analyzer.unit_ids:
        raw_primary_ch_wvfs = template_util.get_raw_unit_primary_ch_wvfs(
            raw_analyzer, analyzer, unit_id)
        filt_primary_ch_wvfs = template_util.get_unit_primary_ch_wvfs(
            analyzer, unit_id)

        centered_raw_wvfs = raw_primary_ch_wvfs - \
            raw_primary_ch_wvfs[:, 30].reshape(-1, 1)

        # Store both 'raw' and 'filt' waveforms in a dictionary per unit ID
        unit_id_waveform_dict[unit_id] = {
            'raw': centered_raw_wvfs,
            'filt': filt_primary_ch_wvfs
        }
    return unit_id_waveform_dict


def remove_bad_waveforms_C(raw_analyzer, analyzer, unit_id, confidence_threshold=0.95):
    """
    Classify waveforms as good (spike) or bad (not_spike) using inference results,
    with an option to filter bad indices by confidence.

    Parameters:
    - unit_id_waveform_dict (dict): Dictionary containing waveforms for each unit.
    - unit_id (int): The unit ID to process.
    - confidence_threshold (float): Maximum confidence for bad indices to be included.

    Returns:
    - good_indices (np.ndarray): Indices of waveforms classified as spikes.
    - bad_indices (np.ndarray): Indices of waveforms classified as not_spikes, filtered by confidence.
    - confidences (pd.Series): Confidence scores for all waveforms.
    """
    # Initialize the WaveformClassifier (lazy import — needs xgboost)
    from batch_process.util.waveform_classifier import WaveformClassifier
    unit_id_waveform_dict = get_unit_id_waveform_dict(raw_analyzer, analyzer)
    wfc = WaveformClassifier(unit_id_waveform_dict)

    # Perform inference if not already done
    # Use re_infer=True if you want to recalculate
    df = wfc.infer(re_infer=False)

    # Filter results for the specified unit ID
    unit_df = df[df['Unit ID'] == unit_id]

    if unit_df.empty:
        print(f"No inference results found for Unit ID {unit_id}.")
        return np.array([]), np.array([]), pd.Series()

    # Separate indices based on predicted labels
    good_indices = unit_df[unit_df['Predicted Label']
                           == "spike"]["Waveform Index"].to_numpy()

    # Filter bad indices based on confidence threshold
    bad_indices = unit_df[
        (unit_df['Predicted Label'] == "not_spike") &
        (unit_df['Confidence'] >= confidence_threshold)
    ]["Waveform Index"].to_numpy()

    # Extract confidence scores
    confidences = unit_df["Confidence"]

    return good_indices, bad_indices, confidences


def validate_remove_bad_waveforms(analyzer, unit_id, primary_ch_idx_dense, good_indices, bad_indices):
    """
    Validate and plot good and bad waveforms on separate axes with shared y-axis limits.
    """
    wvf_ext = analyzer.get_extension("waveforms")
    wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)
    num_wvfs = wvfs.shape[0]

    # Ensure indices are within bounds
    good_indices = [idx for idx in good_indices if idx < num_wvfs]
    bad_indices = [idx for idx in bad_indices if idx < num_wvfs]

    # Truncate indices to avoid plotting too many waveforms
    good_indices_trunc = good_indices[: np.min([len(good_indices), 100])]
    bad_indices_trunc = bad_indices[: np.min([len(bad_indices), 100])]

    templates = analyzer.get_extension("templates")
    template = templates.get_unit_template(unit_id)

    # Extract waveforms for the primary channel
    good_waveforms = wvfs[good_indices_trunc, :, primary_ch_idx_dense].T
    bad_waveforms = wvfs[bad_indices_trunc, :, primary_ch_idx_dense].T

    # Calculate y-axis limits based on both good and bad waveforms
    all_waveforms = np.hstack([good_waveforms, bad_waveforms])
    y_min, y_max = np.min(all_waveforms), np.max(all_waveforms)

    # Create subplots
    fig, axes = plt.subplots(2, 1, figsize=(
        8, 6), sharey=True, constrained_layout=True)

    # Plot good waveforms
    axes[0].plot(good_waveforms, "g", linewidth=0.5, alpha=0.2)
    axes[0].set_title("Good Waveforms", fontsize=10)
    axes[0].set_ylabel("Amplitude (µV)", fontsize=8)
    axes[0].set_ylim([y_min, y_max])  # Set y-axis limits

    # Plot bad waveforms
    axes[1].plot(bad_waveforms, "r", linewidth=0.5, alpha=0.2)
    axes[1].set_title("Bad Waveforms", fontsize=10)
    axes[1].set_xlabel("Time (samples)", fontsize=8)
    axes[1].set_ylabel("Amplitude (µV)", fontsize=8)
    axes[1].set_ylim([y_min, y_max])  # Set y-axis limits

    plt.suptitle(
        f"Unit {unit_id}: Validation of Waveform Classification", fontsize=12)
    plt.show()


def load_raw_analyzer(data_folder, analyzer):
    # create raw analyzer with raw recording and preprocessed sorting object
    from util.load_data import load_data
    dataloader = load_data(
        data_folder, make_folder=True, save_folder_name="batch_sort", first_N_files=4, server_mount_drive="S:"
    )
    rec_raw = dataloader.recording

    raw_analyzer_path = Path(data_folder) / "stage3_raw_analyzer.zarr"

    # if os.path.exists(raw_analyzer_path):
    #     raw_analyzer = si.load_sorting_analyzer(raw_analyzer_path)

    # else:
    raw_analyzer = si.create_sorting_analyzer(
        sorting=analyzer.sorting,
        recording=rec_raw,
        # format="zarr",
        # folder=raw_analyzer_path,  # Todo
        sparse=False,
        overwrite=True,
        max_spikes_per_unit=None,
    )

    if len(raw_analyzer.extensions) < 4:
        raw_analyzer.compute('random_spikes', method='all')
        # raw_analyzer.compute('waveforms')
        # raw_analyzer.compute('templates', operators=['average', 'median'])
        # raw_analyzer.compute('unit_locations')

    return raw_analyzer
# %% Visualize

# ppt_inserter = PPTImageInserter(
#     grid_dims=(2, 2), spacing=(0.05, 0.05), title_font_size=16
# )

# analyzer.sparsity = None

# analyzer.compute("random_spikes", method="all")
# analyzer.compute("waveforms")
# analyzer.compute("templates")

# for unit_id in unit_ids:
#     remove_bad_waveforms_A(analyzer, unit_id, template_ch_dict, plot_flag=True)

#     img_path = "img.png"
#     plt.savefig(str(img_path))
#     ppt_inserter.add_image(str(img_path))
#     plt.close()

# ppt_inserter.save("test.pptx")
