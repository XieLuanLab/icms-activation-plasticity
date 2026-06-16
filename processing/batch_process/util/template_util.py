import batch_process.util.plotting as plot_util
import spikeinterface.full as si
import numpy as np
import os
import shutil
import matplotlib.pyplot as plt


def save_sparse_analyzer(analyzer, method="memory", radius_um=60, job_kwargs=None):
    if job_kwargs is None:
        job_kwargs = {}
    # save sparse analyzer to disk to avoid recomputing extensions
    original_path = analyzer.folder
    stem = original_path.stem
    suffix = original_path.suffix
    new_filename = f"{stem}_sparse{suffix}"
    new_path = original_path.with_name(new_filename)

    if os.path.isdir(new_path):
        shutil.rmtree(new_path)
    sparse_analyzer = analyzer.copy()
    # sparsify
    sparse_analyzer.sparsity = si.compute_sparsity(
        analyzer, method="radius", radius_um=radius_um)

    if method == "memory":
        sparse_analyzer.save_as()
    else:
        sparse_analyzer.save_as(format=method, folder=new_path)
        sparse_analyzer.folder = new_path

    extensions_to_compute = [
        "random_spikes",
        "waveforms",
        "templates",
        "template_similarity",
        "correlograms",
        "spike_amplitudes",
        "unit_locations",
    ]
    # do I need to compute all random spikes?
    # extension_params = {"random_spikes": {"method": "all"},
    #                     "unit_locations": {"method": "center_of_mass"}}
    extension_params = {"unit_locations": {"method": "center_of_mass"}}
    sparse_analyzer.compute(extensions_to_compute,
                            extension_params=extension_params, **job_kwargs)

    return sparse_analyzer


def get_dense_primary_ch_index(analyzer, unit_id):
    template_extremum_ch = si.get_template_extremum_channel(
        analyzer, outputs="id")
    all_ch_ids = analyzer.channel_ids
    primary_ch = template_extremum_ch[unit_id]
    return np.where(all_ch_ids == str(primary_ch))[0][0]


def plot_mean_waveform_with_std_shading(waveforms, color='k', ax=None):
    # Calculate mean and standard deviation
    mean_template = np.mean(waveforms, axis=0)
    std_template = np.std(waveforms, axis=0)

    # Generate x-values (assuming sampling rate of 30 Hz)
    x_values = np.arange(mean_template.shape[0]) / 30  # Adjust as needed

    # Use the provided axis, or get the current one if None is given
    if ax is None:
        ax = plt.gca()

    # Plot the mean waveform
    ax.plot(x_values, mean_template, label='Mean', color=color)

    # Plot the shaded standard deviation
    ax.fill_between(x_values, mean_template - std_template, mean_template + std_template,
                    color=color, alpha=0.2)

    # Optionally, you can customize other aspects of the plot here
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude')


def get_template_ch(analyzer, all_wvfs=True, job_kwargs=None):
    if job_kwargs is None:
        job_kwargs = {}
    original_path = analyzer.folder
    stem = original_path.stem
    suffix = original_path.suffix
    new_filename = f"{stem}_sparse{suffix}"
    new_path = original_path.with_name(new_filename)

    try:
        sparse_analyzer = si.load_sorting_analyzer(new_path)
    except AssertionError:
        print(f"Sparse analyzer not found at {new_path}. Create a sparse analyzer using template.util.save_sparse_analyzer_to_disk().")
        # sparse_analyzer = save_sparse_analyzer(analyzer)

    if sparse_analyzer.sparsity is None:
        print("Error! sparse_analyzer is not actually sparse.")
        return None

    all_ch_ids = sparse_analyzer.channel_ids
    unit_ids = sparse_analyzer.unit_ids

    unit_id_to_ch_ids = sparse_analyzer.sparsity.unit_id_to_channel_ids
    template_extremum_ch = si.get_template_extremum_channel(
        sparse_analyzer, outputs="id")

    template_ch_dict = {}
    for unit_id in unit_ids:
        ch_ids = unit_id_to_ch_ids[unit_id]
        primary_ch = template_extremum_ch[unit_id]
        template_ch_dict[unit_id] = {
            "ch_ids": ch_ids,
            "primary_ch": primary_ch,
            "primary_ch_idx_sparse": np.where(ch_ids == str(primary_ch))[0][0],
            "primary_ch_idx_dense": np.where(all_ch_ids == str(primary_ch))[0][0],
        }
    return template_ch_dict


def get_unit_primary_ch_templates(analyzer):
    # dense analyzer
    if analyzer.sparsity is not None:
        print("Error! sparse_analyzer is not dense.")
        return None
    if analyzer.get_extension('templates') is None:
        analyzer.compute_one_extension('templates')
    templates_ext = analyzer.get_extension("templates")

    template_extremum_indices = si.get_template_extremum_channel(
        analyzer, outputs="index")

    unit_primary_ch_wvfs_dict = {}
    for unit_id in analyzer.unit_ids:
        primary_ch_idx = template_extremum_indices[unit_id]
        primary_ch_wvf = templates_ext.get_unit_template(unit_id)[
            :, primary_ch_idx]
        unit_primary_ch_wvfs_dict[unit_id] = primary_ch_wvf
    return unit_primary_ch_wvfs_dict


def get_unit_primary_ch_wvfs(analyzer, unit_id):
    # Return unit id primary channel waveforms using dense analyzer
    num_ch = analyzer.get_num_channels()

    template_ext = analyzer.get_extension("templates")
    if template_ext is None:
        analyzer.compute("templates")
        template_ext = analyzer.get_extension("templates")
    template = template_ext.get_unit_template(unit_id)
    wvf_ext = analyzer.get_extension("waveforms")
    wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)
    if wvfs.shape[2] == num_ch:
        # dense case
        template_extremum_indices = si.get_template_extremum_channel(
            analyzer, outputs="index")
        primary_ch_index = template_extremum_indices[unit_id]
        return wvfs[:, :, primary_ch_index]
    else:
        raise ValueError(
            "Analyzer must be dense. The waveforms do not cover all channels.")


def get_raw_unit_primary_ch_wvfs(raw_analyzer, analyzer, unit_id):
    # Return unit id primary channel waveforms using dense analyzer
    num_ch = analyzer.get_num_channels()

    template_ext = raw_analyzer.get_extension("templates")
    template = template_ext.get_unit_template(unit_id)
    wvf_ext = raw_analyzer.get_extension("waveforms")
    wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)
    if wvfs.shape[2] == num_ch:
        # dense case
        template_extremum_indices = si.get_template_extremum_channel(
            analyzer, outputs="index")
        primary_ch_index = template_extremum_indices[unit_id]
        return wvfs[:, :, primary_ch_index]
    else:
        raise ValueError(
            "Analyzer must be dense. The waveforms do not cover all channels.")


def get_unit_primary_ch_template(analyzer, unit_id):
    if analyzer.sparsity is not None:
        print("Error! sparse_analyzer is not dense.")
        return None

    templates_ext = analyzer.get_extension("templates")
    template_extremum_indices = si.get_template_extremum_channel(
        analyzer, outputs="index")
    primary_ch_idx = template_extremum_indices[unit_id]
    primary_ch_wvf = templates_ext.get_unit_template(unit_id)[
        :, primary_ch_idx]

    return primary_ch_wvf


def get_unit_dense_wvfs(analyzer, unit_id, raw_analyzer=None):
    """
    Retrieve dense waveforms, sparse channel indices, and the primary channel index for a given unit.

    Parameters:
    analyzer (object): Analyzer object to fetch waveforms and related data.
    unit_id (int): ID of the unit for which to retrieve waveforms.
    raw_analyzer (object, optional): Analyzer object to fetch raw waveforms. If provided, returns raw dense waveforms.
    radius (float, optional): Radius for sparsity calculation (default: 60).
    peak_sign (str, optional): Sign of the waveform peak for sparsity ("neg" or "pos", default: "neg").

    Returns:
    tuple: (dense_wvfs, sparse_ch_indices, dense_primary_ch_idx[, raw_dense_wvfs])
        - dense_wvfs (ndarray): Dense waveforms for the specified unit.
        - sparse_ch_indices (ndarray): Indices of sparse channels for the unit.
        - dense_primary_ch_idx (int): Index of the primary channel for the unit.
        - raw_dense_wvfs (ndarray, optional): Raw dense waveforms for the unit (if `raw_analyzer` is provided).
    """
    if raw_analyzer:
        raw_wvf_ext = raw_analyzer.get_extension("waveforms")
        raw_dense_wvfs = raw_wvf_ext.get_waveforms_one_unit(unit_id)

    wvf_ext = analyzer.get_extension("waveforms")
    dense_wvfs = wvf_ext.get_waveforms_one_unit(unit_id)

    sparsity = si.ChannelSparsity.from_radius(analyzer, 60, peak_sign="neg")
    sparse_ch_indices = sparsity.unit_id_to_channel_indices[unit_id]

    dense_primary_ch_idx = get_dense_primary_ch_index(analyzer, unit_id)

    # Return raw dense waveforms if available
    if raw_analyzer:
        return dense_wvfs, sparse_ch_indices, dense_primary_ch_idx, raw_dense_wvfs

    return dense_wvfs, sparse_ch_indices, dense_primary_ch_idx


# %%


def plot_template(unit_id, analyzer, plot_scalebar=True, y_negative=-80, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    pre_peak_space = 70
    plot_color = "k"
    x_offset = 0

    # templates_ext = analyzer.get_extension("templates")
    templates_ext = analyzer.compute("templates")

    template_mean = templates_ext.get_unit_template(
        unit_id, operator="average")
    template_std = templates_ext.get_unit_template(unit_id, operator="std")

    ch_loc = analyzer.get_channel_locations()
    ch_idx = analyzer.channel_ids_to_indices(analyzer.channel_ids)

    # first index deep, last index shallow
    depth_order = np.argsort(ch_loc[:, 1])
    abs_depth_id = 31 - np.sort(ch_loc[:, 1]) / 60 + 1

    t1 = 0
    t2 = 181
    ordered_template_mean = template_mean[t1:t2, depth_order]
    ordered_template_std = template_std[t1:t2, depth_order]

    primary_ch = np.argmin(np.min(ordered_template_mean, axis=0))
    channels_to_plot = [primary_ch - 1, primary_ch, primary_ch + 1]

    # Handle edge cases
    if primary_ch == 0:  # If primary channel is at the top
        channels_to_plot = [primary_ch, primary_ch + 1, primary_ch + 2]
    # If primary channel is at the bottom
    elif primary_ch == len(ch_idx) - 1:
        channels_to_plot = [primary_ch - 2, primary_ch - 1, primary_ch]

    peak_indices = [np.argmin(ordered_template_mean[:, ch_idx])
                    for ch_idx in channels_to_plot]

    # Calculate the offsets based on the std values before the peak
    offsets = [i * pre_peak_space for i in range(len(channels_to_plot))]

    for i, ch_idx in enumerate(channels_to_plot):
        # Get the mean and std for the current channel
        depth_id = "D" + str(int(abs_depth_id[ch_idx]))
        mean_waveform = ordered_template_mean[:, ch_idx]
        std_waveform = ordered_template_std[:, ch_idx]
        x = np.arange(mean_waveform.shape[0]) + x_offset

        # Plot the mean waveform
        ax.plot(x, mean_waveform + offsets[i], "k")
        ax.text(5, offsets[i] - 20, depth_id, fontsize=8)

        # Add shading for the std
        ax.fill_between(
            x,  # Use the same x-axis values as for the mean waveform
            mean_waveform - std_waveform + offsets[i],
            mean_waveform + std_waveform + offsets[i],
            color="k",
            alpha=0.2,  # Adjust transparency of the shading
        )

    ax.set_ylim([y_negative, offsets[-1] + 50])
    if plot_scalebar:
        h_pos = (0, y_negative + 10)  # Position for horizontal scale bar
        v_pos = (0, y_negative + 10)  # Position for vertical scale bar
        h_length = 30  # Length of horizontal scale bar (in samples)
        # Length of vertical scale bar (in amplitude units)
        v_length = 100

        plot_util.add_scale_bars_wvf(
            ax,
            h_pos,
            v_pos,
            h_length,
            v_length,
            line_width=1,
            h_label="1 ms",
            v_label="100 uV",
            v_label_x_offset=0,
            v_label_y_offset=50,
            h_label_x_offset=20,
            h_label_y_offset=-10,
        )
    ax.axis("off")
# %%


# wvfs = get_unit_primary_ch_wvfs(analyzer, unit_id=23)

# plt.plot(wvfs.T, 'k', linewidth=0.5, alpha=0.1)

# #%%

# # When loading analyzer from stage2, it is dense but should handl eiehter case
# # First handle dense case.

# # Template dense, waveform dense, ch idx from template ext ch function
# unit_id = 6

# wvf_ext = analyzer.get_extension("waveforms")
# wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)[:, :, primary_ch_index]
# plt.plot(wvfs.T)


# #%%
# unit_id = analyzer.unit_ids[0]
# wvfs = get_unit_primary_ch_wvfs(analyzer, unit_id)

# plt.plot(wvfs)
