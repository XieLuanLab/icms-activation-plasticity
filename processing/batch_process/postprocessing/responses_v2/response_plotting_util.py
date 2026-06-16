from matplotlib.patches import Patch
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import pandas as pd
from pathlib import Path
import os

import dill as pickle
import shutil
import re
from collections import defaultdict
from matplotlib.lines import Line2D

# Lazy imports for heavy/optional dependencies
# batch_process.util.plotting (plot_util) — only needed by plot_template
# batch_process.util.classify_cell_type — only needed if cell_type_df accessed
# batch_process.util.file_util — only needed by make_ppt
# batch_process.util.template_util — only needed by plot_template
# batch_process.util.ppt_image_inserter — only needed by make_ppt


def extract_unit_id(image_path):
    """
    Extracts the unit ID from the filename.
    Assumes filenames are in the format 'Unit_<unit_id>_*.png'.
    """
    # Use a regular expression to extract the unit number from the filename
    match = re.search(r"Unit_(\d+)_", image_path.name)
    if match:
        return int(match.group(1))
    return None


def make_ppt(session_path, figures_dir=None):
    import batch_process.util.file_util as file_util
    from batch_process.util.ppt_image_inserter import PPTImageInserter

    date_str = file_util.get_date_str(session_path)
    if figures_dir is None:
        figures_dir = Path(session_path) / "batch_sort" / \
            "figures" / "stim_condition_figures"

    # Collect all image paths from the directory
    # Adjust the pattern if needed
    image_paths = list(figures_dir.glob("*.png"))

    # Group images by unit ID
    images_by_unit = defaultdict(list)
    for image_path in image_paths:
        unit_id = extract_unit_id(image_path)
        if unit_id is not None:
            images_by_unit[unit_id].append(image_path)

    ppt_inserter = PPTImageInserter(grid_dims=(
        2, 2), spacing=(0.05, 0.05), title_font_size=16)

    # Add images to the PPT, grouping them by unit ID
    for unit_id, image_paths in sorted(images_by_unit.items()):
        slide_title = f"{date_str}: unit {unit_id}"
        ppt_inserter.add_slide(slide_title)
        for image_path in image_paths:
            # Pass the image path as a string (required by pptx)
            ppt_inserter.add_image(str(image_path), slide_title=slide_title)

    # Save or return the ppt_inserter object as needed
    ppt_inserter.save(Path(session_path) / "batch_sort" /
                      "figures" / "stim_responses.pptx")

    ppt_inserter.save(Path(session_path).parent / "analysis" /
                      f"stim_responses_{Path(session_path).stem}.pptx")


def save_figure(fig, figures_dir, unit_id, stim_channel):

    # Define the file name and path
    fig_filename = f"Unit_{unit_id}_Stim_Channel_{stim_channel}.png"
    fig_path = figures_dir / fig_filename

    # Save the figure
    fig.savefig(fig_path)
    plt.close(fig)  # Close the figure after saving to free up memory

    print(f"Figure saved to {fig_path}")


def plot_template(session_responses, unit_id, axs, plot_scalebar=True, y_negative=-80):
    import batch_process.util.plotting as plot_util

    if isinstance(axs, matplotlib.axes._axes.Axes):
        ax = axs
    else:
        ax = axs['template']

    pre_peak_space = 70
    plot_color = "k"
    x_offset = 0

    analyzer = session_responses.sorting_analyzer
    templates_ext = analyzer.get_extension("templates")

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
        ax.text(-15, offsets[i] - 20, depth_id, fontsize=8)

        # Add shading for the std
        ax.fill_between(
            x,
            mean_waveform - std_waveform + offsets[i],
            mean_waveform + std_waveform + offsets[i],
            color="k",
            alpha=0.2,
        )

    ax.set_ylim([y_negative, offsets[-1] + 50])
    if plot_scalebar:
        h_pos = (-20, y_negative + 10)  # Position for horizontal scale bar
        v_pos = (-20, y_negative + 10)  # Position for vertical scale bar
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


def plot_primary_channel_template(stim_response, plot_scalebar=True, y_negative=-80, ax=None):
    import batch_process.util.plotting as plot_util

    session_responses = stim_response.session
    unit_id = stim_response.unit_id

    plot_color = "k"
    x_offset = 0

    analyzer = session_responses.sorting_analyzer
    templates_ext = analyzer.get_extension("templates")

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

    # Find primary channel as the channel with the minimum peak amplitude
    primary_ch = np.argmin(np.min(ordered_template_mean, axis=0))

    # Get the peak index and depth ID for the primary channel
    depth_id = "D" + str(int(abs_depth_id[primary_ch]))
    mean_waveform = ordered_template_mean[:, primary_ch]
    std_waveform = ordered_template_std[:, primary_ch]
    x = np.arange(mean_waveform.shape[0]) + x_offset

    # Plot the mean waveform for the primary channel
    ax.plot(x, mean_waveform, plot_color)

    # Add shading for the std
    ax.fill_between(
        x,  # Use the same x-axis values as for the mean waveform
        mean_waveform - std_waveform,
        mean_waveform + std_waveform,
        color=plot_color,
        alpha=0.2,  # Adjust transparency of the shading
    )

    # ax.set_ylim([y_negative, y_negative + 100])
    if plot_scalebar:
        h_pos = (0, y_negative + 10)  # Position for horizontal scale bar
        v_pos = (0, y_negative + 10)  # Position for vertical scale bar
        h_length = 30  # Length of horizontal scale bar (in samples)
        v_length = 100  # Length of vertical scale bar (in amplitude units)
        plot_util.add_scale_bars_wvf(
            ax,
            h_pos,
            v_pos,
            h_length,
            v_length,
            line_width=1,
            h_label="1 ms",
            v_label="100 uV",
        )
    ax.axis("off")


def save_results_to_csv(results, session_path, stim_trial_type):
    # Define the CSV file path
    results_dir = Path(session_path)
    csv_path = results_dir / f"stim_condition_results_{stim_trial_type}.csv"
    pkl_path = results_dir / f"stim_condition_results_{stim_trial_type}.pkl"

    # Convert results list to DataFrame
    df = pd.DataFrame(results)

    # Auto-detect all possible columns
    all_columns = set(col for result in results for col in result.keys())

    # ✅ Ensure missing columns are filled with NaN
    for col in all_columns:
        if col not in df.columns:
            df[col] = np.nan

    # Detect numeric columns dynamically and convert
    numeric_columns = df.select_dtypes(include=['number']).columns.tolist()
    df[numeric_columns] = df[numeric_columns].apply(
        pd.to_numeric, errors='coerce')

    # Save DataFrame to CSV and Pickle
    df.to_csv(csv_path, index=False)
    df.to_pickle(pkl_path)

    print(f"Results saved to {csv_path}")


def get_stim_colormap(max_current=20):
    """
    Generates a colormap for currents, assigning specific colors for currents between 3 to 10 µA,
    and assigning the last color in the palette to 2 µA. Black is assigned for currents outside this range (0, 1, >10 µA).

    Args:
        max_current (int): Maximum current value for assigning colors. Defaults to 20 µA.

    Returns:
        dict: A dictionary mapping currents to specific colors.
    """
    # Define the range of currents that should get a color (4 to 10 µA inclusive)
    valid_currents = list(range(4, 11))

    # Generate a color palette for the valid currents (from 4 to 10 µA)
    palette = sns.color_palette("deep", len(valid_currents))

    # Create a fixed mapping of valid currents (3 to 10 µA) to the color palette
    colormap = {curr: palette[i] for i, curr in enumerate(valid_currents)}

    # Assign black for currents < 4 and above 10 µA
    for invalid_curr in range(0, max_current + 1):
        if invalid_curr < 4 or invalid_curr > 10:
            # Assign black color to invalid currents
            colormap[invalid_curr] = (0, 0, 0)

    return colormap


def initialize_plotting():
    fig = plt.figure(figsize=(8, 6))  # Adjust size as needed
    gs = fig.add_gridspec(4, 5)  # 4 rows, 5 columns

    axs = {
        "pulse_raster": fig.add_subplot(gs[0:2, 0:2]),
        "train_raster": fig.add_subplot(gs[0:2, 2:4]),
        "spike_prob": fig.add_subplot(gs[2:4, 0:2]),
        "train_fr": fig.add_subplot(gs[2:4, 2:4]),
        "template": fig.add_subplot(gs[0:3, 4]),
        "z_score": fig.add_subplot(gs[3, 4])
    }

    return fig, axs


def plot_train_rasters(responses, axs):

    color_map = get_stim_colormap()

    all_rasters = []
    all_colors = []
    all_offsets = []

    linewidth_factor = 40
    offset = 0

    ax = axs['train_raster']
    # Sort by current for nice plotting
    for current, response in sorted(responses):
        raster_array = response.train_response.raster_array
        line_count = len(raster_array)
        line_height = np.ceil(line_count / linewidth_factor)
        offsets = np.arange(line_count) + offset

        all_rasters.extend(raster_array)
        all_colors.extend([color_map[current]] * line_count)
        all_offsets.extend(offsets)

        offset += line_count

    ax.eventplot(
        all_rasters,
        orientation="horizontal",
        colors=all_colors,
        linelengths=line_height,
        linewidths=1,
        lineoffsets=all_offsets,
        label=f'{current} µA'
    )

    handles = [
        Line2D([0], [0], color=color_map[current], lw=3, label=f"{current} µA")
        for current in sorted(set([r[0] for r in responses]))
    ]

    # Derive x-limits from timing params
    tparams = responses[0][1].train_response.timing_params
    stim_dur = tparams['stim_duration_ms']
    x_min = -stim_dur
    x_max = tparams['window_size_ms']

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax.axvline(x=stim_dur, color='gray', linestyle='--', linewidth=1)
    ax.set_xlim([x_min, x_max])
    ax.set_ylabel("Trial")
    ax.legend(handles=handles)


def plot_train_firing_rate(responses, axs):
    color_map = get_stim_colormap()
    # ax = axs['train_fr']
    if isinstance(axs, matplotlib.axes._axes.Axes):
        ax = axs
    else:
        ax = axs['train_fr']

    if isinstance(responses, list):
        tparams = responses[0][1].pulse_response.timing_params
        # list of (current, response) tuples
        sorted_responses = sorted(responses)
    else:
        tparams = responses.train_response.timing_params
        # single response with dummy current
        sorted_responses = [(responses.stim_current, responses)]

    for current, response in sorted_responses:
        tr = response.train_response
        tr_data = tr.ten_ms_smoothed_data
        x = tr_data['bin_centers']
        y = tr_data['firing_rate']
        ax.plot(
            x,
            y,
            color=color_map[current],
            alpha=0.1 if len(tr.raster_array) < 5 else 1,
            label=f"{current} µA",
        )

    # Derive x-limits from timing params
    stim_dur = tparams['stim_duration_ms']
    x_min = -stim_dur
    x_max = tparams['window_size_ms']

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax.axvline(x=stim_dur, color='gray', linestyle='--', linewidth=1)
    ax.set_xlim([x_min, x_max])
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Firing rate (Hz)")


def plot_pulse_raster(responses, axs):
    if isinstance(axs, matplotlib.axes._axes.Axes):
        ax = axs
    else:
        ax = axs['pulse_raster']

    color_map = get_stim_colormap()

    all_rasters = []
    all_colors = []
    all_offsets = []
    legend_handles = []
    used_currents = set()

    offset = 0  # Start offset

    for current, response in sorted(responses):  # Each current condition
        pr = response.pulse_response
        tparams = pr.timing_params
        pulse_raster_array = pr.raster_array
        num_trials = len(pulse_raster_array)

        if num_trials == 0:
            continue  # Skip if empty

        lineoffsets = np.arange(num_trials) + offset

        all_rasters.extend(pulse_raster_array)
        all_colors.extend([color_map[current]] * num_trials)
        all_offsets.extend(lineoffsets)

        if current not in used_currents:
            legend_handles.append(
                Patch(color=color_map[current], label=f"{current} µA"))
            used_currents.add(current)

        offset += num_trials  # Update offset for next current

    # Scale linelength to total number of pulses so raster isn't too tall/short
    total_pulses = len(all_rasters)
    pulse_linelength = max(1, np.ceil(total_pulses / 40))

    # Plot all currents together
    ax.eventplot(
        all_rasters,
        orientation="horizontal",
        linewidths=0.25,
        linelengths=pulse_linelength,
        lineoffsets=all_offsets,
        colors=all_colors,
    )

    # Use the last tparams seen for shading — assumes same timing for all currents
    ax.axvspan(0, tparams["post_stim_blank_ms"], color="lightgray",
               alpha=0.5, linewidth=0, edgecolor='none', zorder=0)
    ax.axvspan(
        tparams["pulse_win_ms"] - tparams["pre_stim_blank_ms"],
        tparams["pulse_win_ms"],
        color="lightgray",
        alpha=0.5, linewidth=0, edgecolor='none', zorder=0
    )
    ax.set_xlim([0, tparams["pulse_win_ms"]])
    ax.set_ylabel("Pulse")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.legend(handles=legend_handles, loc='best')


def plot_probability(responses, axs):

    if isinstance(axs, dict):
        # fallback if key is missing
        ax = axs.get('spike_prob', axs[list(axs.keys())[0]])
    else:
        ax = axs

    color_map = get_stim_colormap()

    if isinstance(responses, list):
        tparams = responses[0][1].pulse_response.timing_params
        # list of (current, response) tuples
        sorted_responses = sorted(responses)
    else:
        tparams = responses.pulse_response.timing_params
        # single response with dummy current
        sorted_responses = [(responses.stim_current, responses)]

    legend_handles = []
    used_currents = set()

    for current, response in sorted_responses:
        pr = response.pulse_response
        tr = response.train_response
        baseline = pr.baseline_metrics
        stim = pr.stim_metrics
        color = color_map[current]

        # Plot stim response
        if current not in used_currents:
            legend_handles.append(
                Patch(color=color_map[current], label=f"{current} µA"))
            used_currents.add(current)

        ax.plot(stim["bin_centers"], stim["mean_prob_spike"],
                color=color, alpha=0.2 if len(tr.raster_array) < 5 else 1, zorder=2)
        ax.fill_between(
            stim["bin_centers"],
            stim["mean_prob_spike"] - stim["std_prob_spike"],
            stim["mean_prob_spike"] + stim["std_prob_spike"],
            color=color, alpha=0.1 if len(tr.raster_array) < 5 else 0.3,
            linewidth=0, edgecolor='none', zorder=1
        )

    # Shade the blanking periods
    ax.axvspan(0, tparams["post_stim_blank_ms"], color="lightgray",
               alpha=0.5, linewidth=0, edgecolor='none', zorder=0)
    ax.axvspan(
        tparams["pulse_win_ms"] - tparams["pre_stim_blank_ms"],
        tparams["pulse_win_ms"],
        color="lightgray", alpha=0.5, linewidth=0, edgecolor='none', zorder=0
    )

    ax.set_xlim([0, tparams["pulse_win_ms"]])
    ax.set_ylabel("P(spike)")
    ax.set_xlabel("Time (ms)")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if isinstance(responses, list):
        handles = [
            Line2D([0], [0], color=color_map[current],
                   lw=3, label=f"{current} µA")
            for current in sorted(set([r[0] for r in responses]))
        ]
        ax.legend(handles=handles)
    ax.legend(handles=legend_handles, loc='best', edgecolor='none')


def plot_z_scores(responses, axs, pre_stim_threshold=0.5):
    ax = axs['z_score']
    color_map = get_stim_colormap()
    currents = [response[0] for response in responses]
    for current, response in sorted(responses):
        tr = response.train_response
        pr = response.pulse_response
        z_score = tr.z_score

        alpha = 0.2 if tr.pre_stim_mean_fr < pre_stim_threshold or len(
            tr.raster_array) < 5 else 1

        is_modulated = tr.paired_p_val < 0.05 and z_score > 0
        hatch = '///' if (is_modulated and pr and pr.is_pulse_locked) else None
        edgecolor = 'black' if is_modulated else None

        ax.bar(
            current, z_score, color=color_map[current], alpha=alpha,
            hatch=hatch, edgecolor=edgecolor, linewidth=1.5 if is_modulated else 0.5)

        if is_modulated and pr:
            label = "PL" if pr.is_pulse_locked else "NPL"
            ax.text(current, z_score, label, ha='center', va='bottom', fontsize=6, fontweight='bold')

    ax.set_xlabel("Current (µA)")
    ax.set_xticks(currents)
    ax.set_xticklabels([str(curr) for curr in currents])
    ax.set_ylabel("Z-Score")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
