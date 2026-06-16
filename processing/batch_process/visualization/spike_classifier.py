import numpy as np
from spikeinterface import full as si
import os
import pickle
from scipy.signal import find_peaks, peak_prominences, peak_widths
import matplotlib.pyplot as plt
import hashlib
import batch_process.util.template_util as template_util

# Module for waveform acceptance criteria


def get_MAD_threshold_primary_ch(analyzer, unit_id):
    template_extremum_ids = si.get_template_extremum_channel(
        analyzer, outputs="id")
    primary_ch = template_extremum_ids[unit_id]
    rec = analyzer.recording
    trace = rec.get_traces(
        start_frame=0, end_frame=30000 * 10, channel_ids=[primary_ch])

    # raw_rec = raw_analyzer.recording
    # raw_trace = raw_rec.get_traces(
    #     start_frame=0, end_frame=30000 * 10, channel_ids=[primary_ch])

    # plt.plot(trace)
    # plt.plot(raw_trace)
    mad = np.median(np.abs(trace - np.median(trace)))
    threshold = 2 * mad / 0.6745

    return threshold


def fit_line_and_get_slope(waveform, start, end):
    x = np.arange(start, end)
    y = waveform[start:end]

    if len(x) < 2:
        return 0

    coeffs = np.polyfit(x, y, deg=1)
    return coeffs[0]  # Return the slope


def detrend_waveform(waveform, samples_to_fit, order=2):
    coefficients = np.polyfit(samples_to_fit, waveform[samples_to_fit], order)

    x = np.arange(len(waveform))
    baseline = np.polyval(coefficients, x)

    detrended_waveform = waveform - baseline  # Subtract baseline
    return detrended_waveform, baseline


def find_waveform_base(waveform, trough_index, search_start=20, search_end=30):
    if search_end is None:
        search_end = len(waveform)

    # Compute the derivative
    derivative = np.diff(waveform)

    derivative_range = derivative[search_start:search_end]
    most_negative_index = np.argmin(derivative_range) + search_start
    # HARD CODED? seems earlier is better than later
    base_index = most_negative_index - 2

    return base_index


def find_width(waveform, trough_index):
    peak_amplitude = waveform[trough_index]
    threshold = 0.5 * peak_amplitude
    left_index = trough_index
    while left_index > 0 and waveform[left_index] < threshold:
        left_index -= 1

    right_index = trough_index
    while right_index < len(waveform) - 1 and waveform[right_index] < threshold:
        right_index += 1

    spike_width = right_index - left_index
    return spike_width, left_index, right_index


def find_up_down_prevalence(detrended_waveform):
    pattern = [1, -1, 1, -1]

    transitions = np.sign(np.diff(detrended_waveform))

    pattern_length = len(pattern)
    matches = []
    for i in range(len(transitions) - pattern_length + 1):
        if (transitions[i:i + pattern_length] == pattern).all():
            matches.append(i)

    prevalence = len(matches) / (len(detrended_waveform) - len(pattern) + 1)
    return prevalence


def get_peak_find_metrics(detrended_waveform):
    sample_range_start = 18  # HARDCODED
    sample_range_end = 50

    detrended_segment = detrended_waveform[sample_range_start:sample_range_end]

    # Find peaks within the segment
    peaks_segment, properties_segment = find_peaks(
        -detrended_segment, prominence=(50, None), width=3
    )

    if len(peaks_segment) > 0:
        # Adjust the peak indices to the full waveform
        peaks = peaks_segment + sample_range_start

        # Calculate prominences within the segment
        prominences = peak_prominences(-detrended_segment, peaks_segment)[0]

        # Calculate widths and heights at rel_height=0.5 within the segment
        prom_width_half, width_height_half, left_ips_half, right_ips_half = peak_widths(
            -detrended_segment, peaks_segment, rel_height=0.5
        )

        # Calculate widths and heights at rel_height=0.8 within the segment
        prom_width_0_8, width_height_0_8, left_ips_0_8, right_ips_0_8 = peak_widths(
            -detrended_segment, peaks_segment, rel_height=0.8
        )

        # Adjust the indices of left_ips and right_ips to match the full waveform
        left_ips_half += sample_range_start
        right_ips_half += sample_range_start
        left_ips_0_8 += sample_range_start
        right_ips_0_8 += sample_range_start

        # Update properties to include both widths and heights
        properties_segment["widths_0.5"] = prom_width_half
        properties_segment["left_ips_0.5"] = left_ips_half
        properties_segment["right_ips_0.5"] = right_ips_half
        properties_segment["width_heights_0.5"] = width_height_half

        properties_segment["widths_0.8"] = prom_width_0_8
        properties_segment["left_ips_0.8"] = left_ips_0_8
        properties_segment["right_ips_0.8"] = right_ips_0_8
        properties_segment["width_heights_0.8"] = width_height_0_8

        # Find the closest peak
        closest_peak_idx = peaks[np.argmin(np.abs(peaks - 30))]
        closest_peak_prominence = prominences[np.argmin(
            np.abs(peaks - 30) - sample_range_start)]
        prom_width_half = prom_width_half[np.argmin(
            np.abs(peaks - 30) - sample_range_start)]
        prom_width_0_8 = prom_width_0_8[np.argmin(
            np.abs(peaks - 30) - sample_range_start)]

        peak_find_metrics = {
            "peaks": peaks,
            "prominences": prominences,
            "properties": properties_segment,
            "closest_peak_idx": closest_peak_idx,
            "closest_peak_prominence": closest_peak_prominence,
            "prom_width_half": prom_width_half,
            "prom_width_0_8": prom_width_0_8,
            "prom_width_diff": prom_width_0_8-prom_width_half
        }

        return peak_find_metrics

    peak_find_metrics = {}
    return peak_find_metrics


def plot_peak_find_metrics(detrended_waveform, peak_find_metrics):

    if len(peak_find_metrics) == 0:
        return

    peaks = peak_find_metrics['peaks']
    prominences = peak_find_metrics['prominences']
    properties = peak_find_metrics['properties']
    closest_peak_idx = peak_find_metrics['closest_peak_idx']
    closest_peak_prominence = peak_find_metrics['closest_peak_prominence']
    prom_width_half = peak_find_metrics['prom_width_half']
    prom_width_0_8 = peak_find_metrics['prom_width_0_8']

    if len(peaks) > 0:
        # Plot the peaks
        plt.plot(peaks, detrended_waveform[peaks], "x", color="C2")

        # Plot prominence lines
        contour_heights = detrended_waveform[peaks] + prominences
        plt.vlines(
            x=peaks,
            ymin=contour_heights,
            ymax=detrended_waveform[peaks],
            color="C2",
        )

        # Plot widths at rel_height=0.5
        plt.hlines(
            y=-properties["width_heights_0.5"],
            xmin=properties["left_ips_0.5"],
            xmax=properties["right_ips_0.5"],
            color="C3",
        )

        # Plot widths at rel_height=0.8
        plt.hlines(
            y=-properties["width_heights_0.8"],
            xmin=properties["left_ips_0.8"],
            xmax=properties["right_ips_0.8"],
            color="C4",
        )

        # Annotate the closest peak
        if closest_peak_idx is not None:
            plt.text(
                0,
                -150,
                f"Prominence: {closest_peak_prominence:.0f}\n"
                f"Width @ 0.5: {prom_width_half:.2f}\n"
                f"Width @ 0.8: {prom_width_0_8:.2f}",
                bbox=dict(facecolor="white", alpha=0.7),
                color="black",
                fontsize=6
            )


def classify_waveform(
    waveform,
    baseline_to_trough_threshold=52,
    trough_index=30,
    pos_peak_max_range=20,
    gap_constant=5,
):
    """
    Classify waveform as spike-like or not based on amplitude and slope features.

    Args:
        waveform (np.array): The input waveform (raw, not detrended).
        trough_index (int): Index of the trough in the waveform.
        pos_peak_max_range (int): Max range (in samples) to search for the peak after the trough.
        gap_constant (int): Gap to separate trough-to-peak and post-peak segments.

    Returns:
        is_spike_like (bool): Whether the waveform is spike-like.
        results (dict): Dictionary of classification results and features.
    """
    # Detrend the waveform
    ignore_post_trough_slope = False
    seg2_start = trough_index + pos_peak_max_range
    seg2_end = seg2_start + 10
    dv = np.diff(waveform[seg2_start:seg2_end])
    max_dv_index = np.argmax(dv)
    if (dv[max_dv_index] > 500):
        # ignore post trough
        ignore_post_trough_slope = True
        seg2_end = seg2_start + max_dv_index

    samples_to_fit = np.concatenate(
        [np.arange(18, 25), np.arange(seg2_start, seg2_end)]
    )
    detrended_waveform, baseline = detrend_waveform(
        waveform, samples_to_fit, order=2
    )

    # Validate peak index range
    peak_index_range = min(trough_index + pos_peak_max_range, len(waveform))
    peak_index = (
        np.argmax(detrended_waveform[trough_index + 1:peak_index_range])
        + trough_index + 1
    )

    # Find base index
    base_index = find_waveform_base(waveform, trough_index)

    # Early pre-peak slope
    early_pre_peak_end = base_index - 1
    early_pre_peak_start = early_pre_peak_end - 10
    early_pre_peak_slope = fit_line_and_get_slope(
        waveform, early_pre_peak_start, early_pre_peak_end
    )

    # Immediate pre-peak slope
    immediate_pre_peak_start = base_index
    immediate_pre_peak_end = trough_index
    immediate_pre_peak_slope = fit_line_and_get_slope(
        waveform, immediate_pre_peak_start, immediate_pre_peak_end
    )

    # Post-trough slope
    if ignore_post_trough_slope:
        post_trough_slope = 999
        post_trough_slope_x2 = 10000
    else:
        # Post trough slope using samples at trough index to trough index + 6 inclusive
        post_trough_slope_x1 = trough_index
        post_trough_slope_x2 = trough_index + 6
        post_trough_slope = fit_line_and_get_slope(
            waveform, trough_index, post_trough_slope_x2 + 1  # Hardcoded range
        )

    # Post-peak slope
    post_peak_index_start = min(peak_index, len(waveform) - 1)
    post_peak_index_end = min(
        post_peak_index_start + gap_constant, len(waveform) - 1)
    post_peak_slope = fit_line_and_get_slope(
        waveform, post_peak_index_start, post_peak_index_end + 1
    )

    # Amplitude calculations
    base_to_trough = detrended_waveform[base_index] - \
        detrended_waveform[trough_index]
    trough_to_peak = detrended_waveform[peak_index] - \
        detrended_waveform[trough_index]
    base_to_peak = detrended_waveform[peak_index] - \
        detrended_waveform[base_index]

    up_down_prevalence = find_up_down_prevalence(detrended_waveform)

    peak_find_metrics = get_peak_find_metrics(detrended_waveform)

    prom_width_half = peak_find_metrics.get('prom_width_half', float('inf'))
    prom_width_0_8 = peak_find_metrics.get('prom_width_0_8', float('inf'))
    closest_peak_prominence = peak_find_metrics.get(
        'closest_peak_prominence', 0)
    num_peaks = len(peak_find_metrics.get('peaks', []))

    # W-ratio: detect W-shape between trough and peak on detrended waveform
    w_ratio_threshold = 0.3
    det_trough_idx = int(np.argmin(detrended_waveform[28:36]) + 28)
    if peak_index > det_trough_idx and trough_to_peak > 0:
        w_seg = detrended_waveform[det_trough_idx:peak_index + 1]
        w_line = np.linspace(w_seg[0], w_seg[-1], len(w_seg))
        w_max_dip = float(np.max(w_line - w_seg))
        w_ratio = w_max_dip / trough_to_peak
    else:
        w_max_dip = 0.0
        w_ratio = 0.0

    # Define thresholds
    amplitude_ratio_max = 2
    min_prominence_threshold = 80
    max_half_width_threshold = 20
    base_to_trough_threshold = 52
    up_down_prevalence_threshold = 0.3

    # Define criteria
    criteria = [
        (immediate_pre_peak_slope < early_pre_peak_slope),
        # (post_trough_slope > np.max([0, post_peak_slope])),
        (post_trough_slope > 0),
        (base_to_trough > base_to_trough_threshold),
        (num_peaks == 1),
        (closest_peak_prominence >= min_prominence_threshold),
        ((prom_width_half > 6) and (prom_width_half < 19)),
        (prom_width_0_8 > 8) and (prom_width_0_8 < 24),
        ((base_to_peak / base_to_trough) < amplitude_ratio_max),
        (up_down_prevalence < up_down_prevalence_threshold),
        (w_ratio < w_ratio_threshold),
    ]

    reasons = [
        "immediate_pre_peak_slope too high",
        "post_trough_slope too low",
        "base_to_trough amplitude too low",
        "More than one peak detected",
        "Closest peak prominence below threshold",
        "Prominence width at half height not in range",
        "Prominence width at 80% not in range",
        "Amplitude ratio too high",
        "Up-down prevalence too high",
        "W-ratio too high (W-shaped waveform)",
    ]

    criteria_info = [
        {"condition": immediate_pre_peak_slope < early_pre_peak_slope,
            "reason": "immediate_pre_peak_slope too high"},
        {"condition": post_trough_slope > 0, "reason": "post_trough_slope too low"},
        {"condition": base_to_trough > base_to_trough_threshold,
            "reason": "base_to_trough amplitude too low"},
        {"condition": num_peaks == 1, "reason": "More than one peak detected"},
        {"condition": closest_peak_prominence >= min_prominence_threshold,
            "reason": "Closest peak prominence below threshold"},
        {"condition": 6 < prom_width_half < 19,
            "reason": "Prominence width at half height not in range"},
        {"condition": 8 < prom_width_0_8 < 24,
            "reason": "Prominence width at 80% not in range"},
        {"condition": (base_to_peak / base_to_trough) <
         amplitude_ratio_max, "reason": "Amplitude ratio too high"},
        {"condition": up_down_prevalence < up_down_prevalence_threshold,
            "reason": "Up-down prevalence too high"},
        {"condition": w_ratio < w_ratio_threshold,
            "reason": "W-ratio too high (W-shaped waveform)"},
    ]

    # Classification logic
    # is_spike_like = all(criteria)

    # # Track failed criteria
    # failed_criteria = [
    #     reasons[i] for i, passed in enumerate(criteria) if not passed
    # ]

    failed_criteria = [info["reason"]
                       for info in criteria_info if not info["condition"]]
    is_spike_like = not failed_criteria

    # Compile results
    results = {
        "detrended_waveform": detrended_waveform,
        "baseline": baseline,
        "base_index": base_index,
        "peak_index": peak_index,
        "base_to_trough": base_to_trough,
        "trough_to_peak": trough_to_peak,


        "trough_x1": trough_index,
        "trough_x2": post_trough_slope_x2,
        "post_trough_slope": post_trough_slope,

        "post_trough_x1": post_peak_index_start,
        "post_trough_x2": post_peak_index_end,
        "post_peak_slope": post_peak_slope,

        "immediate_pre_peak_slope": immediate_pre_peak_slope,
        "early_pre_peak_slope": early_pre_peak_slope,
        "up_down_prevalence": up_down_prevalence,
        "w_ratio": w_ratio,
        "w_max_dip": w_max_dip,
        "criteria": failed_criteria if failed_criteria else [],
        "peak_find_metrics": peak_find_metrics
    }

    return is_spike_like, results


def plot_slopes(waveform, results, ax):

    x1, x2 = results['trough_x1'], results['trough_x2']
    x3, x4 = results['post_trough_x1'], results['post_trough_x2']
    y1, y2 = waveform[x1], waveform[x2]
    y3, y4 = waveform[x3], waveform[x4]

    # add 1 since arange excludes last index
    slope_1 = fit_line_and_get_slope(waveform, x1, x2 + 1)
    slope_2 = fit_line_and_get_slope(waveform, x3, x4 + 1)

    if ax is None:
        fig, ax = plt.subplot()
        ax.plot(waveform)

    ax.scatter([x1, x2], [y1, y2], color='C0', label="Trough Points")
    ax.scatter([x3, x4], [y3, y4], color='C1', label="Post-Trough Points")
    ax.plot([x1, x2], [y1, y2], color='C0',
            linestyle='--', label="Trough Slope")
    ax.plot([x3, x4], [y3, y4], color='C1',
            linestyle='--', label="Post-Trough Slope")

    # Annotate slopes with arrows and text
    ax.annotate(f"Slope: {slope_1:.2f}",
                xy=((x1 + x2) / 2, (y1 + y2) / 2),
                xytext=(-40, 20),
                textcoords='offset points',
                arrowprops=dict(facecolor='black', arrowstyle='->'),
                fontsize=10, color='C0')

    ax.annotate(f"Slope: {slope_2:.2f}",
                xy=((x3 + x4) / 2, (y3 + y4) / 2),
                xytext=(-40, 20),
                textcoords='offset points',
                arrowprops=dict(facecolor='black', arrowstyle='->'),
                fontsize=10, color='C1')


def get_results_string(r):
    def safe_format(value, fmt="{:.2f}"):
        """Helper function to format numbers while handling non-numeric values gracefully."""
        try:
            return fmt.format(value)
        except (ValueError, TypeError):
            return str(value)  # Return as-is for non-numeric values

    post_trough = safe_format(r.get('post_trough_slope', 'N/A'))
    post_peak = safe_format(r.get('post_peak_slope', 'N/A'))
    imm_pre_peak = safe_format(r.get('immediate_pre_peak_slope', 'N/A'))
    early_pre_peak = safe_format(r.get('early_pre_peak_slope', 'N/A'))
    b2t = safe_format(r.get('base_to_trough', 'N/A'), fmt='{:.0f}')
    base_idx = r.get('base_index', 'N/A')
    ud_prev = safe_format(r.get('up_down_prevalence', 0), fmt='{:.2%}')
    w_diff = safe_format(r.get('peak_find_metrics', {}).get('prom_width_diff', 'N/A'))
    criteria = r.get('criteria', 'N/A')
    results_string = (
        f"Post-trough slope: {post_trough}, "
        f"Post-peak slope: {post_peak},\n"
        f"Immediate Pre-Peak Slope: {imm_pre_peak}, "
        f"Early Pre-Peak Slope: {early_pre_peak},\n"
        f"Base to trough: {b2t}uV, "
        f"Base index: {base_idx}, "
        f"Up-down prevalence: {ud_prev}, "
        f"Width diff: {w_diff},\n"
        f"{criteria}"
    )
    return results_string


def rejected_by(criteria, target_reason):
    reasons = [
        "immediate_pre_peak_slope too high",
        "post_trough_slope too low",
        "base_to_trough amplitude too low",
        "More than one peak detected",
        "Closest peak prominence below threshold",
        "Prominence width at half height not in range",
        "Prominence width at 80% not in range",
        "Amplitude ratio too high",
        "Up-down prevalence too high",
    ]

    reasons = [
        "immediate_pre_peak_slope too high",
        "post_trough_slope too low",
        "base_to_trough amplitude too low",
        "More than one peak detected",
        "Closest peak prominence below threshold",
        "Prominence width at half height not in range",
        "Prominence width at 80% not in range",
        "Amplitude ratio too high",
        "Up-down prevalence too high",
    ]
    for condition, reason in zip(criteria, reasons):
        if not condition and reason == target_reason:
            return True
    return False


def compute_waveform_hash(waveform):
    return hashlib.sha256(waveform.tobytes()).hexdigest()


def load_labeled_waveform_dataset():
    file_name = "batch_process/visualization/labeled_waveform_data.pkl"
    if not os.path.exists(file_name):
        print(f"File {file_name} not found. Returning an empty dictionary.")
        return {}

    with open(file_name, "rb") as file:
        data_loaded = pickle.load(file)
    return data_loaded  # Return the entire dataset as a dictionary


def add_waveforms(labeled_waveforms_dict):
    """
    Add labeled waveforms to the dataset with hash-based deduplication.

    Args:
        labeled_waveforms_dict (dict): A dictionary where each key is a unit_id,
                                       and the value is another dictionary with keys
                                       'good_waveforms' and 'bad_waveforms'.
    """
    # Validate input format
    for unit_id, waveforms in labeled_waveforms_dict.items():
        if 'good_waveforms' not in waveforms:
            print(f"Missing 'good_waveforms' key for unit_id {unit_id}.")
            return
        if 'bad_waveforms' not in waveforms:
            print(f"Missing 'bad_waveforms' key for unit_id {unit_id}.")
            return

    # Load existing dataset
    unit_id_waveforms_dict = load_labeled_waveform_dataset()

    # Update or add waveforms for each unit_id
    for unit_id, waveforms in labeled_waveforms_dict.items():
        new_good_waveforms = waveforms['good_waveforms']
        new_bad_waveforms = waveforms['bad_waveforms']

        # Initialize unit if not present
        if unit_id not in unit_id_waveforms_dict:
            unit_id_waveforms_dict[unit_id] = {
                "good_waveforms": np.array([]).reshape(0, new_good_waveforms.shape[1]),
                "bad_waveforms": np.array([]).reshape(0, new_bad_waveforms.shape[1]),
                "good_waveform_hashes": set(),
                "bad_waveform_hashes": set(),
            }

        # Append new waveforms with hash checks
        existing_good_waveforms = unit_id_waveforms_dict[unit_id]["good_waveforms"]
        existing_bad_waveforms = unit_id_waveforms_dict[unit_id]["bad_waveforms"]
        good_hashes = unit_id_waveforms_dict[unit_id]["good_waveform_hashes"]
        bad_hashes = unit_id_waveforms_dict[unit_id]["bad_waveform_hashes"]

        # Deduplicate good waveforms
        unique_good_waveforms = []
        for waveform in new_good_waveforms:
            hash_value = compute_waveform_hash(waveform)
            if hash_value not in good_hashes:
                unique_good_waveforms.append(waveform)
                good_hashes.add(hash_value)

        # Deduplicate bad waveforms
        unique_bad_waveforms = []
        for waveform in new_bad_waveforms:
            hash_value = compute_waveform_hash(waveform)
            if hash_value not in bad_hashes:
                unique_bad_waveforms.append(waveform)
                bad_hashes.add(hash_value)

        # Update the dataset
        if unique_good_waveforms:
            existing_good_waveforms = np.concatenate(
                [existing_good_waveforms, unique_good_waveforms])
        if unique_bad_waveforms:
            existing_bad_waveforms = np.concatenate(
                [existing_bad_waveforms, unique_bad_waveforms])

        unit_id_waveforms_dict[unit_id] = {
            "good_waveforms": existing_good_waveforms,
            "bad_waveforms": existing_bad_waveforms,
            "good_waveform_hashes": good_hashes,
            "bad_waveform_hashes": bad_hashes,
        }

    # Save updated dataset
    file_name = "batch_process/visualization/labeled_waveform_data.pkl"
    with open(file_name, "wb") as file:
        pickle.dump(unit_id_waveforms_dict, file)
    print(f"Updated data saved to {file_name} for unit_ids: {list(labeled_waveforms_dict.keys())}")


def center_waveforms(raw_waveforms):
    return raw_waveforms - raw_waveforms[:, 30].reshape(-1, 1)


def extract_waveform_features(raw_analyzer, analyzer, unit_id, good_indices):
    raw_primary_ch_wvfs = template_util.get_raw_unit_primary_ch_wvfs(
        raw_analyzer, analyzer, unit_id)
    centered_raw_wvfs = center_waveforms(raw_primary_ch_wvfs)
    return centered_raw_wvfs[good_indices], raw_primary_ch_wvfs[good_indices]


def plot_umap(umap_x, umap_y, labels, title, ax, n=5):
    unique_labels = np.unique(labels)
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    for i, label in enumerate(unique_labels):
        color = colors[i % len(colors)]
        mask = (labels == label)
        # downsample
        x_sampled = umap_x[mask][::n]
        y_sampled = umap_y[mask][::n]
        ax.scatter(x_sampled, y_sampled, color=color,
                   s=3, label=f'Cluster {label}')

    ax.set_title(title, fontsize=10)
    ax.set_xlabel('UMAP Dim 1', fontsize=8)  # Adjust axis label font size
    ax.set_ylabel('UMAP Dim 2', fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(fontsize=6)


def plot_mean_waveforms(waveforms, labels, ax, n=5):
    """
    Plot the mean waveform for each cluster using specific colors ('C0', 'C1', 'C2').
    """
    unique_labels = np.unique(labels)
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    for i, label in enumerate(unique_labels):
        color = colors[i % len(colors)]
        cluster_wvfs = waveforms[labels == label][::n]
        mean_wvf = np.mean(cluster_wvfs, axis=0)
        std_wvf = np.std(cluster_wvfs, axis=0)

        # Plot mean waveform with shaded error region
        ax.plot(mean_wvf, color=color, label=f'Cluster {label}')
        ax.fill_between(np.arange(len(mean_wvf)),
                        mean_wvf - std_wvf,
                        mean_wvf + std_wvf,
                        color=color, alpha=0.3)

    ax.set_title("Mean Waveforms", fontsize=10)
    ax.set_xlabel("Samples", fontsize=8)
    ax.set_ylabel("Amplitude (µV)", fontsize=8)
    ax.legend(fontsize=6)
    ax.set_ylim([-300, 300])


def plot_mean_subcluster_waveforms(analyzer, unit_id, max_waveforms=20, ax=None, raw_analyzer=None):
    # unit ids correspond to child ids after stage3
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    if raw_analyzer is not None:
        raw_wvfs = template_util.get_raw_unit_primary_ch_wvfs(
            raw_analyzer, analyzer, unit_id)
        num_waveforms_to_plot = min(len(raw_wvfs), max_waveforms)
        wvfs_to_plot = raw_wvfs - raw_wvfs[:, 30][:, np.newaxis]
        wvfs_to_plot = wvfs_to_plot[:num_waveforms_to_plot]
    else:
        wvfs = template_util.get_unit_primary_ch_wvfs(analyzer, unit_id)
        num_waveforms_to_plot = min(len(wvfs), max_waveforms)
        wvfs_to_plot = wvfs[:num_waveforms_to_plot]

    t = np.arange(wvfs_to_plot.shape[1])/30
    mean_wvf = np.mean(wvfs_to_plot, axis=0)
    ax.plot(t, wvfs_to_plot.T, 'k', linewidth=0.5)
    ax.plot(t, mean_wvf, color='r', linewidth=2)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude (µV)")

    return ax
