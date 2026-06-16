"""
Automatic unit curation for ICMS recordings.

Adapted from batch_process/util/curate_util.py for use with the new
standalone pipeline. Removes stim artifacts and noise units using:
  1) Global artifact check - units with large signals on ALL channels
  2) ACG stim-locked check - units with autocorrelogram peaks at stim rate
  3) Waveform shape classification - rejects non-spike waveforms
"""

import numpy as np
import spikeinterface as si
from scipy.spatial.distance import cdist


def _get_sparse_channel_indices(analyzer, radius_um=60):
    """
    For each unit, get channel indices within `radius_um` of the peak channel.
    Pure numpy — no si.compute_sparsity (avoids Spyder isinstance issues).
    """
    ch_locs = analyzer.get_channel_locations()
    extremum = si.get_template_extremum_channel(analyzer, outputs="index")
    result = {}
    for uid in analyzer.unit_ids:
        dists = cdist(ch_locs[extremum[uid]].reshape(1, -1), ch_locs).flatten()
        result[uid] = np.where(dists <= radius_um)[0]
    return result


# =============================================================================
# Core curation checks
# =============================================================================

def _sparse_template_peak_timing(template):
    """Get the sample index of the negative peak on each channel of a sparse template."""
    return np.argmin(template, axis=0)


def _get_timing_offsets(analyzer):
    """
    For each unit, compute the mean absolute timing offset of the negative peak
    across sparse channels relative to the primary (first) channel.

    Large offsets indicate artifact (simultaneous peak across all channels)
    rather than a propagating spike.
    """
    sparse_map = _get_sparse_channel_indices(analyzer, radius_um=60)
    unit_ids = analyzer.unit_ids
    templates_ext = analyzer.get_extension("templates")

    ts_offsets = {}
    for unit_id in unit_ids:
        template = templates_ext.get_unit_template(unit_id)
        sparse_ch_inds = sparse_map[unit_id]
        sparse_template = template[:, sparse_ch_inds]
        ts = _sparse_template_peak_timing(sparse_template)
        primary_ts = ts[0]
        offsets = ts - primary_ts
        ts_offsets[unit_id] = np.mean(np.abs(offsets[1:])) if len(offsets) > 1 else 0
    return ts_offsets


def _exclude_artifact_unit(unit_id, analyzer, threshold=0.3):
    """
    Check if a unit is a stim artifact by measuring whether channels far from
    the primary channel still have large amplitudes.

    For real spikes, only nearby channels should have signal.
    For artifacts, ALL channels have ~equal amplitude.

    Returns True if the unit should be EXCLUDED (is an artifact).
    """
    num_channels = analyzer.get_num_channels()
    sparse_map = _get_sparse_channel_indices(analyzer, radius_um=120)
    sparse_ch_indices = sparse_map[unit_id]

    templates_ext = analyzer.get_extension("templates")
    dense_template = templates_ext.get_unit_template(unit_id)
    assert dense_template.shape[1] == num_channels

    # Find the peak sample (negative peak on primary channel)
    min_idx = np.argmin(np.min(dense_template, axis=1))
    range_offset = 10

    # Get max absolute amplitude around the peak on all channels
    peak_region = dense_template[max(0, min_idx - range_offset):min_idx + range_offset + 1, :]
    max_peaks = np.max(np.abs(peak_region), axis=0)

    # Channels NOT in the sparse radius
    all_ch_indices = np.arange(num_channels)
    non_sparse_ch_indices = np.setdiff1d(all_ch_indices, sparse_ch_indices)

    if len(non_sparse_ch_indices) == 0:
        return False  # All channels are "nearby" — can't tell

    max_peaks_non_sparse = max_peaks[non_sparse_ch_indices]
    template_peak = np.abs(np.min(dense_template[min_idx, :]))

    if template_peak == 0:
        return True

    # Normalized ratio: how large are distant channels vs the peak?
    normalized = max_peaks_non_sparse / template_peak
    top5 = np.sort(normalized)[-5:]
    mean_norm_ratio = np.mean(top5)

    return mean_norm_ratio > threshold


def _classify_waveform(waveform):
    """
    Simple classifier: checks if the waveform shape at specific indices
    looks like a real neural spike vs noise/artifact.

    Based on empirically determined boundaries from the old pipeline.
    """
    point1 = waveform[6]
    point2 = waveform[74]
    return -60 <= point1 <= 100 and -21 <= point2 <= 35


def _get_primary_ch_templates(analyzer):
    """Get each unit's template on its primary (extremum) channel."""
    templates_ext = analyzer.get_extension("templates")
    extremum_indices = si.get_template_extremum_channel(analyzer, outputs="index")

    result = {}
    for unit_id in analyzer.unit_ids:
        primary_idx = extremum_indices[unit_id]
        result[unit_id] = templates_ext.get_unit_template(unit_id)[:, primary_idx]
    return result


def _sort_units_by_depth(analyzer, unit_ids):
    """Sort unit IDs by channel depth (deepest first)."""
    if len(unit_ids) == 0:
        return np.array([], dtype=int)

    ch_loc = analyzer.get_channel_locations()[:, 1]
    ch_ids = analyzer.channel_ids
    ch_idx = analyzer.channel_ids_to_indices(ch_ids)
    ch_idx_to_loc = {idx: loc for idx, loc in zip(ch_idx, ch_loc)}

    extremum_idx = si.get_template_extremum_channel(analyzer, outputs="index")

    unit_locs = {uid: ch_idx_to_loc.get(extremum_idx[uid], 0) for uid in unit_ids}
    sorted_ids = sorted(unit_locs.keys(), key=lambda u: unit_locs[u], reverse=True)
    return np.array(sorted_ids)


# =============================================================================
# Main curation function
# =============================================================================

def curate_units(analyzer, min_spikes=100):
    """
    Curate sorted units to remove stim artifacts and noise.

    Checks:
      1. Global artifact: large amplitude on ALL channels → artifact
      2. ACG stim-locked: autocorrelogram peak at ~33ms (30 Hz stim) + large
         timing offset across channels → artifact
      3. Waveform shape: primary channel template doesn't look like a spike → noise
      4. Min spike count: fewer than `min_spikes` → reject

    Parameters
    ----------
    analyzer : SortingAnalyzer
        Must have 'templates' and 'correlograms' extensions computed.
        Must be DENSE (sparsity=None).
    min_spikes : int
        Minimum number of spikes to keep a unit.

    Returns
    -------
    good_ids : np.ndarray
        Unit IDs that pass all checks, sorted by depth.
    bad_ids : np.ndarray
        Unit IDs that fail one or more checks, sorted by depth.
    labels : dict
        Reason each unit was rejected: 'artifact', 'noise', 'low_count', or 'good'.
    """
    unit_ids = analyzer.unit_ids
    labels = {}

    # Ensure required extensions exist
    if analyzer.get_extension("templates") is None:
        print("  Computing templates...")
        analyzer.compute("templates")
    if analyzer.get_extension("correlograms") is None:
        print("  Computing correlograms...")
        analyzer.compute("correlograms", window_ms=100)

    # Get autocorrelograms
    ccgs, bins = analyzer.get_extension("correlograms").get_data()
    bins = bins[len(bins) // 2:]  # positive half only (symmetric)
    ccgs = ccgs[:, :, len(ccgs[0, 0]) // 2:]

    # Stim-locked bins (around 33 ms for 30 Hz stim)
    acg_range = np.where((bins > 31) & (bins < 35))[0]

    # Timing offsets
    ts_offsets = _get_timing_offsets(analyzer)
    ts_offset_threshold = 10  # samples

    # Spike counts
    spike_counts = analyzer.sorting.count_num_spikes_per_unit()

    # Primary channel templates for waveform classification
    primary_templates = _get_primary_ch_templates(analyzer)

    # ---- Pass 1: Artifact rejection ----
    pass1_good = []
    for i, unit_id in enumerate(unit_ids):
        # Check spike count first
        if spike_counts.get(unit_id, 0) < min_spikes:
            labels[unit_id] = "low_count"
            continue

        # Global artifact check
        if _exclude_artifact_unit(unit_id, analyzer):
            labels[unit_id] = "artifact"
            continue

        # ACG + timing offset check
        acg = ccgs[i, i, :]
        sort_idx = np.argsort(acg)[::-1]
        if np.any(np.isin(sort_idx[:5], acg_range)) and (ts_offsets.get(unit_id, 0) > ts_offset_threshold):
            labels[unit_id] = "artifact"
            continue

        pass1_good.append(unit_id)

    # ---- Pass 2: Waveform shape classification ----
    final_good = []
    for unit_id in pass1_good:
        template = primary_templates[unit_id]
        if _classify_waveform(template):
            labels[unit_id] = "good"
            final_good.append(unit_id)
        else:
            labels[unit_id] = "noise"

    good_ids = _sort_units_by_depth(analyzer, final_good)
    bad_unit_list = [uid for uid in unit_ids if uid not in final_good]
    bad_ids = _sort_units_by_depth(analyzer, bad_unit_list)

    # Print summary
    n_artifact = sum(1 for v in labels.values() if v == "artifact")
    n_noise = sum(1 for v in labels.values() if v == "noise")
    n_low = sum(1 for v in labels.values() if v == "low_count")
    n_good = len(good_ids)
    print(f"\n  Curation Summary:")
    print(f"    {len(unit_ids)} total units")
    print(f"    {n_good} good units")
    print(f"    {n_artifact} artifacts")
    print(f"    {n_noise} noise (bad waveform shape)")
    print(f"    {n_low} low spike count (<{min_spikes})")

    return good_ids, bad_ids, labels
