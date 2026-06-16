# Waveform Curation Process Reference

This documents the waveform curation pipeline used in `stage2_curate.py`.
For each unit, spikes go through 4 sequential steps. A spike must survive
all steps to be kept in `analyzer_curated.zarr`.

---

## Step 1: Remove Bad Waveforms

Two functions run in sequence. A spike rejected by either is removed.

### remove_bad_waveforms_A (`curate_util.py:242`)
Detects blanking artifacts and corrupted waveforms on the primary channel.

| Check | Criteria | Rationale |
|-------|----------|-----------|
| Consecutive zeros (range 1) | 3+ zeros in samples [24:31] | Data dropout / blanking |
| Consecutive zeros (range 2) | 3+ zeros in samples [32:38] | Data dropout / blanking |
| Pre-peak amplitude | Any sample in [29:32] > -40 uV | Abnormal baseline |

### remove_bad_waveforms_B (`curate_util.py:304`)
Removes amplitude outliers using a hybrid threshold.

| Parameter | Default | Description |
|-----------|---------|-------------|
| hard_threshold | -50 uV | Absolute amplitude cutoff |
| num_stds | 10 | Statistical range = median +/- 10*SD |

A spike is kept only if:
- amplitude >= median - 10*SD
- amplitude <= median + 10*SD
- amplitude <= -50 uV (hard threshold)

---

## Step 2: Classify Spike-Like vs Non-Spike-Like

Uses raw (unfiltered) waveforms centered at sample 30.

### MAD Threshold (`spike_classifier.py:13`)
- Extracts 10 seconds of raw recording from the primary channel
- Computes MAD (Median Absolute Deviation)
- Converts to threshold: `2 * MAD / 0.6745`

### classify_waveform (`spike_classifier.py:226`)
A waveform must pass ALL 9 criteria to be classified as "spike-like":

| # | Criterion | Threshold | Description |
|---|-----------|-----------|-------------|
| 1 | Slope check | immediate < early | Pre-peak slope must steepen approaching trough |
| 2 | Post-trough slope | > 0 | Must rise after trough (samples [30:36]) |
| 3 | Base-to-trough amplitude | > 60 uV | Minimum spike amplitude |
| 4 | Peak count | == 1 | Only one peak (prominence >= 50, width >= 3) in samples [18:50] |
| 5 | Peak prominence | >= 80 uV | Minimum repolarization peak prominence |
| 6 | Peak width (50%) | 4 - 19 samples | Width at half-height |
| 7 | Peak width (80%) | 6 - 22 samples | Width at 80% height |
| 8 | Amplitude ratio | < 2 | base_to_peak / base_to_trough |
| 9 | Oscillation | < 0.15 | Prevalence of alternating [+,-,+,-] pattern |

Non-spike-like waveforms are labeled 0 and excluded from subclustering.

---

## Step 3: Subcluster with UMAP + Isosplit6

Operates on spike-like waveforms only, using samples [20:42] of the
primary channel (the peak region).

### cluster_waveforms (`subcluster_util.py:273`)
Iterative dimensionality reduction + clustering:

| Iteration | n_neighbors | min_dist | Action |
|-----------|-------------|----------|--------|
| 1 | 50 | 0.1 | If <= 4 clusters, stop |
| 2 | 60 | 0.2 | If <= 4 clusters, stop |
| 3 | 70 | 0.3 | If <= 4 clusters, stop |
| 4 | 80 | 0.4 | If <= 4 clusters, stop |
| 5 | 90 | 0.5 | Final attempt |

Each iteration increases smoothing (larger n_neighbors, larger min_dist)
to reduce over-splitting. Max 4 subclusters allowed.

---

## Step 4: Evaluate Subclusters

Each subcluster is independently evaluated. Checks run in order;
first failure rejects the subcluster.

### accept_subcluster (`subcluster_util.py:226`)

| # | Check | Threshold | Reject label |
|---|-------|-----------|--------------|
| 1 | Artifact (MNR) | MNR > 0.25 | "artifact" |
| 2 | Peak location | abs max not in [29, 31] | "max_idx_out_of_range" |
| 3 | Trough-to-peak duration | > 30 samples | "large_T2P" |
| 4 | Pre-trough slope | slope < -5 | "negative_slope" |
| 5 | SNR | < 2.0 | "low_snr" |
| 6 | Spike count | <= 50 | (too few spikes) |

**MNR (Mean Normalized Relative peak):** Measures how much the spike
extends beyond the primary channel's local neighborhood (120 um radius).
High MNR = artifact spreading across many channels.

---

## What to Do If Curation Looks Wrong

### Too many spikes rejected
- Check `curation_diagnostic.pptx` — middle panel (accepted) vs right panel (rejected)
- If good spikes are in the rejected panel:
  - `classify_waveform` may be too strict (amplitude, width thresholds)
  - `remove_bad_waveforms_B` hard_threshold may be wrong for this animal
  - Location: `spike_classifier.py:226` and `curate_util.py:304`

### Good unit entirely dropped
- Check which step killed it in the stage2 output log
- Common causes:
  - `low_snr`: SNR < 2.0 threshold in `accept_subcluster`
  - `large_T2P`: trough-to-peak > 30 samples (slow waveform?)
  - Too few spike-like waveforms (< 20) after classification

### Over-splitting
- UMAP produced too many subclusters for the same neuron
- Solution: add the split IDs as a merge pair in `curation_config.json`
  e.g., `"merge_pairs": [[101001, 101002]]`

### Under-splitting (mixed unit survived)
- Isosplit6 didn't separate distinct waveform shapes
- May need to lower `n_neighbors` or `min_dist` starting values
- Location: `subcluster_util.py:273`

---

## File Locations

| File | Functions |
|------|-----------|
| `batch_process/util/curate_util.py` | remove_bad_waveforms_A, _B |
| `batch_process/visualization/spike_classifier.py` | get_MAD_threshold_primary_ch, classify_waveform, extract_waveform_features |
| `batch_process/util/subcluster_util.py` | cluster_waveforms, accept_subcluster |
| `control/stage2_curate.py` | curate_unit_waveforms (orchestrates all 4 steps) |
