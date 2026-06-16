"""
Merge Suggestion Utility
========================
Computes pairwise merge candidates among good units using three criteria:
  1. Template similarity (cosine on primary channel)
  2. Cross-correlogram refractory dip
  3. Channel proximity

Usage:
    candidates = compute_merge_candidates(analyzer, good_ids)
    print_merge_report(candidates)
"""

import numpy as np
from itertools import combinations


def _cosine_similarity(a, b):
    """Cosine similarity between two 1D vectors."""
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    if norm == 0:
        return 0.0
    return float(np.dot(a, b) / norm)


def get_template_similarity(analyzer, unit_a, unit_b):
    """
    Cosine similarity between the templates of two units on their
    respective primary channels (channel with max PTP).

    Returns
    -------
    sim : float
        Cosine similarity in [0, 1].  Higher = more similar.
    """
    templates_ext = analyzer.get_extension("templates")
    tpl_a = templates_ext.get_unit_template(unit_a)  # (n_samples, n_channels)
    tpl_b = templates_ext.get_unit_template(unit_b)

    # Primary channel for each unit = channel with largest peak-to-peak
    primary_ch_a = int(np.argmax(np.ptp(tpl_a, axis=0)))
    primary_ch_b = int(np.argmax(np.ptp(tpl_b, axis=0)))

    # Compare unit A's template on A's primary channel
    # vs unit B's template on the SAME channel (A's primary)
    # AND vice versa; take the max similarity
    sim_on_a = _cosine_similarity(tpl_a[:, primary_ch_a], tpl_b[:, primary_ch_a])
    sim_on_b = _cosine_similarity(tpl_a[:, primary_ch_b], tpl_b[:, primary_ch_b])

    return max(sim_on_a, sim_on_b)


def get_ccg_refractory_score(analyzer, unit_a, unit_b,
                              refractory_ms=1.5, bin_ms=0.2, window_ms=50.0):
    """
    Compute the cross-correlogram between two units and check for a
    refractory period dip in the central bins.

    Returns
    -------
    ratio : float
        Ratio of mean count in refractory zone to peak CCG bin (cross-contamination estimate).
        Lower = stronger refractory dip = more likely same neuron.
        Returns 1.0 if CCG is flat/empty.
    violations : int
        The raw number of spikes in the CCG within the refractory period.
    """
    correlograms_ext = analyzer.get_extension("correlograms")
    if correlograms_ext is None:
        return 1.0, 0

    ccg_data = correlograms_ext.get_data()
    # get_data() returns (ccgs_array, bins_array) tuple
    if isinstance(ccg_data, tuple):
        ccgs, time_bins = ccg_data
    else:
        ccgs = ccg_data
        time_bins = None

    actual_bin_ms = correlograms_ext.params.get("bin_ms", bin_ms)

    unit_ids = list(analyzer.unit_ids)
    idx_a = unit_ids.index(unit_a)
    idx_b = unit_ids.index(unit_b)

    ccg = ccgs[idx_a, idx_b, :]  # 1D array of bin counts
    n_bins = len(ccg)
    center = n_bins // 2

    # Number of bins that fall within the refractory zone
    refractory_bins = int(np.ceil(refractory_ms / actual_bin_ms))

    # Refractory zone: center ± refractory_bins
    ref_start = max(0, center - refractory_bins)
    ref_end = min(n_bins, center + refractory_bins + 1)

    refractory_sum = int(np.sum(ccg[ref_start:ref_end])) if ref_end > ref_start else 0
    refractory_count = float(np.mean(ccg[ref_start:ref_end])) if ref_end > ref_start else 0.0

    # Peak: max bin count outside the refractory zone
    mask = np.ones(n_bins, dtype=bool)
    mask[ref_start:ref_end] = False
    outside = ccg[mask]
    peak_count = float(np.max(outside)) if len(outside) > 0 and np.max(outside) > 0 else 1.0

    ratio = refractory_count / peak_count
    return float(ratio), refractory_sum

def get_amplitude_ratio(analyzer, unit_a, unit_b):
    """
    Compute the ratio of peak-to-peak amplitudes between two units
    on their shared best channel. Returns a value between 0 and 1.
    """
    templates_ext = analyzer.get_extension("templates")
    tpl_a = templates_ext.get_unit_template(unit_a)
    tpl_b = templates_ext.get_unit_template(unit_b)

    primary_ch_a = int(np.argmax(np.ptp(tpl_a, axis=0)))
    primary_ch_b = int(np.argmax(np.ptp(tpl_b, axis=0)))
    
    # Use the channel with the larger combined signal
    combined_a = np.ptp(tpl_a[:, primary_ch_a]) + np.ptp(tpl_b[:, primary_ch_a])
    combined_b = np.ptp(tpl_a[:, primary_ch_b]) + np.ptp(tpl_b[:, primary_ch_b])
    shared_ch = primary_ch_a if combined_a >= combined_b else primary_ch_b
    
    ptp_a = np.ptp(tpl_a[:, shared_ch])
    ptp_b = np.ptp(tpl_b[:, shared_ch])
    
    if max(ptp_a, ptp_b) == 0:
        return 0.0
        
    return float(min(ptp_a, ptp_b) / max(ptp_a, ptp_b))


def get_channel_distance(analyzer, unit_a, unit_b):
    """
    Distance between primary channels of two units, measured in
    electrode positions (integer index difference).

    Returns
    -------
    dist : int
        Absolute index difference between primary channels.
    """
    templates_ext = analyzer.get_extension("templates")
    tpl_a = templates_ext.get_unit_template(unit_a)
    tpl_b = templates_ext.get_unit_template(unit_b)

    primary_ch_a = int(np.argmax(np.ptp(tpl_a, axis=0)))
    primary_ch_b = int(np.argmax(np.ptp(tpl_b, axis=0)))

    return abs(primary_ch_a - primary_ch_b)


def compute_merge_candidates(analyzer, good_ids, slay_merges=None):
    """
    Compute pairwise merge stats for all good_ids using custom biological logic,
    and include any pairs flagged by SLAy.

    Parameters
    ----------
    analyzer : SortingAnalyzer
        Must have 'templates' and 'correlograms' extensions computed.
    good_ids : list
        List of unit IDs to consider for merging.
    slay_merges : list of tuples
        Optional list of (unit_a, unit_b) returned natively by SLAy.

    Returns
    -------
    candidates : list of dict
        Each dict has: unit_a, unit_b, template_sim, ccg_contam_ratio, ccg_violations,
        channel_distance, amp_ratio, recommendation.
    """
    candidates = []
    
    # Custom thresholds
    sim_thresh = 0.90
    ccg_contam_thresh = 0.20
    ch_dist_thresh = 2
    amp_ratio_thresh = 0.50
    
    # Convert slay_merges to a set of sorted tuples for easy lookup
    slay_set = set()
    if slay_merges:
        for uA, uB in slay_merges:
            slay_set.add(tuple(sorted([uA, uB])))

    for uA, uB in combinations(good_ids, 2):
        tpl_sim = get_template_similarity(analyzer, uA, uB)
        ccg_ratio, ccg_viol = get_ccg_refractory_score(analyzer, uA, uB)
        ch_dist = get_channel_distance(analyzer, uA, uB)
        amp_ratio = get_amplitude_ratio(analyzer, uA, uB)
        
        is_slay = tuple(sorted([uA, uB])) in slay_set
        
        # Check Custom Biological Criteria
        pass_sim = tpl_sim >= sim_thresh
        pass_ccg = ccg_ratio <= ccg_contam_thresh
        pass_ch = ch_dist <= ch_dist_thresh
        pass_amp = amp_ratio >= amp_ratio_thresh
        
        # Determine strict biological recommendation
        recs = []
        if pass_sim and pass_ch and pass_amp:
            # We explicitly ignore CCG pass/fail here because of ICMS artifacts!
            # If they look identical, are on the exact same channels, and same size: MERGE
            recs.append("MERGE (Sim/Ch/Amp)")
        elif pass_sim and pass_ch:
            recs.append("REVIEW (Sim/Ch)")
            
        if is_slay:
            recs.append("SLAy")
            
        if not recs:
            continue
            
        recommendation = " + ".join(recs)

        candidates.append({
            "unit_a": uA,
            "unit_b": uB,
            "template_sim": round(tpl_sim, 3),
            "ccg_contam_ratio": round(ccg_ratio, 3),
            "ccg_violations": int(ccg_viol),
            "channel_distance": ch_dist,
            "amp_ratio": round(amp_ratio, 3),
            "recommendation": recommendation
        })

    return candidates


def print_merge_report(candidates, session_name=""):
    """Pretty-print merge candidates."""
    header = f"MERGE CANDIDATE REPORT"
    if session_name:
        header += f" -- {session_name}"

    print(f"\n{'=' * 72}")
    print(header)
    print(f"{'=' * 72}")

    merge = candidates

    if not merge:
        print("\n  No merge candidates found by custom logic or SLAy.\n")
        print(f"{'=' * 72}")
        return

    def _print_row(c):
        print(f"  | {c['unit_a']:>5} | {c['unit_b']:>5} | "
              f"{c['template_sim']:>10.3f} | "
              f"{c['ccg_contam_ratio']:>10.3f} | "
              f"{c['ccg_violations']:>7d} | "
              f"{c['channel_distance']:>8} | "
              f"{c['amp_ratio']:>10.3f} | "
              f"{c['recommendation']} |")

    print(f"\n  Merge Candidates:")
    print(f"  {'-' * 110}")
    print(f"  | {'UnitA':>5} | {'UnitB':>5} | {'TplSim':>10} | "
          f"{'CCGContam':>10} | {'CCGViol':>7} | {'ChDist':>8} | {'AmpRatio':>10} | {'Recommend'} |")
    print(f"  {'-' * 110}")
    for c in sorted(merge, key=lambda x: str(x['recommendation']), reverse=True):
        _print_row(c)
    print(f"  {'-' * 110}")

    print(f"\n  --> Copy confirmed pairs to MERGE_PAIRS and re-run stage2.")
    print(f"      e.g. MERGE_PAIRS = [({merge[0]['unit_a']}, {merge[0]['unit_b']})]"
          if merge else "")
    print(f"{'=' * 72}\n")


def _extract_wide_snippets(recording, spike_train, channel_idx,
                           n_snippets=100, half_width_ms=5.0, fs=30000.0):
    """Extract wide raw snippets centered on spike times for one channel."""
    rng = np.random.default_rng(42)
    half_samples = int(half_width_ms * fs / 1000)
    n_total = recording.get_num_samples()

    # Filter valid spike indices
    valid = spike_train[(spike_train > half_samples) &
                        (spike_train < n_total - half_samples)]
    if len(valid) == 0:
        return None, None

    n_pick = min(n_snippets, len(valid))
    chosen = valid[rng.choice(len(valid), n_pick, replace=False)]

    snippets = np.zeros((n_pick, 2 * half_samples))
    for i, t in enumerate(chosen):
        seg = recording.get_traces(
            start_frame=int(t - half_samples),
            end_frame=int(t + half_samples),
            return_scaled=True,
        )
        snippets[i] = seg[:, channel_idx]

    time_ms = np.arange(2 * half_samples) / fs * 1000 - half_width_ms
    return snippets, time_ms


def generate_merge_comparison_figure(analyzer, unit_a, unit_b,
                                     raw_recording=None, context_ms=15.0):
    """
    Generate a 2-row comparison figure for a merge candidate pair.

    Row 1: [Unit A template] [Unit B template] [Overlay] [Cross-correlogram]
    Row 2: [Unit A raw waveforms (wide)] [Unit B raw waveforms (wide)] [Overlay raw]

    Parameters
    ----------
    analyzer : SortingAnalyzer
    unit_a, unit_b : unit IDs
    raw_recording : BaseRecording, optional
        Unprocessed recording for raw waveform display. If None, uses
        analyzer.recording (preprocessed).
    context_ms : float
        Half-width of context window in ms (default ±15 ms).

    Returns
    -------
    fig : matplotlib.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    templates_ext = analyzer.get_extension("templates")
    tpl_a = templates_ext.get_unit_template(unit_a)
    tpl_b = templates_ext.get_unit_template(unit_b)

    # Primary channels
    primary_a = int(np.argmax(np.ptp(tpl_a, axis=0)))
    primary_b = int(np.argmax(np.ptp(tpl_b, axis=0)))
    ptp_a = np.ptp(tpl_a[:, primary_a])
    ptp_b = np.ptp(tpl_b[:, primary_b])

    # For overlay, use the channel that is closest to both primaries
    if primary_a == primary_b:
        shared_ch = primary_a
    else:
        combined_a = np.ptp(tpl_a[:, primary_a]) + np.ptp(tpl_b[:, primary_a])
        combined_b = np.ptp(tpl_a[:, primary_b]) + np.ptp(tpl_b[:, primary_b])
        shared_ch = primary_a if combined_a >= combined_b else primary_b

    # Spike counts
    spike_train_a = analyzer.sorting.get_unit_spike_train(unit_a)
    spike_train_b = analyzer.sorting.get_unit_spike_train(unit_b)
    n_spikes_a = spike_train_a.size
    n_spikes_b = spike_train_b.size

    fig = plt.figure(figsize=(16, 8))
    gs = GridSpec(2, 4, width_ratios=[1, 1, 1, 1.2], hspace=0.35)

    # ===== ROW 1: Templates + CCG =====

    # --- Panel 1: Unit A template ---
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(tpl_a[:, primary_a], 'b-', linewidth=1.5)
    ax0.set_title(f"Unit {unit_a} (n={n_spikes_a})\nCh {primary_a} | PTP {ptp_a:.0f} uV", fontsize=9)
    ax0.set_ylabel("Amplitude (uV)", fontsize=8)

    # --- Panel 2: Unit B template ---
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)
    ax1.plot(tpl_b[:, primary_b], 'r-', linewidth=1.5)
    ax1.set_title(f"Unit {unit_b} (n={n_spikes_b})\nCh {primary_b} | PTP {ptp_b:.0f} uV", fontsize=9)

    # --- Panel 3: Overlay on shared channel ---
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(tpl_a[:, shared_ch], 'b-', linewidth=1.5, alpha=0.8, label=f"Unit {unit_a}")
    ax2.plot(tpl_b[:, shared_ch], 'r-', linewidth=1.5, alpha=0.8, label=f"Unit {unit_b}")
    sim = _cosine_similarity(tpl_a[:, shared_ch], tpl_b[:, shared_ch])
    ch_dist = abs(primary_a - primary_b)
    amp_ratio = get_amplitude_ratio(analyzer, unit_a, unit_b)
    ax2.set_title(f"Overlay (Ch {shared_ch})\nSim: {sim:.3f} | ChDist: {ch_dist} | AmpRatio: {amp_ratio:.3f}", fontsize=9)
    ax2.legend(fontsize=7, loc="lower right")

    # --- Panel 4: Cross-correlogram ---
    ax3 = fig.add_subplot(gs[0, 3])
    correlograms_ext = analyzer.get_extension("correlograms")
    if correlograms_ext is not None:
        ccg_data = correlograms_ext.get_data()
        if isinstance(ccg_data, tuple):
            ccgs, _ = ccg_data
        else:
            ccgs = ccg_data
        unit_ids = list(analyzer.unit_ids)
        idx_a = unit_ids.index(unit_a)
        idx_b = unit_ids.index(unit_b)
        ccg = ccgs[idx_a, idx_b, :]
        bin_ms = correlograms_ext.params.get("bin_ms", 0.5)
        n_bins = len(ccg)
        time_axis = (np.arange(n_bins) - n_bins // 2) * bin_ms
        ax3.bar(time_axis, ccg, width=bin_ms * 0.9, color="k", alpha=0.7)
        ax3.axvspan(-1.5, 1.5, color="red", alpha=0.15, label="Refractory (1.5ms)")
        _, ccg_viol = get_ccg_refractory_score(analyzer, unit_a, unit_b)
        ax3.set_title(f"Cross-correlogram\n{unit_a} x {unit_b} | Violations: {int(ccg_viol)}", fontsize=9)
        ax3.set_xlabel("Lag (ms)", fontsize=8)
        ax3.legend(fontsize=7)
    else:
        ax3.text(0.5, 0.5, "No CCG data", transform=ax3.transAxes, ha="center")

    # ===== ROW 2: RAW waveforms with wide context =====
    rec = raw_recording if raw_recording is not None else analyzer.recording
    rec_label = "RAW" if raw_recording is not None else "preprocessed"
    has_rec = rec is not None

    if has_rec:
        snippets_a, time_a = _extract_wide_snippets(
            rec, spike_train_a, primary_a, n_snippets=25,
            half_width_ms=context_ms)
        snippets_b, time_b = _extract_wide_snippets(
            rec, spike_train_b, primary_b, n_snippets=25,
            half_width_ms=context_ms)
        snippets_a_shared, _ = _extract_wide_snippets(
            rec, spike_train_a, shared_ch, n_snippets=25,
            half_width_ms=context_ms)
        snippets_b_shared, _ = _extract_wide_snippets(
            rec, spike_train_b, shared_ch, n_snippets=25,
            half_width_ms=context_ms)
        # Median-subtract each snippet to zero-center (remove DC offset)
        for snips in [snippets_a, snippets_b, snippets_a_shared, snippets_b_shared]:
            if snips is not None:
                snips -= np.median(snips, axis=1, keepdims=True)
    else:
        snippets_a = snippets_b = None
        snippets_a_shared = snippets_b_shared = None

    # Fixed y-limits for raw waveform panels
    _ylim = (-500, 200)

    # --- Panel 5: Unit A raw waveforms ---
    ax4 = fig.add_subplot(gs[1, 0])
    if snippets_a is not None:
        ax4.plot(time_a, snippets_a.T, color='b', alpha=0.08, linewidth=0.5)
        ax4.plot(time_a, np.median(snippets_a, axis=0), 'b-', linewidth=1.5)
        ax4.set_title(f"Unit {unit_a} {rec_label} (Ch {primary_a})\n25 wvf, +/-{context_ms:.0f} ms", fontsize=9)
        if _ylim:
            ax4.set_ylim(_ylim)
    else:
        ax4.text(0.5, 0.5, "No recording", transform=ax4.transAxes, ha="center")
    ax4.set_xlabel("Time (ms)", fontsize=8)
    ax4.set_ylabel("Amplitude (uV)", fontsize=8)

    # --- Panel 6: Unit B raw waveforms ---
    ax5 = fig.add_subplot(gs[1, 1], sharey=ax4)
    if snippets_b is not None:
        ax5.plot(time_b, snippets_b.T, color='r', alpha=0.08, linewidth=0.5)
        ax5.plot(time_b, np.median(snippets_b, axis=0), 'r-', linewidth=1.5)
        ax5.set_title(f"Unit {unit_b} {rec_label} (Ch {primary_b})\n25 wvf, +/-{context_ms:.0f} ms", fontsize=9)
    else:
        ax5.text(0.5, 0.5, "No recording", transform=ax5.transAxes, ha="center")
    ax5.set_xlabel("Time (ms)", fontsize=8)

    # --- Panel 7: Raw overlay on shared channel ---
    ax6 = fig.add_subplot(gs[1, 2], sharey=ax4)
    if snippets_a_shared is not None and snippets_b_shared is not None:
        ax6.plot(time_a, snippets_a_shared.T, color='b', alpha=0.05, linewidth=0.5)
        ax6.plot(time_a, snippets_b_shared.T, color='r', alpha=0.05, linewidth=0.5)
        ax6.plot(time_a, np.median(snippets_a_shared, axis=0), 'b-', linewidth=1.5, label=f"{unit_a}")
        ax6.plot(time_a, np.median(snippets_b_shared, axis=0), 'r-', linewidth=1.5, label=f"{unit_b}")
        ax6.set_title(f"{rec_label} overlay (Ch {shared_ch})\n+/-{context_ms:.0f} ms context", fontsize=9)
        ax6.legend(fontsize=7, loc="lower right")
    else:
        ax6.text(0.5, 0.5, "No recording", transform=ax6.transAxes, ha="center")
    ax6.set_xlabel("Time (ms)", fontsize=8)

    # --- Panel 8: empty ---
    ax7 = fig.add_subplot(gs[1, 3])
    ax7.axis("off")

    plt.tight_layout()
    return fig
