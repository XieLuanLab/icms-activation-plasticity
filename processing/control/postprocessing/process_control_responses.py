"""
Process and plot control session responses.

Workflow:
1. Create session responses (or load from pickle)
2. Extract modulation metrics → CSV
3. Generate per-unit × per-channel figures → PNG
4. Compile figures into PowerPoint

Usage:
    python process_control_responses.py <session_dir>
    python process_control_responses.py <session_dir> --load  # load existing pickle
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for batch processing

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
import dill as pickle
import shutil
import re

import batch_process.postprocessing.responses_v2.response_plotting_util as rpu


def get_output_dirs(session_dir):
    """Get standard output directories for a control session."""
    spike_dir = Path(session_dir) / "spike_sorting"
    figures_dir = spike_dir / "figures" / "stim_condition_figures"
    return spike_dir, figures_dir


def process_stim_responses(session_responses):
    """Extract modulation metrics from all units and save to CSV.

    Mirrors batch_process/postprocessing/responses/process_stim_responses.py
    but adapted for control sessions.
    """
    results = []

    for unit_id in session_responses.unit_ids:
        print(f"  Extracting metrics for unit {unit_id}...")
        unit_response = session_responses.get_unit_response(unit_id)
        stim_channels = sorted(
            {ch for ch, _ in unit_response.stim_responses.keys()})

        for stim_channel in stim_channels:
            currents = sorted(
                [curr for ch, curr in unit_response.stim_responses.keys() if ch == stim_channel])

            for stim_current in currents:
                if stim_current == 0:
                    continue
                scr = unit_response.get_stim_response(stim_channel, stim_current)
                pr = scr.pulse_response
                tr = scr.train_response

                result = {
                    "unit_id": unit_id,
                    "stim_channel": stim_channel,
                    "stim_current": stim_current,
                }

                if pr and tr:
                    result.update({
                        "is_pulse_locked": pr.is_pulse_locked,
                        "max_spike_prob": np.max(pr.stim_metrics["mean_prob_spike"]),
                        "pli": pr.pli,
                        "latency": pr.stim_metrics["mean_latency"],
                        "latency_std": pr.stim_metrics["std_latency"],
                        "modulated": tr.paired_p_val < 0.05,
                        "mod_p_val": tr.paired_p_val,
                        "z_score": tr.z_score,
                        "train_mean_fr": tr.stim_mean_fr,
                        "pre_stim_mean_fr": tr.pre_stim_mean_fr,
                        "pre_stim_std_fr": tr.pre_stim_sigma_fr,
                        "num_spikes": len(pr.rel_spike_timestamps),
                        "num_trials": len(tr.raster_array),
                        "baseline_too_slow": tr.pre_stim_mean_fr < 0.5,
                        "t_to_max_10ms_smoothed": tr.ten_ms_smoothed_data.get('time_to_max_fr', np.nan),
                        "t_to_max_2ms_smoothed": tr.two_ms_smoothed_data.get('time_to_max_fr', np.nan),
                    })

                results.append(result)

    return results


def plot_unit_stim_responses(session_responses, figures_dir):
    """Generate per-unit × per-channel figures using shared rpu plotting."""
    if figures_dir.exists():
        shutil.rmtree(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    for unit_id in session_responses.unit_ids:
        print(f"  Plotting unit {unit_id}...")
        ur = session_responses.get_unit_response(unit_id)
        stim_conditions = ur.show_stim_conditions()
        responses_by_channel = defaultdict(list)

        for stim_ch, current in stim_conditions:
            sr = ur.get_stim_response(stim_ch, current)
            responses_by_channel[stim_ch].append((current, sr))

        for stim_ch, responses in responses_by_channel.items():
            fig, axs = rpu.initialize_plotting()
            rpu.plot_train_rasters(responses, axs)
            rpu.plot_train_firing_rate(responses, axs)
            rpu.plot_pulse_raster(responses, axs)
            rpu.plot_probability(responses, axs)
            rpu.plot_z_scores(responses, axs)
            rpu.plot_template(session_responses, unit_id, axs)

            plt.suptitle(f"Unit {unit_id} - Stim Channel {stim_ch}")
            plt.tight_layout()

            rpu.save_figure(fig, figures_dir, unit_id, stim_ch)


def make_control_ppt(session_dir, figures_dir):
    """Create PowerPoint from saved figures."""
    from batch_process.util.ppt_image_inserter import PPTImageInserter

    image_paths = sorted(figures_dir.glob("*.png"))
    if not image_paths:
        print("  No figures found, skipping PPT creation.")
        return

    # Group by unit ID
    images_by_unit = defaultdict(list)
    for img in image_paths:
        match = re.search(r"Unit_(\d+)_", img.name)
        if match:
            images_by_unit[int(match.group(1))].append(img)

    session_name = Path(session_dir).name
    ppt = PPTImageInserter(grid_dims=(2, 2), spacing=(0.05, 0.05), title_font_size=16)

    for unit_id, imgs in sorted(images_by_unit.items()):
        title = f"{session_name}: unit {unit_id}"
        ppt.add_slide(title)
        for img in imgs:
            ppt.add_image(str(img), slide_title=title)

    spike_dir = Path(session_dir) / "spike_sorting"
    ppt_path = spike_dir / "figures" / f"stim_responses_{session_name}.pptx"
    ppt.save(ppt_path)
    print(f"  PPT saved to {ppt_path}")


def run(session_dir, load_pickle=False):
    """Full pipeline: create responses → metrics CSV → plots → PPT."""
    spike_dir, figures_dir = get_output_dirs(session_dir)

    if load_pickle:
        pkl_path = spike_dir / "session_responses.pkl"
        print(f"Loading from {pkl_path}...")
        with open(pkl_path, "rb") as f:
            session_responses = pickle.load(f)
        # Ensure sorting analyzer is loaded for plotting
        session_responses.load_sorting_analyzer()
    else:
        from control.postprocessing.create_session_responses import (
            create_session_responses, save_session_responses,
        )
        print(f"Creating session responses for {session_dir}...")
        session_responses = create_session_responses(session_dir)
        save_session_responses(session_responses, spike_dir)

    # Extract metrics and save CSV
    print("Extracting modulation metrics...")
    results = process_stim_responses(session_responses)
    rpu.save_results_to_csv(results, spike_dir, "control")

    # Generate figures
    print("Generating figures...")
    plot_unit_stim_responses(session_responses, figures_dir)

    # Create PPT
    print("Creating PowerPoint...")
    make_control_ppt(session_dir, figures_dir)

    print("Done!")
    return session_responses


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python process_control_responses.py <session_dir> [--load]")
        sys.exit(1)

    session_dir = sys.argv[1]
    load_flag = "--load" in sys.argv
    run(session_dir, load_pickle=load_flag)
