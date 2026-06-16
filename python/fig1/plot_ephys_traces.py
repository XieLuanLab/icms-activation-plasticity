"""Fig 1F: Raw and artifact-removed electrophysiological traces during ICMS.

Shows 6 channels of raw (left) and preprocessed (right) traces during
100 Hz stimulation. Gray shading marks stimulation period.

Usage:
    python python/fig1/plot_ephys_traces.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style

apply_global_style()

FIG_DIR = OUTPUT_DIR / 'fig1'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    npz_path = DATA_DIR / 'fig1' / 'fig1f_traces.npz'
    if not npz_path.exists():
        print(f"Fig 1F trace data ({npz_path}) is not in the public dataset; "
              "available on request.")
        return
    d = np.load(npz_path)
    raw = d['raw_traces']         # (n_samples, 6)
    preproc = d['preproc_traces']  # (n_samples, 6)
    fs = int(d['fs'])
    stim_onset = float(d['stim_onset_sample'])  # samples from trace start

    n_samples, n_ch = raw.shape
    time_ms = np.arange(n_samples) / (fs / 1000)
    stim_start_ms = stim_onset / (fs / 1000)
    stim_end_ms = stim_start_ms + 700  # 700 ms stim train

    ch_spacing = 6000  # uV offset between channels (doubled from 3000)

    # Trim to -200 ms pre-stim to +200 ms post-stim
    pre_ms = 200
    post_ms = 200
    xlim_start = stim_start_ms - pre_ms
    xlim_end = stim_end_ms + post_ms

    fig, (ax_raw, ax_filt) = plt.subplots(1, 2, figsize=(6, 5), sharey=True)

    for ax, traces, title in [(ax_raw, raw, 'Raw'),
                               (ax_filt, preproc, 'Artifact-removed')]:
        for i in range(n_ch):
            offset = i * ch_spacing
            ax.plot(time_ms, traces[:, i] + offset, 'k', linewidth=0.3)

        # Stim shading
        ax.axvspan(stim_start_ms, stim_end_ms, color='gray', alpha=0.15, zorder=0)

        ax.set_xlim([xlim_start, xlim_end])
        ax.set_title(title, fontsize=8)
        ax.set_xlabel('Time (ms)')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)

    # Scale bars on left panel
    sb_x = xlim_start
    sb_y = -ch_spacing * 0.3
    bar_ms = 100
    ax_raw.plot([sb_x, sb_x + bar_ms], [sb_y, sb_y], 'k', linewidth=1)
    ax_raw.text(sb_x + bar_ms / 2, sb_y - ch_spacing * 0.1, '100 ms',
                ha='center', fontsize=5)
    ax_raw.plot([sb_x, sb_x], [sb_y, sb_y + 500], 'k', linewidth=1)
    ax_raw.text(sb_x - 10, sb_y + 250, r'500 $\mu$V', ha='right', fontsize=5,
                rotation=90, va='center')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'ephys_traces.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "ephys_traces.svg"}')


if __name__ == '__main__':
    main()
