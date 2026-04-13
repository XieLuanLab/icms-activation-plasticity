"""Raw neural trace with stim artifacts and detected spike overlay."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from utils.plotting import apply_global_style, PALETTE
from utils.config import DATA_DIR, OUTPUT_DIR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')


FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()


def main():
    d = np.load(DATA_DIR / 'fig4' / 'raw_trace_example.npz')
    raw_trace = d['raw_trace']
    spike_frames = d['spike_frames']
    stim_ts = d['stim_ts']
    fs = int(d['fs'])

    time_ms = np.arange(len(raw_trace)) / (fs / 1000)

    fig, ax = plt.subplots(figsize=(3, 2))

    stim_ms = stim_ts / (fs / 1000)
    blank_before = 0.5  # ms before pulse
    blank_after = 1.5   # ms after pulse
    x_start = stim_ms[0] - blank_before
    x_end = stim_ms[-1] + blank_after

    # Raw trace
    ax.plot(time_ms, raw_trace, 'k', linewidth=0.5)

    for i, s in enumerate(stim_ms):
        ax.axvspan(s - blank_before, s + blank_after,
                   color='red', alpha=0.08, zorder=0,
                   label='Blanked' if i == 0 else '')
        ax.axvline(s, ymin=0, ymax=0.85, color=PALETTE[3], linewidth=0.5,
                   label='Stim onset' if i == 0 else '')

    spike_ms = spike_frames / (fs / 1000)
    visible = (spike_ms >= x_start) & (spike_ms <= x_end)
    spike_ms_vis = spike_ms[visible]
    spike_y_vis = raw_trace[spike_frames[visible].astype(int)]
    ax.scatter(spike_ms_vis, spike_y_vis - 180, color=PALETTE[0], marker='^', s=12,
               zorder=3, label='Detected spike')

    ax.set_xlim(x_start, x_end)
    ax.set_ylim([raw_trace.min() - 100, raw_trace.max() + 200])
    ax.axis('off')

    y_lo = ax.get_ylim()[0]
    sb_x = x_start
    sb_y = y_lo - 50
    ax.plot([sb_x, sb_x + 5], [sb_y, sb_y], 'k', linewidth=1, clip_on=False)
    ax.text(sb_x + 2.5, sb_y - 150, '5 ms', ha='center', fontsize=5)
    ax.plot([sb_x, sb_x], [sb_y, sb_y + 200], 'k', linewidth=1, clip_on=False)
    ax.text(sb_x - 0.5, sb_y + 100, r'200 $\mu$V', ha='right', fontsize=5,
            rotation=90, va='center')

    # ax.legend(fontsize=5, loc='lower center', frameon=False)

    fig.savefig(FIG_DIR / 'raw_trace.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "raw_trace.svg"}')


if __name__ == '__main__':
    main()
