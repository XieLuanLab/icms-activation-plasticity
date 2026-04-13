"""Raw vs filtered individual spike waveforms."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from utils.plotting import apply_global_style
from utils.config import DATA_DIR, OUTPUT_DIR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')


FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()


def main():
    d = np.load(DATA_DIR / 'fig4' / 'raw_filt_spikes.npz')
    raw_spikes = d['raw_spikes']
    filt_spikes = d['filt_spikes']
    fs = int(d['fs'])

    # Filter to good spikes
    good = np.where(
        (np.argmax(raw_spikes, axis=1) < 30) &
        (np.min(raw_spikes, axis=1) > -4000) &
        (np.min(filt_spikes, axis=1) > -220) &
        (np.all(filt_spikes[:, 75:] > -100, axis=1)) &
        (raw_spikes[:, 23] < -1700) &
        (raw_spikes[:, 37] > -2900) &
        (filt_spikes[:, 119] < 100)
    )[0]

    good_raw = raw_spikes[good][:, 12:]  # trim leading samples
    good_filt = filt_spikes[good][:, 12:]

    fig, axes = plt.subplots(2, 1, figsize=(1.2, 2.2))

    # Raw spikes
    axes[0].plot(good_raw.T, 'k', linewidth=0.5, alpha=0.25)
    axes[0].set_ylim([-3300, -1200])
    axes[0].set_title('Raw spikes', fontsize=6)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].axis('off')

    # Scale bars for raw
    axes[0].plot([10, 40], [-3150, -3150], 'k', linewidth=1)  # 1 ms
    axes[0].text(25, -3420, '1 ms', ha='center', fontsize=5)
    axes[0].plot([10, 10], [-3150, -2950], 'k', linewidth=1)  # 200 uV
    axes[0].text(5, -3050, '200 \u00b5V', ha='right',
                 fontsize=5, rotation=90, va='center')

    # Filtered spikes
    axes[1].plot(good_filt.T, 'k', linewidth=0.5, alpha=0.25)
    axes[1].set_ylim([-350, 200])
    axes[1].set_title('Filtered spikes', fontsize=6)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].axis('off')

    # Scale bars for filtered
    axes[1].plot([10, 40], [-320, -320], 'k', linewidth=1)
    axes[1].text(25, -380, '1 ms', ha='center', fontsize=5)
    axes[1].plot([10, 10], [-320, -220], 'k', linewidth=1)
    axes[1].text(5, -270, '100 \u00b5V', ha='right',
                 fontsize=5, rotation=90, va='center')

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'raw_filt_spikes.svg',
                format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "raw_filt_spikes.svg"}')


if __name__ == '__main__':
    main()
