# T2P histogram of all units vs modulated units, pyramidal/interneuron split.
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils.config import RAW_DF_PATH, OUTPUT_DIR
from utils.plotting import apply_global_style, PALETTE
from utils.filters import filter_modulated
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

matplotlib.use('Agg')


FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    apply_global_style()
    threshold = 0.45
    bins = np.linspace(0.2, 0.9, 23)

    raw_df = pd.read_pickle(RAW_DF_PATH)
    fig, axes = plt.subplots(2, 1, figsize=(
        2.3, 1.8), sharex=True, sharey=True)

    # All units
    df_unique = raw_df.drop_duplicates(
        subset=['animal_id', 'session', 'unit_id'])
    t2p = df_unique['t2p_ms']
    t2p = t2p[t2p > 0]

    axes[0].hist(t2p[t2p < threshold], bins=bins, color=PALETTE[7])
    axes[0].hist(t2p[t2p >= threshold], bins=bins, color=PALETTE[3])
    axes[0].axvline(threshold, c='k', linestyle='--')
    axes[0].set_title("All Units")
    axes[0].set_ylabel("Count")
    axes[0].text(threshold * 1.01, 70,
                 f'{threshold:.2f} ms', rotation=270, va='top')

    # Modulated units
    df_mod = filter_modulated(raw_df)
    df_unique_mod = df_mod.drop_duplicates(
        subset=['animal_id', 'session', 'unit_id'])
    t2p_mod = df_unique_mod['t2p_ms']
    t2p_mod = t2p_mod[t2p_mod > 0]

    axes[1].hist(t2p_mod[t2p_mod < threshold], bins=bins, color=PALETTE[7])
    axes[1].hist(t2p_mod[t2p_mod >= threshold], bins=bins, color=PALETTE[3])
    axes[1].axvline(threshold, c='k', linestyle='--')
    axes[1].set_title("Modulated Units")
    axes[1].set_xlabel("Trough-to-Peak Time (ms)")
    axes[1].set_ylabel("Count")
    axes[1].text(threshold * 1.01, 70,
                 f'{threshold:.2f} ms', rotation=270, va='top')

    axes[0].set_ylim([0, 80])
    axes[1].set_ylim([0, 80])
    axes[1].set_xlim([0.2, 0.9])

    N_py_all = len(t2p[t2p >= threshold])
    N_in_all = len(t2p[t2p < threshold])
    N_py_mod = len(t2p_mod[t2p_mod >= threshold])
    N_in_mod = len(t2p_mod[t2p_mod < threshold])
    print(f"Pyramidal (all): {N_py_all}, Interneuron (all): {N_in_all}")
    print(f"Pyramidal (mod): {N_py_mod}, Interneuron (mod): {N_in_mod}")

    plt.tight_layout()
    fig.savefig(FIG_DIR / 't2p_histogram.svg',
                format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "t2p_histogram.svg"}')


if __name__ == '__main__':
    main()
