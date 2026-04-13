"""NET probe diagram with unit templates colored by stim-evoked firing rate."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import json
from matplotlib.patches import Polygon, Circle

from utils.config import DATA_DIR, OUTPUT_DIR
from utils.plotting import apply_global_style

FIG_DIR = OUTPUT_DIR / 'fig4'
FIG_DIR.mkdir(parents=True, exist_ok=True)

apply_global_style()

PROBE_JSON = DATA_DIR / 'net32Ch.json'


def main():
    # Load probe geometry
    with open(PROBE_JSON) as f:
        probe_data = json.load(f)

    probe = probe_data['probes'][0]
    contact_positions = np.array(probe['contact_positions'])
    contour = np.array(probe['probe_planar_contour'])
    contour[:, 0] *= 3
    contour[:, 1] *= 1.2

    # Load unit template data
    d = np.load(DATA_DIR / 'fig4' / 'probe_unit_templates.npz', allow_pickle=True)
    stim_ch = int(d['stim_ch'])
    stim_cu = int(d['stim_current'])
    unit_ids = list(d['unit_ids'])

    # Build unit dict sorted by y (depth)
    units = []
    for uid in unit_ids:
        units.append({
            'uid': uid,
            'template': d[f'template_{uid}'],
            'y': float(d[f'y_{uid}']),
            'fr': float(d[f'fr_{uid}']),
        })
    units.sort(key=lambda u: u['y'], reverse=True)

    # Plot
    fig, ax = plt.subplots(figsize=(2, 5))

    # Probe body
    polygon = Polygon(contour, closed=True, fill=True, color='#f6f0e0')
    ax.add_patch(polygon)

    # Contact sites
    for pos in contact_positions:
        circle = Circle((0, pos[1]), radius=14, color='#eed79d', fill=True)
        ax.add_patch(circle)

    # Stim channel highlight
    stim_y = 32 * 60 - stim_ch * 60
    stim_circle = Circle((0, stim_y), radius=14, color='red', fill=True, alpha=0.7)
    ax.add_patch(stim_circle)

    # Colormap for firing rate
    cnorm = plt.Normalize(vmin=0, vmax=50)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=cnorm)
    sm.set_array([])

    # Overlay unit templates
    for i, unit in enumerate(units):
        x_offset = [-80, -40, 0][i % 3]
        wvf = unit['template']
        normalized_fr = np.clip(unit['fr'] / 50.0, 0, 1)
        c = cm.viridis(normalized_fr)
        wvf_x = np.arange(len(wvf))
        ax.plot(x_offset + wvf_x, unit['y'] + wvf, c=c, linewidth=0.8)

    # Scale bars
    ax.plot([130, 100], [-100, -100], 'k', linewidth=1)  # 1 ms = 30 samples
    ax.text(115, -130, '1 ms', ha='center', fontsize=5)
    ax.plot([130, 130], [-100, 0], 'k', linewidth=1)  # 100 uV
    ax.text(145, -50, '100 \u00b5V', ha='left', fontsize=5, rotation=90, va='center')

    ax.set_ylim(-200, 1600)
    ax.set_xlim(-200, 200)
    ax.set_aspect(1)
    ax.axis('off')

    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=-0.05, shrink=0.5)
    cbar.set_label('Stim-Evoked FR (Hz)', fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'probe_with_units.svg', format='svg', bbox_inches='tight')
    plt.close()
    print(f'Saved {FIG_DIR / "probe_with_units.svg"}')


if __name__ == '__main__':
    main()
