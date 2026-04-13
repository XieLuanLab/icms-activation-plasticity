"""Global configuration for figure generation scripts.

All paths are relative to the repository root. No hardcoded absolute paths.
"""
from pathlib import Path

# Repository root: two levels up from python/utils/
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / 'data'
OUTPUT_DIR = REPO_ROOT / 'output'

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Animal IDs
EXPERIMENTAL_ANIMALS = ['ICMS92', 'ICMS93', 'ICMS98', 'ICMS100', 'ICMS101']
ALL_ANIMALS = ['ICMS83'] + EXPERIMENTAL_ANIMALS  # ICMS83 + experimental

# Display names for manuscript
MOUSE_NAMES = {
    'ICMS92': 'Mouse 1',
    'ICMS93': 'Mouse 2',
    'ICMS98': 'Mouse 3',
    'ICMS100': 'Mouse 4',
    'ICMS101': 'Mouse 5',
    'ICMS83': 'Mouse 6',
}

# ICMS83 sessions
ICMS83_SESSIONS = [
    '20-Jul-2023', '24-Jul-2023', '26-Jul-2023', '28-Jul-2023',
    '31-Jul-2023', '02-Aug-2023', '04-Aug-2023', '07-Aug-2023', '10-Aug-2023',
]

# Analysis parameters
EARLY_WEEKS = [0, 1]
LATE_WEEKS = [2, 3, 4]
Z_CLIP = 15
TMAX_CLIP = 500
PSPIKE_CLIP = 0.15
CURRENTS = [4, 5, 6]

# Data file paths
RAW_DF_PATH = DATA_DIR / 'raw_df_700ms.pkl'
RAW_DF_120MS_PATH = DATA_DIR / 'raw_df_sym120ms.pkl'
CELL_TYPE_PKL = DATA_DIR / 'cell_type.pkl'
CONTROL_CSV = DATA_DIR / 'control_results.csv'
DF_MERGED_PATH = DATA_DIR / 'df_merged.pkl'
POP_COUPLING_DIR = DATA_DIR / 'pop_coupling'
TRACK_DIR = DATA_DIR / 'unit_tracking'
TRACKED_TEMPLATES_PATH = DATA_DIR / 'tracked_templates.pkl'

# Figure S-Ephys data
STIM_NONSTIM_PATH = DATA_DIR / 'fig4' / 'stim_nonstim_waveforms.npz'
FILTERED_TRACES_PATH = DATA_DIR / 'filtered_traces.npz'
ICMS83_THRESHOLDS_CSV = DATA_DIR / 'ICMS83_thresholds.csv'

# Control animals
CONTROL_ANIMALS = ['ICMS45', 'ICMS48', 'ICMS54', 'ICMS56']  # ephys modulation
CONTROL_PC_ANIMALS = ['ICMS43', 'ICMS45', 'ICMS48', 'ICMS54', 'ICMS56']  # pop coupling
CONTROL_DATA_DIR = DATA_DIR / 'control'

# Control animal sessions (for relative week computation in pop coupling)
# Date format in pop coupling pkls is '%d-%b-%Y'
CONTROL_SESSIONS = {
    'ICMS43': ['21-Jan-2022', '05-Feb-2022', '12-Feb-2022'],
    'ICMS45': ['05-Feb-2022', '12-Feb-2022', '26-Feb-2022', '04-Mar-2022'],
    'ICMS48': ['08-Apr-2022', '15-Apr-2022', '21-Apr-2022', '25-Apr-2022'],
    'ICMS54': ['20-May-2022', '01-Jun-2022'],
    'ICMS56': ['21-Jun-2022', '30-Jun-2022', '06-Jul-2022', '11-Jul-2022',
               '15-Jul-2022', '23-Jul-2022', '29-Jul-2022'],
}

# Pop coupling animal lists
ANIMALS_PC = EXPERIMENTAL_ANIMALS  # for behavioral pop coupling
ANIMALS_PC_ALL = ALL_ANIMALS       # including ICMS83

# PSTH and cluster cache paths (Figure 6)
PSTH_NPZ_PATH = DATA_DIR / 'hit_miss_psth.npz'
CLUSTERS_SEM_PATH = DATA_DIR / 'clusters_sem.json'
CLUSTERS_ZSCORE_PATH = DATA_DIR / 'clusters_zscore.json'

# Session lists per animal (for unit tracking / relative week computation)
ANIMALS_SESSIONS = {
    'ICMS83': ['20-Jul-2023', '24-Jul-2023', '26-Jul-2023', '28-Jul-2023',
               '31-Jul-2023', '02-Aug-2023', '04-Aug-2023', '07-Aug-2023',
               '10-Aug-2023'],
    'ICMS92': ['30-Aug-2023', '01-Sep-2023', '06-Sep-2023', '08-Sep-2023',
               '12-Sep-2023', '14-Sep-2023', '19-Sep-2023', '21-Sep-2023',
               '25-Sep-2023', '27-Sep-2023'],
    'ICMS93': ['30-Aug-2023', '06-Sep-2023', '12-Sep-2023', '14-Sep-2023',
               '20-Sep-2023', '22-Sep-2023', '26-Sep-2023', '29-Sep-2023',
               '04-Oct-2023', '06-Oct-2023'],
    'ICMS98': ['20-Oct-2023', '24-Oct-2023', '26-Oct-2023', '31-Oct-2023',
               '02-Nov-2023', '07-Nov-2023', '14-Nov-2023', '17-Nov-2023',
               '20-Nov-2023', '22-Nov-2023'],
    'ICMS100': ['26-Oct-2023', '02-Nov-2023', '03-Nov-2023', '09-Nov-2023',
                '16-Nov-2023', '21-Nov-2023', '05-Dec-2023'],
    'ICMS101': ['27-Oct-2023', '01-Nov-2023', '03-Nov-2023', '16-Nov-2023',
                '21-Nov-2023', '22-Nov-2023', '29-Nov-2023', '01-Dec-2023',
                '05-Dec-2023'],
}
