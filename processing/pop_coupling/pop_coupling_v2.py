from statsmodels.stats.multitest import multipletests
from scipy.stats import sem, mannwhitneyu
from scipy.stats import mannwhitneyu
from pop_coupling.shared_plotting import apply_global_style, PALETTE
import seaborn as sns
from collections import defaultdict
import traceback
import pop_coupling.fig_utils as fig_utils
import pandas as pd
import os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import pickle
import batch_process.util.file_util as file_util
from scipy.ndimage import gaussian_filter1d
import util.load_data as load_data


def compute_stpr_and_pc(pop_raster, neuron_id, bin_size=0.001, win=0.50, gauss_sigma=0.01):
    """
    pop_raster: dict[id -> np.array of spike times in seconds]
    neuron_id: id of neuron i
    Returns:
        lags (s), stpr (array length 2*win_bins+1), pc (float)
    """
    T = max((spikes.max() if len(spikes) else 0.0)
            for spikes in pop_raster.values())
    edges = np.arange(0.0, T + bin_size, bin_size)
    nbins = edges.size - 1

    pop_counts = np.zeros(nbins, dtype=float)
    for j, spikes in pop_raster.items():
        if j == neuron_id or len(spikes) == 0:
            continue
        # histogram into 1 ms bins
        pop_counts += np.histogram(spikes, bins=edges)[0]

    pop_rate = pop_counts / bin_size
    sigma_bins = gauss_sigma / bin_size  # convert seconds → bins
    pop_rate_smooth = gaussian_filter1d(
        pop_rate, sigma=sigma_bins, mode='nearest')
    # pop_rate_smooth = pop_rate_smooth - np.mean(pop_rate_smooth, axis=0)

    spikes_i = pop_raster[neuron_id]
    stpr_windows = []
    win_bins = int(round(win / bin_size))
    # convert spike times to bin indices
    spike_bins = np.floor(spikes_i / bin_size).astype(int)

    stpr_windows = []
    for sb in spike_bins:
        a = sb - win_bins
        b = sb + win_bins
        if a < 0 or (b + 1) > nbins:
            continue
        stpr_windows.append(pop_rate_smooth[a:b+1])

    # stpr = np.vstack(stpr_windows).mean(axis=0)  # average across spikes
    # stpr -= stpr.mean()  # subtract mean
    # lags = np.arange(-win_bins, win_bins+1) * bin_size  # seconds
    # pc = stpr[win_bins]  # population coupling value or stPR at zero lag

    # average population rate across entire recording
    global_mean = pop_rate_smooth.mean()
    stpr = np.vstack(stpr_windows).mean(axis=0)
    stpr -= global_mean
    pc = stpr[win_bins]

    return lags, stpr, pc


def compute_stpr_and_pc(pop_raster, neuron_id, bin_size=0.001, win=0.50, gauss_sigma=0.01,
                        n_shuffles=1, rng=None):
    """
    Returns: lags (s), stpr (2*win_bins+1), pc (float), n_spikes_used (int), pc_norm (float or None)
    pc_norm is pc normalized by median(|pc_shuffled|) if n_shuffles>0, else None.
    """
    # --- bin edges ---
    T = max((spikes.max() if len(spikes) else 0.0)
            for spikes in pop_raster.values())
    if T <= 0:
        return None, None, np.nan, 0, None
    edges = np.arange(0.0, T + bin_size, bin_size)
    nbins = edges.size - 1

    # --- population rate (exclude target) ---
    pop_counts = np.zeros(nbins, dtype=float)
    for j, spikes in pop_raster.items():
        if j == neuron_id or len(spikes) == 0:
            continue
        pop_counts += np.histogram(spikes, bins=edges)[0]
    pop_rate = pop_counts / bin_size
    sigma_bins = max(1e-9, gauss_sigma / bin_size)  # seconds -> bins
    pop_rate_smooth = gaussian_filter1d(
        pop_rate, sigma=sigma_bins, mode='nearest')

    # mean-center (equivalent to sum_j(f_j - mu_j))
    pop_mc = pop_rate_smooth - np.mean(pop_rate_smooth)

    # --- spike-triggered average around neuron i's spikes ---
    spikes_i = np.asarray(pop_raster.get(neuron_id, []), float)
    if spikes_i.size == 0:
        win_bins = int(round(win / bin_size))
        lags = np.arange(-win_bins, win_bins+1) * bin_size
        return lags, np.full(lags.size, np.nan), np.nan, 0, None

    win_bins = int(round(win / bin_size))
    lags = np.arange(-win_bins, win_bins+1) * bin_size
    spike_bins = np.floor(spikes_i / bin_size).astype(int)

    windows = []
    for sb in spike_bins:
        a, b = sb - win_bins, sb + win_bins
        if a < 0 or b >= nbins:
            continue
        windows.append(pop_mc[a:b+1])

    if len(windows) == 0:
        return lags, np.full(lags.size, np.nan), np.nan, 0, None

    # average over spikes => /||f_i||
    stpr = np.vstack(windows).mean(axis=0)
    pc = float(stpr[win_bins])                 # value at zero lag
    n_used = len(windows)

    # --- optional shuffle normalization ---
    pc_norm = None
    if n_shuffles and n_used > 0:
        if rng is None:
            rng = np.random.default_rng()
        # circularly shift the mean-centered population rate by random offsets
        # (preserves spectrum & rate, breaks alignment with spikes)
        shuf_vals = []
        for _ in range(int(n_shuffles)):
            # avoid trivial near-zero shifts
            shift = rng.integers(low=win_bins+1, high=nbins-win_bins-1)
            pop_mc_sh = np.roll(pop_mc, shift)
            # same windows, but on shuffled trace
            vals = [pop_mc_sh[sb - win_bins: sb + win_bins + 1][win_bins]
                    for sb in spike_bins
                    if (sb - win_bins) >= 0 and (sb + win_bins) < nbins]
            if len(vals):
                shuf_vals.append(np.mean(vals))
        if len(shuf_vals):
            den = np.median(np.abs(shuf_vals))
            pc_norm = np.nan if den <= 1e-12 else float(pc / den)

    return lags, stpr, pc, n_used, pc_norm


def lookup_pl_or_npl_units(df, session_path, ch, cur, is_pl=True):

    animal_string = file_util.get_animal_id(session_path)
    date_string = file_util.get_date_str(session_path)

    mask = (
        (df['animal_id'] == animal_string) &
        (df['session'] == date_string) &
        (df['stim_channel'] == ch) &
        (df['stim_current'] == cur) &
        (df['is_pulse_locked'] == is_pl)
    )

    return df.loc[mask, 'unit_id'].tolist()


def get_stim_condition_segments_dict(trial_df):
    currents, channels = get_all_valid_stim_conditions(trial_df)
    out = {}
    for current in currents:
        for channel in channels:
            sub_df = trial_df[(trial_df['current'] == current) &
                              (trial_df['channel'] == channel)]
            if sub_df.empty:
                continue
            segments = []
            for _, trial in sub_df.iterrows():
                ts = trial.stim_timestamps
                if ts is None:
                    continue
                # handle scalar (int/float) vs array/list
                if np.isscalar(ts):
                    if isinstance(ts, float) and np.isnan(ts):
                        continue
                    start_ts = end_ts = float(ts)
                    segments.append((start_ts, end_ts))
                else:
                    ts_arr = np.asarray(ts)
                    if ts_arr.size == 0 or np.all(~np.isfinite(ts_arr)):
                        continue
                    segments.append((float(ts_arr[0]), float(ts_arr[-1])))
            if segments:
                out[(current, channel)] = segments
    return out


def segment_spike_train(spike_train, segments):

    segmented_spikes = []

    # Use np.searchsorted to efficiently find the spikes within each segment
    for start, end in segments:
        # Find the indices in spike_train where spikes are within the segment (start, end)
        start_idx = np.searchsorted(spike_train, start, side='right')
        end_idx = np.searchsorted(spike_train, end, side='left')
        segmented_spikes.append(spike_train[start_idx:end_idx])

    return segmented_spikes


def subdict(d, keys):
    missing = [k for k in keys if k not in d]
    if missing:
        print(f"[warn] missing unit_ids: {missing}")
    return {k: d[k] for k in keys if k in d}


def get_all_valid_stim_conditions(df):
    currents = df['current'].unique()
    channels = df['channel'].unique()
    currents = currents[(currents >= 4) & (currents < 7)]
    channels = channels[channels > 0]
    return currents, channels


def get_session_path(animal_id, session):
    root = "C:\\data"
    if animal_id in ['ICMS92', 'ICMS93']:
        subfolder = 'Behavior'
    else:
        subfolder = ''
    return os.path.join(root, animal_id, subfolder, session)


def concat_segments_rebased(spike_train, segments):
    """Return 1D array of spike times where each segment is re-based to start at 0
    and segments are concatenated end-to-end (no gaps)."""
    out = []
    offset = 0.0
    for a, b in segments:
        # spikes within [a,b)
        seg = spike_train[(spike_train >= a) & (spike_train < b)]
        if seg.size:
            out.append((seg - a) + offset)
        offset += (b - a)
    return np.concatenate(out) if out else np.array([], dtype=float)


def get_pop_coupling_session_pl_npl(session_path, session_data, df_modulated_only, df_pl, df_not_pl):
    '''
    session_data: value from spike data dictionary with session key 
    '''
    unit_spike_train_dict = session_data['unit_spike_train_dict']
    trial_df = load_data.get_dataframe(session_path, make_folder=False)
    condition_stim_trains_dict = get_stim_condition_segments_dict(trial_df)

    condition_results = {}
    FS = 30000.0

    for cur, ch in condition_stim_trains_dict.keys():
        stim_segments = condition_stim_trains_dict[(cur, ch)]
        stim_segments_s = [(a/FS, b/FS) for (a, b) in stim_segments]

        pl_unit_ids = lookup_pl_or_npl_units(
            df_modulated_only, session_path, ch, cur, is_pl=True)
        npl_unit_ids = lookup_pl_or_npl_units(
            df_modulated_only, session_path, ch, cur, is_pl=False)
        mod_ids = sorted(set(pl_unit_ids + npl_unit_ids))

        pl_unit_spike_train_dict = subdict(unit_spike_train_dict, pl_unit_ids)
        npl_unit_spike_train_dict = subdict(
            unit_spike_train_dict, npl_unit_ids)
        mod_unit_spike_train_dict = subdict(unit_spike_train_dict, mod_ids)

        modulated_pop_raster, pl_pop_raster, npl_pop_raster = {}, {}, {}
        for unit_id, spike_train in mod_unit_spike_train_dict.items():
            # stim_segmented_spikes = segment_spike_train(
            #     spike_train, stim_segments)
            # spikes_s = np.concatenate(stim_segmented_spikes) / 30000
            spikes_s = concat_segments_rebased(
                spike_train.astype(float)/FS, stim_segments_s)
            modulated_pop_raster[unit_id] = spikes_s
            if unit_id in pl_unit_ids:
                pl_pop_raster[unit_id] = spikes_s
            elif unit_id in npl_unit_ids:
                npl_pop_raster[unit_id] = spikes_s

        pl_stpr_dict, pl_pc_dict, pl_pc_norm_dict = {}, {}, {}
        npl_stpr_dict, npl_pc_dict, npl_pc_norm_dict = {}, {}, {}

        rng = np.random.default_rng(0)
        for unit_id in modulated_pop_raster.keys():
            lags, stpr, pc, n_used, pc_norm = compute_stpr_and_pc(
                modulated_pop_raster, unit_id, win=0.1, n_shuffles=200, rng=rng
            )

            MIN_SPIKES = 50  # choose based on your data stability
            if n_used < MIN_SPIKES:
                stpr = None
                pc = np.nan
                pc_norm = np.nan

            if unit_id in pl_unit_ids:
                pl_stpr_dict[unit_id] = stpr
                pl_pc_dict[unit_id] = pc
                pl_pc_norm_dict[unit_id] = pc_norm
            elif unit_id in npl_unit_ids:
                npl_stpr_dict[unit_id] = stpr
                npl_pc_dict[unit_id] = pc
                npl_pc_norm_dict[unit_id] = pc_norm

        condition_results[(cur, ch)] = {
            "pl_stpr_dict": pl_stpr_dict,
            "pl_pc_dict": pl_pc_dict,               # Hz
            "pl_pc_norm_dict": pl_pc_norm_dict,     # unitless
            "npl_stpr_dict": npl_stpr_dict,
            "npl_pc_dict": npl_pc_dict,             # Hz
            "npl_pc_norm_dict": npl_pc_norm_dict,   # unitless
        }

    return condition_results


# %%

raw_df = pd.read_pickle('raw_df_700ms.pkl')
df_modulated_only = fig_utils.get_modulated_only_responses(raw_df)
df_pl = fig_utils.get_filtered_responses(raw_df)
df_not_pl = fig_utils.get_filtered_responses(raw_df, use_filt_B=True)

animal_ids = ['ICMS92', 'ICMS93', 'ICMS98', 'ICMS100', 'ICMS101']
# animal_ids = ['ICMS93', 'ICMS98', 'ICMS100', 'ICMS101']

for animal_id in animal_ids:
    pickle_file_name = f'{animal_id}_spike_data.pkl'
    with open(pickle_file_name, 'rb') as f:
        data = pickle.load(f)
    print(f'Processing {animal_id}...')
    session_names = list(data.keys())
    session_dict = OrderedDict()
    # session_names = ['20-Oct-2023']
    for session in session_names:
        print(f'\tProcessing {session}...')
        try:
            session_path = get_session_path(animal_id, session)
            condition_results = get_pop_coupling_session_pl_npl(
                session_path, data[session], df_modulated_only, df_pl, df_not_pl)
            session_dict[session] = condition_results

        except Exception as e:
            print(f"[WARN] {session}: {e}")
            traceback.print_exc()

    out_pkl = f"pop_coupling/{animal_id}_pop_coupling_v5.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(session_dict, f)


# %%


def _merge_unit_stprs(dicts_list):
    """Merge list of {unit_id -> stpr(np.ndarray)} into a single dict."""
    merged = {}
    for d in dicts_list:
        # unit_ids assumed unique across channels; if not, last wins
        merged.update(d)
    return merged


def plot_session_by_current(session_data, lags, animal_id, session_name,
                            ymax=50, show_mean=False):
    # Collect available currents and map current -> list of channels
    currents = sorted({cur for (cur, ch) in session_data.keys()})
    cur_to_channels = defaultdict(list)
    for (cur, ch) in session_data.keys():
        cur_to_channels[cur].append(ch)

    n_rows = len(currents)
    fig, axes = plt.subplots(n_rows, 2, figsize=(
        4, 2*n_rows), sharex=True, sharey=True)
    if n_rows == 1:
        axes = np.atleast_2d(axes)  # ensure 2D indexing

    for r, cur in enumerate(currents):
        # Gather PL/NPL dicts across all channels for this current
        pl_dicts = [session_data[(cur, ch)]["pl_stpr_dict"]
                    for ch in cur_to_channels[cur] if "pl_stpr_dict" in session_data[(cur, ch)]]
        npl_dicts = [session_data[(cur, ch)]["npl_stpr_dict"]
                     for ch in cur_to_channels[cur] if "npl_stpr_dict" in session_data[(cur, ch)]]

        pl_all = _merge_unit_stprs(pl_dicts)
        npl_all = _merge_unit_stprs(npl_dicts)

        # Color palettes: deeper sequential shades, avoid very light end
        pl_colors = plt.cm.Blues(np.linspace(0.35, 0.9, max(1, len(pl_all))))
        npl_colors = plt.cm.Oranges(
            np.linspace(0.35, 0.9, max(1, len(npl_all))))

        # --- PL subplot (col 0)
        ax = axes[r, 0]
        for color, (uid, stpr) in zip(pl_colors, pl_all.items()):
            ax.plot(lags, stpr, color=color, lw=1)
        if show_mean and len(pl_all) > 0:
            try:
                M = np.vstack([v for v in pl_all.values()
                              if v is not None and np.all(np.isfinite(v))])
                ax.plot(lags, M.mean(axis=0),
                        color='navy', lw=2.5, label='Mean')
            except ValueError:
                pass
        ax.set_title(f"PL ({cur} uA)")
        ax.set_ylabel("Firing rate (Hz)")
        ax.set_ylim(top=ymax)

        # --- NPL subplot (col 1)
        ax = axes[r, 1]
        for color, (uid, stpr) in zip(npl_colors, npl_all.items()):
            ax.plot(lags, stpr, color=color, lw=1)
        if show_mean and len(npl_all) > 0:
            try:
                M = np.vstack([v for v in npl_all.values()
                              if v is not None and np.all(np.isfinite(v))])
                ax.plot(lags, M.mean(axis=0),
                        color='darkorange', lw=2.5, label='Mean')
            except ValueError:
                pass
        ax.set_title(f"NPL ({cur} uA)")
        ax.set_ylim(top=ymax)

        if r == n_rows - 1:
            axes[r, 0].set_xlabel("Lag (s)")
            axes[r, 1].set_xlabel("Lag (s)")

    fig.suptitle(f"{animal_id}: {session_name}")
    plt.tight_layout()
    plt.show()


animal_id = 'ICMS92'
pickle_file_name = f"pop_coupling/{animal_id}_pop_coupling_v5.pkl"
with open(pickle_file_name, 'rb') as f:
    pc_animal_data = pickle.load(f)
session_names = list(pc_animal_data.keys())

bin_size = 0.001
win = 0.1
win_bins = int(round(win / bin_size))
lags = np.arange(-win_bins, win_bins+1) * bin_size

for session in session_names:
    session_data = pc_animal_data[session]
    plot_session_by_current(
        session_data, lags, animal_id='ICMS92', session_name=session)


# %%
def _collect_session_cur_vals(session_data, cur):
    """Return two flat lists of scalars for this current: (pl_vals, npl_vals)."""
    pl_vals, npl_vals = [], []
    chs = [ch for (c, ch) in session_data.keys() if c == cur]
    for ch in chs:
        d = session_data[(cur, ch)]
        if "pl_pc_dict" in d:
            pl_vals.extend([float(v) for v in d["pl_pc_dict"].values()
                            if v is not None and np.isfinite(v)])
        if "npl_pc_dict" in d:
            npl_vals.extend([float(v) for v in d["npl_pc_dict"].values()
                             if v is not None and np.isfinite(v)])
    return pl_vals, npl_vals


def plot_pc_points_across_sessions(pc_animal_data, session_names,
                                   jitter=0.15, s=12,
                                   pl_color="#4C72B0", npl_color="#DD8452",
                                   alpha=0.8, show_median=True):
    """
    Scatter all individual PC scalars per session/current.
    Overlay group medians as lines.
    - PL: deep blue, NPL: deep orange
    """
    currents_all = sorted({
        cur for session in session_names
        for (cur, ch) in pc_animal_data[session].keys()
    })
    n_cur = len(currents_all)
    if n_cur == 0:
        print("No currents found.")
        return

    fig, axes = plt.subplots(n_cur, 1,
                             figsize=(max(1, 0.7*len(session_names)), 2*n_cur),
                             sharex=True)
    if n_cur == 1:
        axes = np.atleast_1d(axes)

    rng = np.random.default_rng(42)
    x = np.arange(len(session_names))

    for i, cur in enumerate(currents_all):
        ax = axes[i]
        pl_medians, npl_medians = [], []

        for j, session in enumerate(session_names):
            sd = pc_animal_data[session]
            pl_vals, npl_vals = _collect_session_cur_vals(sd, cur)

            # Scatter PL
            if pl_vals:
                xj = x[j] + rng.uniform(-jitter, jitter,
                                        size=len(pl_vals)) - 0.2
                ax.scatter(xj, pl_vals, s=s, color=pl_color,
                           alpha=alpha, label="PL" if (i == 0 and j == 0) else None)
                pl_medians.append(np.median(pl_vals))
            else:
                pl_medians.append(np.nan)

            # Scatter NPL
            if npl_vals:
                xj = x[j] + rng.uniform(-jitter, jitter,
                                        size=len(npl_vals)) + 0.2
                ax.scatter(xj, npl_vals, s=s, color=npl_color,
                           alpha=alpha, label="NPL" if (i == 0 and j == 0) else None)
                npl_medians.append(np.median(npl_vals))
            else:
                npl_medians.append(np.nan)

        # Overlay median lines
        if show_median:
            ax.plot(x - 0.2, pl_medians, color=pl_color,
                    lw=2, marker="o", ms=1)
            ax.plot(x + 0.2, npl_medians, color=npl_color,
                    lw=2, marker="o", ms=1)

        ax.set_ylabel(f"{cur} µA")
        ax.grid(True, axis='y', alpha=0.25)

    axes[-1].set_xticks(x)
    session_names_clean = [name.replace("-2023", "") for name in session_names]
    axes[-1].set_xticklabels(session_names_clean, rotation=65, ha='right')
    fig.supylabel("Population coupling")
    fig.suptitle(f"PC across Sessions {animal_id}", y=0.995)
    axes[0].legend(frameon=False, loc="best")
    plt.tight_layout()
    plt.show()


plot_pc_points_across_sessions(pc_animal_data, session_names)


# %% Aggregate

'''
Compare early vs late weeks for each current? 
Group currents first
Then group by time (need early and late weeks for each animal)
Group by pulse-locking status
'''
animal_ids = ['ICMS92', 'ICMS93', 'ICM98', 'ICMS100', 'ICMS101']

def is_early_or_late_session(df_modulated_only, animal_id, session):
    sel = df_modulated_only[(df_modulated_only['session'] == session) &
                            (df_modulated_only['animal_id'] == animal_id)]
    vals = sel['rel_week'].dropna().unique()
    if len(vals) == 0:
        return 'nan'
    rel_week = float(vals[0])  # take the (presumably) single value
    return 'early' if rel_week < 2 else 'late'


pl_early_dict = {4: [], 5: [], 6: []}  # current key, pop coupling value
pl_late_dict = {4: [], 5: [], 6: []}
npl_early_dict = {4: [], 5: [], 6: []}
npl_late_dict = {4: [], 5: [], 6: []}

for animal_id in animal_ids:
    pickle_file_name = f"pop_coupling/{animal_id}_pop_coupling_v4.pkl"
    with open(pickle_file_name, 'rb') as f:
        pc_animal_data = pickle.load(f)
    session_names = list(pc_animal_data.keys())
    for session in session_names:
        time_period = is_early_or_late_session(
            df_modulated_only, animal_id, session)
        if time_period == 'nan':
            continue
        session_data = pc_animal_data[session]
        conditions = session_data.keys()
        for cur, ch in conditions:
            condition_data = session_data[(cur, ch)]

            pl_pcs = condition_data['pl_pc_norm_dict'].values()
            npl_pcs = condition_data['npl_pc_norm_dict'].values()

            # pl_pcs = condition_data['pl_pc_dict'].values()
            # npl_pcs = condition_data['npl_pc_dict'].values()

            if time_period == 'early':
                pl_early_dict[cur].extend(pl_pcs)
                npl_early_dict[cur].extend(npl_pcs)
            elif time_period == 'late':
                pl_late_dict[cur].extend(pl_pcs)
                npl_late_dict[cur].extend(npl_pcs)
            else:
                print('error!')

# %% Plot

# Helper: p-value -> stars


def p_to_stars(p):
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def add_sig_bracket(ax, x1, x2, y, h=0.9, text="ns", lw=1):
    """Draw a bracket from x1 to x2 at height y, with little vertical caps of height h."""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color="k", lw=lw, clip_on=False)
    if text == "ns":
        text_offset = 0.2
    else:
        text_offset = -0.1
    ax.text((x1+x2)/2, y+text_offset, text,
            ha="center", va="bottom", fontsize=9)


def clean_NaNs(pl_early_dict, pl_late_dict, npl_early_dict, npl_late_dict, cur):
    pl_e = np.asarray(pl_early_dict[cur], float)
    pl_l = np.asarray(pl_late_dict[cur],  float)
    npl_e = np.asarray(npl_early_dict[cur], float)
    npl_l = np.asarray(npl_late_dict[cur],  float)
    pl_e, pl_l = pl_e[np.isfinite(pl_e)], pl_l[np.isfinite(pl_l)]
    npl_e, npl_l = npl_e[np.isfinite(npl_e)], npl_l[np.isfinite(npl_l)]
    return pl_e, pl_l, npl_e, npl_l


tests = []  # store all tests

for cur in range(4, 7):
    xpos = cur - 4

    # Pull arrays and clean NaN
    pl_e, pl_l, npl_e, npl_l = clean_NaNs(
        pl_early_dict, pl_late_dict, npl_early_dict, npl_late_dict, cur)

    # y baseline above this current's data
    all_vals_cur = np.concatenate([pl_e, pl_l, npl_e, npl_l]) if (
        pl_e.size or pl_l.size or npl_e.size or npl_l.size) else np.array([0.0])
    y_base = np.nanmax(all_vals_cur) + 0.05 * max(1.0, np.nanstd(all_vals_cur))

    # --- Early vs Late within each group
    if pl_e.size and pl_l.size:
        stat, p = mannwhitneyu(pl_e, pl_l, alternative="two-sided")
        tests.append(dict(cur=cur, label="PL early vs late",
                          x1=xpos-0.2, x2=xpos-0.1, y=y_base, p=p, prio=0))
    if npl_e.size and npl_l.size:
        stat, p = mannwhitneyu(npl_e, npl_l, alternative="two-sided")
        tests.append(dict(cur=cur, label="NPL early vs late",
                          x1=xpos+0.1, x2=xpos+0.2, y=y_base*1.05, p=p, prio=1))

    # --- PL vs NPL within each time period
    if pl_e.size and npl_e.size:
        stat, p = mannwhitneyu(pl_e, npl_e, alternative="two-sided")
        tests.append(dict(cur=cur, label="Early PL vs NPL",
                          x1=xpos-0.2, x2=xpos+0.1, y=y_base*1.3, p=p, prio=2))
    if pl_l.size and npl_l.size:
        stat, p = mannwhitneyu(pl_l, npl_l, alternative="two-sided")
        tests.append(dict(cur=cur, label="Late PL vs NPL",
                          x1=xpos-0.1, x2=xpos+0.2, y=y_base*1.55, p=p, prio=3))

# --- Apply multiple-comparison correction (FDR/BH)
pvals = np.array([t["p"] for t in tests], float)
try:
    from statsmodels.stats.multitest import multipletests
    _, pvals_corr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
except Exception:
    pvals_corr = pvals.copy()

# --- Global baseline (keep your existing GLOBAL_Y0 computation) ---
SPACING = 3
MARGIN = 8.0   # positive; above data
GLOBAL_Y0 = 14

tests_by_cur = defaultdict(list)
for t in tests:
    tests_by_cur[t["cur"]].append(t)

for cur, items in tests_by_cur.items():
    items.sort(key=lambda d: d.get("prio", 99))
    for i, t in enumerate(items):
        t["y"] = GLOBAL_Y0 + i * SPACING

apply_global_style()
pl_c = PALETTE[0]
npl_c = PALETTE[1]

fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))

for cur in range(4, 7):
    xpos = cur - 4  # base position: 0, 1, 2

    pl_e, pl_l, npl_e, npl_l = clean_NaNs(
        pl_early_dict, pl_late_dict, npl_early_dict, npl_late_dict, cur)

    pl_early_median = np.median(pl_e)
    pl_late_median = np.median(pl_l)
    npl_early_median = np.median(npl_e)
    npl_late_median = np.median(npl_l)

    pl_early_iqr = np.nanpercentile(
        pl_e, 75) - np.nanpercentile(pl_e, 25)
    pl_late_iqr = np.nanpercentile(
        pl_l, 75) - np.nanpercentile(pl_l, 25)
    npl_early_iqr = np.nanpercentile(
        npl_e, 75) - np.nanpercentile(npl_e, 25)
    npl_late_iqr = np.nanpercentile(
        npl_l, 75) - np.nanpercentile(npl_l, 25)

    ax.errorbar(xpos - 0.2, pl_early_median, yerr=pl_early_iqr,
                fmt='o', color=pl_c, alpha=0.6, capsize=3,
                label='PL early' if cur == 4 else "")
    ax.errorbar(xpos - 0.1, pl_late_median, yerr=pl_late_iqr,
                fmt='o', color=pl_c, capsize=3,
                label='PL late' if cur == 4 else "")

    ax.errorbar(xpos + 0.1, npl_early_median, yerr=npl_early_iqr,
                fmt='o', color=npl_c, alpha=0.6, capsize=3,
                label='NPL early' if cur == 4 else "")
    ax.errorbar(xpos + 0.2, npl_late_median, yerr=npl_late_iqr,
                fmt='o', color=npl_c, capsize=3,
                label='NPL late' if cur == 4 else "")


ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['4 µA', '5 µA', '6 µA'])
ax.set_xlim(-0.5, 2.5)
# ax.set_ylim([-8, 50])


# --- Draw the brackets
for t, p_adj in zip(tests, pvals_corr):
    stars = p_to_stars(p_adj)
    add_sig_bracket(ax, t["x1"], t["x2"], y=t["y"] * 0.45, h=0.1, text=stars)

ax.legend(
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False
)
plt.ylabel('Population coupling')
plt.tight_layout()
plt.show()

# %%

# --- helpers ---


def concat_all(d, currents=(4, 5, 6)):
    out = []
    for c in currents:
        out.extend(list(d.get(c, [])))
    a = np.asarray(out, float)
    return a[np.isfinite(a)]


def agg(a, mode='mean_sem'):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan, 0.0, 0
    if mode == 'mean_sem':
        return float(np.nanmean(a)), float(sem(a, nan_policy='omit')), a.size
    # median ± IQR
    med = float(np.nanmedian(a))
    iqr = float(np.nanpercentile(a, 75) - np.nanpercentile(a, 25))
    return med, iqr, a.size


def p_to_stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def add_sig_bracket(ax, x1, x2, y, h=0.06, text="ns", lw=1):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color="k", lw=lw, clip_on=False)
    ax.text((x1+x2)/2, y+h, text, ha="center", va="bottom", fontsize=9)

# --- build pooled arrays (use pc_norm dicts ideally) ---
# Expecting you already have these dicts keyed by current: 4,5,6
# pl_early_dict, pl_late_dict, npl_early_dict, npl_late_dict


pl_e_all = concat_all(pl_early_dict)
pl_l_all = concat_all(pl_late_dict)
npl_e_all = concat_all(npl_early_dict)
npl_l_all = concat_all(npl_late_dict)

# choose summary mode: 'mean_sem' or 'median_iqr'
MODE = 'median_iqr'

# --- plot ---
# pl_c = "#1f77b4"
# npl_c = "#ff7f0e"

fig, ax = plt.subplots(figsize=(2.5, 2.5))
x = np.array([0, 1, 2, 3])
labels = ["PL early", "PL late", "NPL early", "NPL late"]
cols = [pl_c,       pl_c,       npl_c,       npl_c]

groups = [pl_e_all, pl_l_all, npl_e_all, npl_l_all]
stats = [agg(g, MODE) for g in groups]  # (center, err, n)

# scatter means with errors
for xi, (g, (m, e, n)), c, lab in zip(x, zip(groups, stats), cols, labels):
    if np.isfinite(m):
        ax.errorbar(xi, m, yerr=e, fmt='o', color=c, capsize=3,
                    label=f"{lab} (n={n})")

# --- significance tests (Mann–Whitney; FDR across all 4)
tests = []


def add_test(name, x1, x2, a, b):
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size and b.size:
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        tests.append(dict(name=name, x1=x1, x2=x2, p=p))


add_test("PL early vs late",  0, 1, pl_e_all,  pl_l_all)
add_test("NPL early vs late", 2, 3, npl_e_all, npl_l_all)
add_test("Early PL vs NPL",  0, 2, pl_e_all,  npl_e_all)
add_test("Late PL vs NPL",  1, 3, pl_l_all,  npl_l_all)

pvals = np.array([t["p"] for t in tests], float)
_, pvals_corr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
for t, p_adj in zip(tests, pvals_corr):
    t["p_adj"] = float(p_adj)

# place brackets in a compact band above the data
finite_vals = np.concatenate([g[np.isfinite(g)] for g in groups if g.size])
bracket_vertical_offsets = []
if finite_vals.size:
    ymax = float(np.nanmax(finite_vals))
    ymax = 5
    pad = 0.06 * max(1.0, np.nanstd(finite_vals))
else:
    ymax, pad = 0.0, 0.2
band = 0.35
step = band / max(1, len(tests))
step = 1.8
for i, t in enumerate(tests):
    y = ymax + pad + i*step
    add_sig_bracket(ax, t["x1"], t["x2"], y=y, h=0.4,
                    text=p_to_stars(t["p_adj"]))

# cosmetics
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=15)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim([-5, 8])
ax.set_ylabel("Population coupling (unitless)" if 'norm' in str(
    pl_early_dict).lower() else "Population coupling")
# ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
plt.tight_layout()
plt.show()


# %% Example PL and NPL stPR


def plot_session_by_current(session_data, lags, animal_id, session_name, current,
                            pl_id=None, npl_id=None, ymax=60, show_mean=False,
                            gray_alpha=0.6):
    """
    Create TWO figures:
      1) PL STPR figure (all curves light gray; pl_id highlighted deep blue if present)
      2) NPL STPR figure (all curves light gray; npl_id highlighted deep orange if present)

    Parameters
    ----------
    session_data : dict
        Keys: (current_uA, channel). Values: dict with optional keys
            "pl_stpr_dict" : { unit_id -> stpr_vector }
            "npl_stpr_dict" : { unit_id -> stpr_vector }
    lags : 1D array-like
        Time lags for STPR x-axis.
    animal_id, session_name : str
        For figure titles (not otherwise used).
    current : int or float
        Which current to plot (matches the first element in session_data keys).
    pl_id, npl_id : hashable or None
        Unit IDs to highlight in their respective figures.
    ymax : float
        Y-axis upper limit for both figures.
    show_mean : bool
        If True, overlay mean curve in each figure.
    gray_alpha : float
        Alpha for background curves.
    """

    # Collect channels that have the requested current
    cur_to_channels = defaultdict(list)
    for (cur, ch) in session_data.keys():
        cur_to_channels[cur].append(ch)

    if current not in cur_to_channels:
        raise ValueError(f"No entries found for current={current} in session_data.")

    # Merge unit->stpr dicts across channels for this current
    def _merge_unit_stprs(dict_list):
        merged = {}
        for d in dict_list:
            for k, v in d.items():
                merged[k] = v
        return merged

    pl_dicts = [
        session_data[(current, ch)]["pl_stpr_dict"]
        for ch in cur_to_channels[current]
        if "pl_stpr_dict" in session_data[(current, ch)]
    ]
    npl_dicts = [
        session_data[(current, ch)]["npl_stpr_dict"]
        for ch in cur_to_channels[current]
        if "npl_stpr_dict" in session_data[(current, ch)]
    ]

    pl_all = _merge_unit_stprs(pl_dicts) if pl_dicts else {}
    npl_all = _merge_unit_stprs(npl_dicts) if npl_dicts else {}

    # ---- Figure 1: PL ----
    fig_pl, ax_pl = plt.subplots(1, 1, figsize=(3.0, 2.5))
    # background: all PL in light gray
    for uid, stpr in pl_all.items():
        if stpr is None or not np.all(np.isfinite(stpr)):
            continue
        ax_pl.plot(lags, stpr, color='0.7', alpha=gray_alpha, lw=1.0)

    # highlight selected pl_id
    if pl_id is not None and pl_id in pl_all and pl_all[pl_id] is not None:
        if np.all(np.isfinite(pl_all[pl_id])):
            ax_pl.plot(lags, pl_all[pl_id], color=plt.cm.Blues(
                0.95), lw=2.2, label=f'PL {pl_id}')

    # optional mean
    if show_mean and len(pl_all) > 0:
        try:
            M = np.vstack([v for v in pl_all.values()
                          if v is not None and np.all(np.isfinite(v))])
            if M.size > 0:
                ax_pl.plot(lags, M.mean(axis=0), color=plt.cm.Blues(
                    0.65), lw=2.0, label='Mean')
        except ValueError:
            pass

    ax_pl.set_title(f"PL STPR • {animal_id} {session_name} • {current} µA")
    ax_pl.set_ylabel("Firing rate (Hz)")
    ax_pl.set_xlabel("Lag (s)")
    ax_pl.set_ylim(top=ymax)
    if show_mean or (pl_id is not None and pl_id in pl_all):
        ax_pl.legend(frameon=False, fontsize=8)
    fig_pl.tight_layout()

    # ---- Figure 2: NPL ----
    fig_npl, ax_npl = plt.subplots(1, 1, figsize=(3.0, 2.5))
    # background: all NPL in light gray
    for uid, stpr in npl_all.items():
        if stpr is None or not np.all(np.isfinite(stpr)):
            continue
        ax_npl.plot(lags, stpr, color='0.7', alpha=gray_alpha, lw=1.0)

    # highlight selected npl_id
    if npl_id is not None and npl_id in npl_all and npl_all[npl_id] is not None:
        if np.all(np.isfinite(npl_all[npl_id])):
            ax_npl.plot(lags, npl_all[npl_id], color=plt.cm.Oranges(
                0.95), lw=2.2, label=f'NPL {npl_id}')

    # optional mean
    if show_mean and len(npl_all) > 0:
        try:
            M = np.vstack([v for v in npl_all.values()
                          if v is not None and np.all(np.isfinite(v))])
            if M.size > 0:
                ax_npl.plot(lags, M.mean(axis=0), color=plt.cm.Oranges(
                    0.7), lw=2.0, label='Mean')
        except ValueError:
            pass

    ax_npl.set_title(f"NPL STPR • {animal_id} {session_name} • {current} µA")
    ax_npl.set_ylabel("Firing rate (Hz)")
    ax_npl.set_xlabel("Lag (s)")
    ax_npl.set_ylim(top=ymax)
    if show_mean or (npl_id is not None and npl_id in npl_all):
        ax_npl.legend(frameon=False, fontsize=8)
    fig_npl.tight_layout()

    return (fig_pl, ax_pl), (fig_npl, ax_npl)


animal_id = 'ICMS92'
pickle_file_name = f"pop_coupling/{animal_id}_pop_coupling.pkl"
with open(pickle_file_name, 'rb') as f:
    pc_animal_data = pickle.load(f)
session_names = list(pc_animal_data.keys())

bin_size = 0.001
win = 0.1
win_bins = int(round(win / bin_size))
lags = np.arange(-win_bins, win_bins+1) * bin_size

session_name = [session_names[-1]]
for session in session_name:
    session_data = pc_animal_data[session]
    plot_session_by_current(
        session_data, lags, animal_id='ICMS92', session_name=session)

# %% Rasters


def plot_session_by_current(session_data, lags, animal_id, session_name, current,
                            pl_id=None, npl_id=None, ymax=50, show_mean=False,
                            gray_alpha=0.6):
    """
    Create TWO figures:
      1) PL STPR figure (all curves light gray; pl_id highlighted deep blue if present)
      2) NPL STPR figure (all curves light gray; npl_id highlighted deep orange if present)

    Parameters
    ----------
    session_data : dict
        Keys: (current_uA, channel). Values: dict with optional keys
            "pl_stpr_dict" : { unit_id -> stpr_vector }
            "npl_stpr_dict" : { unit_id -> stpr_vector }
    lags : 1D array-like
        Time lags for STPR x-axis.
    animal_id, session_name : str
        For figure titles (not otherwise used).
    current : int or float
        Which current to plot (matches the first element in session_data keys).
    pl_id, npl_id : hashable or None
        Unit IDs to highlight in their respective figures.
    ymax : float
        Y-axis upper limit for both figures.
    show_mean : bool
        If True, overlay mean curve in each figure.
    gray_alpha : float
        Alpha for background curves.
    """

    # Collect channels that have the requested current
    cur_to_channels = defaultdict(list)
    for (cur, ch) in session_data.keys():
        cur_to_channels[cur].append(ch)

    if current not in cur_to_channels:
        raise ValueError(f"No entries found for current={current} in session_data.")

    def _merge_unit_stprs(dict_list):
        merged = {}
        for d in dict_list:
            for k, v in d.items():
                merged[k] = v
        return merged

    pl_dicts = [
        session_data[(current, ch)]["pl_stpr_dict"]
        for ch in cur_to_channels[current]
        if "pl_stpr_dict" in session_data[(current, ch)]
    ]
    npl_dicts = [
        session_data[(current, ch)]["npl_stpr_dict"]
        for ch in cur_to_channels[current]
        if "npl_stpr_dict" in session_data[(current, ch)]
    ]

    pl_all = _merge_unit_stprs(pl_dicts) if pl_dicts else {}
    npl_all = _merge_unit_stprs(npl_dicts) if npl_dicts else {}

    # ---- Figure 1: PL ----
    fig_pl, ax_pl = plt.subplots(1, 1, figsize=(3.0, 2.5))
    # background: all PL in light gray
    for uid, stpr in pl_all.items():
        if stpr is None or not np.all(np.isfinite(stpr)):
            continue
        ax_pl.plot(lags, stpr, color='0.7', alpha=gray_alpha, lw=1.0)

    # highlight selected pl_id
    if pl_id is not None and pl_id in pl_all and pl_all[pl_id] is not None:
        if np.all(np.isfinite(pl_all[pl_id])):
            ax_pl.plot(lags, pl_all[pl_id], color=PALETTE[0],
                       lw=2.2, label=f'PL {pl_id}')

    # optional mean
    if show_mean and len(pl_all) > 0:
        try:
            M = np.vstack([v for v in pl_all.values()
                          if v is not None and np.all(np.isfinite(v))])
            if M.size > 0:
                ax_pl.plot(lags, M.mean(axis=0), color=plt.cm.Blues(
                    0.65), lw=2.0, label='Mean')
        except ValueError:
            pass

    ax_pl.set_title(f"Pulse-locked units stPR")
    ax_pl.set_ylabel("Population\nfiring rate (Hz)")
    ax_pl.set_xlabel("Lag (s)")
    ax_pl.set_ylim(top=ymax)

    ax_pl.set_ylim(-20, ymax)
    ax_pl.autoscale(False)
    fig_pl.tight_layout()

    # ---- Figure 2: NPL ----
    fig_npl, ax_npl = plt.subplots(1, 1, figsize=(3, 2.5))

    for uid, stpr in npl_all.items():
        if stpr is None or not np.all(np.isfinite(stpr)):
            continue
        ax_npl.plot(lags, stpr, color='0.7', alpha=gray_alpha, lw=1.0)

    if npl_id is not None and npl_id in npl_all and npl_all[npl_id] is not None:
        if np.all(np.isfinite(npl_all[npl_id])):
            ax_npl.plot(
                lags, npl_all[npl_id], color=PALETTE[1], lw=2.2, label=f'NPL {npl_id}')

    # optional mean
    if show_mean and len(npl_all) > 0:
        try:
            M = np.vstack([v for v in npl_all.values()
                          if v is not None and np.all(np.isfinite(v))])
            if M.size > 0:
                ax_npl.plot(lags, M.mean(axis=0), color=plt.cm.Oranges(
                    0.7), lw=2.0, label='Mean')
        except ValueError:
            pass

    ax_npl.set_title(f"Non-pulse-locked units stPR")
    ax_npl.set_ylabel("Population\nfiring rate (Hz)")
    ax_npl.set_xlabel("Lag (s)")

    ax_npl.set_ylim(-20, ymax)
    ax_npl.autoscale(False)
    fig_npl.tight_layout()

    return (fig_pl, ax_pl), (fig_npl, ax_npl)


animal_id = 'ICMS92'
pickle_file_name = f"pop_coupling/{animal_id}_pop_coupling_v4.pkl"
with open(pickle_file_name, 'rb') as f:
    pc_animal_data = pickle.load(f)
session_names = list(pc_animal_data.keys())

bin_size = 0.001
win = 0.1
win_bins = int(round(win / bin_size))
lags = np.arange(-win_bins, win_bins+1) * bin_size


for session in session_names[3:4]:
    session_data = pc_animal_data[session]
    (fig_pl, ax_pl), (fig_npl, ax_npl) = plot_session_by_current(
        session_data, lags, animal_id='ICMS92', session_name=session,
        current=4, pl_id=4, npl_id=11)
