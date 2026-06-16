from scipy.stats import mannwhitneyu
from batch_process.util.plotting import add_sig_bracket, p_to_stars
import batch_process.util.template_util as template_util
import seaborn as sns
import batch_process.postprocessing.responses_v2.response_plotting_util as rpu
from batch_process.util import curate_util
import spikeinterface.full as si
import numpy as np
import pandas as pd
import batch_process.postprocessing.longitudinal_data_utils as utils
import pickle
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from batch_process.util.plotting import add_scale_bars_wvf
import sys
from pathlib import Path

project_root = Path(
    "C:/Users/robin/Desktop/Python/research/spikeinterface/icms-plasticity")
sys.path.insert(0, str(project_root))

# helper functions


def filter_responses(df, max_z_score=100, min_spikes=50):
    # Filter stim responses using the correct comparison for boolean and numerical conditions
    for ch in df['stim_channel'].unique():
        ch_df = df[df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            df = df[df['stim_channel'] != ch]

    df = df[
        (df["rel_week"] < 5) &
        (df["baseline_too_slow"] == False) &
        (df['modulated'] == True) &
        (df['mod_p_val'] < 0.05) &
        (df["z_score"] > 0) &
        # (df["z_score"] < 15) &
        (df["is_pulse_locked"] == True) &
        (df["num_trials"] > 5) &
        # (df["pli"] >= 0.02) &
        # (df['jitter']**2 < 0.1) &
        # (
        #     # If z_score > 0, must be pulse-locked
        #     ((df["z_score"] > 0) & (df["is_pulse_locked"] == True)) |
        #     # Allow negative or zero z_scores without pulse-locked condition
        #     (df["z_score"] <= 0)
        # ) &
        (df['num_spikes'] > 50) &
        (df['stim_current'] < 7) &
        (df['stim_current'] > 3)]

    return df


def filter_responses_B(df, max_z_score=100, min_spikes=50):
    # Filter stim responses using the correct comparison for boolean and numerical conditions
    for ch in df['stim_channel'].unique():
        ch_df = df[df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            df = df[df['stim_channel'] != ch]

    df = df[
        (df["rel_week"] < 5) &
        (df["baseline_too_slow"] == False) &
        (df['modulated'] == True) &
        (df['mod_p_val'] < 0.05) &
        (df["z_score"] > 0) &
        (df["z_score"] < max_z_score) &
        (df["is_pulse_locked"] == False) &
        (df["num_trials"] > 5) &
        # (df["pli"] < 0.02) &
        # (df['jitter']**2 >= 0.1) &
        # # (
        #     # If z_score > 0, must be pulse-locked
        #     ((df["z_score"] > 0) & (df["is_pulse_locked"] == True)) |
        #     # Allow negative or zero z_scores without pulse-locked condition
        #     (df["z_score"] <= 0)
        # ) &
        (df['num_spikes'] > min_spikes) &
        (df['stim_current'] < 7) &
        (df['stim_current'] > 3)]

    return df


def filter_responses_just_modulated(df, max_z_score=100, min_spikes=50):
    # Filter stim responses using the correct comparison for boolean and numerical conditions
    for ch in df['stim_channel'].unique():
        ch_df = df[df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            df = df[df['stim_channel'] != ch]

    df = df[
        (df["rel_week"] < 5) &
        (df["baseline_too_slow"] == False) &
        (df['modulated'] == True) &
        (df['mod_p_val'] < 0.05) &
        (df["num_trials"] > 5) &
        (df["z_score"] > 0) &
        (df["z_score"] < max_z_score) &
        (df['num_spikes'] > min_spikes) &
        (df['stim_current'] < 7) &
        (df['stim_current'] > 3)]

    return df


def _icms83_depth(probe_height, outside_length=60):
    """Compute cortical depth for ICMS83 (flipped probe orientation).

    ICMS83 probe is vertically flipped relative to experimental animals:
    probe y=0 is the shallow end (surface), y=1860 is deep.
    So probe_height IS the depth from surface on the probe.
    """
    angle_from_vertical_deg = 30
    projection = np.cos(np.deg2rad(angle_from_vertical_deg))
    return (probe_height - outside_length) * projection


# ICMS83 stim depth index -> probe y-position (from analyzer channel locations)
_ICMS83_STIM_Y = {3: 120, 4: 180, 5: 240}


def load_icms83_data():
    """Load ICMS83 session CSVs and compute rel_day/rel_week + depth."""
    from datetime import datetime
    icms83_root = Path("E:/ICMS83")
    sessions = [
        "20-Jul-2023", "24-Jul-2023", "26-Jul-2023", "28-Jul-2023",
        "31-Jul-2023", "02-Aug-2023", "04-Aug-2023", "07-Aug-2023",
        "10-Aug-2023",
    ]
    all_dfs = []
    for session_name in sessions:
        csv_path = icms83_root / session_name / "batch_sort" / "stim_condition_results_all.csv"
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        df["session"] = session_name
        df["session_date"] = datetime.strptime(session_name, "%d-%b-%Y")
        all_dfs.append(df)

    if not all_dfs:
        return pd.DataFrame()

    master = pd.concat(all_dfs, ignore_index=True).sort_values("session_date")
    first_date = master["session_date"].min()
    master["rel_day"] = (master["session_date"] - first_date).dt.days
    master["rel_week"] = np.ceil(master["rel_day"] / 7).astype(int)
    master.loc[master["rel_day"] == 0, "rel_week"] = 0
    master["animal_id"] = "ICMS83"

    # Compute depth (flipped probe: y IS depth from surface)
    if "unit_location" in master.columns:
        unit_y = master["unit_location"].apply(
            lambda s: np.fromstring(str(s).strip("[]"), sep=" ")[1]
            if pd.notna(s) else np.nan)
        master["unit_depth_um"] = _icms83_depth(unit_y)

    # Stim channel depth (depth index -> probe y -> cortical depth)
    master["stim_ch_depth_um"] = master["stim_channel"].map(_ICMS83_STIM_Y)
    master["stim_ch_depth_um"] = _icms83_depth(master["stim_ch_depth_um"])

    return master


def get_raw_df(stim_trial_type='all', include_icms83=True):
    animal_ids = ["ICMS92", "ICMS93", "ICMS98",
                  "ICMS100", "ICMS101"]

    base_path = Path("C:/data/")

    # Load session data for each animal
    dfs = [utils.load_session_data(animal_id, base_path, stim_trial_type)
           for animal_id in animal_ids]

    # Concatenate the DataFrames for all animals
    df = pd.concat(dfs, ignore_index=True)
    df = utils.add_data_to_df(df)

    if include_icms83:
        icms83_df = load_icms83_data()
        if not icms83_df.empty:
            # Align columns — ICMS83 CSVs already have cell_type, unit_location
            # add_data_to_df columns that ICMS83 won't have
            for col in ['t2p_ms', 'stim_ch_depth_um', 'unit_depth_um', 'detection_threshold']:
                if col not in icms83_df.columns:
                    icms83_df[col] = np.nan
            # Keep only columns that exist in both
            common_cols = sorted(set(df.columns) & set(icms83_df.columns))
            # Align datetime resolution to avoid concat error
            if 'session_date' in common_cols:
                df['session_date'] = df['session_date'].astype('datetime64[us]')
                icms83_df['session_date'] = icms83_df['session_date'].astype('datetime64[us]')
            df = pd.concat([df[common_cols], icms83_df[common_cols]], ignore_index=True)
            print(f"Added ICMS83: {len(icms83_df)} rows")

    return df


def get_modulated_only_responses(raw_df, max_z_score=100, min_spikes=50):
    df = filter_responses_just_modulated(raw_df, max_z_score, min_spikes)
    return df


def get_filtered_responses(raw_df,  max_z_score=100, min_spikes=50, use_filt_B=False):
    if use_filt_B:
        df = filter_responses_B(raw_df, max_z_score, min_spikes)
    else:
        df = filter_responses(raw_df, max_z_score, min_spikes)
    return df


def plot_modulation_grouped_by_current(df, var, ax, use_median=True, alt='less'):

    ax = utils.plot_single_metric_aggregated(
        df, var=var, ax=ax, use_median=use_median, animal_id='Aggregated', comparison_alternative=alt)

    ax.set_xlabel('Weeks Relative to First Session')


def compare_two_groups(df, var, ax, time_group1=[1, 2], time_group2=[4, 5], alt='less'):
    utils.plot_single_metric_aggregated_compare_two_time_groups(
        df, var=var, time_group1=time_group1, time_group2=time_group2,
        use_median=True, animal_id='Aggregated', alt=alt, ax=ax)

    # ax.set_xlabel('Weeks Relative to First Session')


def plot_modulation_at_threshold_current(df, raw_df, ax):

    ax1, ax2 = utils.plot_single_metric_at_threshold_current(
        df, raw_df, var="z_score", ax=ax)
    ax1.set_ylabel('Z-Score')
    ax1.set_xlabel('Weeks Relative to First Session')
    ax1.set_title('ICMS Train Modulation Score\nat Threshold Current')
    ax2.spines.right.set_visible(True)


def plot_unit_counts(raw_df, ax):
    ax = utils.plot_animal_unit_counts_by_week(raw_df, ax=ax)
    ax.set_xticks(np.arange(5))
    ax.set_ylim([0, 40])
    ax.set_xlabel("Weeks")
    ax.set_ylabel('Unit Count')


def load_data_for_raw_trace(data_folder):
    data_folder = Path(data_folder)

    # Load session responses
    pkl_path = data_folder / "batch_sort/session_responses_all.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)

    # Load analyzers
    analyzer = si.load_sorting_analyzer(
        folder=data_folder / "batch_sort/merge/hmerge_analyzer.zarr")
    raw_analyzer = curate_util.load_raw_analyzer(data_folder, analyzer)

    return session_responses, raw_analyzer, analyzer


def load_data_for_waveforms(plt_data):
    raw_analyzer = plt_data["raw_analyzer"]
    analyzer = plt_data["analyzer"]
    stim_ts = plt_data["stim_ts"]
    t1 = 0
    t2 = int(len(stim_ts)) - 1
    rec_raw = raw_analyzer.recording
    rec_filt = analyzer.recording
    raw_traces = rec_raw.get_traces(
        start_frame=stim_ts[t1], end_frame=stim_ts[t2], return_scaled=True)
    filt_traces = rec_filt.get_traces(
        start_frame=stim_ts[t1], end_frame=stim_ts[t2], return_scaled=True)
    return raw_traces, filt_traces


def plot_all_raw_and_filt_waveforms(raw_traces, filt_traces, plt_data, axes=None):
    analyzer = plt_data["analyzer"]
    session_responses = plt_data['session_responses']
    unit_id = plt_data["unit_id"]
    stim_ts = plt_data["stim_ts"]
    scr = plt_data['scr']

    stim_spikes = scr.stim_spikes
    stim_spikes_norm = stim_spikes - stim_ts[0]
    primary_ch_index = template_util.get_dense_primary_ch_index(
        analyzer, unit_id)

    filt_trace = filt_traces[:, primary_ch_index]
    raw_trace = raw_traces[:, primary_ch_index]

    raw_spikes = []
    filt_spikes = []

    for spike_ts in stim_spikes_norm:
        spike_ts = int(spike_ts)
        raw_spike = raw_trace[spike_ts - 70: spike_ts + 60]
        filt_spike = filt_trace[spike_ts - 70: spike_ts + 60]
        if np.any(raw_spike[0:10] > -1200):
            if len(raw_spike) > 0:
                raw_spikes.append(raw_spike)
            if len(filt_spike) > 0:
                filt_spikes.append(filt_spike)

    raw_spikes = np.array(raw_spikes)
    filt_spikes = np.array(filt_spikes)

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(4, 2))

    # Plot raw spikes
    axes[0].plot(raw_spikes.T, 'k', linewidth=0.5, alpha=0.3)
    axes[0].set_ylim([-3200, -1300])
    axes[0].set_title("Raw Spikes", fontsize=9)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')

    # Plot filtered spikes
    axes[1].plot(filt_spikes.T, 'k', linewidth=0.5, alpha=0.3)
    axes[1].set_ylim([-300, 200])
    axes[1].set_title("Preprocessed Spikes", fontsize=9)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].set_xlabel('')
    axes[1].set_ylabel('')

    return raw_spikes, filt_spikes, axes


def plot_good_raw_and_filt_waveforms(raw_spikes, filt_spikes, axes=None):

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(4, 2))

    good_spike_indices = np.where(
        (np.argmax(raw_spikes, axis=1) < 30) &
        (np.min(raw_spikes, axis=1) > -4000) &
        (np.min(filt_spikes, axis=1) > -220) &
        (np.all(filt_spikes[:, 75:] > -100, axis=1)) &
        (raw_spikes[:, 23] < -1700) &
        (raw_spikes[:, 37] > -2900) &
        (filt_spikes[:, 119] < 100)
    )

    good_raw_spikes = raw_spikes[good_spike_indices]
    good_filt_spikes = filt_spikes[good_spike_indices]

    # Plot the good raw spikes in the first subplot
    axes[0].plot(good_raw_spikes[:, 12:].T, 'k', linewidth=0.5, alpha=0.25)
    axes[0].set_ylim([-3300, -1200])
    axes[0].set_title("Raw Spikes", fontsize=6, pad=-1000)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')

    add_scale_bars_wvf(
        axes[0],
        h_pos=[10, -3150],
        v_pos=[10, -3150],
        h_length=30,
        v_length=200,
        line_width=1,
        h_label="1 ms",
        v_label="200 uV",
        v_label_x_offset=-0.5,
        v_label_y_offset=150,
        h_label_x_offset=10,
        h_label_y_offset=-45,
        font_size=5
    )

    # Plot the good filtered spikes in the second subplot
    axes[1].plot(good_filt_spikes[:, 12:].T, 'k', linewidth=0.5, alpha=0.25)
    axes[1].set_ylim([-350, 200])
    axes[1].set_title("Filtered Spikes", fontsize=6, pad=-1000)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].set_xlabel('')
    axes[1].set_ylabel('')

    add_scale_bars_wvf(
        axes[1],
        h_pos=[10, -300],
        v_pos=[10, -300],
        h_length=30,
        v_length=100,
        line_width=1,
        h_label="1 ms",
        v_label="100 uV",
        v_label_x_offset=-0.5,
        v_label_y_offset=50,
        h_label_x_offset=10,
        h_label_y_offset=-10,
        font_size=5
    )

    for ax in axes.flat:
        ax.axis('off')

    return axes


def plot_fig4_raw_trace_with_spikes(data, ax):
    colors = sns.color_palette('deep')

    raw_analyzer = data["raw_analyzer"]
    analyzer = data["analyzer"]
    unit_id = data["unit_id"]
    stim_ts = data["stim_ts"]
    start_pulse = data["start_pulse"]
    end_pulse = data["end_pulse"]
    colors = sns.color_palette("deep")

    ax = utils.plot_primary_ch_trace_with_unit_spikes(
        raw_analyzer, analyzer, unit_id, stim_ts, start_pulse, end_pulse, ch_id=None, ax=ax)

    # raw_analyzer, analyzer, unit_id,
    #                                            stim_ts, start_pulse, end_pulse,
    #                                            ch_id=None, marker_y_offset=150, ax=None

    start_frame = stim_ts[start_pulse]
    stim_ts_ms = (stim_ts[start_pulse:end_pulse] - start_frame) / 30
    for i, stim_onset in enumerate(stim_ts_ms):
        ax.vlines(
            x=stim_onset,
            ymin=-2900,
            ymax=-1300,
            color=colors[3],
            linestyle='-',
            linewidth=1,
            label='Stim onset' if i == 0 else ""
        )

        ax.fill_betweenx(
            y=[-2900, -1300],  # Specify the vertical range
            x1=stim_onset - 0.5,
            x2=stim_onset + 1.5,
            color=colors[3],
            alpha=0.4,  # 0.2
            linewidth=0,
            label='Blanked' if i == 0 else ""
        )

    add_scale_bars_wvf(
        ax,
        h_pos=[-2, -2950],
        v_pos=[-2, -2950],
        h_length=5,
        v_length=200,
        line_width=1,
        h_label="5 ms",
        v_label="200 uV",
        v_label_x_offset=-0.5,
        v_label_y_offset=100,
        h_label_x_offset=2,
        h_label_y_offset=-20,
        font_size=4
    )

    # add_scale_bars_wvf(
    #     ax,
    #     h_pos=[-1, -2950],
    #     v_pos=[-1, -2950],
    #     h_length=5,
    #     v_length=200,
    #     line_width=1,
    #     h_label="5 ms",
    #     v_label="200 uV",
    #     v_label_x_offset=-0.2,
    #     v_label_y_offset=100,
    #     h_label_x_offset=2,
    #     h_label_y_offset=-20,
    #     font_size=4
    # )
    # ax.scatter(x=24.42, y=-2520, color=colors[0], marker='^', s=12)
    ax.scatter(x=14.42, y=-2520, color=colors[0], marker='^', s=12)

    ax.set_ylim([-3000, -1300])
    plt.legend(loc='upper center',
               bbox_to_anchor=(0.5, 0.05),
               ncol=3,
               frameon=False)
    plt.tight_layout(rect=[-0.1, 0, 1, 1])
    return ax, data


def _get_param(name, local_value, param_dict):
    """Helper to resolve parameter from user input or dict, or raise error."""
    value = local_value if local_value is not None else param_dict.get(name)
    if value is None:
        raise ValueError(
            f"`{name}` must be provided either as an argument or in raw_trace_plt_params.")
    return value


def plot_fig4_raster_and_fr(axes):
    data_folder = 'C:\\data\\ICMS93\\behavior\\22-Sep-2023'
    data_folder = 'C:\\data\\ICMS98\\14-Nov-2023'

    pkl_path = Path(data_folder) / "batch_sort/session_responses_all.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)

    unit_response = session_responses.get_unit_response(15)
    sc1 = unit_response.stim_responses[(6, 4)]
    sc2 = unit_response.stim_responses[(6, 5)]
    sc3 = unit_response.stim_responses[(6, 6)]

    responses = [(4, sc1), (5, sc2), (6, sc3)]
    utils.plot_train_rasters(responses, axes[0])
    utils.plot_train_firing_rate(responses, axes[1])

    color_map = utils.get_stim_colormap()
    handles = [
        Line2D([0], [0], color=color_map[current], label=f"{current} µA")
        for current in [4, 5, 6]
    ]
    axes[1].legend(handles=handles)
    axes[1].legend(frameon=False)
    axes[1].legend(loc="center left", bbox_to_anchor=(
        0.6, 0.75), handlelength=1.5)

    plt.tight_layout()


def plot_units_depth_distribution(df_raw, ax):

    df = df_raw[
        (df_raw['unit_depth_um'] > 0) &
        (df_raw['unit_depth_um'] < 1400)
        # (df_raw['mod_p_val'] < 0.005) &
        # (df_raw['z_score'] > 0)
    ].copy()

    df_unique = df.drop_duplicates(subset=['animal_id', 'session', 'unit_id'])

    cortical_depths = df_unique['unit_depth_um']
    bin_size = 100
    bins = np.arange(0, np.max(cortical_depths) + bin_size, bin_size)

    plt.hist(cortical_depths, orientation='horizontal',
             color='k', bins=bins,  histtype="stepfilled")
    plt.gca().invert_yaxis()
    plt.xlabel('Number of units')
    plt.ylabel('Depth (µm)')
    # plt.title('Unit cortical depth')

    plt.tight_layout()
    plt.show()


def plot_units_t2p_time(df_raw, ax):

    df = df_raw[
        (df_raw['unit_depth_um'] > 0) &
        (df_raw['unit_depth_um'] < 1400)
        # (df_raw['mod_p_val'] < 0.005) &
        # (df_raw['z_score'] > 0)
    ].copy()

    df_unique = df.drop_duplicates(subset=['animal_id', 'session', 'unit_id'])
    # all_r2 = df_unique['r2']
    all_t2p = df_unique['t2p_ms']
    all_t2p = all_t2p[all_t2p > 0]
    threshold = 0.45
    in_t2p = all_t2p[all_t2p < threshold]
    py_t2p = all_t2p[all_t2p >= threshold]
    num_bins = 23
    # num_bins = 17
    bins = np.linspace(np.min(all_t2p), np.max(all_t2p), num_bins)
    colors = sns.color_palette('deep')
    ax.hist(in_t2p, bins=bins, histtype="stepfilled",
            color=colors[7], rwidth=1.0, edgecolor='none')
    ax.hist(py_t2p, bins=bins, histtype="stepfilled",
            color=colors[3], rwidth=1.0, edgecolor='none')

    ymax = plt.ylim()[1]
    ax.text(threshold * 1.01, ymax * 0.8,
            f'{threshold:.2f} ms', rotation=270, va='top', fontsize=5)
    ax.axvline(threshold, c='k', linestyle='--', linewidth=1)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Number of units')
    # ax.title('Unit trough-to-peak time')
    plt.tight_layout()


def plot_cell_type_feature_over_weeks(df, var, ax, alt='less', cell_type='pyramidal', star_size=8, stars_y_offset=0, stars_spacing=0.07):
    cell_df = df[df['cell_type'] == cell_type]
    ax, stars_arr, all_y_vals = utils.plot_single_metric_two_group_comparison(
        cell_df, var=var, alt=alt, ax=ax)
    ax = utils.draw_two_group_error_bars(
        stars_arr, all_y_vals, stars_y_offset=stars_y_offset, star_size=star_size, y_scale_factor=1.05, stars_spacing=stars_spacing, ax=ax)
    ax.legend(bbox_to_anchor=(0.2, 0.8), borderaxespad=0.)
    plt.legend().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()


def plot_fig5_PL_raw_trace_with_spikes(raw_trace_plt_params, start_pulse=None,
                                       end_pulse=None, ch=None, cu=None,
                                       unit_id=None, plot_ch_id=None, ax=None, marker_y_offset=150):

    session_responses = raw_trace_plt_params['session_responses']
    raw_analyzer = raw_trace_plt_params['raw_analyzer']
    analyzer = raw_trace_plt_params['analyzer']

    start_pulse = _get_param('start_pulse', start_pulse, raw_trace_plt_params)
    end_pulse = _get_param('end_pulse', end_pulse, raw_trace_plt_params)
    ch = _get_param('ch', ch, raw_trace_plt_params)
    cu = _get_param('cu', cu, raw_trace_plt_params)
    unit_id = _get_param('unit_id', unit_id, raw_trace_plt_params)

    ur = session_responses.get_unit_response(unit_id)
    scr = ur.get_stim_response(ch, cu)
    stim_ts = scr.stim_timestamps

    ax = utils.plot_primary_ch_trace_with_unit_spikes(
        raw_analyzer, analyzer, unit_id, stim_ts, start_pulse, end_pulse, plot_ch_id, marker_y_offset, ax)

    start_frame = stim_ts[start_pulse]
    stim_ts_ms = (stim_ts[start_pulse:end_pulse] - start_frame) / 30
    for i, stim_onset in enumerate(stim_ts_ms):
        ax.vlines(
            x=stim_onset,
            ymin=-3900,
            ymax=-1800,
            color='C3',
            linestyle='-',
            linewidth=0.5,
            label='Stim onset' if i == 0 else ""
        )

        ax.fill_betweenx(
            y=[-3900, -1800],  # Specify the vertical range
            x1=stim_onset - 0.5,
            x2=stim_onset + 1.5,
            color='C3',
            alpha=0.2,
            linewidth=0,
            label='Blanked' if i == 0 else ""
        )

    add_scale_bars_wvf(
        ax,
        h_pos=[-2, -4150],
        v_pos=[-2, -4150],
        h_length=5,
        v_length=200,
        line_width=1,
        h_label="5 ms",
        v_label="200 uV",
        v_label_x_offset=-0.5,
        v_label_y_offset=100,
        h_label_x_offset=2,
        h_label_y_offset=-20,
        font_size=4
    )
    ax.set_ylim([-5000, -3000])
    plt.legend(loc='upper center',
               bbox_to_anchor=(0.5, 0.25),
               ncol=3,
               frameon=False)
    plt.tight_layout(rect=[-0.1, 0, 1, 1])

    raw_trace_plt_params['start_pulse'] = start_pulse
    raw_trace_plt_params['end_pulse'] = end_pulse
    raw_trace_plt_params['ch'] = ch
    raw_trace_plt_params['cu'] = cu
    raw_trace_plt_params['unit_id'] = unit_id
    raw_trace_plt_params['stim_ts'] = stim_ts

    return raw_trace_plt_params


def plot_fig5_NPL_raw_trace_with_spikes(raw_trace_plt_params, start_pulse=None,
                                        end_pulse=None, ch=None, cu=None,
                                        unit_id=None, plot_ch_id=None, ax=None, marker_y_offset=150):

    session_responses = raw_trace_plt_params['session_responses']
    raw_analyzer = raw_trace_plt_params['raw_analyzer']
    analyzer = raw_trace_plt_params['analyzer']

    start_pulse = _get_param('start_pulse', start_pulse, raw_trace_plt_params)
    end_pulse = _get_param('end_pulse', end_pulse, raw_trace_plt_params)
    ch = _get_param('ch', ch, raw_trace_plt_params)
    cu = _get_param('cu', cu, raw_trace_plt_params)
    unit_id = _get_param('unit_id', unit_id, raw_trace_plt_params)

    ur = session_responses.get_unit_response(unit_id)
    scr = ur.get_stim_response(ch, cu)
    stim_ts = scr.stim_timestamps

    ax = utils.plot_primary_ch_trace_with_unit_spikes(
        raw_analyzer, analyzer, unit_id, stim_ts, start_pulse, end_pulse, plot_ch_id, marker_y_offset, ax)

    start_frame = stim_ts[start_pulse]
    stim_ts_ms = (stim_ts[start_pulse:end_pulse] - start_frame) / 30
    for i, stim_onset in enumerate(stim_ts_ms):
        ax.vlines(
            x=stim_onset,
            ymin=-1500,
            ymax=-300,
            color='C3',
            linestyle='-',
            linewidth=0.5,
            label='Stim onset' if i == 0 else ""
        )

        ax.fill_betweenx(
            y=[-1500, -300],  # Specify the vertical range
            x1=stim_onset - 0.5,
            x2=stim_onset + 1.5,
            color='C3',
            alpha=0.2,
            linewidth=0,
            label='Blanked' if i == 0 else ""
        )

    add_scale_bars_wvf(
        ax,
        h_pos=[-2, -1550],
        v_pos=[-2, -1550],
        h_length=5,
        v_length=200,
        line_width=1,
        h_label="5 ms",
        v_label="200 uV",
        v_label_x_offset=-0.5,
        v_label_y_offset=100,
        h_label_x_offset=2,
        h_label_y_offset=-20,
        font_size=4
    )
    ax.set_ylim([-1600, -450])
    plt.legend(loc='upper center',
               bbox_to_anchor=(0.5, 0),
               ncol=3,
               frameon=False)
    plt.tight_layout(rect=[-0.1, 0, 1, 1])

    raw_trace_plt_params['start_pulse'] = start_pulse
    raw_trace_plt_params['end_pulse'] = end_pulse
    raw_trace_plt_params['ch'] = ch
    raw_trace_plt_params['cu'] = cu
    raw_trace_plt_params['unit_id'] = unit_id
    raw_trace_plt_params['stim_ts'] = stim_ts

    return raw_trace_plt_params


def plot_num_responses_vs_time(df, current, ax):

    ax1, ax2 = utils.plot_num_responses_vs_time(df, current, ax)
    ax1.set_ylim([0, 30])
    ax2.set_ylim([0.3, 0.7])
    ax2.spines.right.set_visible(True)


def plot_num_responses_vs_time_at_threshold(df, ax):

    ax1, ax2 = utils.plot_num_responses_vs_time_at_threshold(df, ax)


def plot_responses_2(df, metric, ax):

    ax = utils.plot_responses_2(df, metric, ax)


def plot_pulse_response(pulse_locked=True, axs=None):
    if pulse_locked:
        session_responses = utils.load_pl_session_responses()
        unit_id = 18
        unit_response = session_responses.get_unit_response(unit_id)
        sc1 = unit_response.stim_responses[(11, 4)]
        sc2 = unit_response.stim_responses[(11, 5)]
        sc3 = unit_response.stim_responses[(11, 6)]
    else:
        session_responses = utils.load_npl_session_responses()
        unit_id = 17
        unit_response = session_responses.get_unit_response(unit_id)
        sc1 = unit_response.stim_responses[(4, 4)]
        sc2 = unit_response.stim_responses[(4, 5)]
        sc3 = unit_response.stim_responses[(4, 6)]

    responses = [(4, sc1), (5, sc2), (6, sc3)]
    rpu.plot_pulse_raster(responses, axs)
    rpu.plot_probability(responses, axs)
    return axs


def load_data_for_raw_trace_fig_pulse_locked():
    data_folder = 'C:\\data\\ICMS92\\behavior\\01-Sep-2023'
    # get all
    pkl_path = Path(data_folder) / "batch_sort/session_responses_all.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)

    data_folder = Path(data_folder)
    analyzer = si.load_sorting_analyzer(
        folder=data_folder / "batch_sort/merge/hmerge_analyzer.zarr")
    raw_analyzer = curate_util.load_raw_analyzer(data_folder, analyzer)
    start_pulse = 230
    end_pulse = 238

    unit_id = 2
    session_responses.load_sorting_analyzer()
    rec = session_responses.sorting_analyzer.recording
    ur = session_responses.get_unit_response(unit_id)
    scr = ur.get_stim_response(9, 4)
    stim_ts = scr.stim_timestamps

    return raw_analyzer, analyzer, unit_id, stim_ts, start_pulse, end_pulse


def plot_raw_vs_filt_waveforms(data, ax):
    stim_wvfs = data['stim_wvfs']
    non_stim_wvfs = data['non_stim_wvfs']
    ch_locs = data['ch_locs']
    sparse_ch_indices = data['sparse_ch_indices']

    mean_non_stim = non_stim_wvfs.mean(axis=0)
    mean_stim = stim_wvfs.mean(axis=0)
    std_non_stim = non_stim_wvfs.std(axis=0)
    std_stim = stim_wvfs.std(axis=0)

    utils.plot_stim_vs_nonstim_template(mean_non_stim, mean_stim,
                                        std_non_stim, std_stim,
                                        ch_locs, sparse_ch_indices, ax=ax)


def plot_spike_sorting_validation(ax):

    df = pd.read_pickle('spike_sorting_validation.pkl')
    print(df['type'].value_counts())

    within_dists = df[df['type'] == 'within']['distance']
    between_dists = df[df['type'] == 'between']['distance']
    stat, pval = mannwhitneyu(within_dists, between_dists, alternative='less')

    summary = df.groupby('type')['distance'].agg(
        median='median',
        q1=lambda x: np.percentile(x, 25),
        q3=lambda x: np.percentile(x, 75)
    ).reset_index()
    summary['yerr'] = (summary['q3'] - summary['q1']) / 2
    summary['type'] = pd.Categorical(summary['type'], categories=[
                                     'within', 'between'], ordered=True)
    summary = summary.sort_values('type')
    x = np.arange(len(summary))
    y = summary['median'].values
    yerr = summary['yerr'].values

    ax.errorbar(x, y, yerr=yerr, fmt='o', color='k', capsize=3, markeredgecolor='none',
                elinewidth=0.32,
                capthick=0.32,
                markersize=3,)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['within unit\n(stim vs. nonstim)',
                       'between units\n(stim vs. stim)'], rotation=0)
    ax.set_ylabel('Centroid Distance')
    y_max = (np.array(yerr) + np.array(y)).max() * 1.1
    ax.set_xlabel('')
    ax.set_xticks(x)
    ax.set_xlim([-0.25, len(x) - 0.75])

    # add_sig_line(0, 1, y_max, pval)
    add_sig_bracket(ax, 0, 1, y_max, h=15, text=p_to_stars(pval), lw=0.32)

    plt.ylim(0, y_max * 1.2)

    plt.tight_layout()
