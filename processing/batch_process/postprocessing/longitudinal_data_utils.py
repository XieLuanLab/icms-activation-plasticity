from batch_process.util.plotting import add_scale_bars_wvf
from batch_process.util import template_util
from matplotlib.lines import Line2D
from scipy.stats import ttest_ind, mannwhitneyu, shapiro
from scipy.stats import mannwhitneyu
import os
import numpy as np
import matplotlib.pyplot as plt
import batch_process.util.file_util as file_util
from pathlib import Path
import pandas as pd
from datetime import datetime
import batch_process.postprocessing.responses_v2.response_plotting as plot_util
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import linregress
import pickle
from scipy.stats import norm
import batch_process.postprocessing.responses_v2.response_plotting_util as rpu


def process_data_folder(folder, stim_trial_type='hit'):
    """Process a single data folder and return its DataFrame."""
    print(f"Processing data folder {folder}...")
    save_folder = os.path.join(folder, "batch_sort")
    if not save_folder:
        return None

    animalID = file_util.get_animal_id(save_folder)
    date_str = file_util.get_date_str(save_folder)

    df_path = Path(save_folder) / \
        f"stim_condition_results_{stim_trial_type}.csv"
    if not os.path.exists(df_path):
        return None

    df = pd.read_csv(df_path)
    df["date_str"] = date_str
    df["date"] = datetime.strptime(date_str, "%d-%b-%Y")
    return df


def correct_df_unit_location(unit_probe_height, animal_id):
    """
    Converts unit height along probe to cortical depth (vertical),
    accounting for insertion angle and number of contacts outside cortex.

    Inputs:
    - unit_probe_height: height along the probe (0 at bottom, 1860 at top)
    - animal_id: used to look up how many microns are outside brain

    Output:
    - cortical depth in µm (0 = surface, > 0 = deeper into cortex)
    """
    probe_length = 31 * 60
    unit_probe_depth = probe_length - unit_probe_height
    angle_from_vertical_deg = 30
    projection = np.cos(np.deg2rad(angle_from_vertical_deg))

    # Determine translation based on animal_id
    match animal_id:
        case "ICMS92":  # 5 out
            outside_length = 300
        case "ICMS93":  # 2.5 out
            outside_length = 180
        case "ICMS98":  # 1 out
            outside_length = 60
        case "ICMS100":  # 300 deep along 30 deg slope
            outside_length = -300
        case "ICMS101":  # 7 out
            outside_length = 360
        case _:
            print("Error: invalid animal_id")
            return

    cortical_depth = (unit_probe_depth - outside_length) * projection
    return cortical_depth


def add_data_to_df(input_df):
    """
    Adds additional data columns to dataframe with all animals and sessions.

    Args: 
        pd.DataFrame: Combined DataFrame with data for all animals and sessions.
    Returns:
        pd.DataFrame: Combined DataFrame with additional data columns. 
    """
    def get_response_df(input_df, animal_id, session_id, unit_id):
        return input_df[(input_df['animal_id'] == animal_id) &
                        (input_df['session'] == session_id) &
                        (input_df['unit_id'] == unit_id)]

    with open(Path(__file__).resolve().parents[2] / "cell_type.pkl", 'rb') as f:
        cell_type_dict = pickle.load(f)

    df_new = []

    threshold = 0.45
    for animal_id in cell_type_dict:
        session_unit_dict = cell_type_dict[animal_id]
        for session_id in session_unit_dict:
            unit_dict = session_unit_dict[session_id]

            for unit_id in unit_dict:
                template = unit_dict[unit_id]['template']
                t2p = unit_dict[unit_id]['t2p']
                cell_type = 'interneuron' if t2p < threshold else 'pyramidal'
                if (np.argmax(template) < 80) & (np.max(template[65:]) < 25) & (np.min(template[0:25]) > -100) & (t2p <= 1):
                    # height on probe
                    probe_location = unit_dict[unit_id]['unit_location']
                    cortical_depth = correct_df_unit_location(
                        probe_location, animal_id)
                    # r2 = unit_dict[unit_id]['r2']
                    df = get_response_df(
                        input_df, animal_id, session_id, unit_id).copy()
                    # df['acg_r2'] = r2
                    df['cell_type'] = cell_type
                    df['t2p_ms'] = t2p
                    # df['probe_location'] = 1860 - probe_location
                    stim_ch_height = (32 - df['stim_channel']) * 60

                    df['stim_ch_depth_um'] = correct_df_unit_location(
                        stim_ch_height, animal_id)
                    df['unit_depth_um'] = cortical_depth
                    df_new.append(df)

    df_all = pd.concat(df_new, ignore_index=True)
    df_all = df_all.drop(columns=['unit_location'])  # drop multiple columns
    df_all = df_all[(df_all['stim_ch_depth_um'] > 0)
                    & (df_all['unit_depth_um'] > 0)]

    return df_all


def load_session_data(animal_id, base_path, stim_trial_type='all'):
    """
    Loads the session data for a given animal ID.

    Args:
    Returns:
        pd.DataFrame: Combined DataFrame with session data for the animal.
    """
    if stim_trial_type != '':
        stim_trial_type = f'_{stim_trial_type}'
    if animal_id in {"ICMS92", "ICMS93"}:
        session_files = base_path.glob(
            f"{animal_id}/Behavior/*/batch_sort/stim_condition_results{stim_trial_type}.pkl")
    else:
        session_files = base_path.glob(
            f"{animal_id}/*/batch_sort/stim_condition_results{stim_trial_type}.pkl")

    session_files = list(session_files)
    # animal_id (str): The animal ID to filter by.
    # base_path (Path): The base directory where session data is stored.

    df_list = []
    for session_file in session_files:
        with open(session_file, 'rb') as file:
            df = pickle.load(file)
        session_name = session_file.parent.parent.name
        # Extract session name from the path
        df["session"] = session_file.parent.parent.name
        session_date = pd.to_datetime(session_name, format="%d-%b-%Y")
        df["session_date"] = session_date
        df_list.append(df)

    # Combine all session DataFrames
    df = pd.concat(df_list, ignore_index=True)

    # Order by session date
    df = df.sort_values(
        by="session_date").reset_index(drop=True)

    sessions = df['session'].unique()
    relative_days = file_util.convert_dates_to_relative_days(sessions)
    session_to_day = dict(zip(sessions, relative_days))

    df['rel_day'] = df['session'].map(session_to_day)
    df['rel_week'] = np.ceil(df['rel_day'] / 7).astype(int)

    # df = get_stim_train_segment_firing_rate(df)
    df = add_threshold_to_responses(df, animal_id)
    return df


def add_threshold_to_responses(df, animal_id):
    csv_path = Path(__file__).resolve().parent / "behavioral_data" / "thresholds.csv"
    csv_df = pd.read_csv(csv_path)
    for session in csv_df['date'].unique():
        session_rows = df[df['session'] == session]

        csv_rows = csv_df[(csv_df['date'] == session) &
                          (csv_df['animal_id'] == animal_id)]

        if len(csv_rows) == 0:
            continue

        chA = csv_rows['chA'].iloc[0]
        chB = csv_rows['chB'].iloc[0]
        chC = csv_rows['chC'].iloc[0]

        chA_threshold = csv_rows['chA_thresholds'].iloc[0]
        chB_threshold = csv_rows['chB_thresholds'].iloc[0]
        chC_threshold = csv_rows['chC_thresholds'].iloc[0]

        for i, row in session_rows.iterrows():
            stim_channel = row['stim_channel']
            if stim_channel == chA:
                df.loc[i, 'detection_threshold'] = np.round(chA_threshold)
            elif stim_channel == chB:
                df.loc[i, 'detection_threshold'] = np.round(chB_threshold)
            elif stim_channel == chC:
                df.loc[i, 'detection_threshold'] = np.round(chC_threshold)
    return df


def add_sig_modulated_to_responses(df):
    significantly_modulated = []
    for index, row in df.iterrows():
        z_score = row['z_score']
        p_value = 2 * (1 - norm.cdf(abs(z_score)))
        if p_value < 0.05:
            significantly_modulated.append(True)
        else:
            significantly_modulated.append(False)
    df['significantly_modulated'] = significantly_modulated
    return df


def add_time_to_peak_to_responses(df):
    time_to_peak_list = []
    for index, row in df.iterrows():
        fr_times = row['fr_times']  # Assuming this is a list of time points
        frs = row['firing_rate']  # Assuming this is a list of firing rates
        if isinstance(fr_times, np.ndarray) and isinstance(frs, np.ndarray) and len(fr_times) == len(frs):
            valid_bins = np.where((np.array(fr_times) >= 0)
                                  & (np.array(fr_times) <= 700))[0]
            if len(valid_bins) > 0:
                peak_index = np.argmax(np.array(frs)[valid_bins])
                time_to_peak = np.array(fr_times)[valid_bins][peak_index]
            else:
                time_to_peak = np.nan  # No valid times in the range
        else:
            time_to_peak = np.nan  # If data is not valid
        time_to_peak_list.append(time_to_peak)
    df['time_to_peak'] = time_to_peak_list
    return df


def add_max_fr_pulse_index(df):
    def compute_max_pulse(row):
        train_binned_fr = row['train_binned_fr']
        fr_times = row['fr_times']

        # Apply filtering
        valid_indices = (fr_times > 0) & (fr_times < 700)
        filtered_binned_fr = train_binned_fr[valid_indices]

        if len(filtered_binned_fr) == 0:  # Handle empty case
            return np.nan

        return np.argmax(filtered_binned_fr)

    # Apply function row-wise and create new column
    df['max_pulse_number'] = df.apply(compute_max_pulse, axis=1)

    return df


def get_stim_train_segment_firing_rate(df, fraction_range=[0, 0.5]):
    if df.empty or 'train_smoothed_bins_10ms' not in df.columns:
        print("DataFrame is empty or 'train_smoothed_bins_10ms' column is missing.")
        return df
    fr_times = df['train_smoothed_bins_10ms'].iloc[0]
    N = len(df)  # Loop over all rows in the DataFrame
    first_half_avg_fr = []
    second_half_avg_fr = []
    for i in range(N):
        stim_fr = df['train_smoothed_fr_10ms'].iloc[i]

        if fr_times is None or len(fr_times) == 0 or stim_fr is None or len(stim_fr) == 0:
            first_half_avg_fr.append(np.nan)
            second_half_avg_fr.append(np.nan)
            continue  # Skip this iteration if no valid data
        # Get the firing rate for the time window between 0 and 700 ms
        train_seg = stim_fr[(fr_times > 0) & (fr_times < 700)]
        stim_times = fr_times[(fr_times > 0) & (fr_times < 700)]
        # Divide the time segment into two halves
        bin_range = (len(stim_times) * np.array(fraction_range)).astype(int)
        first_half = train_seg[bin_range[0]:bin_range[1]]
        second_half = train_seg[bin_range[1]:]
        # Calculate the average firing rate for each half
        avg_first_half = np.mean(first_half) if len(first_half) > 0 else np.nan
        avg_second_half = np.mean(second_half) if len(
            second_half) > 0 else np.nan
        # Append the results to their respective lists
        first_half_avg_fr.append(avg_first_half)
        second_half_avg_fr.append(avg_second_half)
    # Add the new columns to the DataFrame
    df['first_half_avg_fr'] = first_half_avg_fr
    df['second_half_avg_fr'] = second_half_avg_fr
    return df


def plot_longitudinal_data(df, var1="z_score", var2="pulse_mean_fr"):
    """
    Plots longitudinal data of Z-scores and Pulse Mean FR across sessions.

    Args:
        animal_id (str): The animal ID to plot.
        base_path (Path): The base directory where session data is stored.
    """

    # Calculate median Z-scores and Pulse Mean FR
    median_z_scores = df.groupby(["rel_day", "stim_channel", "stim_current"])[
        var1].median().reset_index()
    median_pulse_fr = (
        df.groupby(["rel_day", "stim_channel", "stim_current"])[
            var2].median().reset_index()
    )

    # Get unique channels and currents
    stim_channels = sorted(median_z_scores["stim_channel"].unique())
    stim_currents = sorted(median_z_scores["stim_current"].unique())

    # Create subplots
    num_channels = len(stim_channels)
    fig, axes = plt.subplots(num_channels, 2, figsize=(
        12, 3 * num_channels), sharex=True)

    if num_channels == 1:
        axes = [axes]  # Ensure axes is iterable if only one subplot

    for i, stim_channel in enumerate(stim_channels):
        # Plot Z-Scores
        ax_z = axes[i, 0] if num_channels > 1 else axes[0]
        channel_data_z = median_z_scores[median_z_scores["stim_channel"]
                                         == stim_channel]

        for stim_current in stim_currents:
            current_data_z = channel_data_z[channel_data_z["stim_current"]
                                            == stim_current]
            ax_z.plot(
                current_data_z["rel_day"], current_data_z[var1], marker="o", label=f"{stim_current} µA"
            )

        ax_z.set_title(f"Channel {stim_channel} - {var1}")
        ax_z.set_ylabel(f"Median {var1}")
        ax_z.set_xlabel("Days Relative to First Session")

        # Plot Pulse Mean FR
        ax_fr = axes[i, 1] if num_channels > 1 else axes[1]
        channel_data_fr = median_pulse_fr[median_pulse_fr["stim_channel"]
                                          == stim_channel]

        for stim_current in stim_currents:
            current_data_fr = channel_data_fr[channel_data_fr["stim_current"]
                                              == stim_current]
            ax_fr.plot(
                current_data_fr["rel_day"],
                current_data_fr[var2],
                marker="o",
                label=f"{stim_current} µA",
            )

        ax_fr.set_title(f"Channel {stim_channel} - {var2}")
        ax_fr.set_ylabel(f"Median {var2}")
        ax_fr.set_xlabel("Days Relative to First Session")

    # Create a single legend for all subplots
    handles, labels = ax_z.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(
        0.85, 0.5), title="Stimulation Current")

    plt.suptitle(animal_id)
    # Adjust the right padding to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()


def plot_aggregated_longitudinal_data(df):
    """
    Plots aggregated longitudinal data of Z-scores and Pulse Mean FR across sessions.

    Args:
        df (DataFrame): The aggregated data for all channels and all sessions.
    """

    print("taking abs value of z score")
    # Compute median of absolute z-scores
    aggregated_z_scores = df.groupby(["rel_day", "stim_current"])[
        "z_score"].apply(lambda x: x.abs().mean()).reset_index()

    aggregated_pulse_fr = df.groupby(["rel_day", "stim_current"])[
        "pulse_mean_fr"].mean().reset_index()

    # Get unique stimulation currents
    stim_currents = [2, 3, 4, 5, 6, 7]

    # Create a figure with 2 subplots: one for Z-scores and one for Pulse Mean FR
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

    # Plot Z-Scores
    ax_z = axes[0]
    for stim_current in stim_currents:
        current_data_z = aggregated_z_scores[aggregated_z_scores["stim_current"] == stim_current]

        # Plot line
        ax_z.plot(
            current_data_z["rel_day"],
            current_data_z["z_score"],
            label=f"{stim_current} µA"
        )

        # Plot individual points
        ax_z.scatter(
            current_data_z["rel_day"],
            current_data_z["z_score"],
            marker="o", s=50, zorder=5  # 's' for marker size, 'zorder' brings points to the front
        )

    ax_z.set_title("Aggregated Z-Scores Across Channels")
    ax_z.set_ylabel("Median Z-Score")
    ax_z.set_xlabel("Days Relative to First Session")

    # Plot Pulse Mean FR
    ax_fr = axes[1]
    for stim_current in stim_currents:
        current_data_fr = aggregated_pulse_fr[aggregated_pulse_fr["stim_current"] == stim_current]

        # Plot line
        ax_fr.plot(
            current_data_fr["rel_day"],
            current_data_fr["pulse_mean_fr"],
            label=f"{stim_current} µA"
        )

        # Plot individual points
        ax_fr.scatter(
            current_data_fr["rel_day"],
            current_data_fr["pulse_mean_fr"],
            marker="o", s=50, zorder=5  # 's' for marker size, 'zorder' brings points to the front
        )

    ax_fr.set_title("Aggregated Pulse Mean FR Across Channels")
    ax_fr.set_ylabel("Median Pulse FR (Hz)")
    ax_fr.set_xlabel("Days Relative to First Session")

    # Create a single legend for both subplots
    handles, labels = ax_z.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right",
               title="Stimulation Current")

    # Adjust the layout and add a title
    plt.suptitle("Aggregated Longitudinal Data")
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()


def plot_longitudinal_all_animals():
    animal_ids = ["ICMS92", "ICMS93", "ICMS98",
                  "ICMS100", "ICMS101"]  # Example animal IDs
    base_path = Path("C:/data/")  # Replace with your actual base path

    # Load session data for each animal
    dfs = [load_session_data(animal_id, base_path) for animal_id in animal_ids]

    # Concatenate the DataFrames for all animals
    df = pd.concat(dfs, ignore_index=True)

    # Plot the aggregated longitudinal data
    plot_aggregated_longitudinal_data(df)


def plot_prestim_vs_stim_firing_rates_for_cell_type(df, animal_id, y_value='train_mean_fr', margin=0.05, stim_currents=[3, 4, 5, 6]):
    # Filter dataframe to only include the specified stim currents
    df = df[df['stim_current'].isin(stim_currents)]

    # Calculate global min and max for consistent axis limits, with margin adjustment
    global_min = min(df["pre_stim_mean_fr"].min(), df[y_value].min())
    global_max = max(df["pre_stim_mean_fr"].max(), df[y_value].max())

    # Adjust global_min and global_max to add a margin
    global_min = global_min - (global_max - global_min) * margin
    global_max = global_max + (global_max - global_min) * margin

    grouped = df.groupby("cell_type")
    filenames = []

    for i, (cell_type, group) in enumerate(grouped):
        # Drop NaNs
        group = group.dropna(subset=["pre_stim_mean_fr", y_value])

        # Create a jointplot
        g = sns.jointplot(
            data=group,
            x="pre_stim_mean_fr",
            y=y_value,
            hue="stim_current",
            palette="viridis",
            kind="scatter",
            joint_kws={'s': 20, 'alpha': 0.8}
        )

        # Set consistent axis limits
        g.ax_joint.set_xlim(global_min, 40)
        g.ax_joint.set_ylim(global_min, 80)

        # Plot the diagonal line (x=y)
        g.ax_joint.plot([global_min, global_max], [global_min,
                        global_max], linestyle="--", color="gray")

        # Set titles and labels
        g.set_axis_labels("Pre-Stim Mean FR", y_value)
        g.fig.suptitle(f"{cell_type}", y=0.9)

        # Save the figure
        filename = f"celltype_{cell_type}.png"
        g.savefig(filename)
        filenames.append(filename)
        plt.close(g.fig)

    # Combine the saved plots into a single figure
    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, len(filenames))

    for i, filename in enumerate(filenames):
        img = plt.imread(filename)
        ax = fig.add_subplot(gs[i])
        ax.imshow(img)
        ax.axis("off")  # Hide the axes

    plt.suptitle(animal_id)
    plt.tight_layout()
    plt.savefig(f"{animal_id}_combined.png")
    plt.show()


def plot_prestim_vs_stim_firing_rates_for_cell_type_all_animals(y_value="train_mean_fr"):
    animal_id1 = "ICMS92"  # Replace with your animal ID
    animal_id2 = "ICMS93"  # Replace with your second animal ID
    animal_id3 = "ICMS98"  # Replace with your third animal ID
    animal_id4 = "ICMS100"  # Replace with your third animal ID
    animal_id5 = "ICMS101"  # Replace with your third animal ID

    base_path = Path("C:/data/")  # Replace with your base path

    df1 = load_session_data(animal_id1, base_path)
    df2 = load_session_data(animal_id2, base_path)
    df3 = load_session_data(animal_id3, base_path)
    df4 = load_session_data(animal_id4, base_path)
    df5 = load_session_data(animal_id5, base_path)

    # Concatenate the DataFrames
    df = pd.concat([df1, df2, df3, df4, df5])
    df = df.reset_index(drop=True)
    df = filter_data(df)

    plot_prestim_vs_stim_firing_rates_for_cell_type(df, "Aggregated", y_value)


def plot_prestim_vs_stim_firing_rates_for_cell_type_and_other_var(df, animal_id, var="stim_current", margin=0.05):
    # Calculate global min and max for consistent axis limits, with margin adjustment
    global_min = min(df["pre_stim_mean_fr"].min(), df["train_mean_fr"].min())
    global_max = max(df["pre_stim_mean_fr"].max(), df["train_mean_fr"].max())

    # Adjust global_min to add a margin
    global_min = global_min - (global_max - global_min) * margin
    global_max = global_max + (global_max - global_min) * margin

    # Filter the dataframe to include only stim currents 3, 4, 5, 6
    if var == "stim_current":
        df_filtered = df[df["stim_current"].isin([3, 4, 5, 6])]
    elif var == "detection_threshold":
        df_filtered = df[df["detection_threshold"].isin([2, 3, 4, 5])]

    # Get unique stim currents (should only be 3, 4, 5, 6) and cell types
    unique_var_values = sorted(df_filtered[var].unique())
    unique_cell_types = df_filtered["cell_type"].unique()

    # Number of rows = number of unique stim currents (3, 4, 5, 6)
    num_rows = len(unique_var_values)
    num_cols = len(unique_cell_types)

    # Create a grid of subplots dynamically based on the available stim currents and cell types
    fig, axes = plt.subplots(
        nrows=num_rows, ncols=num_cols, figsize=(4 * num_cols, 2 * num_rows))
    fig.tight_layout(pad=5.0)  # Add some padding between plots

    # Group the dataframe by both `cell_type` and `stim_current`
    grouped = df_filtered.groupby(["cell_type", var])

    for (cell_type, var_i), group in grouped:
        # Drop NaN values
        group = group.dropna(subset=["pre_stim_mean_fr", "train_mean_fr"])

        # Get the row and column index for the subplot
        row_idx = unique_var_values.index(var_i)
        col_idx = list(unique_cell_types).index(cell_type)

        ax = axes[row_idx, col_idx]  # Select the correct subplot

        # Create scatter plot directly on the subplot `ax`
        sns.scatterplot(
            data=group,
            x="pre_stim_mean_fr",
            y="train_mean_fr",
            ax=ax,
            palette="viridis",
            s=10,  # Size of the scatter points
            alpha=0.8,  # Transparency for better visibility
            legend=False
        )

        # Set consistent axis limits
        ax.set_xlim(global_min, 40)
        ax.set_ylim(global_min, 100)

        # Linear regression and plot on the same subplot
        if len(group) >= 2:
            slope, intercept, r_value, p_value, std_err = linregress(
                group["pre_stim_mean_fr"], group["train_mean_fr"])
            x_vals = np.array(ax.get_xlim())
            y_vals = intercept + slope * x_vals
            ax.plot(x_vals, y_vals, linestyle="--", color="red")

            # Add slope and R² to the title of each subplot with regular font size
            ax.set_title(
                f"R² = {r_value**2:.2f}, Slope = {slope:.2f}",
                fontsize=10
            )

        # Plot the diagonal line (x=y)
        ax.plot([global_min, global_max], [global_min, global_max],
                linestyle="--", color="gray")

        # Only set the ylabel (current) for the first column in each row
        if col_idx == 0:
            ax.set_ylabel(f"{int(var_i)} uA")

        # Only set the title (cell type) for the top row, but place it above the entire column with larger font size
        if row_idx == 0:
            ax.set_title(
                f"{cell_type}\nR² = {r_value**2:.2f}, Slope = {slope:.2f}", fontsize=10)

        # Remove repetitive x and y labels except for the bottom-left plot
        if col_idx > 0:
            ax.set_ylabel("")
        if row_idx < num_rows - 1:
            ax.set_xlabel("")
        else:
            if col_idx == num_cols - 1:
                ax.set_xlabel("Pre-Stim Mean FR")
                ax.set_ylabel("Stim Mean FR")
            else:
                ax.set_xlabel("")

    # Save and show the figure
    plt.suptitle(animal_id)
    plt.tight_layout()
    plt.savefig(f"{animal_id}_combined.png")
    plt.show()


def plot_prestim_vs_stim_firing_rates_for_cell_type_and_var_all_animals(var):
    animal_id1 = "ICMS92"  # Replace with your animal ID
    animal_id2 = "ICMS93"  # Replace with your second animal ID
    animal_id3 = "ICMS98"  # Replace with your third animal ID
    animal_id4 = "ICMS100"  # Replace with your third animal ID
    animal_id5 = "ICMS101"  # Replace with your third animal ID

    base_path = Path("C:/data/")  # Replace with your base path

    df1 = load_session_data(animal_id1, base_path)
    df2 = load_session_data(animal_id2, base_path)
    df3 = load_session_data(animal_id3, base_path)
    df4 = load_session_data(animal_id4, base_path)
    df5 = load_session_data(animal_id5, base_path)

    # Concatenate the DataFrames
    df = pd.concat([df1, df2, df3, df4, df5])
    df = df.reset_index(drop=True)

    plot_prestim_vs_stim_firing_rates_for_cell_type_and_other_var(
        df, "Aggregated", var)


def fix_list_format(value):
    """
    Fixes the formatting of a list by splitting the string on whitespace
    and converting it into a Python list.
    """
    if isinstance(value, str):
        # Remove brackets and split the string by whitespace
        value = value.strip("[]").split()
        # Convert the elements to floats
        try:
            value = [float(v) for v in value]
        except ValueError:
            value = np.nan  # If conversion fails, return NaN
    return value


def plot_prestim_vs_stim_firing_rates_by_unit_location(df, animal_id, margin=0.05):
    df["unit_location"] = df["unit_location"].apply(fix_list_format)
    df["unit_location_extracted"] = df["unit_location"].apply(
        lambda x: x[1] if len(x) > 1 else np.nan)
    df["unit_location_extracted"] = pd.to_numeric(
        df["unit_location_extracted"], errors="coerce")
    # Define the bins for unit locations as numeric ranges
    bins = [0, 100, 400, 550, 800, 1000, 1500]
    bin_labels = ["0-100", "100-400", "400-550",
                  "550-800", "800-1000", "1000-1500"]
    df["unit_location_bin"] = pd.cut(
        df["unit_location_extracted"], bins=bins, labels=bin_labels, include_lowest=True)
    global_min = min(df["pre_stim_mean_fr"].min(), df["train_mean_fr"].min())
    global_max = max(df["pre_stim_mean_fr"].max(), df["train_mean_fr"].max())
    global_min = global_min - (global_max - global_min) * margin
    global_max = global_max + (global_max - global_min) * margin
    grouped = df.groupby("unit_location_bin")
    filenames = []
    for i, (location_bin, group) in enumerate(grouped):
        # Create a jointplot with scatter and density plots
        g = sns.jointplot(
            data=group,
            x="pre_stim_mean_fr",
            y="train_mean_fr",
            hue="stim_current",
            palette="viridis",
            kind="scatter",  # or 'kde' if you want contour plots
        )

        # Set the same x and y limits for all plots
        g.ax_joint.set_xlim(global_min, 50)
        g.ax_joint.set_ylim(global_min, 100)

        # Plot the diagonal line (x=y)
        g.ax_joint.plot([global_min, global_max], [global_min,
                        global_max], linestyle="--", color="gray")

        # Set titles and labels
        g.set_axis_labels("Pre-Stim Mean FR", "Train Mean FR")
        g.fig.suptitle(f"Location {location_bin} µm", y=0.9)

        # Save the figure
        filename = f"location_{location_bin}.png"
        g.savefig(filename)
        filenames.append(filename)
        plt.close(g.fig)

    # Combine the individual plots into one figure
    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(2, int(len(filenames) / 2))

    for i, filename in enumerate(filenames):
        img = plt.imread(filename)
        ax = fig.add_subplot(gs[i])
        ax.imshow(img)
        ax.axis("off")  # Hide the axes

    plt.suptitle(f"Animal ID: {animal_id}")
    plt.tight_layout()
    plt.savefig(f"{animal_id}_combined_locations.png")
    plt.show()


def plot_prestim_vs_stim_firing_rates_by_unit_location_all_animals():
    animal_id1 = "ICMS92"  # Replace with your animal ID
    animal_id2 = "ICMS93"  # Replace with your second animal ID
    animal_id3 = "ICMS98"  # Replace with your third animal ID
    animal_id4 = "ICMS100"  # Replace with your third animal ID
    animal_id5 = "ICMS101"  # Replace with your third animal ID

    base_path = Path("C:/data/")  # Replace with your base path

    df1 = load_session_data(animal_id1, base_path)
    df2 = load_session_data(animal_id2, base_path)
    df3 = load_session_data(animal_id3, base_path)
    df4 = load_session_data(animal_id4, base_path)
    df5 = load_session_data(animal_id5, base_path)

    # Concatenate the DataFrames
    df = pd.concat([df1, df2, df3, df4, df5])
    df = df.reset_index(drop=True)

    plot_prestim_vs_stim_firing_rates_by_unit_location(
        df, "Aggregated")


def plot_prestim_vs_stim_firing_rates_for_stim_current_and_threshold(df, animal_id, margin=0.05):
    """
    Plots the pre-stim vs stim firing rates for cell types while comparing stim_current and detection_threshold.
    Args:
        df: DataFrame with firing rate data.
        animal_id: ID of the animal.
        margin: Margin to add to axis limits.
    """
    # Calculate global min and max for consistent axis limits, with margin adjustment
    global_min = min(df["pre_stim_mean_fr"].min(), df["train_mean_fr"].min())
    global_max = max(df["pre_stim_mean_fr"].max(), df["train_mean_fr"].max())

    # Adjust global_min to add a margin
    global_min = global_min - (global_max - global_min) * margin
    global_max = global_max + (global_max - global_min) * margin

    # Filter the dataframe to include only stim currents and thresholds
    df_filtered = df[df["stim_current"].isin([3, 4, 5, 6])]
    df_filtered = df_filtered[df_filtered["detection_threshold"].isin([
                                                                      3, 4, 5, 6])]

    # Get unique stim currents and detection thresholds
    unique_stim_currents = sorted(df_filtered["stim_current"].unique())
    unique_detection_thresholds = sorted(
        df_filtered["detection_threshold"].unique())

    # Number of rows = number of unique stim currents
    # Number of columns = number of unique detection thresholds
    num_rows = len(unique_stim_currents)
    num_cols = len(unique_detection_thresholds)

    # Create a grid of subplots dynamically based on the available stim currents and detection thresholds
    fig, axes = plt.subplots(
        nrows=num_rows, ncols=num_cols, figsize=(2 * num_cols, 1.5 * num_rows))
    fig.tight_layout(pad=1.0)  # Add some padding between plots

    # Group the dataframe by `stim_current` and `detection_threshold`
    grouped = df_filtered.groupby(["stim_current", "detection_threshold"])

    for (stim_current, detection_threshold), group in grouped:
        # Drop NaN values
        group = group.dropna(subset=["pre_stim_mean_fr", "train_mean_fr"])

        # Get the row and column index for the subplot
        row_idx = unique_stim_currents.index(stim_current)
        col_idx = unique_detection_thresholds.index(detection_threshold)

        ax = axes[row_idx, col_idx]  # Select the correct subplot

        # Create scatter plot directly on the subplot `ax`
        sns.scatterplot(
            data=group,
            x="pre_stim_mean_fr",
            y="train_mean_fr",
            ax=ax,
            s=10,  # Size of the scatter points
            alpha=0.8,  # Transparency for better visibility
            legend=False
        )

        # Set consistent axis limits
        ax.set_xlim(global_min, 25)
        ax.set_ylim(global_min, 80)

        # Plot the diagonal line (x=y)
        ax.plot([global_min, global_max], [global_min, global_max],
                linestyle="--", color="gray")

        # Set stim current label (y-axis) for every row
        if col_idx == 0:
            ax.set_ylabel(f"Current: {stim_current}")
            ax.set_xlabel("")
        elif (col_idx == num_cols - 1) and (row_idx == num_rows - 1):
            ax.set_ylabel("Stim Mean FR")
            ax.set_xlabel("Pre-Stim Mean FR")
        else:
            ax.set_ylabel("")
            ax.set_xlabel("")

    # Add detection thresholds as column titles above the plots
    for col_idx, detection_threshold in enumerate(unique_detection_thresholds):
        axes[0, col_idx].annotate(f"Threshold: {int(detection_threshold)}",
                                  xy=(0.5, 1.2), xycoords='axes fraction', ha='center', va='center',
                                  fontsize=11)

    # Set super title and save the figure
    # plt.suptitle(f"Animal ID: {animal_id}", fontsize=16)
    # Adjust layout to make space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{animal_id}_combined.png")
    plt.show()


def plot_prestim_vs_stim_firing_rates_by_stim_current_and_threshold_all_animals():
    animal_id1 = "ICMS92"  # Replace with your animal ID
    animal_id2 = "ICMS93"  # Replace with your second animal ID
    animal_id3 = "ICMS98"  # Replace with your third animal ID
    animal_id4 = "ICMS100"  # Replace with your third animal ID
    animal_id5 = "ICMS101"  # Replace with your third animal ID

    base_path = Path("C:/data/")  # Replace with your base path

    df1 = load_session_data(animal_id1, base_path)
    df2 = load_session_data(animal_id2, base_path)
    df3 = load_session_data(animal_id3, base_path)
    df4 = load_session_data(animal_id4, base_path)
    df5 = load_session_data(animal_id5, base_path)

    # Concatenate the DataFrames
    df = pd.concat([df1, df2, df3, df4, df5])
    df = df.reset_index(drop=True)
    df = filter_data(df)

    plot_prestim_vs_stim_firing_rates_for_stim_current_and_threshold(
        df, "Aggregated")


def get_stim_colormap():
    valid_currents = list(range(4, 7))
    palette = sns.color_palette("deep", len(valid_currents))
    colormap = {curr: palette[i] for i, curr in enumerate(valid_currents)}
    return colormap


def plot_train_rasters(responses, ax):
    color_map = get_stim_colormap()

    all_rasters = []
    all_colors = []
    all_offsets = []

    linewidth_factor = 40
    offset = 0

    # Plot same number of trials for all current amplitudes
    num_trials_arr = []
    for current, response in sorted(responses):
        raster_array = response.train_response.raster_array
        num_trials_arr.append(len(raster_array))
    num_trials = np.min(num_trials_arr)

    for current, response in sorted(responses):
        raster_array = response.train_response.raster_array[0:num_trials]
        line_count = num_trials
        line_height = np.ceil(line_count / linewidth_factor)
        offsets = np.arange(line_count) + offset

        all_rasters.extend(raster_array)
        all_colors.extend([color_map[current]] * line_count)
        all_offsets.extend(offsets)

        offset += line_count

    ax.eventplot(
        all_rasters,
        orientation="horizontal",
        colors=all_colors,
        linelengths=line_height,
        linewidths=0.4,
        lineoffsets=all_offsets,
        label=f'{current} µA'
    )

    ax.plot([0, 700], [offset + 2, offset + 2], color='black', linewidth=2)
    ax.set_xticks([0, 700])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.axvline(x=0, y=100, color='gray', linestyle='--', linewidth=1)
    # ax.axvline(x=700, color='gray', linestyle='--', linewidth=1)

    # ax.plot([0, 0], [-1, offset], color='gray',
    #         linestyle='--', linewidth=1)
    # ax.plot([700, 700], [-1, offset], color='gray',
    #         linestyle='--', linewidth=1)

    ax.set_xlim([-700, 1400])
    ax.set_ylabel("Trial")


def plot_train_firing_rate(responses, ax):
    color_map = get_stim_colormap()

    for current, response in sorted(responses):
        tr = response.train_response
        tr_data = tr.ten_ms_smoothed_data
        x = tr_data['bin_centers']
        y = tr_data['firing_rate']
        ax.plot(
            x,
            y,
            color=color_map[current],
            label=f"{current} µA",
        )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([0, 700])

    # ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    # ax.axvline(x=700, color='gray', linestyle='--', linewidth=1)
    ax.set_xlim([-700, 1400])
    ax.plot([0, 700], [100 + 2, 100 + 2], color='black', linewidth=2)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Firing rate (Hz)")


def plot_EI_ratio(df):
    """
    Plots the E/I ratio (Excitatory/Inhibitory) across weeks and stimulation currents.

    Args:
        df (DataFrame): A dataframe containing the columns 'rel_week', 'stim_current', 'cell_type'
                        where 'cell_type' can be 'pyramidal' for excitatory or some other type for inhibitory.
    """
    # Get unique weeks and currents
    df['rel_week'] = (df['rel_day'] // 7).astype(int)
    weeks = sorted(df['rel_week'].unique())
    stim_currents = [2, 3, 4, 5, 6, 7]

    # Dictionary to store E/I ratios for each week and current
    EI_ratios = {current: [] for current in stim_currents}

    # Loop over each week
    for week in weeks:
        df_week = df[df['rel_week'] == week]

        # Loop over each current
        for current in stim_currents:
            df_current = df_week[df_week['stim_current'] == current]

            # Initialize counts
            E_count = 0  # Excitatory (Pyramidal)
            I_count = 0  # Inhibitory (Non-pyramidal)

            # Classify cells and count
            for _, row in df_current.iterrows():
                if row['cell_type'] == 'Pyramidal Cell':
                    E_count += 1
                else:
                    I_count += 1

            # Calculate E/I ratio (avoid division by zero)
            if I_count > 0:
                EI_ratio = E_count / I_count
            else:
                EI_ratio = np.nan  # Set to NaN if no inhibitory cells found

            # Append the ratio to the dictionary
            EI_ratios[current].append(EI_ratio)

    # Plotting the E/I ratio
    plt.figure(figsize=(10, 6))

    for current in stim_currents:
        plt.plot(weeks, EI_ratios[current], marker='o', label=f'{current} µA')

    plt.xlabel('Weeks Relative to First Session')
    plt.ylabel('E/I Ratio (Excitatory/Inhibitory)')
    plt.title('E/I Ratio Across Weeks and Stimulation Currents')
    plt.legend(title='Stimulation Current')
    plt.tight_layout()
    plt.show()


def plot_waveforms_by_cell_type_and_week(df):
    """
    Plots unique unit waveforms by cell type and week with appropriate axis labels.

    Args:
        df (DataFrame): DataFrame containing unit information, cell types, templates, and weeks.
    """
    # Filter the data for each cell type
    pyr_df = df[df['cell_type'] == "Pyramidal Cell"]
    wide_int_df = df[df['cell_type'] == "Wide Interneuron"]
    narrow_int_df = df[df['cell_type'] == "Narrow Interneuron"]

    # List of dataframes and corresponding labels for cell types
    cell_type_dfs = [pyr_df, wide_int_df, narrow_int_df]
    cell_type_labels = ["Pyramidal Cell",
                        "Wide Interneuron", "Narrow Interneuron"]

    # Unique weeks (0 to 5)
    weeks = [0, 1, 2, 3, 4, 5]

    # Create a 3x6 grid for subplots
    fig, axes = plt.subplots(3, 6, figsize=(10, 6), sharey=True)

    for row_idx, (cell_df, cell_label) in enumerate(zip(cell_type_dfs, cell_type_labels)):
        for col_idx, week in enumerate(weeks):
            ax = axes[row_idx, col_idx]

            # Filter for the current week and get unique units by "unit_id"
            week_df = cell_df[cell_df['rel_week']
                              == week].drop_duplicates(subset='unit_id')

            # Plot the waveforms for each unique unit in the current week
            for _, row in week_df.iterrows():
                # Assuming 'template' contains the waveform
                waveform = row['template']
                # Time axis in ms (30 samples per ms)
                time_axis = np.linspace(0, len(waveform) / 30, len(waveform))
                ax.plot(time_axis, waveform, alpha=0.7)

            # Set titles and labels
            if col_idx == 0:
                # Manually place the cell type label to the left
                ax.text(-0.15, 0.5, f"{cell_label}", fontsize=8, rotation=90,
                        transform=ax.transAxes, va='center')

            if row_idx == 0:
                ax.set_title(f"Week {week}", fontsize=12)  # Week title on top

            # Set y-axis ticks for amplitude only on the last column
            if col_idx == 5:
                ax.set_yticks([-200, 0, 50])  # Example y-ticks in µV
                ax.yaxis.set_label_position("right")
                ax.set_ylabel("Amplitude (µV)", fontsize=10)
                ax.yaxis.tick_right()

            # Set x-axis labels only for the bottom-left subplot
            if row_idx == 2 and col_idx == 0:
                ax.set_xlabel('Time (ms)')
                ax.set_xticks([0, len(waveform) / 60, len(waveform) / 30])
                ax.set_xticklabels(
                    [0, len(waveform) / 60, len(waveform) // 30])
            else:
                ax.set_xticks([])  # Hide x-ticks for all other subplots

            # Remove y-tick labels for non-leftmost subplots
            if col_idx != 5:
                ax.set_yticklabels([])

    # Add a super title and adjust the layout
    plt.suptitle("Unit Waveforms by Cell Type and Week", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


def compute_centroid_distances(df):
    from scipy.spatial.distance import euclidean

    """
    Computes centroid distances and categorizes them into "within-group" and "between-group,"
    while removing empty distances or self-comparisons.
    """
    # Filter and prepare data
    df_filtered = df[df["stim_current"].isin([3, 4, 5, 6, 7])]
    df_filtered = df_filtered[df_filtered["detection_threshold"].isin([
                                                                      3, 4, 5, 6, 7])]

    unique_conditions = list(df_filtered.groupby(
        ["stim_current", "detection_threshold"]).groups.keys())

    distances = []
    labels = []

    for i, condition_a in enumerate(unique_conditions):
        for j, condition_b in enumerate(unique_conditions):
            # Skip self-comparisons (when comparing the same condition)
            if i < j:  # Only calculate for distinct pairs and avoid duplicates
                df_a = df_filtered[(df_filtered["stim_current"] == condition_a[0]) &
                                   (df_filtered["detection_threshold"] == condition_a[1])]
                df_b = df_filtered[(df_filtered["stim_current"] == condition_b[0]) &
                                   (df_filtered["detection_threshold"] == condition_b[1])]

                # Check if either DataFrame is empty
                if df_a.empty or df_b.empty:
                    continue  # Skip if there's no data for either condition

                # Calculate centroids
                centroid_a = [df_a["pre_stim_mean_fr"].mean(),
                              df_a["train_mean_fr"].mean()]
                centroid_b = [df_b["pre_stim_mean_fr"].mean(),
                              df_b["train_mean_fr"].mean()]

                # Skip if centroids contain NaN values (i.e., if no valid data for centroid calculation)
                if np.isnan(centroid_a).any() or np.isnan(centroid_b).any():
                    continue  # Skip invalid centroids

                # Compute Euclidean distance between the centroids
                centroid_distance = euclidean(centroid_a, centroid_b)

                # Check if this is a "within-group" or "between-group" comparison
                if condition_a[0] == condition_a[1] and condition_b[0] == condition_b[1]:
                    category = 'Within-Group'
                else:
                    category = 'Between-Group'

                distances.append(centroid_distance)
                labels.append(category)

    return distances, labels


def plot_similarity_comparison(distances, labels):
    """
    Plots a swarm plot comparing the within-group and between-group similarities.
    """
    data = pd.DataFrame({'Distance': distances, 'Category': labels})

    plt.figure(figsize=(10, 6))
    sns.swarmplot(x='Category', y='Distance',
                  data=data, palette="Set2", size=5)
    plt.title(
        "Procrustes Distance Comparison Between Within-Group and Between-Group")
    plt.ylabel("Procrustes Distance")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Category', y='Distance', data=data,
                   inner=None, palette="Set2", alpha=0.5)
    sns.swarmplot(x='Category', y='Distance', data=data,
                  color='k', alpha=0.7, size=5)
    plt.title(
        "Procrustes Distance Comparison Between Within-Group and Between-Group")
    plt.ylabel("Procrustes Distance")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.show()


def plot_centroid_comparison(distances, labels):
    """
    Plots a swarm plot comparing within-group and between-group centroid distances.
    """
    data = pd.DataFrame({'Distance': distances, 'Category': labels})

    plt.figure(figsize=(10, 6))
    sns.swarmplot(x='Category', y='Distance',
                  data=data, palette="Set2", size=5)
    plt.title("Centroid Distance Comparison Between Within-Group and Between-Group")
    plt.ylabel("Centroid Distance (Euclidean)")
    plt.xlabel("Group")
    plt.tight_layout()
    plt.show()


def plot_aggregated_weekly_data(df, var1="z_score", var2="pulse_mean_fr", animal_id="Aggregated", use_median=False):
    """
    Plots aggregated longitudinal data of two selected variables (e.g., Z-Scores and Pulse Mean FR) across weeks
    with error bars to represent either the standard error of the mean (SEM) or interquartile range (IQR).

    Args:
        df (DataFrame): The aggregated data for all channels and all sessions.
        var1 (str): The first variable to plot (default is 'z_score').
        var2 (str): The second variable to plot (default is 'pulse_mean_fr').
        animal_id (str): The identifier for the animal (default is 'Aggregated').
        use_median (bool): If True, plots median with IQR. If False, plots mean with SEM.
    """

    # Create a new column for 'rel_week' by dividing 'rel_day' by 7 and flooring it
    df['rel_week'] = (df['rel_day'] // 7).astype(int)  # TODO
    # df['rel_week'] = df['rel_day']

    # Absolute values of var1 and var2
    df['abs_var1'] = df[var1].abs()
    df['abs_var2'] = df[var2].abs()

    if use_median:
        # Compute aggregated data using median and IQR
        aggregated_data_var1 = df.groupby(['rel_week', 'stim_current']).agg(
            median=('abs_var1', 'median'),
            q1=('abs_var1', lambda x: x.quantile(0.25)),
            q3=('abs_var1', lambda x: x.quantile(0.75))
        ).reset_index()

        aggregated_data_var2 = df.groupby(['rel_week', 'stim_current']).agg(
            median=('abs_var2', 'median'),
            q1=('abs_var2', lambda x: x.quantile(0.25)),
            q3=('abs_var2', lambda x: x.quantile(0.75))
        ).reset_index()
    else:
        # Compute aggregated data using mean and SEM
        aggregated_data_var1 = df.groupby(['rel_week', 'stim_current'])[
            'abs_var1'].agg(['mean', 'sem']).reset_index()
        aggregated_data_var2 = df.groupby(['rel_week', 'stim_current'])[
            'abs_var2'].agg(['mean', 'sem']).reset_index()

    # Get unique stimulation currents
    stim_currents = sorted(df['stim_current'].unique())

    # Create a figure with 2 subplots: one for var1 and one for var2
    fig, axes = plt.subplots(1, 2, figsize=(7, 4), sharex=True)

    # Plot for var1 with bars (error bars)
    ax1 = axes[0]
    for stim_current in stim_currents:
        current_data_var1 = aggregated_data_var1[aggregated_data_var1['stim_current'] == stim_current]
        if use_median:
            # Plot with median and IQR as error bars
            ax1.errorbar(
                current_data_var1['rel_week'],
                current_data_var1['median'],
                yerr=[current_data_var1['median'] - current_data_var1['q1'],
                      current_data_var1['q3'] - current_data_var1['median']],
                fmt='-o',
                label=f"{stim_current} µA"
            )
        else:
            # Plot with mean and SEM as error bars
            ax1.errorbar(
                current_data_var1['rel_week'],
                current_data_var1['mean'],
                yerr=current_data_var1['sem'],
                fmt='-o',
                label=f"{stim_current} µA"
            )

    # Replace underscores with spaces for better formatting
    var1_label = var1.replace('_', ' ').capitalize()
    ax1.set_title(f"{var1_label.capitalize()}")
    ax1.set_ylabel(f"{'Median' if use_median else 'Mean'} {var1_label}")
    ax1.set_xlabel("Week")
    ax1.set_xticks(sorted(df['rel_week'].unique()))

    # Plot for var2 with bars (error bars)
    ax2 = axes[1]
    for stim_current in stim_currents:
        current_data_var2 = aggregated_data_var2[aggregated_data_var2['stim_current'] == stim_current]
        if use_median:
            # Plot with median and IQR as error bars
            ax2.errorbar(
                current_data_var2['rel_week'],
                current_data_var2['median'],
                yerr=[current_data_var2['median'] - current_data_var2['q1'],
                      current_data_var2['q3'] - current_data_var2['median']],
                fmt='-o',
                label=f"{stim_current} µA"
            )
        else:
            # Plot with mean and SEM as error bars
            ax2.errorbar(
                current_data_var2['rel_week'],
                current_data_var2['mean'],
                yerr=current_data_var2['sem'],
                fmt='-o',
                label=f"{stim_current} µA",
                capsize=3, elinewidth=0.5, capthick=0.5
            )

    # Replace underscores with spaces for better formatting
    var2_label = var2.replace('_', ' ').capitalize()

    ax2.set_title(f"Pulse Window FR" if var2 ==
                  "pulse_mean_fr" else f"{var2_label}")
    ax2.set_ylabel(f"{'Median' if use_median else 'Mean'} {var2_label}")
    ax2.set_xlabel("Week")
    ax2.set_xticks(sorted(df['rel_week'].unique()))

    # Create a single legend for both subplots
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles, labels, loc="upper right",
               title="Stimulation Current", fontsize=6)

    # Adjust the layout and show the plot
    plt.suptitle(f"{animal_id} Longitudinal Data")
    # plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.tight_layout()
    plt.show()
    return axes


def plot_single_metric_two_group_comparison(df, var="z_score", alt='less', ax=None, y_clip=None):
    """
    Plot metric over weeks with early vs late comparison.

    Args:
        y_clip: If set, clip scatter points and set y-axis max to this value.
            Stats (median, IQR, p-values) are computed on ALL data.
            Only the plot presentation is clipped.
    """
    from scipy.stats import ttest_ind, mannwhitneyu, shapiro

    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {stim: i * offset - ((len(stim_currents) - 1) / 2)
                  * offset for i, stim in enumerate(stim_currents)}
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {stim: colors[i] for i, stim in enumerate(stim_currents)}

    aggregated_data_var = df.groupby(['rel_week', 'stim_current']).agg(
        median=(var, 'median'),
        q1=(var, lambda x: x.quantile(0.25)),
        q3=(var, lambda x: x.quantile(0.75))
    ).reset_index()

    aggregated_data_var['week_group_offset'] = aggregated_data_var['rel_week'] + \
        aggregated_data_var['stim_current'].map(offset_map)

    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 2.5))

    all_y_vals = []
    stars_arr = {}
    for stim_current in stim_currents:
        current_data_var = aggregated_data_var[aggregated_data_var['stim_current'] == stim_current]

        # Individual data points with jitter
        cur_df = df[df['stim_current'] == stim_current]
        for week in cur_df['rel_week'].unique():
            week_vals = cur_df[cur_df['rel_week'] == week][var].dropna()
            x_base = week + offset_map[stim_current]
            jitter = np.random.uniform(-0.05, 0.05, size=len(week_vals))
            ax.scatter(
                x_base + jitter, week_vals,
                s=1.5, alpha=0.15, color=color_map[stim_current],
                edgecolors='none', rasterized=True, zorder=1)

        # plot median + IQR on top
        ax.errorbar(
            current_data_var['week_group_offset'],
            current_data_var['median'],
            yerr=[current_data_var['median'] - current_data_var['q1'],
                  current_data_var['q3'] - current_data_var['median']],
            fmt='-o',
            label=f"{stim_current} \u00b5A",
            color=color_map[stim_current],
            ms=3,
            capsize=3, elinewidth=0.5, capthick=0.5, zorder=2
        )
        all_y_vals.extend(current_data_var['q3'].values)

        # Compare early vs late groups (0–1 vs 2–4)
        early_vals = df[(df['rel_week'].isin([0, 1])) & (
            df['stim_current'] == stim_current)][var].dropna()
        late_vals = df[(df['rel_week'].isin([2, 3, 4])) & (
            df['stim_current'] == stim_current)][var].dropna()

        # We, pe = shapiro(early_vals)
        # Wl, pl = shapiro(late_vals)
        # print(f"Shapiro-Wilk: W={We:.3f}, p={pe:.3g}")
        # print(f"Shapiro-Wilk: W={Wl:.3f}, p={pe:.3g}")

        stat, p_val = mannwhitneyu(
            early_vals, late_vals, alternative=alt)
        print(f"{stim_current}µA: n_early={len(early_vals)}, n_late={len(late_vals)}")
        print(f"{stim_current}µA p = {p_val:.1e}")

        if p_val <= 0.001:
            stars = '***'
        elif p_val <= 0.01:
            stars = '**'
        elif p_val <= 0.05:
            stars = '*'
        else:
            stars = 'NS'

        stars_arr[stim_current] = stars

    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels([0, 1, 2, 3, 4])
    ax.set_ylabel(var.replace('_', ' ').capitalize())

    from matplotlib.patches import Patch
    handles, labels = ax.get_legend_handles_labels()
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    handles = [Patch(facecolor=c, edgecolor='none') for c in colors]
    # ax.legend(handles, labels, title="Current", loc='upper left')
    ax.legend(handles, labels, frameon=False,
              loc='center left',          # anchor legend relative to the axes
              bbox_to_anchor=(1.02, 0.2),  # just outside the right edge
              borderaxespad=0.0
              )
    ax.set_xlabel('Weeks of training')

    if y_clip is not None:
        ax.set_ylim(top=y_clip * 1.05)

    return ax, stars_arr, all_y_vals


def draw_two_group_error_bars(stars_arr, all_y_vals, star_size=8, y_scale_factor=1.05,
                              bar_y_offset=0, stars_y_offset=0, stars_spacing=0.05, ax=None):
    stim_currents = [4, 5, 6]
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {stim: colors[i] for i, stim in enumerate(stim_currents)}

    if not all_y_vals:
        return ax

    # Use second-highest q3 value if max is an outlier (>2x the median q3)
    sorted_vals = sorted(all_y_vals)
    max_val = sorted_vals[-1]
    median_val = sorted_vals[len(sorted_vals) // 2]
    if len(sorted_vals) > 2 and max_val > 2 * median_val:
        y_bar = sorted_vals[-2] * y_scale_factor + bar_y_offset
    else:
        y_bar = max_val * y_scale_factor + bar_y_offset

    # Get axis data range for dynamic spacing
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min if y_max > y_min else max(all_y_vals) - min(all_y_vals)
    if y_range <= 0:
        y_range = 1

    # Draw bracket
    bracket_h = y_range * 0.03
    ax.plot([0, 1], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([2, 4], [y_bar, y_bar], color='black', lw=0.5)
    ax.plot([0.5, 0.5], [y_bar, y_bar + bracket_h], color='black', lw=0.5)
    ax.plot([3, 3], [y_bar, y_bar + bracket_h], color='black', lw=0.5)
    ax.plot([0.5, 3], [y_bar + bracket_h, y_bar + bracket_h], color='black', lw=0.5)

    # Dynamic spacing: scale with y_range
    base_y = y_bar + bracket_h + y_range * 0.02 + stars_y_offset
    line_spacing = y_range * 0.06  # space between lines

    top_y = base_y
    for i, curr in enumerate([4, 5, 6]):
        if stars_arr.get(curr):
            star_text = stars_arr[curr]
            is_ns = star_text.upper() in ['NS', 'N.S.']
            fontsize = 6 if is_ns else star_size
            extra = y_range * 0.01 if is_ns else 0
            text_y = base_y + i * line_spacing + extra
            ax.text(
                1.5,
                text_y,
                star_text,
                ha='center',
                fontsize=fontsize,
                color=color_map[curr]
            )
            top_y = max(top_y, text_y)

    # Extend ylim if stars would overlap with title
    current_ylim = ax.get_ylim()
    needed_top = top_y + y_range * 0.1  # 10% padding above top star
    if needed_top > current_ylim[1]:
        ax.set_ylim(current_ylim[0], needed_top)

    return ax


def plot_num_responses(df, animal_id="Aggregated", count_pulse_locked=True, normalize=False, normalize_within_current=False):
    """
    Plots the number of responses per week, grouped by stimulation current.

    Args:
        df (pd.DataFrame): Dataframe containing 'rel_day', 'stim_current', and 'is_pulse_locked'.
        animal_id (str, optional): Name of the animal (default: "Aggregated").
        count_pulse_locked (bool, optional): If True, counts pulse-locked responses. If False, counts non-pulse-locked responses.
        normalize (bool, optional): If True, normalizes by total number of responses (PL + NPL) per week.
        normalize_within_current (bool, optional): If True, normalizes within each stimulation current.

    Returns:
        None
    """

    # Ensure the dataframe contains both pulse-locked and non-pulse-locked responses
    unique_pulse_locked = df['is_pulse_locked'].unique()
    if not (True in unique_pulse_locked and False in unique_pulse_locked):
        print("Warning: Dataframe does not contain both pulse-locked and non-pulse-locked responses.")

    # Create a new column for 'rel_week' by dividing 'rel_day' by 7 and flooring it
    df['rel_week'] = (df['rel_day'] // 7).astype(int)

    # Filter based on pulse-locked or non-pulse-locked responses
    df_filtered = df[df['is_pulse_locked'] == count_pulse_locked]

    # Count responses per week and stimulation current
    aggregated_data = df_filtered.groupby(
        ['rel_week', 'stim_current']).size().reset_index(name='count')

    # Get unique weeks and unique stimulation currents
    all_weeks = sorted(df['rel_week'].unique())
    stim_currents = sorted(df['stim_current'].unique())

    # Normalize if required
    if normalize:
        if normalize_within_current:
            # Normalize within each stimulation current
            total_responses_per_current_week = df.groupby(
                ['rel_week', 'stim_current']).size().to_dict()
            aggregated_data['proportion'] = aggregated_data.apply(
                lambda row: row['count'] / total_responses_per_current_week.get(
                    (row['rel_week'], row['stim_current']), 1),
                axis=1
            )
        else:
            # Normalize across all responses (PL + NPL)
            total_responses_per_week = df.groupby(
                'rel_week').size().to_dict()
            aggregated_data['proportion'] = aggregated_data['rel_week'].map(
                total_responses_per_week)
            aggregated_data['proportion'] = aggregated_data['count'] / \
                aggregated_data['proportion']

    # Create a figure for the plot
    plt.figure(figsize=(6, 4))

    # Plot for each stimulation current
    for stim_current in stim_currents:
        current_data = aggregated_data[aggregated_data['stim_current']
                                       == stim_current]
        y_values = current_data['proportion'] if normalize else current_data['count']

        plt.plot(
            current_data['rel_week'],
            y_values,
            marker='o', linestyle='-', label=f"{stim_current} µA"
        )

    # Formatting
    response_type = "Pulse-Locked" if count_pulse_locked else "Non-Pulse-Locked"
    ylabel = "Proportion of Responses" if normalize else "Number of Responses"
    if normalize and normalize_within_current:
        ylabel = "Proportion (within Stimulation Current)"

    plt.title(f"{animal_id} - {response_type} Responses Over Weeks")
    plt.xlabel("Week")
    plt.ylabel(ylabel)
    plt.xticks(all_weeks)
    plt.legend(title="Stimulation Current", fontsize=8)

    # Show the plot
    plt.tight_layout()
    plt.show()


def plot_aggregated_longitudinal_data(df, var1="z_score", var2="pulse_mean_fr", animal_id="Aggregated", use_median=False, group_by_week=True):
    """
    Plots aggregated longitudinal data of two selected variables (e.g., Z-Scores and Pulse Mean FR) across sessions or weeks
    with error bars to represent either the standard error of the mean (SEM) or interquartile range (IQR).
    """

    # Determine whether to group by weeks or sessions
    if group_by_week:
        df['time_group'] = (df['rel_day'] // 7).astype(int)  # Weeks
        xlabel = "Week"
    else:
        unique_sessions_sorted = df[['session', 'rel_day']].drop_duplicates(
        ).sort_values('rel_day')
        session_index_mapping = {session: i for i, session in enumerate(
            unique_sessions_sorted['session'].unique(), start=1)}
        df['time_group'] = df['session'].map(session_index_mapping)
        xlabel = "Session Index"

    # Take absolute values of var1 and var2
    df['abs_var1'] = df[var1].abs()
    df['abs_var2'] = df[var2].abs()

    # Define a small horizontal shift to separate different currents
    jitter_strength = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {stim: i * jitter_strength - ((len(stim_currents) - 1) / 2) * jitter_strength
                  for i, stim in enumerate(stim_currents)}

    # Aggregate data
    if use_median:
        aggregated_data_var1 = df.groupby(['time_group', 'stim_current']).agg(
            median=('abs_var1', 'median'),
            q1=('abs_var1', lambda x: x.quantile(0.25)),
            q3=('abs_var1', lambda x: x.quantile(0.75))
        ).reset_index()
        aggregated_data_var2 = df.groupby(['time_group', 'stim_current']).agg(
            median=('abs_var2', 'median'),
            q1=('abs_var2', lambda x: x.quantile(0.25)),
            q3=('abs_var2', lambda x: x.quantile(0.75))
        ).reset_index()
    else:
        aggregated_data_var1 = df.groupby(['time_group', 'stim_current'])[
            'abs_var1'].agg(['mean', 'sem', 'std']).reset_index()
        aggregated_data_var2 = df.groupby(['time_group', 'stim_current'])[
            'abs_var2'].agg(['mean', 'sem', 'std']).reset_index()

    # Apply x-offset to prevent overlap
    aggregated_data_var1['time_group_offset'] = aggregated_data_var1['time_group'] + \
        aggregated_data_var1['stim_current'].map(offset_map)
    aggregated_data_var2['time_group_offset'] = aggregated_data_var2['time_group'] + \
        aggregated_data_var2['stim_current'].map(offset_map)

    # Create a figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharex=True)

    # Plot for var1 with error bars
    ax1 = axes[0]
    for stim_current in stim_currents:
        current_data_var1 = aggregated_data_var1[aggregated_data_var1['stim_current'] == stim_current]
        if use_median:
            ax1.errorbar(
                current_data_var1['time_group_offset'],
                current_data_var1['median'],
                yerr=[current_data_var1['median'] - current_data_var1['q1'],
                      current_data_var1['q3'] - current_data_var1['median']],
                fmt='-o',
                label=f"{stim_current} µA"
            )
        else:
            ax1.errorbar(
                current_data_var1['time_group_offset'],
                current_data_var1['mean'],
                yerr=current_data_var1['std'],
                fmt='-o',
                label=f"{stim_current} µA"
            )

    # Format first subplot
    var1_label = var1.replace('_', ' ').capitalize()
    ax1.set_title(f"{var1_label}")
    ax1.set_ylabel(f"{'Median' if use_median else 'Mean'} {var1_label}")
    ax1.set_xlabel(xlabel)
    ax1.set_xticks(sorted(df['time_group'].unique()))

    # Plot for var2 with error bars
    ax2 = axes[1]
    for stim_current in stim_currents:
        current_data_var2 = aggregated_data_var2[aggregated_data_var2['stim_current'] == stim_current]
        if use_median:
            ax2.errorbar(
                current_data_var2['time_group_offset'],
                current_data_var2['median'],
                yerr=[current_data_var2['median'] - current_data_var2['q1'],
                      current_data_var2['q3'] - current_data_var2['median']],
                fmt='-o',
                label=f"{stim_current} µA"
            )
        else:
            ax2.errorbar(
                current_data_var2['time_group_offset'],
                current_data_var2['mean'],
                yerr=current_data_var2['std'],
                fmt='-o',
                label=f"{stim_current} µA"
            )

    # Format second subplot
    var2_label = var2.replace('_', ' ').capitalize()
    ax2.set_title(f"{var2_label}")
    ax2.set_ylabel(f"{'Median' if use_median else 'Mean'} {var2_label}")
    ax2.set_xlabel(xlabel)
    ax2.set_xticks(sorted(df['time_group'].unique()))

    # Create a single legend for both subplots
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles, labels, loc="upper right",
               title="Stimulation Current")

    # Adjust layout and show the plot
    plt.suptitle(
        f"{animal_id} Longitudinal Data({'Weekly' if group_by_week else 'Session-wise'})")
    plt.tight_layout()
    plt.show()
    return axes


def plot_single_metric_aggregated(df, var="z_score", animal_id="Aggregated", comparison_alternative='less', use_median=True, ax=None):

    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {stim: i * offset - ((len(stim_currents) - 1) / 2)
                  * offset for i, stim in enumerate(stim_currents)}

    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {stim: colors[i] for i, stim in enumerate(stim_currents)}

    if use_median:
        aggregated_data_var = df.groupby(['rel_week', 'stim_current']).agg(
            median=(var, 'median'),
            q1=(var, lambda x: x.quantile(0.25)),
            q3=(var, lambda x: x.quantile(0.75))
        ).reset_index()
    else:
        aggregated_data_var = df.groupby(['rel_week', 'stim_current']).agg(
            mean=(var, 'mean'),
            # Compute SEM manually
            sem=(var, lambda x: x.std() / np.sqrt(len(x))),
            std=(var, 'std')
        ).reset_index()

    aggregated_data_var['week_group_offset'] = aggregated_data_var['rel_week'] + \
        aggregated_data_var['stim_current'].map(offset_map)

    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 2))

    for stim_current in stim_currents:

        current_data_var = aggregated_data_var[aggregated_data_var['stim_current'] == stim_current]
        # Mann-Whitney U test between Week 0 and Week 4
        week_0_var = df[(df['rel_week'] == 0) & (
            df['stim_current'] == stim_current)][var]
        week_last_var = df[(df['rel_week'] == 4) & (
            df['stim_current'] == stim_current)][var]

        stat, p_val = mannwhitneyu(
            week_0_var, week_last_var, alternative=comparison_alternative)
        print(p_val)

        if p_val <= 0.05:
            # Determine number of stars based on significance level
            if p_val <= 0.001:
                stars = "***"
            elif p_val <= 0.01:
                stars = "**"
            else:
                stars = "*"
            subset = aggregated_data_var[(aggregated_data_var['rel_week'] == 4) &
                                         (aggregated_data_var['stim_current'] == stim_current)]

            last_y = subset['median'].iloc[0] if use_median else subset['mean'].iloc[0]
            upper_error_bar_len = (
                subset['q3'].iloc[0] - last_y) if use_median else subset['sem'].iloc[0]

            x_pos = current_data_var['week_group_offset'].iloc[-1]
            y_pos = last_y + upper_error_bar_len * 1.03  # Adjust y-position

            ax.text(x_pos, y_pos, stars, color='k', fontsize=10, ha='center')

        if use_median:
            ax.errorbar(
                current_data_var['week_group_offset'],
                current_data_var['median'],
                yerr=[current_data_var['median'] - current_data_var['q1'],
                      current_data_var['q3'] - current_data_var['median']],
                fmt='-o',
                label=f"{stim_current} µA",
                color=color_map[stim_current],
                ms=4
            )
        else:
            ax.errorbar(
                current_data_var['week_group_offset'],
                current_data_var['mean'],
                yerr=current_data_var['sem'],
                fmt='-o',
                label=f"{stim_current} µA",
                color=color_map[stim_current],
                ms=4
            )

    var_label = var.replace('_', ' ').capitalize()
    ax.set_xticks(sorted(df['rel_week'].unique()))

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title="Current", loc='upper left')

    plt.tight_layout()
    plt.show()
    return ax


def plot_single_metric_aggregated_compare_two_time_groups(df, var="z_score", animal_id="Aggregated",
                                                          time_group1=[0, 1], time_group2=[3, 4],
                                                          use_median=True, alt='less', ax=None):

    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {stim: i * offset - ((len(stim_currents) - 1) / 2)
                  * offset for i, stim in enumerate(stim_currents)}

    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {stim: colors[i] for i, stim in enumerate(stim_currents)}

    # Aggregate data by time group
    df['time_group'] = np.where(df['rel_week'].isin(time_group1), "Time Group 1",
                                np.where(df['rel_week'].isin(time_group2), "Time Group 2", np.nan))
    df = df.dropna(subset=["time_group"])  # Keep only relevant time groups

    if use_median:
        aggregated_data_var = df.groupby(['time_group', 'stim_current']).agg(
            median=(var, 'median'),
            q1=(var, lambda x: x.quantile(0.25)),
            q3=(var, lambda x: x.quantile(0.75))
        ).reset_index()
    else:
        aggregated_data_var = df.groupby(['time_group', 'stim_current'])[var].agg(
            ['mean', 'sem', 'std']).reset_index()

    aggregated_data_var['group_offset'] = aggregated_data_var['stim_current'].map(
        offset_map)

    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 2))

    # Perform Mann-Whitney U test
    for stim_current in stim_currents:
        group1_var = df[(df['time_group'] == "Time Group 1") &
                        (df['stim_current'] == stim_current)][var]
        group2_var = df[(df['time_group'] == "Time Group 2") &
                        (df['stim_current'] == stim_current)][var]

        if len(group1_var) > 0 and len(group2_var) > 0:
            stat, p_val = mannwhitneyu(
                group1_var, group2_var, alternative=alt)

            # Determine number of stars based on significance level
            if p_val <= 0.001:
                stars = "***"
            elif p_val <= 0.01:
                stars = "**"
            elif p_val <= 0.05:
                stars = "*"
            else:
                stars = None  # No significance, no stars

            if stars:
                subset = aggregated_data_var[(aggregated_data_var['time_group'] == "Time Group 2") &
                                             (aggregated_data_var['stim_current'] == stim_current)]
                if not subset.empty:
                    last_median = subset['median'].iloc[0] if use_median else subset['mean'].iloc[0]
                    upper_error_bar_len = subset['q3'].iloc[0] - \
                        last_median if use_median else subset['std'].iloc[0]

                    # Position on "Time Group 2"
                    x_pos = 1 + offset_map[stim_current]
                    y_pos = last_median + upper_error_bar_len * 1.025  # Adjust y-position

                    ax.text(x_pos, y_pos, stars, color='k',
                            fontsize=10, ha='center')

        # Plot error bars and mean/median
        for time_group, x_pos in zip(["Time Group 1", "Time Group 2"], [0, 1]):
            subset = aggregated_data_var[(aggregated_data_var['time_group'] == time_group) &
                                         (aggregated_data_var['stim_current'] == stim_current)]
            if not subset.empty:
                if use_median:
                    ax.errorbar(
                        x_pos +
                        offset_map[stim_current], subset['median'].iloc[0],
                        yerr=[[subset['median'].iloc[0] - subset['q1'].iloc[0]],
                              [subset['q3'].iloc[0] - subset['median'].iloc[0]]],
                        fmt='-o', label=f"{stim_current} µA" if time_group == "Time Group 1" else "",
                        color=color_map[stim_current], ms=4
                    )
                else:
                    ax.errorbar(
                        x_pos +
                        offset_map[stim_current], subset['mean'].iloc[0],
                        yerr=subset['std'].iloc[0],
                        fmt='-o', label=f"{stim_current} µA" if time_group == "Time Group 1" else "",
                        color=color_map[stim_current], ms=4
                    )

    # ax.set_xticks([0, 1])
    ax.set_ylabel('Z-Score')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper left', fontsize=6)

    plt.tight_layout()
    plt.show()
    return ax


def plot_single_metric_aggregated_two_group_comparison(df, var="z_score", animal_id="Aggregated", comparison_alternative='less', use_median=True, ax=None):
    offset = 0.15
    stim_currents = sorted(df['stim_current'].unique())
    offset_map = {stim: i * offset - ((len(stim_currents) - 1) / 2)
                  * offset for i, stim in enumerate(stim_currents)}
    colors = sns.color_palette("deep", n_colors=len(stim_currents))
    color_map = {stim: colors[i] for i, stim in enumerate(stim_currents)}

    if use_median:
        aggregated_data_var = df.groupby(['rel_week', 'stim_current']).agg(
            median=(var, 'median'),
            q1=(var, lambda x: x.quantile(0.25)),
            q3=(var, lambda x: x.quantile(0.75))
        ).reset_index()
    else:
        aggregated_data_var = df.groupby(['rel_week', 'stim_current']).agg(
            mean=(var, 'mean'),
            sem=(var, lambda x: x.std() / np.sqrt(len(x))),
            std=(var, 'std')
        ).reset_index()

    aggregated_data_var['week_group_offset'] = aggregated_data_var['rel_week'] + \
        aggregated_data_var['stim_current'].map(offset_map)

    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 2.5))

    for stim_current in stim_currents:
        current_data_var = aggregated_data_var[aggregated_data_var['stim_current'] == stim_current]

        # plot individual weeks
        if use_median:
            ax.errorbar(
                current_data_var['week_group_offset'],
                current_data_var['median'],
                yerr=[current_data_var['median'] - current_data_var['q1'],
                      current_data_var['q3'] - current_data_var['median']],
                fmt='-o',
                label=f"{stim_current} µA",
                color=color_map[stim_current],
                ms=4
            )
        else:
            ax.errorbar(
                current_data_var['week_group_offset'],
                current_data_var['mean'],
                yerr=current_data_var['sem'],
                fmt='-o',
                label=f"{stim_current} µA",
                color=color_map[stim_current],
                ms=4
            )

        # Compare early vs late groups (0–1 vs 2–4)
        early_vals = df[(df['rel_week'].isin([0, 1])) & (
            df['stim_current'] == stim_current)][var]
        late_vals = df[(df['rel_week'].isin([2, 3, 4])) & (
            df['stim_current'] == stim_current)][var]

        if len(early_vals) > 0 and len(late_vals) > 0:
            stat, p_val = mannwhitneyu(
                early_vals, late_vals, alternative=comparison_alternative)
            print(f"{stim_current}µA p = {p_val:.4f}")

            if p_val <= 0.05:
                if p_val <= 0.001:
                    stars = '***'
                elif p_val <= 0.01:
                    stars = '**'
                else:
                    stars = '*'

                # Find max y-position across weeks 0–4 for this stim_current
                subset = aggregated_data_var[
                    (aggregated_data_var['rel_week'].isin([0, 1, 2, 3, 4])) &
                    (aggregated_data_var['stim_current'] == stim_current)
                ]
                y_max = subset['q3' if use_median else 'mean'].max()
                y_bar = y_max * 1.05  # space above top bar

                # Draw line above early (avg x of 0 and 1) and late (avg x of 2–4)
                early_xs = [week + offset_map[stim_current] for week in [0, 1]]
                late_xs = [week + offset_map[stim_current]
                           for week in [2, 3, 4]]
                early_center = np.mean(early_xs)
                late_center = np.mean(late_xs)

                ax.plot([early_center, late_center], [
                        y_bar, y_bar], color='black', lw=1)
                ax.text((early_center + late_center) / 2, y_bar * 1.01,
                        stars, ha='center', fontsize=10, color='k')

    ax.set_xticks(sorted(df['rel_week'].unique()))
    ax.set_xticklabels([f"Week {w}" for w in sorted(df['rel_week'].unique())])
    ax.set_ylabel(var.replace('_', ' ').capitalize())

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title="Current", loc='upper left')

    plt.tight_layout()
    plt.show()
    return ax


def plot_single_metric_at_threshold_current(df, raw_df, var="z_score", ax=None):

    thresh_map = (raw_df
                  .drop_duplicates(subset=['animal_id', 'session', 'stim_channel'])
                  [['animal_id', 'session', 'stim_channel', 'rel_week', 'detection_threshold']])

    df_ = df.merge(
        thresh_map,
        on=['animal_id', 'session', 'stim_channel'],
        how='inner',
        suffixes=('', '_chan')
    ).copy()

    df_['thr_dist'] = (df_['stim_current'] - df_['detection_threshold']).abs()
    df_ = df_.sort_values('thr_dist')
    nearest = (df_.groupby(['animal_id', 'session', 'stim_channel'], as_index=False)
                  .head(1))

    df_filt = nearest.dropna(subset=[var, 'rel_week']).copy()

    # df_filt = df[np.isclose(df['stim_current'], df['detection_threshold'])]

    # Mann-Whitney U Test on raw data
    df_phase = df_filt.copy()
    df_phase = df_phase[df_phase['rel_week'].between(0, 4)]
    df_phase['phase'] = np.where(df_phase['rel_week'] <= 1, 'early', 'late')

    early_vals = df_phase.loc[df_phase['phase'] == 'early', var].dropna()
    late_vals = df_phase.loc[df_phase['phase'] == 'late',  var].dropna()

    p_val_el = 1.0
    if len(early_vals) > 0 and len(late_vals) > 0:
        _, p_val_el = mannwhitneyu(
            early_vals, late_vals, alternative='two-sided')
        print(f"Mann-Whitney (weeks 0–1 vs 2–4) p = {p_val_el:.4f}")
    else:
        print("Not enough data for early/late comparison.")

    # --- Optional: draw a simple significance bar across early vs late on ax1 ---
    # Place the bar roughly from x=0.5 (between weeks 0 and 1) to x=3 (between weeks 3 and 4)
    if p_val_el <= 0.05:
        x1_bar, x2_bar = 0.5, 3.0
        # pick a y just above current left-axis data
        y_top = (df_var['q3'].max() if len(df_var) else 0)
        pad = 0.05 * (df_var['q3'].max() -
                      df_var['q1'].min()) if len(df_var) else 0.1
        y_bar = y_top + pad
        # draw bar
        ax1.plot([x1_bar, x1_bar, x2_bar, x2_bar], [y_bar, y_bar +
                 0.02, y_bar+0.02, y_bar], color=color_zscore, lw=1)
        ax1.text((x1_bar+x2_bar)/2.0, y_bar+0.03, '*', ha='center',
                 va='bottom', color=color_zscore, fontsize=10)

    # Aggregate z-score (first variable)
    df_var = df_filt.groupby(['rel_week']).agg(
        median=(var, 'median'),
        q1=(var, lambda x: x.quantile(0.25)),
        q3=(var, lambda x: x.quantile(0.75))
    ).reset_index()

    df_threshold_unique = raw_df.groupby(['animal_id', 'session', 'stim_channel'])[
        'detection_threshold'].median().reset_index()

    df_threshold_unique = df_threshold_unique.merge(raw_df[['animal_id', 'session', 'rel_week']].drop_duplicates(),
                                                    on=['animal_id', 'session'], how='left')
    df_threshold_unique = df_threshold_unique[df_threshold_unique['rel_week'].between(
        0, 4)]

    df_threshold = df_threshold_unique.groupby(['rel_week']).agg(
        median=('detection_threshold', 'median'),
        q1=('detection_threshold', lambda x: x.quantile(0.25)),
        q3=('detection_threshold', lambda x: x.quantile(0.75))
    ).reset_index()

    # Apply offset to avoid overlap
    offset = 0.1
    df_var['rel_week_offset'] = df_var['rel_week'] - \
        offset  # Shift left
    # Shift right
    df_threshold['rel_week_offset'] = df_threshold['rel_week'] + offset

    # Initialize figure
    if ax is None:
        fig, ax1 = plt.subplots(figsize=(4, 3))
    else:
        ax1 = ax

    ax2 = ax1.twinx()  # Create second y-axis

    # Define colors
    colors = sns.color_palette("deep")
    color_zscore = colors[4]  # Deep color 0
    color_threshold = colors[1]  # Deep color 1
    color_threshold = 'k'  # Deep color 1

    # Plot z-score on left y-axis (with slight left shift)
    ax1.errorbar(
        df_var['rel_week_offset'], df_var['median'],
        yerr=[df_var['median'] - df_var['q1'],
              df_var['q3'] - df_var['median']],
        fmt='-o', ms=4, color=color_zscore, label="Z-score"
    )
    ax1.set_ylabel(var.replace("_", " ").capitalize(), color=color_zscore)
    ax1.tick_params(axis='y', labelcolor=color_zscore)

    # Plot threshold on right y-axis (with slight right shift)
    ax2.errorbar(
        df_threshold['rel_week_offset'], df_threshold['median'],
        yerr=[df_threshold['median'] - df_threshold['q1'],
              df_threshold['q3'] - df_threshold['median']],
        fmt='-o', ms=4, color=color_threshold, label="Threshold (µA)"
    )
    ax2.set_ylabel("Threshold (µA)", color=color_threshold)
    ax2.tick_params(axis='y', labelcolor=color_threshold)

    # Mark significance if p-value ≤ 0.05
    if p_val <= 0.05:
        last_median = df_var[df_var['rel_week'] == 4]['median'].iloc[0]
        upper_error_bar_len = df_var[df_var['rel_week']
                                     == 4]['q3'].iloc[0] - last_median

        # Place significance marker on left y-axis
        x_pos = 4 - offset  # Last week
        y_pos = last_median + upper_error_bar_len + 0.1
        ax1.text(x_pos, y_pos, '*', color=color_zscore,
                 fontsize=10, ha='center')

    # Formatting
    ax1.set_xlabel("Weeks")
    ax1.set_title("Modulation at Threshold Current")
    # Ensure correct x-ticks
    ax1.set_xticks([0, 1, 2, 3, 4])
    # Keep weeks readable
    ax1.set_xticklabels([0, 1, 2, 3, 4])

    # Legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

    plt.tight_layout()
    return ax1, ax2


def plot_animal_unit_counts_by_week(raw_df, ax=None):
    weekly_counts = {'All': {}, 'Modulated': {}}
    animal_ids = ["ICMS92", "ICMS93", "ICMS98", "ICMS100", "ICMS101"]
    base_path = Path("C:/data/")

    df = raw_df[raw_df['rel_week'] < 5]
    for animal_id in animal_ids:
        filt_df = df[df['animal_id'] == animal_id]
        modulated_df = filt_df[filt_df['modulated'] == True]
        relative_weeks = np.arange(5)

        for relative_week in relative_weeks:
            units = filt_df[filt_df['rel_week'] == relative_week]
            mod_units = modulated_df[modulated_df['rel_week']
                                     == relative_week]

            # Append counts to the respective week list in weekly_counts
            weekly_counts['All'].setdefault(relative_week, []).append(
                units['unit_id'].nunique())
            weekly_counts['Modulated'].setdefault(relative_week, []).append(
                mod_units['unit_id'].nunique())

    # Plot each count type with different colors and a slight x-offset for clarity
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 3))

    categories = ['All', 'Modulated']
    # Small offset for each category to avoid overlap
    offsets = [-0.25, 0, 0.25]

    # Store median values for each category and week
    median_values_by_category = {category: [] for category in categories}

    for i, category in enumerate(categories):
        color = f'C{i}'
        offset = offsets[i]
        weeks = sorted(weekly_counts[category].keys())

        # Loop through each week to plot individual counts per animal and calculate medians
        for week in weeks:
            counts = weekly_counts[category][week]
            week_index = np.full(len(counts), week) + offset
            ax.scatter(week_index, counts, color=color, alpha=0.7)

            # Plot the median for the current week at the center of the week and store it
            # if len(counts) > 0:
            #     median_value = np.median(counts)
            #     median_values_by_category[category].append(median_value)
            #     plt.plot(week, median_value, marker='o',
            #              color=color, markersize=6)
            # else:
            #     median_values_by_category[category].append(np.nan)

        # # Connect the median points for each category with lines
        # ax.plot(
        #     weeks, median_values_by_category[category], color=color, linewidth=2)

    # Customize legend to display only scatter points
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=category,
                              markerfacecolor=f'C{i}', markersize=6) for i, category in enumerate(categories)]
    ax.legend(handles=legend_elements, fontsize='small')

    # Plot settings

    plt.tight_layout()
    plt.show()
    return ax


def plot_animal_unit_counts_by_week(raw_df, ax=None):
    """
    Plots the number of total units per animal over weeks.
    - Each animal has a distinct color (using seaborn "deep" palette).
    - Skips plotting markers and lines where count = 0.
    - Directly connects non-zero weeks, skipping over zero-count weeks.

    Args:
        raw_df (pd.DataFrame): DataFrame containing unit counts per week per animal.
        ax (matplotlib.axes.Axes, optional): Existing axes to plot on. Defaults to None.

    Returns:
        ax (matplotlib.axes.Axes): The modified axis object.
    """

    animal_ids = ["ICMS92", "ICMS93", "ICMS98", "ICMS100", "ICMS101"]
    # Use seaborn deep palette
    colors = sns.color_palette("deep", len(animal_ids))

    # Filter data to include only weeks < 5
    df = raw_df[raw_df['rel_week'] < 5]

    # Prepare dictionary for storing unit counts
    weekly_counts = {animal_id: {} for animal_id in animal_ids}

    # Collect unit counts for each animal and week
    for animal_id in animal_ids:
        filt_df = df[df['animal_id'] == animal_id]

        for week in range(5):  # Ensure 0-4 weeks
            units_count = filt_df[filt_df['rel_week']
                                  == week]['unit_id'].nunique()
            weekly_counts[animal_id][week] = units_count

    week_means = []
    for week in range(5):
        week_counts = [weekly_counts[animal_id][week]
                       for animal_id in animal_ids]
        week_means.append(np.median(week_counts))

    # Plot median across animals
    plt.plot(range(5), week_means, marker='o',
             markersize=2, color='k', linewidth=1)

    # Create figure if no axis is provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3))

    for i, animal_id in enumerate(animal_ids):
        weeks = []
        counts = []
        for week in range(5):
            count = weekly_counts[animal_id][week]
            if count > 0:
                weeks.append(week)
                counts.append(count)

        if weeks:
            ax.scatter(weeks, counts,
                       color='k', alpha=0.2, edgecolors='none', s=8, label=animal_id)
            ax.plot(weeks, counts, color='k', alpha=0.2)

    # Formatting
    ax.set_xticks(range(5))
    ax.set_xticklabels([f"{w}" for w in range(5)])
    ax.set_xlabel("Week")
    ax.set_ylabel("Total Units")
    # ax.set_title("Total Units per Animal Over Weeks")

    # Legend with animal names
    # legend_elements = [Line2D([0], [0], marker='o', color='w', label=animal_id,
    #                           markerfacecolor=colors[i], markersize=4) for i, animal_id in enumerate(animal_ids)]
    # ax.legend(handles=legend_elements, fontsize=3, loc="lower right")

    plt.tight_layout()
    plt.show()
    return ax


def plot_raw_and_filtered_spikes(raw_analyzer, analyzer, unit_id, ch_id=None, ax=None):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 2))


def plot_primary_ch_trace_with_unit_spikes(raw_analyzer, analyzer, unit_id,
                                           stim_ts, start_pulse, end_pulse,
                                           ch_id=None, marker_y_offset=150, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 2))

    start_frame = stim_ts[start_pulse]
    end_frame = stim_ts[end_pulse] - 290

    spike_train = analyzer.sorting.get_unit_spike_train(unit_id=unit_id)
    stim_spikes = spike_train[
        (spike_train >= start_frame) & (spike_train < end_frame)
    ]
    primary_ch_index = template_util.get_dense_primary_ch_index(
        analyzer, unit_id)

    if ch_id is None:
        ch_id = raw_analyzer.channel_ids[primary_ch_index]
    trace = raw_analyzer.recording.get_traces(
        start_frame=start_frame, end_frame=end_frame, channel_ids=[
            ch_id], return_scaled=True
    )

    time_axis = np.arange(trace.shape[0]) / 30  # Convert to ms

    ax.plot(time_axis, trace, linewidth=0.5, color='k')

    colors = sns.color_palette("deep")
    spikes_x = (stim_spikes - start_frame).astype(int)
    spikes_y = trace[spikes_x] - marker_y_offset

    ax.scatter(
        x=spikes_x / 30,
        y=spikes_y,
        color=colors[0],
        marker='^',
        s=12,
        label='Detected Spike'
    )

    ax.set_ylim([-3000, -1300])
    return ax


def plot_num_responses_vs_time(df_mod_only, current, ax1):
    if ax1 is None:
        fig, ax1 = plt.subplots(figsize=(4, 3))
    colors = sns.color_palette("deep")
    num_pl_responses = df_mod_only[df_mod_only['is_pulse_locked'] == True] \
        .groupby(['rel_week', 'stim_current']).size()
    num_mod_responses = df_mod_only.groupby(
        ['rel_week', 'stim_current']).size()
    # ratio of count
    df_pl = num_pl_responses.reset_index(name='pl_count')
    df_mod = num_mod_responses.reset_index(name='mod_count')

    df_combined = df_pl.merge(
        df_mod, on=['rel_week', 'stim_current'], how='right')
    df_combined['pl_ratio'] = df_combined['pl_count'] / \
        df_combined['mod_count']
    df_combined['pl_count_normalized'] = df_combined['pl_count'] / 5

    df_at_current = df_combined[df_combined['stim_current'] == current]

    ax2 = ax1.twinx()

    ax1.plot(df_at_current['rel_week'], df_at_current['pl_count_normalized'],
             marker='o', linestyle='-', color=colors[0], label="PL Response Count / Animal")

    ax1.set_ylabel("PL Response Count\nPer Animal", color=colors[0])
    ax1.tick_params(axis='y', labelcolor=colors[0])

    # Plot ratio of PL responses to total modulated responses (Right Y-axis)
    ax2.plot(df_at_current['rel_week'], df_at_current['pl_ratio'], marker='o',
             linestyle='-', color=colors[1], label="PL / Modulated Ratio")

    ax2.set_ylabel("PL / Modulated Ratio", color=colors[1])
    ax2.tick_params(axis='y', labelcolor=colors[1])

    # Formatting
    ax1.set_xlabel("Weeks Relative")
    ax1.set_title(f"{current}µA")

    # Legends (combining both axes)
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    # ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

    return ax1, ax2


def plot_pl_ratio_comparison(df_mod_only, current, ax):
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(3, 2))
    # Count number of pulse-locked responses per rel_week
    df_pl_only = df_mod_only[
        (df_mod_only['is_pulse_locked'] == True) &
        (df_mod_only['stim_current'] == current)
    ]

    df_mod_only_current = df_mod_only[df_mod_only['stim_current'] == current]

    # Group by week, animal, and session, then count
    pl_counts = (
        df_pl_only.groupby(['rel_week', 'animal_id', 'session'])
        .size()
        .reset_index(name='pl_count')
    )

    mod_counts = (
        df_mod_only_current.groupby(['rel_week', 'animal_id', 'session'])
        .size()
        .reset_index(name='mod_count')
    )

    sns.stripplot(
        data=pl_counts,
        x='rel_week',
        y='pl_count',
        ax=ax,
        jitter=True,
        alpha=0.7,
        color='tab:blue',
        label='Pulse-Locked'
    )

    # Plot modulated counts
    sns.stripplot(
        data=mod_counts,
        x='rel_week',
        y='mod_count',
        ax=ax,
        jitter=True,
        alpha=0.7,
        color='tab:orange',
        marker='D',
        label='Modulated'
    )

    ax.set_xlabel('Week')
    ax.set_ylabel('Number of Pulse-Locked Responses')
    ax.set_title(f'Number of Pulse-Locked Responses per Week at {current} uA')
    ax.legend().set_visible(False)


def plot_num_responses_vs_time_at_threshold(df_mod_only, ax1):
    if ax1 is None:
        fig, ax1 = plt.subplots(figsize=(4, 3))
    colors = sns.color_palette("deep")
    df_mod_only = df_mod_only[np.isclose(
        df_mod_only['stim_current'], df_mod_only['detection_threshold'])]

    num_pl_responses = df_mod_only[df_mod_only['is_pulse_locked'] == True] \
        .groupby(['rel_week']).size()
    num_mod_responses = df_mod_only.groupby(
        ['rel_week']).size()
    # ratio of count
    df_pl = num_pl_responses.reset_index(name='pl_count')
    df_mod = num_mod_responses.reset_index(name='mod_count')

    df_combined = df_pl.merge(
        df_mod, on=['rel_week'], how='right')
    df_combined['pl_ratio'] = df_combined['pl_count'] / \
        df_combined['mod_count']
    df_combined['pl_count_normalized'] = df_combined['pl_count'] / 5

    ax2 = ax1.twinx()

    ax1.plot(df_combined['rel_week'], df_combined['pl_count_normalized'],
             marker='o', linestyle='-', color=colors[0], label="PL Response Count / Animal")

    ax1.set_ylabel("PL Response Count\nPer Animal", color=colors[0])
    ax1.tick_params(axis='y', labelcolor=colors[0])

    # Plot ratio of PL responses to total modulated responses (Right Y-axis)
    ax2.plot(df_combined['rel_week'], df_combined['pl_ratio'], marker='o',
             linestyle='-', color=colors[1], label="PL / Modulated Ratio")

    ax2.set_ylabel("PL / Modulated Ratio", color=colors[1])
    ax2.tick_params(axis='y', labelcolor=colors[1])

    # Formatting
    ax1.set_xlabel("Weeks Relative")
    ax1.set_title(f"At threshold")

    # Legends (combining both axes)
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    # ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

    return ax1, ax2


def load_pl_session_responses():
    data_folder = 'C:\\data\\ICMS92\\behavior\\01-Sep-2023'
    pkl_path = Path(data_folder) / "batch_sort/session_responses_hit.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)
    return session_responses


def load_npl_session_responses():
    data_folder = 'C:\\data\\ICMS98\\24-Oct-2023'
    pkl_path = Path(data_folder) / "batch_sort/session_responses_hit.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)
    return session_responses


def plot_responses_2(df_mod_only, metric='pl', ax=None):
    if metric not in ['pl', 'npl', 'ratio']:
        raise ValueError("Invalid metric. Choose 'pl', 'npl', or 'ratio'.")
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 3))
    colors = sns.color_palette("deep")
    if metric == 'pl':
        y = df_mod_only[df_mod_only['is_pulse_locked'] == True] \
            .groupby(['rel_week']).size()
    elif metric == 'npl':
        y = df_mod_only[df_mod_only['is_pulse_locked'] == False] \
            .groupby(['rel_week']).size()
    else:
        y1 = df_mod_only[df_mod_only['is_pulse_locked'] == True] \
            .groupby(['rel_week']).size()
        y2 = df_mod_only[df_mod_only['is_pulse_locked'] == False] \
            .groupby(['rel_week']).size()
        y = y1 / y2
    ax.plot(np.arange(5), y)
    return ax


def plot_stim_vs_nonstim_template(mean_non_stim, mean_stim,
                                  std_non_stim, std_stim,
                                  ch_locs, sparse_ch_indices,
                                  plot_scalebar=True, y_negative=-250, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    pre_peak_space = 150
    x_offset = 0

    # Extract sparse channel locations and sort by depth
    sparse_ch_locs = ch_locs[sparse_ch_indices]
    depth_order = np.argsort(sparse_ch_locs[:, 1])
    abs_depth_id = 31 - sparse_ch_locs[depth_order, 1] / 60 + 1

    offsets = [i * pre_peak_space for i in range(len(depth_order))]

    for i, local_idx in enumerate(depth_order):
        depth_id = f"D{int(abs_depth_id[i])}"
        x = np.arange(mean_non_stim.shape[0]) + x_offset
        x_samples = np.arange(0, 90)
        x = x[x_samples]
        mean_non_stim = mean_non_stim[x_samples]
        mean_stim = mean_stim[x_samples]
        std_non_stim = std_non_stim[x_samples]
        std_stim = std_stim[x_samples]

        # Plot non-stim waveform and shading
        mean_ns = mean_non_stim[:, local_idx]
        std_ns = std_non_stim[:, local_idx]
        ax.plot(x, mean_ns + offsets[i], color='k', linewidth=0.75, alpha=0.4)
        ax.plot(-5, 0, color='k',
                label='Non-stim' if i == 0 else "", linewidth=0.75, alpha=0.6)
        ax.fill_between(x, mean_ns - std_ns + offsets[i], mean_ns + std_ns + offsets[i],
                        color='k', alpha=0.2, linewidth=0, edgecolor='none')

        # Plot stim waveform and shading
        mean_s = mean_stim[:, local_idx]
        std_s = std_stim[:, local_idx]
        ax.plot(x, mean_s + offsets[i], color='r',
                linestyle='-', linewidth=0.75, alpha=0.4)
        ax.plot(-5, 0, color='r', linestyle='-',
                label='Stim' if i == 0 else "", linewidth=0.75, alpha=0.6)
        ax.fill_between(x, mean_s - std_s + offsets[i], mean_s + std_s + offsets[i],
                        color='r', alpha=0.2, linewidth=0, edgecolor='none')

        # ax.text(5, offsets[i] - 50, depth_id, fontsize=6)

    ax.set_ylim([y_negative, offsets[-1] + 70])
    ax.axis("off")
    ax.legend(loc='lower right', frameon=True)

    if plot_scalebar:
        h_pos = (0, y_negative + 10)
        v_pos = (0, y_negative + 10)
        h_length = 30
        v_length = 100

        ax.plot([h_pos[0], h_pos[0] + h_length],
                [h_pos[1], h_pos[1]], 'k', lw=1)
        ax.plot([v_pos[0], v_pos[0]], [
                v_pos[1], v_pos[1] + v_length], 'k', lw=1)
        ax.text(h_pos[0], h_pos[1] - 50, '1 ms')
        ax.text(v_pos[0] - 8, v_pos[1], '100 µV', rotation=90)

    return ax


def load_stim_vs_nonstim_example():
    from spikeinterface import full as si

    with open('all_animal_folder_list.pkl', 'rb') as f:
        animal_folder_list = pickle.load(f)
    data_folder = animal_folder_list[0][1]
    analyzer = si.load_sorting_analyzer(
        Path(data_folder) / "batch_sort/merge/hmerge_analyzer.zarr")

    pkl_path = Path(data_folder) / "batch_sort/session_responses_all.pkl"
    with open(pkl_path, "rb") as file:
        session_responses = pickle.load(file)

    unit_id = 9
    wvf_ext = analyzer.load_extension('waveforms')
    ch_locs = analyzer.get_channel_locations()
    wvfs = wvf_ext.get_waveforms_one_unit(unit_id=unit_id)
    dense_wvfs, sparse_ch_indices, dense_primary_ch_idx = template_util.get_unit_dense_wvfs(
        analyzer, unit_id)

    sparse_ch_locs = ch_locs[sparse_ch_indices]
    depth_order = np.argsort(sparse_ch_locs[:, 1])
    abs_depth_id = 31 - sparse_ch_locs[depth_order, 1] / 60 + 1

    ur = session_responses.unit_responses[unit_id]
    non_stim_spike_indices = ur.non_stim_spike_indices
    stim_spike_indices = ur.stim_spike_indices
    non_stim_wvfs = wvfs[non_stim_spike_indices][:, :, sparse_ch_indices]
    stim_wvfs = wvfs[stim_spike_indices][:, :, sparse_ch_indices]

    data = {}
    data['stim_wvfs'] = stim_wvfs
    data['non_stim_wvfs'] = non_stim_wvfs
    data['sparse_ch_indices'] = sparse_ch_indices
    data['ch_locs'] = ch_locs

    return data
