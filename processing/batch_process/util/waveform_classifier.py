from matplotlib.gridspec import GridSpec
from matplotlib.widgets import CheckButtons
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
import batch_process.util.template_util as template_util
import batch_process.util.curate_util as curate_util
from spikeinterface import full as si
from pathlib import Path
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import recall_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os
import random
import matplotlib.pyplot as plt


class WaveformClassifier:
    def __init__(self, unit_id_waveform_dict, label_file="labeled_waveforms.csv", inference_file="inference_results.csv"):
        self.unit_id_waveform_dict = unit_id_waveform_dict
        self.label_file = label_file
        self.labeled_df = self._load_labels()
        self.inference_file = inference_file
        self.inference_results = self._load_inference_results()
        self.current_unit_id = None
        self.current_waveform_index = None
        self.current_waveform = None

    def _load_labels(self):
        if os.path.exists(self.label_file):
            print(f"Loaded existing labels from {self.label_file}")
            return pd.read_csv(self.label_file)
        print("No existing label file found. Starting fresh.")
        return pd.DataFrame(columns=["unit_id", "waveform_index", "label"])

    def save_labels(self):
        self.labeled_df.to_csv(self.label_file, index=False)
        print(f"Labels saved to {self.label_file}")

    def _load_inference_results(self):
        if os.path.exists(self.inference_file):
            print(f"Loaded existing inference results from {self.inference_file}")
            return pd.read_csv(self.inference_file)
        print("No existing inference results file found. You will need to run inference.")
        return None

    def pick_random_waveform(self):
        all_unit_ids = list(self.unit_id_waveform_dict.keys())
        while True:
            unit_id = random.choice(all_unit_ids)
            # Use filtered to get length
            waveforms = self.unit_id_waveform_dict[unit_id]['filt']
            all_indices = set(range(len(waveforms)))
            labeled_indices = set(
                self.labeled_df[self.labeled_df["unit_id"]
                                == unit_id]["waveform_index"]
            )
            remaining_indices = list(all_indices - labeled_indices)
            if remaining_indices:
                waveform_index = random.choice(remaining_indices)
                return unit_id, waveform_index  # No need to return waveform here
            print("All waveforms for this unit are labeled. Moving to the next unit.")

    def start_labeling(self):
        self.current_unit_id, self.current_waveform_index = self.pick_random_waveform()

        # Create figure and axes
        self.fig = plt.figure(figsize=(10, 6))
        gs = GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])
        self.ax_raw = self.fig.add_subplot(
            gs[0, 0])  # Top-left for raw waveforms
        # Bottom-left for filtered waveforms
        self.ax_filt = self.fig.add_subplot(gs[1, 0])
        # Right for checkboxes, spans both rows
        self.ax_checkbox = self.fig.add_subplot(gs[:, 1])
        self.ax_checkbox.axis("off")

        # Plot waveforms
        self.plot_waveform(
            unit_id=self.current_unit_id,
            waveform_index=self.current_waveform_index,
            fig=self.fig,
            ax_raw=self.ax_raw,
            ax_filt=self.ax_filt,
            ax_checkbox=self.ax_checkbox
        )

        # Connect key event handler
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)

        # Block until the plot window is closed
        print("Labeling started. Close the GUI to proceed.")
        plt.show(block=True)

        # Save the labels after closing
        self.save_labels()

    def plot_waveform(self, unit_id, waveform_index, fig, ax_raw, ax_filt, ax_checkbox):
        # Retrieve waveform data
        raw_waveform = self.unit_id_waveform_dict[unit_id]['raw'][waveform_index]
        filt_waveform = self.unit_id_waveform_dict[unit_id]['filt'][waveform_index]
        sparse_indices = self.unit_id_waveform_dict[unit_id]['sparse_ch_indices']
        primary_channel_idx = self.unit_id_waveform_dict[unit_id]['dense_primary_ch_index']

        # Clear previous plots
        ax_raw.clear()
        ax_filt.clear()

        # Configure raw waveform plot
        ax_raw.set_title(f"Unit {unit_id}, Raw Waveforms")
        ax_raw.set_xlabel("Time (samples)")
        ax_raw.set_ylabel("Amplitude (µV)")
        ax_raw.set_ylim(-300, 400)  # Set y-axis limits for raw waveform
        ax_raw.grid(True)

        # Configure filtered waveform plot
        ax_filt.set_title(f"Unit {unit_id}, Filtered Waveforms")
        ax_filt.set_xlabel("Time (samples)")
        ax_filt.set_ylabel("Amplitude (µV)")
        ax_filt.grid(True)

        # Persistent references for plotted lines
        self.lines_raw = {"Dense indices": [],
                          "Sparse indices": [], "Primary channel": None}
        self.lines_filt = {"Dense indices": [],
                           "Sparse indices": [], "Primary channel": None}

        # Plot dense indices
        for ch in range(raw_waveform.shape[1]):
            self.lines_raw["Dense indices"].append(
                ax_raw.plot(raw_waveform[:, ch], color='gray',
                            alpha=0.5, label="Dense" if ch == 0 else "")[0]
            )
            self.lines_filt["Dense indices"].append(
                ax_filt.plot(filt_waveform[:, ch], color='gray',
                             alpha=0.5, label="Dense" if ch == 0 else "")[0]
            )

        # Plot sparse indices
        for ch in sparse_indices:
            self.lines_raw["Sparse indices"].append(
                ax_raw.plot(raw_waveform[:, ch], color='blue',
                            label="Sparse" if ch == sparse_indices[0] else "")[0]
            )
            self.lines_filt["Sparse indices"].append(
                ax_filt.plot(filt_waveform[:, ch], color='blue',
                             label="Sparse" if ch == sparse_indices[0] else "")[0]
            )

        # Plot primary channel
        self.lines_raw["Primary channel"], = ax_raw.plot(
            raw_waveform[:, primary_channel_idx], color='red', label="Primary"
        )
        self.lines_filt["Primary channel"], = ax_filt.plot(
            filt_waveform[:, primary_channel_idx], color='red', label="Primary"
        )

        ax_raw.legend(loc="upper right", fontsize=8)

        # Add text instructions
        ax_raw.text(
            0.1, 0.9, "Press 1: Spike, 2: Not Spike, 3: Unsure",
            transform=ax_raw.transAxes, fontsize=10, color='k'
        )

        # Adjust subplot spacing to prevent overlap
        plt.subplots_adjust(hspace=0.4)

        # Refresh canvas
        fig.canvas.draw_idle()

        # Create checkboxes
        check_labels = list(self.lines_raw.keys())
        if not hasattr(self, 'check'):
            self.check = CheckButtons(ax_checkbox, check_labels, [
                                      True] * len(check_labels))

            # Adjust vertical spacing of checkboxes
            self.check.linespacing = 0.1  # Reduce spacing between checkboxes

            def toggle_lines(label):
                if isinstance(self.lines_raw[label], list):
                    visible = not all(line.get_visible()
                                      for line in self.lines_raw[label])
                    for line in self.lines_raw[label]:
                        line.set_visible(visible)
                    for line in self.lines_filt[label]:
                        line.set_visible(visible)
                else:  # Single Line2D object (e.g., Primary channel)
                    visible = not self.lines_raw[label].get_visible()
                    self.lines_raw[label].set_visible(visible)
                    self.lines_filt[label].set_visible(visible)

                # Redraw the canvas
                fig.canvas.draw_idle()

                # Restore focus to the figure
                fig.canvas.manager.window.activateWindow()

            self.check.on_clicked(toggle_lines)

        # Move checkboxes closer to the left
        # [left, bottom, width, height]
        ax_checkbox.set_position([0.7, 0.5, 0.15, 0.3])

    def on_key(self, event):
        label_map = {"1": "spike", "2": "not_spike", "3": "unsure"}
        if event.key in label_map:
            label = label_map[event.key]
            new_entry = {
                "unit_id": self.current_unit_id,
                "waveform_index": self.current_waveform_index,
                "label": label,
            }
            self.labeled_df = pd.concat(
                [self.labeled_df, pd.DataFrame([new_entry])], ignore_index=True
            )

            # Pick the next waveform
            try:
                self.current_unit_id, self.current_waveform_index = self.pick_random_waveform()

                # Update the plot dynamically
                self.plot_waveform(
                    self.current_unit_id,
                    self.current_waveform_index,
                    self.fig,
                    self.ax_raw,
                    self.ax_filt,
                    self.ax_checkbox
                )
            except StopIteration:
                self.save_labels()
                print("All waveforms labeled!")
                plt.close(self.fig)

    def get_train_test_split(self):
        waveforms_list = []
        labels_list = []
        filtered_labeled_df = self.labeled_df[self.labeled_df["label"] != "unsure"]
        for _, row in filtered_labeled_df.iterrows():
            unit_id = row["unit_id"]
            waveform_index = row["waveform_index"]

            # Validate unit_id and waveform_index
            if unit_id not in self.unit_id_waveform_dict:
                print(f"Invalid unit_id: {unit_id}. Skipping.")
                continue
            if waveform_index >= len(self.unit_id_waveform_dict[unit_id]['raw']):
                print(f"Invalid waveform index {waveform_index} for unit_id {unit_id}. Skipping.")
                continue

            # Access raw waveform
            primary_channel_idx = self.unit_id_waveform_dict[unit_id]['dense_primary_ch_index']
            wvf = self.unit_id_waveform_dict[unit_id]['raw'][waveform_index][:,
                                                                             primary_channel_idx]
            label = row["label"]

            waveforms_list.append(wvf)
            labels_list.append(label)

        X = np.array(waveforms_list)
        y = np.array(labels_list)

        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )

    def train_model(self, clf, evaluate=True):
        """
        Train the model using the provided classifier and evaluate its performance.

        Parameters:
        - clf: Classifier instance (e.g., RandomForestClassifier, MLPClassifier).
        - evaluate (bool): If True, print evaluation metrics after training.

        Returns:
        - y_pred: Predicted labels for the test set.
        - metrics: A dictionary containing classification report and confusion matrix.
        """
        from sklearn.preprocessing import LabelEncoder

        # clf = RandomForestClassifier(n_estimators=100, random_state=42)
        # clf = MLPClassifier(hidden_layer_sizes=(50, 25), random_state=42)
        if self.X_train is None:
            self.get_train_test_split()

        # Ensure clf is a valid classifier
        if not hasattr(clf, "fit") or not callable(getattr(clf, "fit")):
            raise ValueError(
                "The provided 'clf' is not a valid classifier. Ensure it has a 'fit' method.")

        use_label_encoder = isinstance(clf, XGBClassifier)
        if use_label_encoder:
            le = LabelEncoder()
            y_train_encoded = le.fit_transform(self.y_train)
            y_test_encoded = le.transform(self.y_test)
        else:
            y_train_encoded, y_test_encoded = self.y_train, self.y_test

        self.clf = clf
        print(f"Training Classifier: {type(self.clf).__name__}")
        # Train the classifier
        self.clf.fit(self.X_train, y_train_encoded)
        # Predict labels for the test set
        y_pred_encoded = self.clf.predict(self.X_test)
        if use_label_encoder:
            y_pred = le.inverse_transform(y_pred_encoded)
        else:
            y_pred = y_pred_encoded

        # Print evaluation metrics
        print("Classification Report:")
        print(classification_report(self.y_test, y_pred))
        print("Confusion Matrix:")
        print(confusion_matrix(self.y_test, y_pred))

    def infer(self, re_infer=False):
        """
        Perform inference on a set of waveforms and store results in a DataFrame.

        Parameters:
        - re_infer (bool): If True, perform inference again even if the inference file exists.

        Returns:
        - inference_results (pd.DataFrame): A DataFrame with predictions and confidence scores for each waveform.
        """
        if self.inference_results is not None and not re_infer:
            print(
                "Inference already performed. Skipping. Use `re_infer=True` to force re-inference.")
            return self.inference_results

        results = []

        for unit_id, waveforms in self.unit_id_waveform_dict.items():
            # Use the primary channel waveform for inference
            primary_channel_idx = waveforms['dense_primary_ch_index']
            # Extract primary channel waveform
            primary_waveforms = waveforms['raw'][:, :, primary_channel_idx]

            # Ensure the input to the classifier is a valid NumPy array
            primary_waveforms = np.array(primary_waveforms)

            # Perform predictions
            predicted_labels = self.clf.predict(primary_waveforms)
            confidence_scores = self.clf.predict_proba(primary_waveforms)

            # Collect results for each waveform
            for idx, (label, scores) in enumerate(zip(predicted_labels, confidence_scores)):
                results.append({
                    "Unit ID": unit_id,
                    "Waveform Index": idx,
                    "Predicted Label": label,
                    "Confidence": max(scores)  # Highest probability score
                })

        # Convert results to a DataFrame
        self.inference_results = pd.DataFrame(results)

        # Save the results to a CSV file
        self.inference_results.to_csv("inference_results.csv", index=False)
        print("Inference results saved to 'inference_results.csv'")

        return self.inference_results

    def plot_waveforms(self, lower_threshold=0.0, upper_threshold=1.0, max_waveforms=100):
        """
        Plot waveforms for predictions with confidence within a specified range.

        Parameters:
        - lower_threshold (float): Minimum confidence value for inclusion.
        - upper_threshold (float): Maximum confidence value for inclusion.
        - max_waveforms (int): Maximum number of waveforms to plot.

        Returns:
        - raw_waveforms (np.array): Array of raw waveforms plotted.
        - filt_waveforms (np.array): Array of filtered waveforms plotted.
        - confidences (np.array): Adjusted confidence scores.
        """
        if self.inference_results is None or self.inference_results.empty:
            print("No inference results available to plot.")
            return None, None, None

        # Adjust confidence scores
        adjusted_results = self.inference_results.copy()
        adjusted_results["Adjusted Confidence"] = adjusted_results.apply(
            lambda row: row["Confidence"] if row["Predicted Label"] == "spike" else 1 -
            row["Confidence"],
            axis=1
        )

        # Filter results based on adjusted confidence
        filtered_results = adjusted_results[
            (adjusted_results["Adjusted Confidence"] >= lower_threshold) &
            (adjusted_results["Adjusted Confidence"] <= upper_threshold)
        ]

        if filtered_results.empty:
            print(f"No waveforms with confidence in range [{lower_threshold}, {upper_threshold}].")
            return None, None, None

        # Randomly sample waveforms if the filtered set is larger than max_waveforms
        if len(filtered_results) > max_waveforms:
            filtered_results = filtered_results.sample(
                n=max_waveforms, random_state=42)

        # Initialize arrays to store waveforms
        raw_waveforms = []
        filt_waveforms = []
        confidences = []

        # Create subplots
        fig, axes = plt.subplots(2, 1, figsize=(
            8, 8), sharex=True, constrained_layout=True)

        for _, row in filtered_results.iterrows():
            unit_id = row["Unit ID"]
            waveform_index = row["Waveform Index"]
            adjusted_confidence = row["Adjusted Confidence"]

            # Validate unit_id and waveform_index
            if unit_id not in self.unit_id_waveform_dict:
                print(f"Unit ID {unit_id} not found. Skipping.")
                continue
            if waveform_index >= len(self.unit_id_waveform_dict[unit_id]['raw']):
                print(f"Invalid waveform index {waveform_index} for Unit ID {unit_id}. Skipping.")
                continue

            # Retrieve the waveforms
            primary_ch_index = self.unit_id_waveform_dict[unit_id]['dense_primary_ch_index']
            raw_waveform = self.unit_id_waveform_dict[unit_id]['raw'][waveform_index,
                                                                      :, primary_ch_index]
            filt_waveform = self.unit_id_waveform_dict[unit_id]['filt'][waveform_index,
                                                                        :, primary_ch_index]

            # Append to arrays
            raw_waveforms.append(raw_waveform)
            filt_waveforms.append(filt_waveform)
            confidences.append(adjusted_confidence)

            # Plot waveforms
            axes[0].plot(raw_waveform, color='k', linewidth=0.5, alpha=0.7)
            axes[1].plot(filt_waveform, color='r', linewidth=0.5, alpha=0.7)

        # Configure the subplots
        axes[0].set_title("Raw Waveforms")
        axes[0].set_ylabel("Amplitude (µV)")
        axes[0].set_ylim([-100, 1000])
        axes[0].grid(True)

        axes[1].set_title("Filtered Waveforms")
        axes[1].set_xlabel("Time (samples)")
        axes[1].set_ylabel("Amplitude (µV)")
        axes[1].grid(True)

        plt.tight_layout()
        plt.show()

        # Return waveforms as arrays
        return np.array(raw_waveforms), np.array(filt_waveforms), np.array(confidences)


# %%
# Initialize the classifier with the waveform data
# %% Load data if not already loaded
if __name__ == "__main__":
    data_folder = Path("C:\\data\\ICMS92\\behavior\\01-Sep-2023")
    analyzer = si.load_sorting_analyzer(
        folder=data_folder / "batch_sort/stage2/stage2_analyzer.zarr")
    raw_analyzer = curate_util.load_raw_analyzer(data_folder, analyzer)
    unit_ids = raw_analyzer.unit_ids

# %%
    unit_id_waveform_dict = {}
    for unit_id in unit_ids:
        filt_dense_wvfs, sparse_ch_indices, dense_primary_ch_idx, raw_dense_wvfs = template_util.get_unit_dense_wvfs(
            analyzer, unit_id, raw_analyzer)

        timesample_30 = raw_dense_wvfs[:, 30, :]
        centered_raw_wvfs = raw_dense_wvfs - timesample_30[:, None, :]

        # Store both 'raw' and 'filt' waveforms in a dictionary per unit ID
        unit_id_waveform_dict[unit_id] = {
            'raw': centered_raw_wvfs,
            'filt': filt_dense_wvfs,
            'sparse_ch_indices': sparse_ch_indices,
            'dense_primary_ch_index': dense_primary_ch_idx
        }

    # %%
    wfc = WaveformClassifier(unit_id_waveform_dict)

    # Start labeling
    wfc.start_labeling()

    wfc.get_train_test_split()

    # %%

    # random forest seems to be the best
    rf_clf = RandomForestClassifier(n_estimators=100, random_state=42)
    wfc.train_model(rf_clf)

    # adaboost_clf = AdaBoostClassifier(random_state=42)
    # wfc.train_model(adaboost_clf)

    # mlp_clf = MLPClassifier(hidden_layer_sizes=(
    #     40, 40), random_state=42, early_stopping=True)
    # wfc.train_model(mlp_clf)

    # xgb_clf = XGBClassifier(
    #     n_estimators=100,        # Number of trees
    #     max_depth=3,             # Maximum tree depth
    #     learning_rate=0.1,       # Learning rate (shrinkage)
    #     random_state=42,         # For reproducibility
    #     use_label_encoder=False  # For newer versions of XGBoost
    # )

    # wfc.train_model(xgb_clf)
    # %%
    inference_results = wfc.infer(re_infer=True)
    # %%
    lower_threshold = 0.3
    upper_threshold = 0.4
    raw_waveforms, filt_waveforms, confidences = wfc.plot_waveforms(
        lower_threshold=lower_threshold, upper_threshold=upper_threshold, max_waveforms=2)
    # %%
    wvf_range = np.arange(0, 100)
    fig, axes = plt.subplots(2, 1, figsize=(
        5, 4), sharex=True, constrained_layout=True)
    t = np.arange(90)/30
    axes[0].plot(t, raw_waveforms[wvf_range].T, 'k', linewidth=0.5)
    axes[0].set_ylim([-200, 800])
    axes[0].set_title("Raw waveforms", fontsize=10)
    axes[1].plot(t, filt_waveforms[wvf_range].T, 'r', linewidth=0.5)
    axes[1].set_title("Filtered waveforms", fontsize=10)
    axes[1].set_xlabel("Time (ms)", fontsize=8)
    axes[1].set_ylabel("Amplitude (uV)", fontsize=8)

    for ax in axes:
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.tick_params(axis='both', which='minor', labelsize=8)

    fig.suptitle(f'Confidence range: [{lower_threshold}, {upper_threshold}]')
    plt.tight_layout()

# %%
