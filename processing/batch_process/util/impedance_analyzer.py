import os
import glob
import re
from datetime import datetime
import numpy as np
import pandas as pd
# from util.test_impedance_analyzer import TestImpedanceAnalyzer


class ImpedanceAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def intan_to_ripple(intan_channels):
        intan_to_ripple_map = np.array([31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11,
                                       9, 7, 5, 3, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32])
        ripple_channels = intan_to_ripple_map[intan_channels]
        return ripple_channels

    @staticmethod
    def ripple_to_depth(unordered_indices):
        ordered_indices = np.array([32, 15, 28, 11, 24, 7, 18, 3, 29, 16, 25, 10,
                                   21, 6, 17, 2, 30, 13, 26, 9, 22, 5, 20, 1, 31, 14, 27, 12, 23, 8, 19, 4])
        unordered_to_ordered_indices = np.where(
            unordered_indices == ordered_indices[:, None])[1]
        return unordered_to_ordered_indices  # 0-index

    @staticmethod
    def ripple_to_depth2(unordered_indices):
        ordered_indices = np.array([32, 15, 28, 11, 24, 7, 18, 3, 29, 16, 25, 10,
                                   21, 6, 17, 2, 30, 13, 26, 9, 22, 5, 20, 1, 31, 14, 27, 12, 23, 8, 19, 4])
        unordered_indices = np.array(unordered_indices)
        unordered_indices = unordered_indices[unordered_indices != 0]
        index_mapping = {value: idx for idx,
                         value in enumerate(ordered_indices)}
        ordered_positions = np.array(
            [index_mapping[val] for val in unordered_indices])
        return ordered_positions + 1

    @staticmethod
    def intan_to_depth(unordered_indices):
        ordered_indices = np.array([31, 8, 29, 10, 27, 12, 24, 14, 1, 23, 3, 20,
                                   5, 18, 7, 16, 30, 9, 28, 11, 26, 13, 25, 15, 0, 22, 2, 21, 4, 19, 6, 17])
        unordered_to_ordered_indices = np.where(
            unordered_indices == ordered_indices[:, None])[1]
        return unordered_to_ordered_indices

    @staticmethod
    def _parse_csv_date(csv_path):
        """Parse date from impedance CSV filename (e.g. '7-24-23.csv' or '9-05-23.csv')."""
        name = os.path.splitext(os.path.basename(csv_path))[0]
        # Try M-D-YY or M-DD-YY format
        m = re.match(r'(\d{1,2})-(\d{1,2})-(\d{2,4})', name)
        if m:
            month, day, year = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if year < 100:
                year += 2000
            return datetime(year, month, day)
        return None

    @staticmethod
    def _parse_session_date(date_str):
        """Parse session date from folder name (e.g. '24-Jul-2023' or '08-Apr-2022')."""
        for fmt in ('%d-%b-%Y', '%m-%d-%Y', '%Y-%m-%d'):
            try:
                return datetime.strptime(date_str, fmt)
            except ValueError:
                continue
        return None

    def _find_closest_csv(self, csv_files, session_date):
        """Find the CSV with date closest to (but not after) the session date.
        Falls back to the closest CSV overall if none are before/on the session date."""
        dated = []
        for f in csv_files:
            d = self._parse_csv_date(f)
            if d is not None:
                dated.append((f, d))

        if not dated:
            # Can't parse dates, fall back to newest by mtime
            csv_files.sort(key=os.path.getmtime, reverse=True)
            return csv_files[0]

        # Prefer CSVs on or before session date (closest preceding measurement)
        before = [(f, d) for f, d in dated if d <= session_date]
        if before:
            return max(before, key=lambda x: x[1])[0]

        # No CSV before session — use closest overall
        return min(dated, key=lambda x: abs((x[1] - session_date).days))[0]

    def get_intan_impedances(self, animal_id, server_mount_drive, reorder=False, to_ripple=False, session_date=None):
        self.animal_id = animal_id
        self.server_mount_drive = server_mount_drive
        file_pattern = os.path.join(
            self.server_mount_drive, self.animal_id, "Impedance", "*.csv")
        csv_files = glob.glob(file_pattern)

        if not csv_files:
            # Try alternate folder name
            file_pattern = os.path.join(
                self.server_mount_drive, self.animal_id, "Impedances", "*.csv")
            csv_files = glob.glob(file_pattern)

        if not csv_files:
            print(f"No CSV files found for animalID {self.animal_id}")
            return None

        if session_date is not None:
            if isinstance(session_date, str):
                session_date = self._parse_session_date(session_date)
            if session_date is not None:
                chosen = self._find_closest_csv(csv_files, session_date)
            else:
                csv_files.sort(key=os.path.getmtime, reverse=True)
                chosen = csv_files[0]
        else:
            csv_files.sort(key=os.path.getmtime, reverse=True)
            chosen = csv_files[0]

        print(f"Using impedance file: {chosen}")
        df = pd.read_csv(chosen)
        impedances = df['Impedance Magnitude at 1000 Hz (ohms)']

        if reorder:
            if to_ripple:
                ripple_ch = self.intan_to_ripple(np.arange(32))
                indices = self.ripple_to_depth(ripple_ch)
            else:
                indices = self.intan_to_depth(np.arange(32))
        elif to_ripple:
            indices = self.intan_to_ripple(np.arange(32))
        else:
            indices = np.arange(32)
        self.impedances = impedances[indices]
        return impedances[indices]

    def get_good_impedances(self, threshold):
        impedances = self.impedances
        impedances = impedances[impedances <= threshold]
        return impedances, impedances.index.tolist()


# if __name__ == "__main__":
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestImpedanceAnalyzer)
    # unittest.TextTestRunner().run(suite)
    # impedance_analyzer = ImpedanceAnalyzer()
    # impedances = impedance_analyzer.get_intan_impedances('ICMS92')
    # filtered_impedances, good_channels = impedance_analyzer.get_good_impedances(
    #     threshold=1e6)
    # good_ripple_channels = impedance_analyzer.intan_to_ripple(good_channels)
    # print(good_ripple_channels)
