"""Data filtering utilities for electrophysiology analyses.

Standard filters for modulated, pulse-locked, and non-pulse-locked units.
"""
import numpy as np
import pandas as pd


def filter_modulated(df, max_z_score=50, min_spikes=50):
    """Filter to modulated units (both PL and NPL)."""
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


def filter_pl(df, max_z_score=50, min_spikes=50):
    """Filter to pulse-locked modulated units."""
    for ch in df['stim_channel'].unique():
        ch_df = df[df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            df = df[df['stim_channel'] != ch]

    df = df[
        (df["rel_week"] < 5) &
        (df["baseline_too_slow"] == False) &
        (df['modulated'] == True) &
        (df['mod_p_val'] < 0.05) &
        (df["is_pulse_locked"] == True) &
        (df["num_trials"] > 5) &
        (df["z_score"] > 0) &
        (df["z_score"] < max_z_score) &
        (df['num_spikes'] > min_spikes) &
        (df['stim_current'] < 7) &
        (df['stim_current'] > 3)]
    return df


def filter_npl(df, max_z_score=50, min_spikes=50):
    """Filter to non-pulse-locked modulated units."""
    for ch in df['stim_channel'].unique():
        ch_df = df[df['stim_channel'] == ch]
        if len(ch_df['rel_day'].unique()) < 5:
            df = df[df['stim_channel'] != ch]

    df = df[
        (df["rel_week"] < 5) &
        (df["baseline_too_slow"] == False) &
        (df['modulated'] == True) &
        (df['mod_p_val'] < 0.05) &
        (df["is_pulse_locked"] == False) &
        (df["num_trials"] > 5) &
        (df["z_score"] > 0) &
        (df["z_score"] < max_z_score) &
        (df['num_spikes'] > min_spikes) &
        (df['stim_current'] < 7) &
        (df['stim_current'] > 3)]
    return df
