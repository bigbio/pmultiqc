"""
Statistical utility functions for PMultiQC

This module contains various statistical and mathematical functions used across
different PMultiQC modules for data analysis and visualization.
"""

import numpy as np
import pandas as pd


def nanmedian(values: np.ndarray, all_nan_fallback: np.float64) -> np.float64:
    """
    Compute the median of the given array, ignoring NaNs; if all values are NaN (or the array is empty), return a fallback.

    Parameters:
    -----------
    values: numpy array with values
    all_nan_fallback: if all *values* are NaN, return this number instead of NaN

    Returns
    -------
    float
        The median of the non-NaN values in *values*. If all entries are NaN (or the array is empty),
        returns *all_nan_fallback*.
    """
    if np.isnan(values).all():
        return all_nan_fallback
    else:
        return np.nanmedian(values)


def qual_uniform(group_df_rt):
    """
    Calculate quality uniformity metric for retention time data.

    Parameters:
    -----------
    group_df_rt: pandas Series
        Group retention time data (group["Retention time"] or group["retention_time"])

    Returns:
    --------
    float
        Quality uniformity score between 0 and 1
    """
    x = group_df_rt / np.nansum(group_df_rt)
    n = group_df_rt.notna().sum()
    y = np.nansum(x) / n
    worst = ((1 - y) ** 0.5) * 1 / n + (y**0.5) * (n - 1) / n
    sc = np.sum(np.abs(x - y) ** 0.5) / n
    result = 1.0 if worst == 0 else float((worst - sc) / worst)
    return result


def calculate_delta_mass_distribution(df, col):
    """
    Calculate delta mass distribution with count and frequency bins.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data
    col: str
        Column name for delta mass calculation

    Returns:
    --------
    dict
        Dictionary containing count and frequency bin data
    """
    count_bin = df[col].value_counts(sort=False, bins=1000)
    count_bin_data = dict()
    for index in count_bin.index:
        count_bin_data[float(index.mid)] = int(count_bin[index])

    frequency_bin = df[col].value_counts(sort=False, bins=1000, normalize=True)
    frequency_bin_data = dict()
    for index in frequency_bin.index:
        frequency_bin_data[float(index.mid)] = float(frequency_bin[index])

    delta_mass = {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }
    return delta_mass


def calculate_modification_percentage(group):
    """
    Calculate modification percentages for a group of data.

    Parameters:
    -----------
    group: pandas DataFrame
        Group data containing modifications

    Returns:
    --------
    pandas DataFrame
        DataFrame with modification percentages
    """
    if "Modifications" in group.columns:
        group.rename(columns={"Modifications": "modifications"}, inplace=True)

    counts = group["modifications"].str.split(",").explode().value_counts()
    percentage_df = (counts / len(group["modifications"]) * 100).reset_index()
    percentage_df.columns = ["modifications", "percentage"]

    # Modified (Total)
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"] = (
        100 - percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"]
    )
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "modifications"] = (
        "Modified (Total)"
    )

    return percentage_df


def calculate_na_statistics(df, cols):
    """
    Calculate NA/Non-NA statistics for specified columns grouped by species.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data with species column
    cols: list
        List of column names to analyze

    Returns:
    --------
    dict
        Dictionary with NA/Non-NA counts per species and column
    """
    data_dict = {
        specie + ": " + col: {"Non-NA": group[col].notna().sum(), "NA": group[col].isna().sum()}
        for specie, group in df.groupby("species")
        for col in cols
    }
    return data_dict


def calculate_line_statistics(df, cols, dict_key, only_one_col=False):
    """
    Calculate line plot statistics with binning for specified columns.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data with species column
    cols: list
        List of column names to analyze
    dict_key: str
        Type of statistics (e.g., "cv")
    only_one_col: bool
        Whether to use only one column for grouping

    Returns:
    --------
    dict
        Dictionary with binned statistics data
    """
    data_union = df[cols].stack().reset_index(drop=True)
    data_range_max = data_union.max()

    if dict_key == "cv":
        data_range_min = 0
    else:
        data_range_min = data_union.min() - (data_union.min() * 0.05)

    bin_step = 1000
    data_bins = np.arange(
        data_range_min, data_range_max + (data_range_max * 0.05), (data_range_max / bin_step)
    )
    bin_mid = (data_bins[:-1] + data_bins[1:]) / 2

    data_dict = dict()
    for specie, group in df.groupby("species"):
        for col in cols:
            counts, _ = np.histogram(group[col].dropna(), bins=data_bins)
            cv_counts = pd.DataFrame({"value": bin_mid, "counts": counts})

            if only_one_col:
                data_dict[specie] = dict(zip(cv_counts["value"], cv_counts["counts"]))
            else:
                data_dict[specie + ": " + col] = dict(zip(cv_counts["value"], cv_counts["counts"]))

    return data_dict


def calculate_box_plot_statistics(df, cols):
    """
    Calculate box plot statistics for specified columns.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data with species column
    cols: list
        List of column names to analyze

    Returns:
    --------
    dict
        Dictionary with box plot data per species and column
    """
    boxplot_data = {
        specie + ": " + col: group[col].dropna().tolist()
        for specie, group in df.groupby("species")
        for col in cols
    }
    return boxplot_data


def calculate_intensity_counts_per_file(df, runs_col):
    """
    Calculate intensity counts per file for specified run columns.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data
    runs_col: str
        Column name pattern for run columns

    Returns:
    --------
    dict
        Dictionary with Non-NA and NA counts per file
    """
    import os
    cols = [col for col in df.columns if runs_col in col]
    non_na_dict = df[cols].notna().sum().to_dict()
    na_dict = df[cols].isna().sum().to_dict()

    plot_data = dict()
    for k, v in non_na_dict.items():
        plot_data[os.path.splitext(k)[0]] = {
            "Non-NA": v,
            "NA": na_dict[k],
        }
    return plot_data


def calculate_log_intensity(df, runs_col):
    """
    Calculate log intensity values for specified run columns.

    Parameters:
    -----------
    df: pandas DataFrame
        Input data
    runs_col: str
        Column name pattern for run columns

    Returns:
    --------
    pandas DataFrame
        DataFrame with log intensity values
    """
    cols = [col for col in df.columns if runs_col in col]
    df_copy = df.copy()
    for col in cols:
        df_copy[col] = np.log10(df_copy[col])
    return df_copy