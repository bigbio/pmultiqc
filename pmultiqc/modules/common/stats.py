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
    Calculate quality score based on uniformity of retention time distribution.

    Parameters:
    -----------
    group_df_rt: group["Retention time"] or group["retention_time"]

    Returns:
    --------
    float: Quality score between 0 and 1, where 1 indicates perfect uniformity.
    """
    n = group_df_rt.notna().sum()
    if n == 0:
        return 0.0

    total_sum = np.nansum(group_df_rt)
    if total_sum == 0:
        return 0.0

    x = group_df_rt / total_sum
    y = np.nansum(x) / n
    worst = ((1 - y) ** 0.5) * 1 / n + (y**0.5) * (n - 1) / n
    sc = np.sum(np.abs(x - y) ** 0.5) / n
    result = 1.0 if worst == 0 else float((worst - sc) / worst)

    return result


def cal_delta_mass_dict(df, col, num_bins: int = 1000):
    """
    Calculate delta mass distribution as counts and frequencies.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing the mass delta column.
    col : str
        Name of the column containing mass delta values.
    num_bins : int, optional
        Number of bins for histogram (default: 1000).

    Returns:
    --------
    dict: Dictionary with 'count' and 'frequency' keys containing binned data.
    """
    # Compute value_counts once and derive frequency from counts
    count_bin = df[col].value_counts(sort=False, bins=num_bins)

    # Build count dictionary
    count_bin_data = {
        float(interval.mid): int(count)
        for interval, count in count_bin.items()
    }

    # Derive frequency from counts (more efficient than calling value_counts twice)
    total_count = count_bin.sum()
    if total_count > 0:
        frequency_bin_data = {
            float(interval.mid): float(count / total_count)
            for interval, count in count_bin.items()
        }
    else:
        frequency_bin_data = {k: 0.0 for k in count_bin_data.keys()}

    return {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }


def cal_hm_charge(df: pd.DataFrame, run_col: str, charge_col: str):

    if run_col not in df.columns or charge_col not in df.columns:
        return {}

    charge = dict()

    for raw_file, group in df[[run_col, charge_col]].groupby(run_col):
        vc = group[charge_col].value_counts(dropna=True, normalize=True)
        charge[raw_file] = vc.get(2, 0)
    charge_median = np.median(list(charge.values())) if charge else 0

    return {
        k: float(max(0.0, 1.0 - abs(v - charge_median)))
        for k, v in charge.items()
    }
