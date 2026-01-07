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
    Parameters:
    -----------
    group_df_rt: group["Retention time"] or group["retention_time"]

    """
    x = group_df_rt / np.nansum(group_df_rt)
    n = group_df_rt.notna().sum()
    y = np.nansum(x) / n
    worst = ((1 - y) ** 0.5) * 1 / n + (y**0.5) * (n - 1) / n
    sc = np.sum(np.abs(x - y) ** 0.5) / n
    result = 1.0 if worst == 0 else float((worst - sc) / worst)

    return result


def cal_delta_mass_dict(df, col):

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
