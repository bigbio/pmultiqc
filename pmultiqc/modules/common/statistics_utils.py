## statistic functions


import numpy as np


def nanmedian(values : np.ndarray, all_nan_fallback : np.float64) -> np.float64:
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
