import numpy as np


def qualUniform(group_df_rt):
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
