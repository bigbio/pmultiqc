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