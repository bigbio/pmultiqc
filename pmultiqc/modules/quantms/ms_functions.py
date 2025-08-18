"""
======================
@author:Chengxin Dai

@time:2023/10/23:17:25

======================
"""

import pandas as pd
import numpy as np

# The time resolution in seconds.
# Larger values will result in smaller data files as outputs
# and will slightly smooth the data. 5 seconds seems to be
# a good value for qc purposes.
SECOND_RESOLUTION = 5


def get_ms_qc_info(ms_info: pd.DataFrame):
    """
    Note that here I am using min and not mean for purely qc reasons.
    Since the diagnostic aspect here is mainly to see major fluctuations
    in the intensity, and usually these are scans with very low intensity
    due to bubbles or ionization issues; thus the mean would hide that.

    @param ms_info:
    @return:
    """
    ms_info = ms_info.rename(
        columns={
            "ms_level": "MSLevel",
            "rt": "Retention_Time",
            "summed_peak_intensities": "Summed_Peak_Intensities",
            "num_peaks": "MS_peaks",
            "acquisition_datetime": "AcquisitionDateTime",
        }
    )

    ms1_info = ms_info[ms_info["MSLevel"] == 1].copy()
    ms2_info = ms_info[ms_info["MSLevel"] == 2].copy()
    ms1_info["rt_normalize"] = (
        ms1_info.sort_values(by="Retention_Time")["Retention_Time"] / SECOND_RESOLUTION
    ).astype(int)
    tic_data = ms1_info.groupby("rt_normalize")[
        ["Retention_Time", "Summed_Peak_Intensities"]
    ].min()
    tic_data = dict(zip(tic_data["Retention_Time"], tic_data["Summed_Peak_Intensities"]))

    bpc_data = dict(
        zip(
            ms1_info.groupby("rt_normalize")["Retention_Time"].min(),
            ms1_info.groupby("rt_normalize")["Summed_Peak_Intensities"].max(),
        )
    )

    ms1_peaks = dict(
        zip(
            ms1_info.groupby("rt_normalize")["Retention_Time"].min(),
            ms1_info.groupby("rt_normalize")["MS_peaks"].mean(),
        )
    )

    general_stats = {
        "AcquisitionDateTime": ms1_info["AcquisitionDateTime"][0],
        "log10(TotalCurrent)": np.log10(ms1_info["Summed_Peak_Intensities"].sum()),
        "log10(ScanCurrent)": calculate_scan_current_log10(ms2_info),
    }

    return tic_data, bpc_data, ms1_peaks, general_stats


def calculate_scan_current_log10(ms2_info):
    """
    Calculate log10(ScanCurrent) with better validation to avoid extremely high values
    as mentioned in the DIANN improvements issue.
    """
    if ms2_info.empty or "Summed_Peak_Intensities" not in ms2_info.columns:
        return 0.0
    
    # Filter out zero, negative, inf, or NaN values
    valid_intensities = ms2_info["Summed_Peak_Intensities"]
    valid_intensities = valid_intensities[
        (valid_intensities > 0) & 
        (~np.isinf(valid_intensities)) & 
        (~np.isnan(valid_intensities))
    ]
    
    if valid_intensities.empty:
        return 0.0
    
    total_intensity = valid_intensities.sum()
    
    # If sum results in extremely large log10 values (>30), use mean-based approach
    if total_intensity == 0:
        return 0.0
    
    log_total = np.log10(total_intensity)
    if log_total > 30:  # This would be log10 > 30, corresponding to 10^30
        # Use average intensity per scan instead of total
        avg_intensity = valid_intensities.mean()
        num_scans = len(valid_intensities)
        # This gives a more reasonable scale: log10(avg_intensity) + log10(num_scans)
        if avg_intensity > 0:
            return np.log10(avg_intensity) + np.log10(num_scans)
        else:
            return 0.0
    
    return log_total
