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
    ms1_info = ms_info[ms_info["ms_level"] == 1].copy()
    ms2_info = ms_info[ms_info["ms_level"] == 2].copy()
    ms1_info["rt_normalize"] = (
        ms1_info.sort_values(by="rt")["rt"] / SECOND_RESOLUTION
    ).astype(int)
    tic_data = ms1_info.groupby("rt_normalize")[
        ["rt", "summed_peak_intensities"]
    ].min()
    tic_data = dict(zip(tic_data["rt"], tic_data["summed_peak_intensities"]))

    bpc_data = dict(
        zip(
            ms1_info.groupby("rt_normalize")["rt"].min(),
            ms1_info.groupby("rt_normalize")["summed_peak_intensities"].max(),
        )
    )

    ms1_peaks = dict(
        zip(
            ms1_info.groupby("rt_normalize")["rt"].min(),
            ms1_info.groupby("rt_normalize")["num_peaks"].mean(),
        )
    )

    general_stats = {
        "AcquisitionDateTime": ms1_info["acquisition_datetime"][0],
        "log10(TotalCurrent)": np.log10(ms1_info["summed_peak_intensities"].sum()),
        "log10(ScanCurrent)": np.log10(ms2_info["summed_peak_intensities"].sum()),
    }

    return tic_data, bpc_data, ms1_peaks, general_stats
