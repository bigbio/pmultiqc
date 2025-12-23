"""
Functions for reading and processing mass spectrometry data files
"""

from __future__ import absolute_import
import os
import re

import numpy as np
import pandas as pd

# Initialise the module logger via central logging
from pmultiqc.modules.common.logging import get_logger


log = get_logger("pmultiqc.modules.common.ms_io")

# The time resolution in seconds. Larger values produce smaller outputs and slight smoothing.
SECOND_RESOLUTION = 5

def get_ms_qc_info(ms_info: pd.DataFrame):
    """
    Compute MS QC summary data structures from MS info DataFrame (or from mzML DataFrame).

    Note: Use min instead of mean to better expose major fluctuations (low intensity scans).

    Returns:
        (tic_data, bpc_data, ms1_peaks, general_stats)
    """
    ms1_info = ms_info[ms_info["ms_level"] == 1].copy()
    ms2_info = ms_info[ms_info["ms_level"] == 2].copy()

    if ms1_info.empty:
        log.warning("MS1 MS info DataFrame is empty. 'MS1 analysis' cannot be generated.")

        return None, None, None, None

    ms1_info["rt_normalize"] = (ms1_info.sort_values(by="rt")["rt"] / SECOND_RESOLUTION).astype(int)

    tic_data = ms1_info.groupby("rt_normalize")[
        ["rt", "summed_peak_intensities"]
    ].min()
    tic_data = dict(
        zip(
            tic_data["rt"],
            tic_data["summed_peak_intensities"],
            strict=True
        )
    )

    bpc_data = dict(
        zip(
            ms1_info.groupby("rt_normalize")["rt"].min(),
            ms1_info.groupby("rt_normalize")["summed_peak_intensities"].max(),
            strict=True
        )
    )

    ms1_peaks = dict(
        zip(
            ms1_info.groupby("rt_normalize")["rt"].min(),
            ms1_info.groupby("rt_normalize")["num_peaks"].mean(),
            strict=True
        )
    )

    total_curr = float(ms1_info["summed_peak_intensities"].sum())
    scan_curr = float(ms2_info["summed_peak_intensities"].sum())
    general_stats = {
        "AcquisitionDateTime": ms1_info["acquisition_datetime"].iloc[0],
        "log10(TotalCurrent)": np.log10(max(total_curr, 1e-12)),
        "log10(ScanCurrent)": np.log10(max(scan_curr, 1e-12)),
    }

    return tic_data, bpc_data, ms1_peaks, general_stats

def add_ms_values_df(
        info_df,
        ms_name,
        ms_with_psm,
        identified_spectrum_scan_id,
        mzml_charge_plot,
        mzml_peak_distribution_plot,
        mzml_peaks_ms2_plot,
        mzml_charge_plot_1,
        mzml_peak_distribution_plot_1,
        mzml_peaks_ms2_plot_1,
        ms_without_psm,
        enable_dia: bool = False,
    ):
        """
        Process MS values from a dataframe row and add them to the appropriate histograms

        Args:
            info_df: MS information
            ms_name: Name of the MS file
            ms_with_psm: List of MS files with PSMs
            identified_spectrum_scan_id: List of identified spectra by MS file
            # identified_spectrum: Dictionary of identified spectra by MS file
            mzml_charge_plot: Histogram for charge distribution of identified spectra
            mzml_peak_distribution_plot: Histogram for peak distribution of identified spectra
            mzml_peaks_ms2_plot: Histogram for peaks per MS2 of identified spectra
            mzml_charge_plot_1: Histogram for charge distribution of unidentified spectra
            mzml_peak_distribution_plot_1: Histogram for peak distribution of unidentified spectra
            mzml_peaks_ms2_plot_1: Histogram for peaks per MS2 of unidentified spectra
            ms_without_psm: List of MS files without PSMs
            enable_dia: Whether DIA mode is enabled
        """

        info_df["precursor_charge"] = info_df["precursor_charge"].astype("Int64")
        info_df["base_peak_intensity"] = info_df["base_peak_intensity"].astype("float")
        info_df["num_peaks"] = info_df["num_peaks"].astype("Int64")

        def data_add_value(histogram_data, a_value):
            if not pd.notna(a_value):
                a_value = None
            histogram_data.add_value(a_value)

        if enable_dia:
            for charge_state in info_df["precursor_charge"]:
                data_add_value(mzml_charge_plot, charge_state)
            for base_peak_inte in info_df["base_peak_intensity"]:
                data_add_value(mzml_peak_distribution_plot, base_peak_inte)
            for num_peaks in info_df["num_peaks"]:
                data_add_value(mzml_peaks_ms2_plot, num_peaks)
            return

        if ms_name not in ms_with_psm:
            ms_without_psm.add(ms_name)
            return

        identified_scans = info_df["scan"].isin(identified_spectrum_scan_id)

        for charge_state in info_df[identified_scans]["precursor_charge"]:
            data_add_value(mzml_charge_plot, charge_state)
        for base_peak_inte in info_df[identified_scans]["base_peak_intensity"]:
            data_add_value(mzml_peak_distribution_plot, base_peak_inte)
        for peak_per_ms2 in info_df[identified_scans]["num_peaks"]:
            data_add_value(mzml_peaks_ms2_plot, peak_per_ms2)

        for charge_state in info_df[~identified_scans]["precursor_charge"]:
            data_add_value(mzml_charge_plot_1, charge_state)
        for base_peak_inte in info_df[~identified_scans]["base_peak_intensity"]:
            data_add_value(mzml_peak_distribution_plot_1, base_peak_inte)
        for peak_per_ms2 in info_df[~identified_scans]["num_peaks"]:
            data_add_value(mzml_peaks_ms2_plot_1, peak_per_ms2)


def spectra_ref_check(spectra_ref):
    match_scan = re.search(r"scan=(\d+)", spectra_ref)
    if match_scan:
        return match_scan.group(1)

    match_spectrum = re.search(r"spectrum=(\d+)", spectra_ref)
    if match_spectrum:
        return match_spectrum.group(1)

    try:
        if int(spectra_ref):
            return spectra_ref

    except ValueError:
        raise ValueError("Please check the 'spectra_ref' field in your mzTab file.")


# Remove output files generated by OpenMS().openms_convert after running
def del_openms_convert_tsv():

    files = ["experimental_design.tsv", "openms.tsv"]

    for file_path in files:
        if os.path.exists(file_path):
            os.remove(file_path)
            log.info(f"{file_path} has been deleted.")