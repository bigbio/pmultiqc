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

    # Combine multiple groupby operations into a single aggregation for better performance
    grouped_agg = ms1_info.groupby("rt_normalize").agg(
        rt_min=("rt", "min"),
        summed_min=("summed_peak_intensities", "min"),
        summed_max=("summed_peak_intensities", "max"),
        num_peaks_mean=("num_peaks", "mean"),
    )

    tic_data = dict(zip(grouped_agg["rt_min"], grouped_agg["summed_min"], strict=True))
    bpc_data = dict(zip(grouped_agg["rt_min"], grouped_agg["summed_max"], strict=True))
    ms1_peaks = dict(zip(grouped_agg["rt_min"], grouped_agg["num_peaks_mean"], strict=True))

    total_curr = float(ms1_info["summed_peak_intensities"].sum())
    scan_curr = float(ms2_info["summed_peak_intensities"].sum())
    general_stats = {
        "AcquisitionDateTime": ms1_info["acquisition_datetime"].iloc[0],
        "log10(TotalCurrent)": float(np.log10(max(total_curr, 1e-12))),
        "log10(ScanCurrent)": float(np.log10(max(scan_curr, 1e-12))),
    }
    current_sum = {
        "total_curr": total_curr,
        "scan_curr": scan_curr
    }

    return tic_data, bpc_data, ms1_peaks, general_stats, current_sum

def add_ms_values(
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

    if enable_dia:
        # Use batch method for better performance
        mzml_charge_plot.add_values_batch(info_df["precursor_charge"].dropna())
        mzml_peak_distribution_plot.add_values_batch(info_df["base_peak_intensity"].dropna())
        mzml_peaks_ms2_plot.add_values_batch(info_df["num_peaks"].dropna())
        return

    if ms_name not in ms_with_psm:
        ms_without_psm.add(ms_name)
        return

    identified_scans = info_df["scan"].isin(identified_spectrum_scan_id)
    identified_df = info_df[identified_scans]
    unidentified_df = info_df[~identified_scans]

    # Use batch method for identified spectra
    mzml_charge_plot.add_values_batch(identified_df["precursor_charge"].dropna())
    mzml_peak_distribution_plot.add_values_batch(identified_df["base_peak_intensity"].dropna())
    mzml_peaks_ms2_plot.add_values_batch(identified_df["num_peaks"].dropna())

    # Use batch method for unidentified spectra
    mzml_charge_plot_1.add_values_batch(unidentified_df["precursor_charge"].dropna())
    mzml_peak_distribution_plot_1.add_values_batch(unidentified_df["base_peak_intensity"].dropna())
    mzml_peaks_ms2_plot_1.add_values_batch(unidentified_df["num_peaks"].dropna())

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


def get_ms_long_trends(df: pd.DataFrame):

    if df.empty or "acquisition_datetime" not in df.columns:
        log.warning("No acquisition_datetime found; skipping long trends.")
        return {}

    ac_time = df["acquisition_datetime"].iloc[0]

    ms1_rts = df.loc[df["ms_level"] == 1, "rt"]
    if not ms1_rts.empty and (m_rt := ms1_rts.median()) > 0:
        rt = float(m_rt / 60)
    else:
        rt = 0

    ms2_pre_int = df.loc[df["ms_level"] == 2, "precursor_intensity"].dropna()
    if not ms2_pre_int.empty and (m_pre_int := ms2_pre_int.median()) > 0:
        prec_intensity = float(np.log2(m_pre_int))
    else:
        prec_intensity = 0

    ms1_sum_int = df.loc[df["ms_level"] == 1, "summed_peak_intensities"].dropna()
    if not ms1_sum_int.empty and (m_sum_int := ms1_sum_int.median()) > 0:
        ms1_intensity = float(np.log2(m_sum_int))
    else:
        ms1_intensity = 0

    return {
        "time": {"acquisition_datetime": ac_time},
        "rt": {ac_time: rt},
        "ms2_prec_intensity": {ac_time: prec_intensity},
        "ms1_summed_intensity": {ac_time: ms1_intensity}
    }

def process_long_trends(df, m_name, trends_data):

    run_trends = get_ms_long_trends(df=df)
    if run_trends:
        trends_data["time"][m_name] = run_trends["time"]
        trends_data["rt"].update(run_trends["rt"])
        trends_data["ms2_prec_intensity"].update(run_trends["ms2_prec_intensity"])
        trends_data["ms1_summed_intensity"].update(run_trends["ms1_summed_intensity"])
