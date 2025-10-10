"""
Functions for reading and processing mass spectrometry data files
"""

from __future__ import absolute_import
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
from pyopenms import MzMLFile, MSExperiment
from pyteomics import mzid

from pmultiqc.modules.common.file_utils import file_prefix

# Initialise the module logger via central logging
from pmultiqc.modules.common.logging import get_logger
log = get_logger("pmultiqc.modules.common.ms_io")

# The time resolution in seconds. Larger values produce smaller outputs and slight smoothing.
SECOND_RESOLUTION = 5

def get_ms_qc_info(ms_info: pd.DataFrame):
    """
    Compute MS QC summary data structures from MS info DataFrame.

    Note: Use min instead of mean to better expose major fluctuations (low intensity scans).

    Returns:
        (tic_data, bpc_data, ms1_peaks, general_stats)
    """
    ms1_info = ms_info[ms_info["ms_level"] == 1].copy()
    ms2_info = ms_info[ms_info["ms_level"] == 2].copy()
    ms1_info["rt_normalize"] = (ms1_info.sort_values(by="rt")["rt"] / SECOND_RESOLUTION).astype(int)

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
    enable_dia=False,
):
    """
    Process MS values from a dataframe row and add them to the appropriate histograms

    Args:
        info_df: Series row containing MS information
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
    # info_df is a Pandas.Seires not a DataFrame
    # "precursor_charge" --> "Charge",
    # "base_peak_intensity" --> "Base_Peak_Intensity",
    # "num_peaks": --> "MS_peaks"

    charge_state = (
        int(info_df["precursor_charge"]) if pd.notna(info_df["precursor_charge"]) else None
    )

    base_peak_intensity = (
        float(info_df["base_peak_intensity"]) if pd.notna(info_df["base_peak_intensity"]) else None
    )

    peak_per_ms2 = int(info_df["num_peaks"]) if pd.notna(info_df["num_peaks"]) else None

    if enable_dia:
        mzml_charge_plot.add_value(charge_state)
        mzml_peak_distribution_plot.add_value(base_peak_intensity)
        mzml_peaks_ms2_plot.add_value(peak_per_ms2)
        return

    if ms_name in ms_with_psm:
        # only "scan" in info_df not "SpectrumID"
        if info_df["scan"] in identified_spectrum_scan_id:
            mzml_charge_plot.add_value(charge_state)
            mzml_peak_distribution_plot.add_value(base_peak_intensity)
            mzml_peaks_ms2_plot.add_value(peak_per_ms2)
        else:
            mzml_charge_plot_1.add_value(charge_state)
            mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
            mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
    else:
        ms_without_psm.append(ms_name)


def read_mzmls(
    ms_paths,
    ms_with_psm,
    identified_spectrum,
    mzml_charge_plot,
    mzml_peak_distribution_plot,
    mzml_peaks_ms2_plot,
    mzml_charge_plot_1,
    mzml_peak_distribution_plot_1,
    mzml_peaks_ms2_plot_1,
    ms_without_psm,
    enable_dia=False,
    enable_mzid=False,
):
    """
    Read mzML files and extract information

    Args:
        ms_paths: List of paths to mzML files
        ms_with_psm: List of MS files with PSMs
        identified_spectrum: Dictionary of identified spectra by MS file
        mzml_charge_plot: Histogram for charge distribution of identified spectra
        mzml_peak_distribution_plot: Histogram for peak distribution of identified spectra
        mzml_peaks_ms2_plot: Histogram for peaks per MS2 of identified spectra
        mzml_charge_plot_1: Histogram for charge distribution of unidentified spectra
        mzml_peak_distribution_plot_1: Histogram for peak distribution of unidentified spectra
        mzml_peaks_ms2_plot_1: Histogram for peaks per MS2 of unidentified spectra
        ms_without_psm: List of MS files without PSMs
        enable_dia: Whether DIA mode is enabled
        enable_mzid: Whether mzid_plugin mode is enabled

    Returns:
        tuple: (mzml_table, heatmap_charge, total_ms2_spectra)
    """
    mzml_table = {}
    heatmap_charge = {}
    total_ms2_spectra = 0

    mzml_ms_dicts = list()
    for m in ms_paths:
        ms1_number = 0
        ms2_number = 0
        log.info("{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
        exp = MSExperiment()
        MzMLFile().load(m, exp)
        log.info("{}: Done parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
        m_name = file_prefix(m)
        log.info(
            "{}: Aggregating mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m_name)
        )

        charge_2 = 0
        for i in exp:
            if i.getMSLevel() == 1:
                ms1_number += 1
            elif i.getMSLevel() == 2:
                ms2_number += 1
                charge_state = i.getPrecursors()[0].getCharge()
                peaks_tuple = i.get_peaks()

                if enable_mzid:
                    # retention_time: minute
                    mzml_ms_dicts.append(
                        {
                            "spectrumID": i.getNativeID(),
                            "intensity": float(peaks_tuple[1].max()),
                            "retention_time": i.getRT() / 60,
                            "filename": m_name,
                        }
                    )

                peak_per_ms2 = len(peaks_tuple[0])
                if i.getMetaValue("base peak intensity"):
                    base_peak_intensity = i.getMetaValue("base peak intensity")
                else:
                    base_peak_intensity = max(peaks_tuple[1]) if len(peaks_tuple[1]) > 0 else None

                if charge_state == 2:
                    charge_2 += 1

                if enable_dia:
                    mzml_charge_plot.add_value(charge_state)
                    mzml_peak_distribution_plot.add_value(base_peak_intensity)
                    mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                    continue

                if m_name in ms_with_psm:
                    if i.getNativeID() in identified_spectrum[m_name]:
                        mzml_charge_plot.add_value(charge_state)
                        mzml_peak_distribution_plot.add_value(base_peak_intensity)
                        mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                    else:
                        mzml_charge_plot_1.add_value(charge_state)
                        mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                        mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
                else:
                    if m_name not in ms_without_psm:
                        ms_without_psm.append(m_name)

        heatmap_charge[m_name] = charge_2 / ms2_number if ms2_number > 0 else 0
        total_ms2_spectra += ms2_number
        mzml_table[m_name] = {"MS1_Num": ms1_number}
        mzml_table[m_name]["MS2_Num"] = ms2_number
        log.info(
            "{}: Done aggregating mzML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), m_name
            )
        )

    if enable_mzid:
        return mzml_table, heatmap_charge, total_ms2_spectra, pd.DataFrame(mzml_ms_dicts)
    else:
        return mzml_table, heatmap_charge, total_ms2_spectra


## Removed read_ms_info function; use MsInfoReader class in ms/msinfo.py instead.



def read_mzids(file_paths):

    def parse_location(location):
        if "\\" in location:
            location = location.replace("\\", "/")
        return os.path.basename(location)

    def process_modification(modification):
        if not isinstance(modification, list):
            modifications = None
        else:
            modifi_list = list()
            for i in modification:
                if i.get("name", None) is not None:
                    modifi_list.append(str(i.get("location")) + "-" + i.get("name", None))
                elif i.get("cross-link receiver", None) is not None:
                    modifi_list.append(str(i.get("location")) + "-CrossLinkReceiver")
            modifications = ";".join(modifi_list)
        return modifications

    mzid_dicts = list()
    for file_path in file_paths:

        log.info(
            "{}: Parsing MzIdentML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), file_path
            )
        )

        mzid_data = mzid.MzIdentML(file_path)

        if len(mzid_data) == 0:
            raise ValueError("Please check your MzIdentML", file_path)

        log.info(
            "{}: Done parsing MzIdentML file {}.".format(
                datetime.now().strftime("%H:%M:%S"), file_path
            )
        )
        m = file_prefix(file_path)
        log.info(
            "{}: Aggregating MzIdentML file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
        )

        for mzid_tmp in mzid_data:
            mzid_tmp_part = {
                k: v for k, v in mzid_tmp.items() if k not in ["SpectrumIdentificationItem"]
            }

            for spectrum_item in mzid_tmp.get("SpectrumIdentificationItem", []):
                spectrum_item_part = {
                    k: v
                    for k, v in spectrum_item.items()
                    if k not in ["PeptideEvidenceRef", "PeptideSequence"]
                }

                rank = spectrum_item_part.get("rank")
                peptide_pass = spectrum_item_part.get("peptide passes threshold", "true") == "true"
                pass_threshold = spectrum_item.get("passThreshold", False)

                if rank != 1 or not peptide_pass or not pass_threshold:
                    continue

                for peptide_ref in spectrum_item.get("PeptideEvidenceRef", []):

                    if "name" in peptide_ref:
                        peptide_ref["PeptideEvidenceRef_name"] = peptide_ref.pop("name")
                    if "location" in peptide_ref:
                        peptide_ref["PeptideEvidenceRef_location"] = peptide_ref.pop("location")
                    if "FileFormat" in peptide_ref:
                        peptide_ref["PeptideEvidenceRef_FileFormat"] = peptide_ref.pop(
                            "FileFormat"
                        )

                    mzid_dict = {
                        **mzid_tmp_part,
                        **spectrum_item_part,
                        **peptide_ref,
                        "mzid_file_name": m,
                    }

                    need_keys = [
                        "SEQUEST:xcorr",
                        "Mascot:score",
                        "PEAKS:peptideScore",
                        "xi:score",
                        "retention time",
                        "location",
                        "Modification",
                        "spectrumID",
                        "isDecoy",
                        "accession",
                        "PeptideSequence",
                        "experimentalMassToCharge",
                        "calculatedMassToCharge",
                        "chargeState",
                        "mzid_file_name",
                        "Andromeda:score",
                    ]
                    mzid_dicts.append({k: v for k, v in mzid_dict.items() if k in need_keys})
        log.info(
            "{}: Done aggregating MzIdentML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), m
            )
        )

    mzid_df = pd.DataFrame(mzid_dicts)

    # Check columns
    check_list = [
        "spectrumID",
        "PeptideSequence",
        "chargeState",
        "accession",
        "Modification",
        "experimentalMassToCharge",
        "calculatedMassToCharge",
    ]
    missing_cols = [col for col in check_list if col not in mzid_df.columns]
    if missing_cols:
        log.warning(f"MzIdentML file is missing required fields: {missing_cols}")

    # Filter out contaminants: "Cont"
    filtered_mzid_df = mzid_df[~mzid_df["accession"].str.lower().str.startswith("cont")].copy()

    search_engines = [
        "SEQUEST:xcorr",
        "Mascot:score",
        "PEAKS:peptideScore",
        "xi:score",
        "Andromeda:score",
    ]
    filtered_mzid_df.rename(
        columns=lambda x: "search_engine_score" if x in search_engines else x, inplace=True
    )

    if "search_engine_score" not in filtered_mzid_df.columns:
        log.warning("Please check the 'search_engine_score' field in the mzIdentML file.")

    if "retention time" in filtered_mzid_df.columns:
        filtered_mzid_df.rename(columns={"retention time": "retention_time"}, inplace=True)

    if "location" in filtered_mzid_df.columns:
        filtered_mzid_df["location"] = filtered_mzid_df["location"].apply(parse_location)

    filtered_mzid_df["Modifications"] = filtered_mzid_df["Modification"].apply(
        lambda x: process_modification(x)
    )

    filtered_mzid_df["accession_group"] = filtered_mzid_df.groupby(
        ["spectrumID", "PeptideSequence"]
    )["accession"].transform(lambda x: ";".join(x.unique()))

    if "isDecoy" not in filtered_mzid_df.columns:
        filtered_mzid_df["isDecoy"] = False

    # location: path of mzML file
    if "location" in filtered_mzid_df.columns:
        filtered_mzid_df["filename"] = filtered_mzid_df.apply(
            lambda x: file_prefix(x.location), axis=1
        )
    else:
        filtered_mzid_df["filename"] = filtered_mzid_df["mzid_file_name"]

    return filtered_mzid_df


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