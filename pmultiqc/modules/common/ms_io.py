#!/usr/bin/env python

"""
Functions for reading and processing mass spectrometry data files
"""

from __future__ import absolute_import
import os
import logging
import pandas as pd
import numpy as np
from datetime import datetime
from collections import OrderedDict
from pyopenms import IdXMLFile, MzMLFile, MSExperiment
import math
import re

from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.quantms.ms_functions import get_ms_qc_info

# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def file_prefix(path):
    """Extract the file prefix from a path"""
    try:
        return os.path.splitext(os.path.basename(path))[0]
    except:
        raise SystemExit(f"Illegal file path: {path}")


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
        int(info_df["precursor_charge"])
        if pd.notna(info_df["precursor_charge"])
        else None
    )

    base_peak_intensity = (
        float(info_df["base_peak_intensity"])
        if pd.notna(info_df["base_peak_intensity"])
        else None
    )

    peak_per_ms2 = (
        int(info_df["num_peaks"])
        if pd.notna(info_df["num_peaks"])
        else None
    )

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

    Returns:
        tuple: (mzml_table, heatmap_charge, total_ms2_spectra)
    """
    mzml_table = {}
    heatmap_charge = {}
    total_ms2_spectra = 0

    exp = MSExperiment()
    for m in ms_paths:
        ms1_number = 0
        ms2_number = 0
        log.info("{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
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

    return mzml_table, heatmap_charge, total_ms2_spectra


def read_ms_info(
        ms_info_path,
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
):
    """
    Read MS info files and extract information

    Args:
        ms_info_path: List of paths to MS info files
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

    Returns:
        tuple: (mzml_table, heatmap_charge, total_ms2_spectra, ms1_tic, ms1_bpc, ms1_peaks, ms1_general_stats)
    """
    mzml_table = {}
    heatmap_charge = {}
    total_ms2_spectra = 0
    ms1_tic = {}
    ms1_bpc = {}
    ms1_peaks = {}
    ms1_general_stats = {}

    for file in ms_info_path:
        log.info(
            "{}: Parsing ms_statistics dataframe {}...".format(
                datetime.now().strftime("%H:%M:%S"), file
            )
        )
        mzml_df = pd.read_parquet(file)
        m = file_prefix(file).replace("_ms_info", "")
        if m not in mzml_table:
            mzml_table[m] = dict.fromkeys(["MS1_Num", "MS2_Num", "Charge_2"], 0)
        charge_group = mzml_df.groupby("precursor_charge").size()
        ms_level_group = mzml_df.groupby("ms_level").size()
        charge_2 = charge_group[2] if 2 in charge_group else 0
        ms1_number = int(ms_level_group[1]) if 1 in ms_level_group else 0
        ms2_number = int(ms_level_group[2]) if 2 in ms_level_group else 0
        total_ms2_spectra += ms2_number
        mzml_table[m].update({"MS1_Num": mzml_table[m]["MS1_Num"] + ms1_number})
        mzml_table[m].update({"MS2_Num": mzml_table[m]["MS2_Num"] + ms2_number})
        mzml_table[m].update({"Charge_2": mzml_table[m]["Charge_2"] + charge_2})

        (
            ms1_tic[os.path.basename(file).replace("_ms_info.parquet", "")],
            ms1_bpc[os.path.basename(file).replace("_ms_info.parquet", "")],
            ms1_peaks[os.path.basename(file).replace("_ms_info.parquet", "")],
            ms1_general_stats[os.path.basename(file).replace("_ms_info.parquet", "")],
        ) = get_ms_qc_info(mzml_df)

        group = mzml_df[mzml_df["ms_level"] == 2]
        del mzml_df
        
        if enable_dia:
            identified_spectrum_scan_id = []
        else:
            if not identified_spectrum:
                raise ValueError("ms_io: The identified_spectrum is missing. Please check your mzTab file!")
            identified_spectrum_scan_id = [
                spectra_ref_check(spectrum_id)
                for spectrum_id in identified_spectrum[m]
            ]

        # Apply add_ms_values to each row
        for _, row in group.iterrows():
            add_ms_values(
                row,
                m,
                ms_with_psm,
                identified_spectrum_scan_id,
                mzml_charge_plot,
                mzml_peak_distribution_plot,
                mzml_peaks_ms2_plot,
                mzml_charge_plot_1,
                mzml_peak_distribution_plot_1,
                mzml_peaks_ms2_plot_1,
                ms_without_psm,
                enable_dia,
            )

        for m in mzml_table.keys():
            if mzml_table[m]["MS2_Num"] > 0:
                heatmap_charge[m] = mzml_table[m]["Charge_2"] / mzml_table[m]["MS2_Num"]
            else:
                heatmap_charge[m] = 0

        log.info(
            "{}: Done aggregating ms_statistics dataframe {}...".format(
                datetime.now().strftime("%H:%M:%S"), file
            )
        )

    return (
        mzml_table,
        heatmap_charge,
        total_ms2_spectra,
        ms1_tic,
        ms1_bpc,
        ms1_peaks,
        ms1_general_stats,
    )


def parse_idxml(
        idx_paths,
        mzml_table,
        xcorr_hist_range,
        hyper_hist_range,
        spec_evalue_hist_range,
        pep_hist_range,
        ml_spec_ident_final,
        mzml_peptide_map,
        remove_decoy=True,
):
    """
    Parse idXML files and extract information

    Args:
        idx_paths: List of paths to idXML files
        mzml_table: Dictionary of mzML information
        xcorr_hist_range: Range for xcorr histogram
        hyper_hist_range: Range for hyper histogram
        spec_evalue_hist_range: Range for spec_evalue histogram
        pep_hist_range: Range for PEP histogram
        ml_spec_ident_final: Dictionary of identified spectra counts
        mzml_peptide_map: Dictionary of peptide maps
        remove_decoy: Whether to remove decoy hits

    Returns:
        tuple: (search_engine, MSGF_label, Comet_label, Sage_label)
    """
    consensus_paths = []
    for raw_id in idx_paths:
        if "consensus" in os.path.split(raw_id)[1]:
            consensus_paths.append(raw_id)

    for raw_id in consensus_paths:
        if raw_id in idx_paths:
            idx_paths.remove(raw_id)

    msgf_label, comet_label, sage_label = False, False, False
    search_engine = {
        "SpecE": OrderedDict(),
        "xcorr": OrderedDict(),
        "hyper": OrderedDict(),
        "PEPs": OrderedDict(),
        "consensus_support": OrderedDict(),
        "data_label": OrderedDict(),
    }
    spec_e_label, xcorr_label, hyper_label, peps_label, consensus_label = [], [], [], [], []

    for raw_id in idx_paths:
        log.info(
            "{}: Parsing search result file {}...".format(
                datetime.now().strftime("%H:%M:%S"), raw_id
            )
        )

        protein_ids = []
        peptide_ids = []
        IdXMLFile().load(raw_id, protein_ids, peptide_ids)
        raw_id_name = file_prefix(raw_id)

        if remove_decoy:
            identified_num = len(
                set(
                    [
                        i.getMetaValue("spectrum_reference")
                        for i in peptide_ids
                        if i.getHits()[0].getMetaValue("target_decoy") == "target"
                    ]
                )
            )
        else:
            identified_num = len(peptide_ids)

        ms_name = file_prefix(protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8"))
        search_engine_name = protein_ids[0].getSearchEngine()

        search_engine["SpecE"][raw_id_name] = OrderedDict()
        search_engine["xcorr"][raw_id_name] = OrderedDict()
        search_engine["hyper"][raw_id_name] = OrderedDict()
        search_engine["PEPs"][raw_id_name] = OrderedDict()

        xcorr_breaks = list(
            np.arange(
                xcorr_hist_range["start"],
                xcorr_hist_range["end"] + xcorr_hist_range["step"],
                xcorr_hist_range["step"],
            ).round(2)
        )

        hyper_breaks = list(
            np.arange(
                hyper_hist_range["start"],
                hyper_hist_range["end"] + hyper_hist_range["step"],
                hyper_hist_range["step"],
            ).round(2)
        )

        spec_e_breaks = list(
            np.arange(
                spec_evalue_hist_range["start"],
                spec_evalue_hist_range["end"] + spec_evalue_hist_range["step"],
                spec_evalue_hist_range["step"],
            ).round(2)
        )
        spec_e_breaks.append(float("inf"))
        spec_e_breaks.sort()

        pep_breaks = list(
            np.concatenate(
                [
                    np.arange(pep_hist_range["start"], pep_hist_range["low_thresh"], pep_hist_range["low_step"]),
                    np.arange(pep_hist_range["low_thresh"], pep_hist_range["high_thresh"], pep_hist_range["high_step"]),
                    np.arange(pep_hist_range["high_thresh"], pep_hist_range["end"] + 0.01, pep_hist_range["low_step"])
                ]
            ).round(2))

        bar_stacks = ["target", "decoy", "target+decoy"]
        cross_corr = Histogram(
            "Comet cross-correlation score",
            plot_category="range",
            stacks=bar_stacks,
            breaks=xcorr_breaks,
        )
        hyper = Histogram(
            "Sage hyperscore", plot_category="range", stacks=bar_stacks, breaks=hyper_breaks
        )
        spectral_e = Histogram(
            "MSGF spectral E-value",
            plot_category="range",
            stacks=bar_stacks,
            breaks=spec_e_breaks,
        )
        posterior_error = Histogram(
            "Posterior error probability",
            plot_category="range",
            stacks=bar_stacks,
            breaks=pep_breaks,
        )

        if search_engine_name == "MSGF+" or "msgf" in raw_id_name:
            mzml_table[ms_name]["MSGF"] = identified_num
            msgf_label = True
            spec_e_label.append({"name": raw_id_name, "ylab": "Counts"})
            peps_label.append({"name": raw_id_name, "ylab": "Counts"})
            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    spec_e = (
                        hit.getMetaValue("SpecEvalue-score")
                        if hit.getMetaValue("SpecEvalue-score")
                        else hit.getMetaValue("MS:1002052")
                    )
                    log_spec_e = -math.log(spec_e, 10)
                    pep = (
                        hit.getMetaValue("MS:1001493")
                        if hit.getMetaValue("MS:1001493")
                        else hit.getScore()
                    )
                    spectral_e.add_value(log_spec_e, stack=hit.getMetaValue("target_decoy"))
                    posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

            spectral_e.to_dict()
            posterior_error.to_dict()
            search_engine["SpecE"][raw_id_name] = spectral_e.dict["data"]
            search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

        elif search_engine_name == "Comet" or "comet" in raw_id_name:
            comet_label = True
            mzml_table[ms_name]["Comet"] = identified_num
            xcorr_label.append({"name": raw_id_name, "ylab": "Counts"})
            peps_label.append({"name": raw_id_name, "ylab": "Counts"})
            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    xcorr = hit.getMetaValue("MS:1002252")
                    pep = (
                        hit.getMetaValue("MS:1001493")
                        if hit.getMetaValue("MS:1001493")
                        else hit.getScore()
                    )
                    cross_corr.add_value(xcorr, stack=hit.getMetaValue("target_decoy"))
                    posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

            cross_corr.to_dict()
            posterior_error.to_dict()
            search_engine["xcorr"][raw_id_name] = cross_corr.dict["data"]
            search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

        elif search_engine_name == "Sage" or "sage" in raw_id_name:
            sage_label = True
            mzml_table[ms_name]["Sage"] = identified_num
            hyper_label.append({"name": raw_id_name, "ylab": "Counts"})
            peps_label.append({"name": raw_id_name, "ylab": "Counts"})
            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    hyper_score = hit.getMetaValue("hyperscore")
                    pep = (
                        hit.getMetaValue("MS:1001493")
                        if hit.getMetaValue("MS:1001493")
                        else hit.getScore()
                    )
                    hyper.add_value(hyper_score, stack=hit.getMetaValue("target_decoy"))
                    posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

            hyper.to_dict()
            posterior_error.to_dict()
            search_engine["hyper"][raw_id_name] = hyper.dict["data"]
            search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

        else:
            mzml_table[ms_name][search_engine_name] = identified_num

        mzml_table[ms_name]["num_quant_psms"] = (
            ml_spec_ident_final[ms_name] if ms_name in ml_spec_ident_final.keys() else 0
        )
        mzml_table[ms_name]["num_quant_peps"] = (
            len(mzml_peptide_map[ms_name]) if ms_name in ml_spec_ident_final.keys() else 0
        )

    for raw_id in consensus_paths:
        log.info(
            "{}: Parsing consensus file {}...".format(
                datetime.now().strftime("%H:%M:%S"), format(raw_id)
            )
        )
        protein_ids = []
        peptide_ids = []
        IdXMLFile().load(raw_id, protein_ids, peptide_ids)
        raw_id_name = file_prefix(raw_id)

        consensus_label.append({"name": raw_id_name, "ylab": "Counts"})

        consensus_support = Histogram(
            "Consensus PSM number", plot_category="frequency", stacks=bar_stacks
        )

        for peptide_id in peptide_ids:
            for hit in peptide_id.getHits():
                support = hit.getMetaValue("consensus_support")
                consensus_support.add_value(support, stack=hit.getMetaValue("target_decoy"))
        consensus_support.to_dict()

        for i in consensus_support.dict["data"].keys():
            search_engine["consensus_support"][f"{raw_id_name} ({i})"] = consensus_support.dict["data"][i]

    search_engine["data_label"] = {
        "score_label": [spec_e_label, xcorr_label, hyper_label],
        "peps_label": peps_label,
        "consensus_label": consensus_label,
    }

    return search_engine, msgf_label, comet_label, sage_label


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