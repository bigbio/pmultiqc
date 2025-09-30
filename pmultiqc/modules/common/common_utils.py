import logging
import os

from sdrf_pipelines.openms.openms import OpenMS
from multiqc import config
from datetime import datetime

import numpy as np
import pandas as pd
import copy
from collections import OrderedDict
from typing import List

from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.ms_functions import get_ms_qc_info
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common import ms_io


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def get_exp_sdrf(find_log_files):

    exp_design = None
    enable_exp = False
    enable_sdrf = False

    for f in find_log_files("pmultiqc/exp_design", filecontents=False):
        exp_design = os.path.join(f["root"], f["fn"])
        enable_exp = True

    if not enable_exp:

        for f in find_log_files("pmultiqc/sdrf", filecontents=False):
            sdrf = os.path.join(f["root"], f["fn"])
            OpenMS().openms_convert(
                sdrf,
                config.kwargs["raw"],
                False,
                True,
                False,
                config.kwargs["condition"],
            )
            # experimental_design.tsv is the default output name
            # experimental_design.tsv will be in the folder where pmultiqc is executed.
            exp_design = "experimental_design.tsv"
            enable_sdrf = True

    return exp_design, enable_exp, enable_sdrf

def read_openms_design(desfile):
    with open(desfile, "r") as f:
        data = f.readlines()
        s_row = False
        f_table = []
        s_table = []
        for row in data:
            if row == "\n":
                continue
            if "MSstats_Condition" in row:
                s_row = True
                s_header = row.replace("\n", "").split("\t")
            elif s_row:
                s_table.append(row.replace("\n", "").split("\t"))
            elif "Spectra_Filepath" in row:
                f_header = row.replace("\n", "").split("\t")
            else:
                f_table.append(row.replace("\n", "").split("\t"))

        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table["Run"] = f_table.apply(
            lambda x: file_prefix(x["Spectra_Filepath"]), axis=1
        )
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table

def get_ms_path(find_log_files):

    ms_paths = []
    for mzml_current_file in find_log_files("pmultiqc/mzML", filecontents=False):
        ms_paths.append(os.path.join(mzml_current_file["root"], mzml_current_file["fn"]))

    ms_info_path = []
    for ms_info in find_log_files("pmultiqc/ms_info", filecontents=False):
        ms_info_path.append(os.path.join(ms_info["root"], ms_info["fn"]))
        ms_info_path.sort()
    
    read_ms_info = False
    if len(ms_info_path) > 0:
        read_ms_info = True
        ms_paths = [
            file_prefix(i).replace("_ms_info", ".mzML") for i in ms_info_path
        ]

    return ms_info_path, read_ms_info, ms_paths

def get_msstats_path(find_log_files):

    msstats_input_path = None
    msstats_input_valid = False

    for msstats_input in find_log_files("pmultiqc/msstats", filecontents=False):
        msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])
        msstats_input_valid = True

    return msstats_input_path, msstats_input_valid

def condition_split(conditions):
    items = conditions.split(';')
    key_value_pairs = [item.split('=') for item in items if '=' in item]

    result_dict = {k.strip(): v.strip() for k, v in key_value_pairs}
    return result_dict

def evidence_rt_count(evidence_data):
    if any(column not in evidence_data.columns for column in ["retention time", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    rt_range = [evidence_data["retention time"].min(), evidence_data["retention time"].max()]

    def hist_compute(rt_list, rt_range):

        if rt_range[0] < 5:
            rt_range_min = 0
        else:
            rt_range_min = rt_range[0] - 5
            
        counts, bin_edges = np.histogram(
            rt_list, bins=np.arange(rt_range_min, rt_range[1] + 5, 1)
        )
        bin_mid = (bin_edges[:-1] + bin_edges[1:]) / 2
        rt_counts = pd.DataFrame({"retention_time": bin_mid, "counts": counts})

        return dict(zip(rt_counts["retention_time"], rt_counts["counts"]))

    rt_count_dict = {}
    for raw_file, group in evidence_data.groupby("raw file"):
        rt_count_dict[str(raw_file)] = hist_compute(group["retention time"], rt_range)

    return rt_count_dict

# re-compute mass error
def recommpute_mass_error(evidence_df):

    required_cols = [
        "mass error [ppm]",
        "uncalibrated mass error [ppm]",
        "raw file",
        "mass",
        "charge",
        "m/z",
        "uncalibrated - calibrated m/z [ppm]",
    ]

    if not all(col in evidence_df.columns for col in required_cols):
        log.info("Evidence is missing one or more required columns in recommpute_mass_error.")
        return None

    df = evidence_df[required_cols].copy()

    decal_df = df.groupby("raw file", as_index=False).agg(
        decal=("uncalibrated mass error [ppm]", lambda x: np.median(np.abs(x)) > 1e3)
    )

    if decal_df["decal"].any():

        log.info("Detected at least one raw file with unusually high uncalibrated mass error.")

        df["theoretical_mz"] = df["mass"] / df["charge"] + 1.00726
        df["mass_error_ppm"] = (df["theoretical_mz"] - df["m/z"]) / df["theoretical_mz"] * 1e6
        df["uncalibrated_mass_error_ppm"] = (
            df["mass_error_ppm"] + df["uncalibrated - calibrated m/z [ppm]"]
        )

        idx_overwrite = df["raw file"].isin(decal_df.loc[decal_df["decal"], "raw file"])
        df.loc[idx_overwrite, "mass error [ppm]"] = df.loc[idx_overwrite, "mass_error_ppm"]
        df.loc[idx_overwrite, "uncalibrated mass error [ppm]"] = df.loc[
            idx_overwrite, "uncalibrated_mass_error_ppm"
        ]

    else:
        log.info("No raw files with unusually high uncalibrated mass error detected.")

    return df[["mass error [ppm]", "uncalibrated mass error [ppm]", "raw file"]]

def evidence_calibrated_mass_error(evidence_data, recommpute=False):

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    if recommpute:
        evd_df = recommpute_mass_error(evidence_data)
    else:
        evd_df = evidence_data.copy()

    if evd_df is None:
        if any(column not in evidence_data.columns for column in ["mass error [ppm]", "raw file"]):
            log.warning(
                "evidence_calibrated_mass_error: Required columns 'mass error [ppm]' or 'raw file' are missing in Evidence DataFrame."
            )
            return None
        else:
            log.warning(
                "Missing required columns. Skipping mass error recomputation and falling back to 'mass error [ppm]' and 'raw file' only."
            )
            evd_df = evidence_data[["mass error [ppm]", "raw file"]].copy()

    evd_df.dropna(subset=["mass error [ppm]"], inplace=True)

    count_bin = evd_df["mass error [ppm]"].value_counts(sort=False, bins=1000)
    count_bin_data = dict()
    for index in count_bin.index:
        count_bin_data[float(index.mid)] = int(count_bin[index])

    frequency_bin = evd_df["mass error [ppm]"].value_counts(sort=False, bins=1000, normalize=True)
    frequency_bin_data = dict()
    for index in frequency_bin.index:
        frequency_bin_data[float(index.mid)] = float(frequency_bin[index])

    result_dict = {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }

    return result_dict

def parse_mzml(
        is_bruker: bool = False,
        read_ms_info: bool = False,
        ms_info_path: List[str] = None,
        ms_with_psm: List[str] = None,
        identified_spectrum: List[str] = None,
        enable_dia: bool = False,
        ms_paths: str = None,
        enable_mzid: bool = False,
    ):

    ms1_tic = dict()
    ms1_bpc = dict()
    ms1_peaks = dict()
    ms1_general_stats = dict()

    if is_bruker and read_ms_info:
        for file in ms_info_path:
            log.info(
                "{}: Parsing ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )
            mzml_df = pd.read_csv(file, sep="\t")
            (
                ms1_tic[os.path.basename(file).replace("_ms_info.tsv", "")],
                ms1_bpc[os.path.basename(file).replace("_ms_info.tsv", "")],
                ms1_peaks[os.path.basename(file).replace("_ms_info.tsv", "")],
                ms1_general_stats[os.path.basename(file).replace("_ms_info.tsv", "")],
            ) = get_ms_qc_info(mzml_df)

            log.info(
                "{}: Done aggregating ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )
        return

    mzml_peak_distribution_plot = Histogram(
        "Peak Intensity",
        plot_category="range",
        breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000],
    )

    mzml_charge_plot = Histogram("Precursor Charge", plot_category="frequency")

    mzml_peaks_ms2_plot = Histogram(
        "#Peaks per MS/MS spectrum",
        plot_category="range",
        breaks=[i for i in range(0, 1001, 100)],
    )

    # New instances are used for dictionary construction.
    mzml_peak_distribution_plot_1 = copy.deepcopy(mzml_peak_distribution_plot)
    mzml_charge_plot_1 = copy.deepcopy(mzml_charge_plot)
    mzml_peaks_ms2_plot_1 = copy.deepcopy(mzml_peaks_ms2_plot)

    ms_without_psm = []

    mzml_table = {}
    heatmap_charge = {}

    # Use the refactored functions from ms_io.py
    if read_ms_info:
        result = ms_io.read_ms_info(
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
            enable_dia,
        )
        (
            mzml_table,
            heatmap_charge,
            total_ms2_spectra,
            ms1_tic,
            ms1_bpc,
            ms1_peaks,
            ms1_general_stats,
        ) = result
    else:
        result = ms_io.read_mzmls(
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
            enable_dia,
            enable_mzid,
        )

        if enable_mzid:
            (mzml_table, heatmap_charge, total_ms2_spectra, mzml_ms_df) = result
        else:
            mzml_table, heatmap_charge, total_ms2_spectra = result
            mzml_ms_df = None

    for i in ms_without_psm:
        log.warning("No PSM found in '{}'!".format(i))

    mzml_peaks_ms2_plot.to_dict()
    mzml_peak_distribution_plot.to_dict()

    ms_info = dict()
    ms_info["charge_distribution"] = dict()
    ms_info["peaks_per_ms2"] = dict()
    ms_info["peak_distribution"] = dict()
    
    # Construct compound dictionaries to apply to drawing functions.
    if enable_dia:
        mzml_charge_plot.to_dict()

        ms_info["charge_distribution"] = {
            "Whole Experiment": mzml_charge_plot.dict["data"]
        }
        ms_info["peaks_per_ms2"] = {
            "Whole Experiment": mzml_peaks_ms2_plot.dict["data"]
        }
        ms_info["peak_distribution"] = {
            "Whole Experiment": mzml_peak_distribution_plot.dict["data"]
        }
    else:
        mzml_peaks_ms2_plot_1.to_dict()
        mzml_peak_distribution_plot_1.to_dict()
        mzml_charge_plot.to_dict()
        mzml_charge_plot_1.to_dict()

        mzml_charge_plot.dict["cats"].update(mzml_charge_plot_1.dict["cats"])
        charge_cats_keys = [int(i) for i in mzml_charge_plot.dict["cats"]]
        charge_cats_keys.sort()
        mzml_charge_plot.dict["cats"] = OrderedDict(
            {str(i): mzml_charge_plot.dict["cats"][str(i)] for i in charge_cats_keys}
        )

        ms_info["charge_distribution"] = {
            "identified_spectra": mzml_charge_plot.dict["data"],
            "unidentified_spectra": mzml_charge_plot_1.dict["data"],
        }
        ms_info["peaks_per_ms2"] = {
            "identified_spectra": mzml_peaks_ms2_plot.dict["data"],
            "unidentified_spectra": mzml_peaks_ms2_plot_1.dict["data"],
        }
        ms_info["peak_distribution"] = {
            "identified_spectra": mzml_peak_distribution_plot.dict["data"],
            "unidentified_spectra": mzml_peak_distribution_plot_1.dict["data"],
        }

    median = np.median(list(heatmap_charge.values()))
    heatmap_charge_score = dict(
        zip(
            heatmap_charge.keys(),
            list(map(lambda v: 1 - np.abs(v - median), heatmap_charge.values())),
        )
    )

    return (
        mzml_table,
        mzml_peaks_ms2_plot,
        mzml_peak_distribution_plot,
        ms_info,
        total_ms2_spectra,
        mzml_ms_df,
        heatmap_charge_score,
        mzml_charge_plot,
        ms1_tic,
        ms1_bpc,
        ms1_peaks,
        ms1_general_stats
    )
