import pandas as pd
from datetime import datetime
import os
import copy
from collections import OrderedDict
import numpy as np

from sdrf_pipelines.openms.openms import OpenMS

from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common import ms_io
from pmultiqc.modules.common.logging import get_logger

log = get_logger("pmultiqc.modules.common.common_utils")


def read_openms_design(desfile):
    with open(desfile, "r") as f:

        data = f.readlines()
        s_row = False
        f_table = []
        s_table = []
        f_header = None
        s_header = None

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

        if f_header is None:
            raise ValueError("Cannot find 'Spectra_Filepath' header in file!")
        if s_header is None:
            raise ValueError("Cannot find 'MSstats_Condition' header in file!")

        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table["Run"] = f_table.apply(
            lambda x: file_prefix(x["Spectra_Filepath"]), axis=1
        )
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table


def condition_split(conditions):
    items = conditions.split(';')
    key_value_pairs = [item.split('=') for item in items if '=' in item]

    result_dict = {k.strip(): v.strip() for k, v in key_value_pairs}
    return result_dict


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


def parse_mzml(
        is_bruker: bool = False,
        read_ms_info: bool = False,
        ms_info_path: list[str] | None = None,
        ms_with_psm: list[str] | None = None,
        identified_spectrum: list[str] | None = None,
        enable_dia: bool = False,
        ms_paths: list[str] | None = None,
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
            ) = ms_io.get_ms_qc_info(mzml_df)

            log.info(
                "{}: Done aggregating ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )
        return None

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

    mzml_ms_df = None

    # Use the class-based MS info reader
    if read_ms_info:
        from pmultiqc.modules.common.ms.msinfo import MsInfoReader
        msinfo_reader = MsInfoReader(
            file_paths=ms_info_path,
            ms_with_psm=ms_with_psm,
            identified_spectrum=identified_spectrum,
            mzml_charge_plot=mzml_charge_plot,
            mzml_peak_distribution_plot=mzml_peak_distribution_plot,
            mzml_peaks_ms2_plot=mzml_peaks_ms2_plot,
            mzml_charge_plot_1=mzml_charge_plot_1,
            mzml_peak_distribution_plot_1=mzml_peak_distribution_plot_1,
            mzml_peaks_ms2_plot_1=mzml_peaks_ms2_plot_1,
            ms_without_psm=ms_without_psm,
            enable_dia=enable_dia,
        )
        msinfo_reader.parse()
        mzml_table = msinfo_reader.mzml_table
        heatmap_charge = msinfo_reader.heatmap_charge
        total_ms2_spectra = msinfo_reader.total_ms2_spectra
        ms1_tic = msinfo_reader.ms1_tic
        ms1_bpc = msinfo_reader.ms1_bpc
        ms1_peaks = msinfo_reader.ms1_peaks
        ms1_general_stats = msinfo_reader.ms1_general_stats
    else:
        from pmultiqc.modules.common.ms.mzml import MzMLReader

        mzml_reader = MzMLReader(
            file_paths=ms_paths,
            ms_with_psm=ms_with_psm,
            identified_spectrum=identified_spectrum,
            mzml_charge_plot=mzml_charge_plot,
            mzml_peak_distribution_plot=mzml_peak_distribution_plot,
            mzml_peaks_ms2_plot=mzml_peaks_ms2_plot,
            mzml_charge_plot_1=mzml_charge_plot_1,
            mzml_peak_distribution_plot_1=mzml_peak_distribution_plot_1,
            mzml_peaks_ms2_plot_1=mzml_peaks_ms2_plot_1,
            ms_without_psm=ms_without_psm,
            enable_dia=enable_dia,
            enable_mzid=enable_mzid
        )

        mzml_reader.parse()

        mzml_table = mzml_reader.mzml_table
        heatmap_charge = mzml_reader.heatmap_charge
        total_ms2_spectra = mzml_reader.total_ms2_spectra
        ms1_tic = mzml_reader.ms1_tic
        ms1_bpc = mzml_reader.ms1_bpc
        ms1_peaks = mzml_reader.ms1_peaks
        ms1_general_stats = mzml_reader.ms1_general_stats

        if enable_mzid:
            mzml_ms_df = mzml_reader.mzml_ms_df

    for i in sorted(set(ms_without_psm)):
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


def mod_group_percentage(group):
    if "Modifications" in group.columns:
        group.rename(columns={"Modifications": "modifications"}, inplace=True)

    counts = group["modifications"].str.split(",").explode().value_counts()
    percentage_df = (counts / len(group["modifications"]) * 100).reset_index()
    percentage_df.columns = ["modifications", "percentage"]

    # Modified (Total)
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"] = (
            100 - percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"]
    )
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "modifications"] = (
        "Modified (Total)"
    )

    return percentage_df


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

        counts, bin_edges = np.histogram(rt_list, bins=np.arange(rt_range_min, rt_range[1] + 5, 1))
        bin_mid = (bin_edges[:-1] + bin_edges[1:]) / 2
        rt_counts = pd.DataFrame({"retention_time": bin_mid, "counts": counts})

        return dict(zip(rt_counts["retention_time"], rt_counts["counts"]))

    rt_count_dict = {}
    for raw_file, group in evidence_data.groupby("raw file"):
        rt_count_dict[str(raw_file)] = hist_compute(group["retention time"], rt_range)

    return rt_count_dict


def evidence_calibrated_mass_error(
    evidence_data,
    recommpute=False,
    filter_outliers_ppm: bool = False
):
    # filter_outliers_ppm (if True): Remove rows with mass error [ppm] greater than 1000 (Default: False)

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

    # Remove rows with mass error [ppm] greater than 1000
    if filter_outliers_ppm:
        evd_df = evd_df[evd_df["mass error [ppm]"].abs() <= 1000].copy()

    if evd_df.empty:
        log.warning("No valid mass error [ppm] values found after filtering.")
        return None

    max_abs_ppm = evd_df["mass error [ppm]"].abs().max()

    if max_abs_ppm < 100:
        num_bins = 1000
    elif max_abs_ppm < 1000:
        num_bins = 10000
    elif max_abs_ppm < 5000:
        num_bins = 50000
    else:
        num_bins = 100000

    bin_series = pd.cut(evd_df["mass error [ppm]"], bins=num_bins)

    count_bin = bin_series.value_counts(sort=False)
    count_bin_data = {
        float(interval.mid): int(count) for interval, count in count_bin.items()
    }

    frequency_bin = bin_series.value_counts(sort=False, normalize=True)
    frequency_bin_data = {
        float(interval.mid): float(freq) for interval, freq in frequency_bin.items()
    }

    result_dict = {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }

    return result_dict

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


def parse_sdrf(
        sdrf_path,
        raw_config=None,
        condition_config=None
):
    OpenMS().openms_convert(
        sdrf_path,
        raw_config,  # config.kwargs["raw"],
        False,
        True,
        False,
        condition_config,  # config.kwargs["condition"],
    )