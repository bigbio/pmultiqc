import logging
import os
import pandas as pd
import numpy as np

from multiqc.plots import (
    bargraph,
    linegraph,
    box,
    scatter
)

from ..common.common_plots import remove_subtitle


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def read_file_by_extension(file_path):

    _, ext = os.path.splitext(file_path)

    if ext == ".csv":
        df = pd.read_csv(file_path)
    elif ext == ".tsv":
        df = pd.read_csv(file_path, sep="\t")
    elif ext == ".txt":
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

    return df

def get_pb_data(file_path):

    log.info("Start parsing the ProteoBench result file...")

    pb_df = read_file_by_extension(file_path)

    if any("abundance_" in col for col in pb_df.columns):
        runs_col = "abundance_"
    elif any("_Condition_" in col for col in pb_df.columns):
        runs_col = "_Condition_"
    else:
        log.warning("'_Condition_' or 'abundance_' not found. Check result_performance.csv!")

    # precursor ion charge
    charge_html = draw_precursor_ion_charge(pb_df)

    # log_Intensity_mean
    log_mean_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="log_intensity_mean"
    )

    # log_Intensity
    log_intensity_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="log_intensity",
        runs_col=runs_col
    )

    # intensity_count_per_file
    num_inten_per_file_html = intensity_count_per_file(
        df=pb_df,
        runs_col=runs_col
    )

    # log_intensity_std
    log_std_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="log_intensity_std"
    )

    # CV
    cv_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="cv"
    )

    # log2_A_vs_B
    log_vs_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="log2_A_vs_B"
    )

    # epsilon
    epsilon_html = draw_logmean_std_cv(
        df=pb_df,
        plot_type="epsilon"
    )

    # log2FC vs logIntensityMean
    logfc_logmean_html = draw_logintensitymean_vs_logfc(pb_df)

    return {
        "charge_html": charge_html,
        "log_mean_html": log_mean_html,
        "log_intensity_html": log_intensity_html,
        "num_inten_per_file_html": num_inten_per_file_html,
        "log_std_html": log_std_html,
        "cv_html": cv_html,
        "log_vs_html": log_vs_html,
        "epsilon_html": epsilon_html,
        "logfc_logmean_html": logfc_logmean_html,
    }

def draw_logmean_std_cv(
        df,
        plot_type,
        runs_col=None
    ):

    enable_bar = False
    enable_line = False
    enable_box = False

    only_one_col = False 

    if plot_type == "log_intensity_mean":

        enable_bar = True
        enable_line = True

        cols = ["log_Intensity_mean_A", "log_Intensity_mean_B"]

        if not all(col in df.columns for col in cols):
            log.warning(f"{' and '.join(cols)} not found. Check result_performance.csv!")

        bar_plot_id = "log_intensity_mean_na"
        bar_plot_title = "Missing Values (for Conditions)"
        line_plot_id = "log_intensity_mean_linegraph"
        line_plot_title = "Distribution of Intensity (for Conditions)"
        line_plot_xlab = "log2(Intensity)"

    elif plot_type == "log_intensity":

        enable_bar = True
        enable_line = True

        df, cols = calculate_log_intensity(df, runs_col)

        bar_plot_id = "log_intensity_na"
        bar_plot_title = "Missing Values (for Runs)"
        line_plot_id = "log_intensity_linegraph"
        line_plot_title = "Distribution of Intensity (for Runs)"
        line_plot_xlab = "log2(Intensity)"

    elif plot_type == "log_intensity_std":

        enable_box = True

        cols = ["log_Intensity_std_A", "log_Intensity_std_B"]

        if not all(col in df.columns for col in cols):
            log.warning(f"{' and '.join(cols)} not found. Check result_performance.csv!")

        box_plot_id = "std_intensity_box"
        box_plot_title = "Standard Deviation of Intensity"
        box_plot_xlab = "Standard Deviation of log2(Intensity)"

    elif plot_type == "cv":

        enable_bar = True
        enable_line = True

        cols = ["CV_A", "CV_B"]

        if not all(col in df.columns for col in cols):
            log.warning(f"{' and '.join(cols)} not found. Check result_performance.csv!")

        bar_plot_id = "cv_na"
        bar_plot_title = "Missing Values"
        line_plot_id = "cv_linegraph"
        line_plot_title = "Distribution of CV"
        line_plot_xlab = "Coefficient of Variation"

    elif plot_type == "log2_A_vs_B":

        enable_line = True
        only_one_col = True

        cols = ["log2_A_vs_B"]

        if not all(col in df.columns for col in cols):
            log.warning(f"{cols[0]} not found. Check result_performance.csv!")

        line_plot_id = "log_vs_linegraph"
        line_plot_title = "Log2 Fold Change (A vs B)"
        line_plot_xlab = "log_Intensity_mean_A - log_Intensity_mean_B"

    elif plot_type == "epsilon":

        enable_line = True
        only_one_col = True

        cols = ["epsilon"]

        if not all(col in df.columns for col in cols):
            log.warning(f"{cols[0]} not found. Check result_performance.csv!")

        line_plot_id = "epsilon_linegraph"
        line_plot_title = "Distribution of Epsilon"
        line_plot_xlab = "log2 FC difference"

    # 1. Missing Values
    if enable_bar:

        bar_data = statistics_na_values(df, cols)

        draw_bar_config = {
            "id": bar_plot_id,
            "cpswitch": True,
            "title": bar_plot_title,
            "tt_decimals": 0,
            "ylab": "Count",
        }

        bar_html = bargraph.plot(
            data=bar_data,
            pconfig=draw_bar_config,
        )
        
        bar_html = remove_subtitle(bar_html)
    else:
        bar_html = None
    
    # 2. Distribution
    if enable_line:

        if only_one_col:
            linegraph_data = statistics_line_values(df, cols, "", only_one_col)

        else:
            linegraph_data = statistics_line_values(df, cols, plot_type, only_one_col)

        draw_line_config = {
            "id": line_plot_id,
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": line_plot_title,
            "ymin": 0,
            "tt_decimals": 0,
            "ylab": "Count",
            "xlab": line_plot_xlab,
            "showlegend": True,
        }

        linegraph_html = linegraph.plot(
            data=linegraph_data,
            pconfig=draw_line_config
        )

        linegraph_html = remove_subtitle(linegraph_html)

    else:
        linegraph_html = None

    # 3. BoxPlot
    if enable_box:
        
        box_data = statistics_box_values(df, cols)

        draw_box_config = {
            "id": box_plot_id,
            "cpswitch": False,
            "title": box_plot_title,
            "tt_decimals": 5,
            "xlab": box_plot_xlab,
        }

        box_html = box.plot(
            list_of_data_by_sample=box_data,
            pconfig=draw_box_config,
        )

        box_html = remove_subtitle(box_html)

    else:
        box_html = None

    return {
        "bar_html": bar_html,
        "linegraph_html": linegraph_html,
        "box_html": box_html,
    }

def statistics_na_values(df, cols):

    data_dict = {
        specie + ": " + col: {
            "Non-NA": group[col].notna().sum(),
            "NA": group[col].isna().sum()
        }
        for specie, group in df.groupby("species")
        for col in cols
    }

    return data_dict

def statistics_line_values(df, cols, dict_key, only_one_col):

    data_union = df[cols].stack().reset_index(drop=True)
    data_range_max = data_union.max()
    
    if dict_key == "cv":
        data_range_min = 0
    else:
        data_range_min = data_union.min() - (data_union.min() * 0.05)

    bin_step = 1000
    data_bins = np.arange(
        data_range_min,
        data_range_max + (data_range_max * 0.05),
        (data_range_max / bin_step)
    )
    bin_mid = (data_bins[: -1] + data_bins[1: ]) / 2

    data_dict = dict()
    for specie, group in df.groupby("species"):
        for col in cols:
            
            counts, _ = np.histogram(
                group[col].dropna(),
                bins=data_bins
            )

            cv_counts = pd.DataFrame({"value": bin_mid, "counts": counts})

            if only_one_col:
                data_dict[specie] = dict(zip(cv_counts["value"], cv_counts["counts"]))
            else:
                data_dict[specie + ": " + col] = dict(zip(cv_counts["value"], cv_counts["counts"]))

    return data_dict

def statistics_box_values(df, cols):

    boxplot_data = {
        specie + ": " + col: group[col].dropna().tolist()
        for specie, group in df.groupby("species")
        for col in cols
    }

    return boxplot_data

def intensity_count_per_file(df, runs_col):

    cols = [col for col in df.columns if runs_col in col]

    non_na_dict = df[cols].notna().sum().to_dict()
    na_dict = df[cols].isna().sum().to_dict()

    plot_data = dict()
    for k, v in non_na_dict.items():

        plot_data[os.path.splitext(k)[0]] = {
            "Non-NA": v,
            "NA": na_dict[k],
        }

    draw_bar_config = {
        "id": "num_detected_features_per_run",
        "cpswitch": True,
        "title": "Number of Detected Features per Run",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=plot_data,
        pconfig=draw_bar_config,
    )

    bar_html = remove_subtitle(bar_html)

    return bar_html

def draw_precursor_ion_charge(df):

    df[["modified_sequence", "Z=charge"]] = df["precursor ion"].str.split("|", expand=True)
    df["charge"] = df["Z=charge"].str.extract(r"Z=(\d+)")

    charge_data = (
        df.groupby("species")["charge"]
        .value_counts()
        .sort_index()
        .unstack(fill_value=0)
        .to_dict(orient="index")
    )

    draw_config = {
        "id": "charge",
        "cpswitch": True,
        "title": "Distribution of Precursor Charges",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=charge_data,
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    return bar_html

# log2FC vs logIntensityMean
def draw_logintensitymean_vs_logfc(df):

    df["log_Intensity_mean"] = df[["log_Intensity_mean_A", "log_Intensity_mean_B"]].mean(axis=1, skipna=False)
    df_sub = df[["species", "log_Intensity_mean", "log2_A_vs_B"]].copy().dropna()
    df_sub.rename(
        columns={
            "log_Intensity_mean": "y",
            "log2_A_vs_B": "x"
        },
        inplace=True
    )
    df_sub["color"] = df_sub["species"].map({"ECOLI": "blue", "HUMAN": "green", "YEAST": "red"})

    plot_data = {
        specie: group[["x", "y", "color"]].to_dict(orient="records")
        for specie, group in df_sub.groupby("species")
    }
    species_order = ["HUMAN", "YEAST", "ECOLI"]
    plot_data = {
        key: plot_data[key] for key in species_order if key in plot_data
    }

    draw_config = {
        "id": "log2fc_vs_logintensitymean",
        "title": "log2FC vs logIntensityMean",
        "xlab": "log2FC(A:B)",
        "ylab": "logIntensityMean",
        "showlegend": True
    }

    scatter_html = scatter.plot(
        data=plot_data,
        pconfig=draw_config
    )

    scatter_html == remove_subtitle(scatter_html)

    return scatter_html

def calculate_log_intensity(df, runs_col):

    cols = [col for col in df.columns if runs_col in col]
    log_cols = [f"log2_{os.path.splitext(col)[0]}" for col in cols]

    df_copy = df.copy()

    for col, new_col in zip(cols, log_cols):
        df_copy[new_col] = np.log2(df_copy[col])

    need_cols = ["species"] + log_cols
    plot_df = df_copy[need_cols].copy()
    plot_df.columns = [col.removeprefix("log2_") for col in plot_df.columns]
    cols = [col for col in plot_df.columns if runs_col in col]

    return plot_df, cols
