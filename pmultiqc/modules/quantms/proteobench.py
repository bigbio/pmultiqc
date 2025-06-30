import pandas as pd
import os
import numpy as np

from multiqc.plots import bargraph, linegraph, box
import logging

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

    # log_Intensity_mean
    log_mean_html = draw_coefficient_variation(pb_df, "log_intensity_mean")

    # intensity_count_per_file
    num_inten_per_file_html = intensity_count_per_file(pb_df)

    # log_intensity_std
    log_std_html = draw_coefficient_variation(pb_df, "log_intensity_std")

    # CV
    cv_html = draw_coefficient_variation(pb_df, "cv")

    # log2_A_vs_B
    log_vs_html = draw_coefficient_variation(pb_df, "log2_A_vs_B")

    # epsilon
    epsilon_html = draw_coefficient_variation(pb_df, "epsilon")

    return {
        "log_mean_html": log_mean_html,
        "num_inten_per_file_html": num_inten_per_file_html,
        "log_std_html": log_std_html,
        "cv_html": cv_html,
        "log_vs_html": log_vs_html,
        "epsilon_html": epsilon_html,
    }


def draw_coefficient_variation(df, plot_type):

    enable_bar = False
    enable_line = False
    enable_box = False

    only_one_col = False 

    if plot_type == "log_intensity_mean":

        enable_bar = True
        enable_line = True

        col_A = "log_Intensity_mean_A"
        col_B = "log_Intensity_mean_B"
        bar_plot_id = "log_intensity_mean_na"
        bar_plot_title = "Missing Values"
        line_plot_id = "log_intensity_mean_linegraph"
        line_plot_title = "Distribution of Intensity"
        line_plot_xlab = "log2(Intensity)"

        show_lineplot_legend = True

    elif plot_type == "log_intensity_std":

        enable_box = True

        col_A = "log_Intensity_std_A"
        col_B = "log_Intensity_std_B"

        box_plot_id = "std_intensity_box"
        box_plot_title = "Standard Deviation of Intensity"
        box_plot_xlab = "Standard Deviation of log2(Intensity)"

    elif plot_type == "cv":

        enable_bar = True
        enable_line = True

        col_A = "CV_A"
        col_B = "CV_B"
        bar_plot_id = "cv_na"
        bar_plot_title = "Missing Values"
        line_plot_id = "cv_linegraph"
        line_plot_title = "Distribution of CV"
        line_plot_xlab = "Coefficient of Variation"

        show_lineplot_legend = True

    elif plot_type == "log2_A_vs_B":

        enable_line = True
        only_one_col = True

        col_A = "log2_A_vs_B"

        line_plot_id = "log_vs_linegraph"
        line_plot_title = "Log2 Fold Change (A vs B)"
        line_plot_xlab = "log2(log_Intensity_mean_A) - log2(log_Intensity_mean_B)"

        show_lineplot_legend = False

    elif plot_type == "epsilon":

        enable_line = True
        only_one_col = True

        col_A = "epsilon"

        line_plot_id = "epsilon_linegraph"
        line_plot_title = "Distribution of Epsilon"
        line_plot_xlab = "log2 FC difference"

        show_lineplot_legend = False

    # 1. Missing Values
    if enable_bar:

        bar_data = statistics_na_values(df, col_A, col_B)

        draw_bar_config = {
            "id": bar_plot_id,
            "cpswitch": True,
            "title": bar_plot_title,
            "tt_decimals": 0,
            "ylab": "Count",
        }

        bar_html = bargraph.plot(
            bar_data,
            pconfig=draw_bar_config,
        )
    else:
        bar_html = None
    
    # 2. Distribution
    if enable_line:

        if only_one_col:
            linegraph_data = statistics_one_line_values(df, col_A)
        
        else:
            linegraph_data = statistics_line_values(df, col_A, col_B, plot_type)

        draw_line_config = {
            "id": line_plot_id,
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": line_plot_title,
            "ymin": 0,
            "tt_decimals": 0,
            "ylab": "Count",
            "xlab": line_plot_xlab,
            "showlegend": show_lineplot_legend,
        }
        linegraph_html = linegraph.plot(linegraph_data, pconfig=draw_line_config)
    else:
        linegraph_html = None

    # 3. BoxPlot
    if enable_box:
        
        box_data = statistics_box_values(df, col_A, col_B)

        draw_box_config = {
            "id": box_plot_id,
            "cpswitch": True,
            "title": box_plot_title,
            "tt_decimals": 0,
            "xlab": box_plot_xlab,
        }

        box_html = box.plot(
            box_data,
            pconfig=draw_box_config,
        )
    else:
        box_html = None

    return {
        "bar_html": bar_html,
        "linegraph_html": linegraph_html,
        "box_html": box_html,
    }


def statistics_na_values(df, col1, col2):

    data_dict = {
        col: {
            "Non-NA": df[col].notna().sum(),
            "NA": df[col].isna().sum()
        }
        for col in  [col1, col2]
    }

    return data_dict


def statistics_line_values(df, col1, col2, dict_key):

    data_union = pd.concat([df[col1], df[col2]]).dropna()
    data_range_max = data_union.max()
    
    if dict_key == "cv":
        data_range_min = 0
    else:
        data_range_min = data_union.min() - (data_union.min() * 0.05)

    bin_step = 200
    data_bins = np.arange(data_range_min, data_range_max + (data_range_max * 0.05), (data_range_max / bin_step))
    bin_mid = (data_bins[: -1] + data_bins[1: ]) / 2

    data_dict = dict()
    for col in [col1, col2]:
        
        counts, _ = np.histogram(
            df[col].dropna(),
            bins=data_bins
        )

        cv_counts = pd.DataFrame({dict_key: bin_mid, "counts": counts})
        data_dict[col] = dict(zip(cv_counts[dict_key], cv_counts["counts"]))

    return data_dict


def statistics_one_line_values(df, col):

    data_list = df[col].dropna()
    data_range_max = data_list.max()

    data_range_min = data_list.min() - (data_list.min() * 0.05)

    bin_step = 1000

    data_bins = np.arange(data_range_min, data_range_max + (data_range_max * 0.05), (data_range_max / bin_step))
    bin_mid = (data_bins[: -1] + data_bins[1: ]) / 2

    data_dict = dict()
    counts, _ = np.histogram(
        df[col].dropna(),
        bins=data_bins
    )

    cv_counts = pd.DataFrame({"value": bin_mid, "counts": counts})
    data_dict[col] = dict(zip(cv_counts["value"], cv_counts["counts"]))

    return data_dict


def statistics_box_values(df, col1, col2):

    boxplot_data = dict()

    for col in [col1, col2]:
        boxplot_data[col] = df[col].dropna().tolist()

    return boxplot_data


def intensity_count_per_file(df):

    cols = [col for col in df.columns if "_Condition_" in col]

    non_na_dict = df[cols].notna().sum().to_dict()
    na_dict = df[cols].isna().sum().to_dict()

    plot_data = dict()
    for k, v in non_na_dict.items():

        plot_data[os.path.splitext(k)[0]] = {
            "Non-NA": v,
            "NA": na_dict[k],
        }

    draw_bar_config = {
        "id": "num_detected_features_per_file",
        "cpswitch": True,
        "title": "Number of Detected Features per File",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        plot_data,
        pconfig=draw_bar_config,
    )

    return bar_html
