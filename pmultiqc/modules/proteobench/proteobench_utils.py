import os
import re

import numpy as np
import pandas as pd
from multiqc.plots import bargraph, linegraph, box, scatter

from ..common.plots.general import remove_subtitle

from pmultiqc.modules.common.logging import get_logger
log = get_logger("pmultiqc.modules.proteobench.proteobench_utils")


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
        runs_col = None
        log.warning(
            "'_Condition_' or 'abundance_' not found. Check result_performance.csv!"
        )

    # precursor ion charge
    charge_html = draw_precursor_ion_charge(pb_df)

    # log_Intensity_mean
    log_mean_html = draw_logmean_std_cv(df=pb_df, plot_type="log_intensity_mean")

    # log_Intensity
    log_intensity_html = draw_logmean_std_cv(
        df=pb_df, plot_type="log_intensity", runs_col=runs_col
    )

    # intensity_count_per_file
    num_inten_per_file_html = intensity_count_per_file(df=pb_df, runs_col=runs_col)

    # log_intensity_std
    log_std_html = draw_logmean_std_cv(df=pb_df, plot_type="log_intensity_std")

    # CV
    cv_html = draw_logmean_std_cv(df=pb_df, plot_type="cv")

    # log2_A_vs_B
    log_vs_html = draw_logmean_std_cv(df=pb_df, plot_type="log2_A_vs_B")

    # epsilon
    epsilon_html = draw_logmean_std_cv(df=pb_df, plot_type="epsilon")

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


def draw_logmean_std_cv(df, plot_type, runs_col=None):

    config = _config_for_plot_type(df, plot_type, runs_col)

    bar_html = _maybe_draw_bar(df, config)
    linegraph_html = _maybe_draw_line(df, config)
    box_html = _maybe_draw_box(df, config)

    return {
        "bar_html": bar_html,
        "linegraph_html": linegraph_html,
        "box_html": box_html,
    }


def _config_for_plot_type(df, plot_type, runs_col):
    config = {
        "enable_bar": False,
        "enable_line": False,
        "enable_box": False,
        "only_one_col": False,
        "cols": [],
        "bar": {"id": None, "title": None},
        "line": {"id": None, "title": None, "xlab": None},
        "box": {"id": None, "title": None, "xlab": None},
        "plot_type": plot_type,
    }

    if plot_type == "log_intensity_mean":
        config["enable_bar"] = True
        config["enable_line"] = True
        config["cols"] = ["log_Intensity_mean_A", "log_Intensity_mean_B"]
        if not all(col in df.columns for col in config["cols"]):
            log.warning(f"{' and '.join(config['cols'])} not found. Check result_performance.csv!")
        config["bar"].update({"id": "log_intensity_mean_na", "title": "Missing Values (for Conditions)"})
        config["line"].update({
            "id": "log_intensity_mean_linegraph",
            "title": "Distribution of Intensity (for Conditions)",
            "xlab": "log2(Intensity)",
        })

    elif plot_type == "log_intensity":
        config["enable_bar"] = True
        config["enable_line"] = True
        new_df, cols = calculate_log_intensity(df, runs_col)
        df = new_df if new_df is not None else df
        config["cols"] = cols or []
        config["bar"].update({"id": "log_intensity_na", "title": "Missing Values (for Runs)"})
        config["line"].update({
            "id": "log_intensity_linegraph",
            "title": "Distribution of Intensity (for Runs)",
            "xlab": "log2(Intensity)",
        })

    elif plot_type == "log_intensity_std":
        config["enable_box"] = True
        config["cols"] = ["log_Intensity_std_A", "log_Intensity_std_B"]
        if not all(col in df.columns for col in config["cols"]):
            log.warning(f"{' and '.join(config['cols'])} not found. Check result_performance.csv!")
        config["box"].update({
            "id": "std_intensity_box",
            "title": "Standard Deviation of Intensity",
            "xlab": "Standard Deviation of log2(Intensity)",
        })

    elif plot_type == "cv":
        config["enable_bar"] = True
        config["enable_line"] = True
        config["cols"] = ["CV_A", "CV_B"]
        if not all(col in df.columns for col in config["cols"]):
            log.warning(f"{' and '.join(config['cols'])} not found. Check result_performance.csv!")
        config["bar"].update({"id": "cv_na", "title": "Missing Values"})
        config["line"].update({
            "id": "cv_linegraph",
            "title": "Distribution of CV",
            "xlab": "Coefficient of Variation",
        })

    elif plot_type == "log2_A_vs_B":
        config["enable_line"] = True
        config["only_one_col"] = True
        config["cols"] = ["log2_A_vs_B"]
        if not all(col in df.columns for col in config["cols"]):
            log.warning(f"{config['cols'][0]} not found. Check result_performance.csv!")
        config["line"].update({
            "id": "log_vs_linegraph",
            "title": "Log2 Fold Change (A vs B)",
            "xlab": "log_Intensity_mean_A - log_Intensity_mean_B",
        })

    elif plot_type == "epsilon":
        config["enable_line"] = True
        config["only_one_col"] = True
        config["cols"] = ["epsilon"]
        if not all(col in df.columns for col in config["cols"]):
            log.warning(f"{config['cols'][0]} not found. Check result_performance.csv!")
        config["line"].update({
            "id": "epsilon_linegraph",
            "title": "Distribution of Epsilon",
            "xlab": "log2 FC difference",
        })

    config["df"] = df
    return config


def _maybe_draw_bar(df, config):
    if not config["enable_bar"]:
        return None
    cols = config["cols"]
    if not cols:
        return None
    bar_data = statistics_na_values(df, cols)
    draw_bar_config = {
        "id": config["bar"]["id"],
        "cpswitch": True,
        "title": config["bar"]["title"],
        "tt_decimals": 0,
        "ylab": "Count",
    }
    bar_html = bargraph.plot(data=bar_data, pconfig=draw_bar_config)
    return remove_subtitle(bar_html)


def _maybe_draw_line(df, config):
    if not config["enable_line"]:
        return None
    cols = config["cols"]
    if not cols:
        return None
    dict_key = "" if config["only_one_col"] else config["plot_type"]
    linegraph_data = statistics_line_values(df, cols, dict_key, config["only_one_col"])
    draw_line_config = {
        "id": config["line"]["id"],
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": config["line"]["title"],
        "ymin": 0,
        "tt_decimals": 0,
        "ylab": "Count",
        "xlab": config["line"]["xlab"],
        "showlegend": True,
    }
    linegraph_html = linegraph.plot(data=linegraph_data, pconfig=draw_line_config)
    return remove_subtitle(linegraph_html)


def _maybe_draw_box(df, config):
    if not config["enable_box"]:
        return None
    cols = config["cols"]
    if not cols:
        return None
    box_data = statistics_box_values(df, cols)
    draw_box_config = {
        "id": config["box"]["id"],
        "cpswitch": False,
        "title": config["box"]["title"],
        "tt_decimals": 5,
        "xlab": config["box"]["xlab"],
    }
    box_html = box.plot(list_of_data_by_sample=box_data, pconfig=draw_box_config)
    return remove_subtitle(box_html)


def statistics_na_values(df, cols):

    data_dict = {
        specie + ": " + col: {"Non-NA": group[col].notna().sum(), "NA": group[col].isna().sum()}
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
        data_range_min, data_range_max + (data_range_max * 0.05), (data_range_max / bin_step)
    )
    bin_mid = (data_bins[:-1] + data_bins[1:]) / 2

    data_dict = dict()
    for specie, group in df.groupby("species"):
        for col in cols:

            counts, _ = np.histogram(group[col].dropna(), bins=data_bins)

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


def intensity_count_per_file(df, runs_col=None):

    if runs_col is None:
        log.warning("runs_col is not set; skipping intensity_count_per_file.")
        return None

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

    df[["modified_sequence", "charge"]] = df["precursor ion"].apply(
        lambda x: pd.Series(parse_precursor(x))
    )

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

    df["log_Intensity_mean"] = df[["log_Intensity_mean_A", "log_Intensity_mean_B"]].mean(
        axis=1, skipna=False
    )
    df_sub = df[["species", "log_Intensity_mean", "log2_A_vs_B"]].copy().dropna()
    df_sub.rename(columns={"log_Intensity_mean": "y", "log2_A_vs_B": "x"}, inplace=True)
    df_sub["color"] = df_sub["species"].map({"ECOLI": "blue", "HUMAN": "green", "YEAST": "red"})

    plot_data = {
        specie: group[["x", "y", "color"]].to_dict(orient="records")
        for specie, group in df_sub.groupby("species")
    }
    species_order = ["HUMAN", "YEAST", "ECOLI"]
    plot_data = {key: plot_data[key] for key in species_order if key in plot_data}

    draw_config = {
        "id": "log2fc_vs_logintensitymean",
        "title": "log2FC vs logIntensityMean",
        "xlab": "log2FC(A:B)",
        "ylab": "logIntensityMean",
        "showlegend": True,
    }

    scatter_html = scatter.plot(data=plot_data, pconfig=draw_config)

    scatter_html = remove_subtitle(scatter_html)

    return scatter_html


def calculate_log_intensity(df, runs_col):

    if runs_col is None:
        return None, None

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


def parse_precursor(precursor):
    if "|" in precursor:
        seq, z = precursor.split("|", 1)
        charge = re.search(r"Z=(\d+)", z)
        return seq, charge.group(1) if charge else None
    elif "/" in precursor:
        seq, charge = precursor.split("/", 1)
        return seq, charge
    else:
        return precursor, None