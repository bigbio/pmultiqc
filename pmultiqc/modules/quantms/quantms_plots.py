import logging
from multiqc.plots import bargraph, box

from ..common.common_plots import (
    draw_ids_rt_count,
    remove_subtitle
)
from ..core.section_groups import add_sub_section
from ..maxquant.maxquant_utils import evidence_rt_count
from ..common.common_plots import draw_oversampling

import re
import pandas as pd
import numpy as np


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Normalised"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Normalised"])
    
    draw_dia_intensity_dis(sub_section, df_sub)
    draw_dia_intensity_std(sub_section, df_sub)

def draw_dia_ms2s(sub_section, df):

    draw_dia_ms2_charge(sub_section, df)

    if "MS2.Scan" in df.columns:
        msms_count_data = calculate_msms_count(df)
        draw_oversampling(
            sub_section,
            msms_count_data["plot_data"],
            msms_count_data["cats"],
            False
        )

def draw_dia_time_mass(sub_section, df):

    df = df[["Run", "RT"]].copy()
    df.rename(
        columns={
            "Run": "raw file",
            "RT": "retention time"
        },
        inplace=True
    )

    ids_over_rt = evidence_rt_count(df)

    draw_ids_rt_count(
        sub_section,
        ids_over_rt,
        "dia"
    )

# Intensity Distribution
def draw_dia_intensity_dis(sub_section, df):

    box_data = {
        run: group["log_intensity"].dropna().tolist()
        for run, group in df.groupby("Run")
    }

    draw_config = {
        "id": "intensity_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Intensity Distribution",
        "tt_decimals": 5,
        "xlab": "log2(Precursor.Normalised)",
    }

    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_config
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=3,
        description="log2(Precursor.Normalised) for each Run.",
        helptext="""
            [DIA-NN: report.tsv] log2(Precursor.Normalised) for each Run.
            """
    )

# Charge-state of Per File
def draw_dia_ms2_charge(sub_section, df):

    df = df.copy()
    df["Precursor.Charge"] = df["Precursor.Charge"].astype("str")

    bar_data = (
        df.groupby("Run")["Precursor.Charge"]
        .value_counts()
        .sort_index()
        .unstack(fill_value=0)
        .to_dict(orient="index")
    )

    draw_config = {
        "id": "charge_state_of_per_file",
        "cpswitch": True,
        "title": "Charge-state of Per File",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description="The distribution of the charge-state of the precursor ion.",
        helptext="""
            [DIA-NN: report.tsv] The distribution of the charge-state of the precursor ion (Precursor.Charge).
            """
    )

# DIA: Standard Deviation of Intensity
def draw_dia_intensity_std(sub_section, df):

    box_data = calculate_dia_intensity_std(df)

    draw_box_config = {
        "id": "dia_std_intensity_box",
        "title": "Standard Deviation of Intensity",
        "cpswitch": False,
        "tt_decimals": 5,
        "xlab": "Standard Deviation of log2(Precursor.Normalised)",
    }
    
    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_box_config,
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=6,
        description="Standard deviation of intensity under different experimental conditions.",
        helptext="""
            [DIA-NN: report.tsv] First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Precursor.Normalised).
            """
    )


def extract_condition_and_replicate(run_name):

    match = re.search(r"(.*?)([A-Za-z]*)(\d+)$", run_name)

    if match:
        condition_base = match.group(1) + match.group(2)
        
        if condition_base.endswith("_"):
            condition_base = condition_base[: -1]

        replicate = int(match.group(3))

        return condition_base, replicate
    else:
        log.warning("Failed to identify condition groups in DIA report.tsv!")

def calculate_dia_intensity_std(df):

    df_sub = df.copy()
    df_sub[["run_condition", "run_replicate"]] = df_sub["Run"].apply(
        lambda x: pd.Series(extract_condition_and_replicate(x))
    )

    grouped_std = (
        df_sub
        .groupby(["run_condition", "Modified.Sequence"])["log_intensity"]
        .std()
        .reset_index(name="log_intensity_std")
    )

    plot_data = {
        condition: group["log_intensity_std"].dropna().tolist()
        for condition, group in grouped_std.groupby("run_condition")
    }

    return plot_data

def calculate_msms_count(df):
    count_df = df.groupby(
        [
            "Run",
            "Stripped.Sequence",
            "Precursor.Charge"
        ]
    )["MS2.Scan"].nunique().reset_index(name="msms_count")

    run_counts = count_df.groupby("Run")["msms_count"].value_counts().reset_index(name="msms_count_run")

    run_counts["msms_count_str"] = run_counts["msms_count"].apply(lambda x: ">=3" if x >= 3 else str(x))
    merged_df = run_counts.groupby(["Run", "msms_count_str"])["msms_count_run"].sum().reset_index()

    result_dict = dict()
    for run_name, group in merged_df.groupby("Run"):
        result_dict[run_name] = dict(zip(group["msms_count_str"], group["msms_count_run"]))

    return {
         "plot_data": result_dict,
         "cats": list(merged_df["msms_count_str"].unique())
    }
