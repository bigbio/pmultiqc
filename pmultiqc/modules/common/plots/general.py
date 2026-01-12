from multiqc.plots import heatmap, table, bargraph
from pmultiqc.modules.core.section_groups import add_sub_section
from multiqc.types import SampleGroup, SampleName
from multiqc.plots.table_object import InputRow
from multiqc import config
from typing import Dict, List
import re
import pandas as pd
import numpy as np

from pmultiqc.modules.common.common_utils import (
    read_openms_design,
    condition_split
)


FLAT_THRESHOLD = 100000

def plot_html_check(plot_html):

    checked_html = remove_subtitle(plot_html)
    checked_html = apply_hoverinfo_config(checked_html)

    return checked_html


def plot_data_check(
    plot_data,
    plot_html,
    log_text,
    function_name
):

    from collections.abc import Mapping
    from pmultiqc.modules.common.logging import get_logger

    log = get_logger(log_text)

    def count_elements(plot_data):

        count = 0
        if isinstance(plot_data, list):
            for item in plot_data:
                count += count_elements(item)

        elif isinstance(plot_data, Mapping):
            for v in plot_data.values():
                count += count_elements(v)

        else:
            count += 1

        return count

    data_counts = count_elements(plot_data)

    log.info(f"{function_name}: Plot data count: {data_counts}")

    if data_counts >= FLAT_THRESHOLD:
        plot_html.flat = True
        log.warning(
            f"{function_name}: Number of plotting data points: {data_counts}, exceeds threshold {FLAT_THRESHOLD}, switching to flat plot"
        )

    return plot_html


def remove_subtitle(plot_html):
    for dataset in plot_html.datasets:
        if "subtitle" in dataset.dconfig:
            dataset.dconfig["subtitle"] = ""
        if dataset.layout and "title" in dataset.layout:
            title_text = dataset.layout["title"].get("text", "")
            if "<br><sup>" in title_text:
                dataset.layout["title"]["text"] = title_text.split("<br><sup>")[0]
    return plot_html


def apply_hoverinfo_config(plot_html):

    if config.kwargs.get("disable_hoverinfo", False):
        for dataset in plot_html.datasets:

            dataset.trace_params["hoverinfo"] = "skip"

            if "hovertemplate" in dataset.trace_params:
                dataset.trace_params["hovertemplate"] = ""

    return plot_html


def draw_heatmap(
    sub_sections,
    hm_colors,
    heatmap_data,
    heatmap_xnames,
    heatmap_ynames,
    report_type
):
    pconfig = {
        "id": "heatmap",
        "title": "HeatMap",
        "min": 0,
        "max": 1,
        "xlab": "Metrics",
        "ylab": "RawName",
        "zlab": "Score",
        "tt_decimals": 4,
        "square": False,
        "colstops": hm_colors,
        "cluster_rows": False,
        "cluster_cols": False,
        "save_data_file": False,
    }
    if report_type == "maxquant":
        hm_html = heatmap.plot(data=heatmap_data, pconfig=pconfig)
        description_text = "This heatmap provides an overview of the performance of MaxQuant."
    elif report_type == "fragpipe":
        hm_html = heatmap.plot(data=heatmap_data, pconfig=pconfig)
        description_text = "This heatmap provides an overview of the performance of FragPipe."
    else:
        hm_html = heatmap.plot(heatmap_data, heatmap_xnames, heatmap_ynames, pconfig)
        description_text = "This heatmap provides an overview of the performance of quantms."

    hm_html = plot_html_check(hm_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=hm_html,
        order=2,
        description=description_text,
        helptext="""
            This plot shows the pipeline performance overview.
        """,
    )


def draw_exp_design(sub_sections, exp_design):
    # Currently this only supports the OpenMS two-table format (default in quantms pipeline)
    sample_df, file_df = read_openms_design(exp_design)

    exp_design_runs = file_df["Run"].unique().tolist()

    is_bruker = False
    if not file_df.empty:
        first_path = str(file_df["Spectra_Filepath"].iloc[0])
        is_bruker = first_path.endswith((".d", ".d.tar"))

    pattern = r'^(\w+=[^=;]+)(;\w+=[^=;]+)*$'
    is_multi_conditions = all(sample_df["MSstats_Condition"].apply(lambda x: bool(re.match(pattern, str(x)))))

    rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

    if is_multi_conditions:
        for sample in sorted(
            sample_df["Sample"].tolist(),
            key=lambda x: (str(x).isdigit(), int(x) if str(x).isdigit() else str(x).lower()),
        ):
            file_df_sample = file_df[file_df["Sample"] == sample].copy()
            sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
            row_data: List[InputRow] = []

            sample_data = {}
            for k, v in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                sample_data["MSstats_Condition_" + str(k)] = v
            sample_data["MSstats_BioReplicate"] = sample_df_slice["MSstats_BioReplicate"].iloc[0]
            sample_data["Fraction_Group"] = ""
            sample_data["Fraction"] = ""
            sample_data["Label"] = ""

            row_data.append(
                InputRow(
                    sample=SampleName(f"Sample {str(sample)}"),
                    data=sample_data,
                )
            )

            for row in file_df_sample.itertuples():
                sample_data = {}
                for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                    sample_data["MSstats_Condition_" + str(k)] = ""
                sample_data["MSstats_BioReplicate"] = ""
                sample_data["Fraction_Group"] = row.Fraction_Group
                sample_data["Fraction"] = row.Fraction
                sample_data["Label"] = row.Label

                row_data.append(
                    InputRow(
                        sample=SampleName(row.Run),
                        data=sample_data,
                    )
                )
            group_name: SampleGroup = SampleGroup(sample)
            rows_by_group[group_name] = row_data
        headers = {"Sample": {
            "title": "Sample [Spectra File]",
            "description": "",
            "scale": False,
        }}
        # Use first row of sample_df for condition keys (safer than relying on loop variable)
        first_condition = sample_df["MSstats_Condition"].iloc[0] if not sample_df.empty else ""
        for k, _ in condition_split(first_condition).items():
            headers["MSstats_Condition_" + str(k)] = {
                "title": "MSstats Condition: " + str(k),
                "description": "",
                "scale": False,
            }
        headers["MSstats_BioReplicate"] = {
            "title": "MSstats BioReplicate",
            "description": "",
            "scale": False,
        }
        headers["Fraction_Group"] = {
            "title": "Fraction Group",
            "description": "",
            "scale": False,
        }
        headers["Fraction"] = {
            "title": "Fraction",
            "description": "Fraction Identifier",
            "scale": False,
        }
        headers["Label"] = {
            "title": "Label",
            "description": "",
            "scale": False,
        }
    else:
        for sample in sorted(
            sample_df["Sample"].tolist(),
            key=lambda x: (str(x).isdigit(), int(x) if str(x).isdigit() else str(x).lower()),
        ):
            file_df_sample = file_df[file_df["Sample"] == sample].copy()
            sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
            row_data: List[InputRow] = []
            row_data.append(
                InputRow(
                    sample=SampleName(f"Sample {str(sample)}"),
                    data={
                        "MSstats_Condition": sample_df_slice["MSstats_Condition"].iloc[0],
                        "MSstats_BioReplicate": sample_df_slice["MSstats_BioReplicate"].iloc[0],
                        "Fraction_Group": "",
                        "Fraction": "",
                        "Label": "",
                    },
                )
            )
            for row in file_df_sample.itertuples():
                row_data.append(
                    InputRow(
                        sample=SampleName(row.Run),
                        data={
                            "MSstats_Condition": "",
                            "MSstats_BioReplicate": "",
                            "Fraction_Group": row.Fraction_Group,
                            "Fraction": row.Fraction,
                            "Label": row.Label,
                        },
                    )
                )
            group_name: SampleGroup = SampleGroup(sample)
            rows_by_group[group_name] = row_data

        headers = {
            "Sample": {
                "title": "Sample [Spectra File]",
                "description": "",
                "scale": False,
            },
            "MSstats_Condition": {
                "title": "MSstats Condition",
                "description": "MSstats Condition",
                "scale": False,
            },
            "MSstats_BioReplicate": {
                "title": "MSstats BioReplicate",
                "description": "MSstats BioReplicate",
                "scale": False,
            },
            "Fraction_Group": {
                "title": "Fraction Group",
                "description": "Fraction Group",
                "scale": False,
            },
            "Fraction": {
                "title": "Fraction",
                "description": "Fraction Identifier",
                "scale": False,
            },
            "Label": {
                "title": "Label",
                "description": "Label",
                "scale": False,
            },
        }

    pconfig = {
        "id": "experimental_design",
        "title": "Experimental Design",
        "save_file": False,
        "raw_data_fn": "multiqc_Experimental_Design_table",
        "no_violin": True,
        "save_data_file": False,
    }

    table_html = table.plot(rows_by_group, headers, pconfig)
    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=1,
        description="""
            This table shows the design of the experiment. I.e., which files and channels correspond to which sample/condition/fraction.
            """,
        helptext="""
            You can see details about it in
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            """
    )

    return sample_df, file_df, exp_design_runs, is_bruker, is_multi_conditions


def stat_pep_intensity(intensities: pd.Series):

    stat_result = np.log2(intensities[intensities >= 1])

    return stat_result.tolist()


def search_engine_score_bins(
    bins_start: int,
    bins_end: int,
    bins_step: int,
    df: pd.DataFrame,
    groupby_col: str,
    score_col: str
):

    bins = list(range(bins_start, bins_end + 1, bins_step)) + [float("inf")]
    labels = [
        f"{i} ~ {i + bins_step}" for i in range(bins_start, bins_end, bins_step)
    ] + [
        f"{bins_end} ~ inf"
    ]

    plot_data = []
    data_labels = []
    for name, group in df.groupby(groupby_col):

        group = group.copy()
        group["score_bin"] = pd.cut(group[score_col], bins=bins, labels=labels, right=False)
        score_dist = group["score_bin"].value_counts().sort_index().reset_index()

        plot_data.append(
            {k: {"count": v} for k, v in zip(score_dist["score_bin"], score_dist["count"], strict=True)}
        )
        data_labels.append(name)

    result = {
        "plot_data": plot_data,
        "data_labels": data_labels,
    }

    return result


# Search Engine Scores
def draw_search_engine_scores(sub_section, plot_data, plot_type):

    if plot_type == "maxquant":
        plot_config = {
            "id": "summary_of_andromeda_scores",
            "title": "Summary of Andromeda Scores",
            "helptext": """
                This statistic is extracted from msms.txt. Andromeda score for the best associated MS/MS spectrum.
                """,
        }
    elif plot_type == "fragpipe":
        plot_config = {
            "id": "summary_of_hyperscore",
            "title": "Summary of Hyperscore",
            "helptext": """
                This statistic is extracted from psm.tsv.<br>
                [Hyperscore] Similarity score between observed and theoretical spectra, higher values indicate greater similarity.
                """,
        }
    else:
        raise ValueError("[draw_search_engine_scores] Please check the plot type.")

    pconfig = {
        "id": plot_config["id"],
        "cpswitch": False,
        "title": plot_config["title"],
        "ylab": "Counts",
        "tt_suffix": "",
        "tt_decimals": 0,
        "data_labels": plot_data["data_labels"],
        "save_data_file": False,
    }

    bar_html = bargraph.plot(data=plot_data["plot_data"], pconfig=pconfig)

    bar_html = plot_html_check(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description="",
        helptext=plot_config["helptext"],
    )
