from multiqc.plots import heatmap, table
from pmultiqc.modules.core.section_groups import add_sub_section
from multiqc.types import SampleGroup, SampleName
from multiqc.plots.table_object import InputRow
from multiqc import config
from typing import Dict, List
import re

from pmultiqc.modules.common.common_utils import (
    read_openms_design,
    condition_split
)


def plot_html_check(plot_html):

    checked_html = remove_subtitle(plot_html)
    checked_html = apply_hoverinfo_config(checked_html)

    return checked_html


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
        is_maxquant
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
    }
    if is_maxquant:
        hm_html = heatmap.plot(data=heatmap_data, pconfig=pconfig)
        description_text = "This heatmap provides an overview of the performance of MaxQuant."
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

            for _, row in file_df_sample.iterrows():
                sample_data = {}
                for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                    sample_data["MSstats_Condition_" + str(k)] = ""
                sample_data["MSstats_BioReplicate"] = ""
                sample_data["Fraction_Group"] = row["Fraction_Group"]
                sample_data["Fraction"] = row["Fraction"]
                sample_data["Label"] = row["Label"]

                row_data.append(
                    InputRow(
                        sample=SampleName(row["Run"]),
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
        for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
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
            for _, row in file_df_sample.iterrows():
                row_data.append(
                    InputRow(
                        sample=SampleName(row["Run"]),
                        data={
                            "MSstats_Condition": "",
                            "MSstats_BioReplicate": "",
                            "Fraction_Group": row["Fraction_Group"],
                            "Fraction": row["Fraction"],
                            "Label": row["Label"],
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