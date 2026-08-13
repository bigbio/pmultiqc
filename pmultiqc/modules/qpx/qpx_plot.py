from __future__ import absolute_import

import pandas as pd

from multiqc.plots import table, bargraph

from pmultiqc.modules.common.plots.general import (
    plot_html_check
)
from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.common_utils import group_charge
from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.qpx.qpx_plot")

def draw_summary_table(
    sub_sections,
    total_ms2_spectra_identified: int = 0,
    total_peptide_count: int = 0,
    total_protein_identified: int = 0,
    total_protein_quantified: int = 0
):
    """Global summary table: identified spectra, peptides and proteins."""
    log.info("Summary table generation...")

    summary_table = {
        "Global Summary": {
            "#Identified MS2 Spectra": total_ms2_spectra_identified,
            "#Peptides Identified": total_peptide_count,
            "#Proteins Identified": total_protein_identified,
            "#Proteins Quantified": total_protein_quantified
        }
    }

    t_headers = {
        "#Identified MS2 Spectra": {
            "title": "#Identified MS2 Spectra",
            "description": "Total number of MS/MS spectra identified",
        },
        "#Peptides Identified": {
            "title": "#Peptides Identified",
        },
        "#Proteins Identified": {
            "title": "#Proteins Identified",
        },
        "#Proteins Quantified": {
            "title": "#Proteins Quantified",
        },
    }

    t_pconfig = {
        "id": "identification_summary_table",
        "title": "Summary Table",
        "save_file": False,
        "raw_data_fn": "summary_table_table",
        "sort_rows": False,
        "only_defined_headers": True,
        "col1_header": "Data Scope",
        "no_violin": True,
        "save_data_file": False,
    }

    table_html = table.plot(summary_table, t_headers, t_pconfig)

    description_str = "This table shows the pipeline summary statistics."
    helptext_str = """
        This table shows the pipeline summary statistics.
        """

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=1,
        description=description_str,
        helptext=helptext_str
    )

# Distribution of Precursor Charges
def draw_whole_exp_charge(sub_section, df):
    """Precursor charge distribution across the whole experiment."""

    df["charge"] = df["charge"].astype("str")

    bar_data = {
        "Whole Experiment": df["charge"].value_counts().sort_index().to_dict()
    }

    bar_data = {str(k): v for k, v in bar_data.items()}

    draw_config = {
        "id": "distribution_of_precursor_charges",
        "cpswitch": True,
        "cpswitch_c_active": False,
        "title": "Distribution of Precursor Charges",
        "tt_decimals": 0,
        "ylab": "Count",
        "save_data_file": False,
    }

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = plot_html_check(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=3,
        description="""
            This is a bar chart representing the distribution of the precursor ion charges for a given whole experiment.
        """,
        helptext="""
            [QPX: *feature.parquet] distribution of the precursor ion charges for a given whole experiment.
            charge: the charge of the feature.
            """,
    )


# Charge-state
def draw_qpx_ms2_charge(sub_section, df=None, sdrf_file_df=None):
    """
    Per-run (and per-sample) precursor charge state distribution.
    """

    if df is None or df.empty or "charge" not in df.columns:
        return 
    if sdrf_file_df is None:
        sdrf_file_df = pd.DataFrame()

    df["charge"] = df["charge"].astype("str")

    if "run" not in df.columns:
        return

    stat_data_by_run = group_charge(df, "run", "charge")

    draw_config = {
        "id": "charge_state",
        "cpswitch": True,
        "cpswitch_c_active": False,
        "title": "Charge-state",
        "tt_decimals": 0,
        "ylab": "Count",
        "save_data_file": False,
    }

    if sdrf_file_df.empty:
        bar_data = stat_data_by_run.to_dict(orient="index")
    else:
        df = df.merge(
            right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
            left_on="run",
            right_on="Run"
        )
        stat_data_by_sample = group_charge(df, "Sample", "charge")
        bar_data = [
            stat_data_by_run.to_dict(orient="index"),
            stat_data_by_sample.to_dict(orient="index")
        ]
        draw_config["data_labels"] = ["by Run", "by Sample"]

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = plot_html_check(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=4,
        description="The distribution of the charge-state of the precursor ion.",
        helptext="""
            [QPX: *feature.parquet] The distribution of the charge-state of feature (charge).
            """,
    )
