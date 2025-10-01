"""
Contaminant-related plotting functions for PMultiQC

This module contains functions for plotting contaminant analysis results.
"""

from multiqc.plots import bargraph
from pmultiqc.modules.core.section_groups import add_sub_section


def draw_potential_contaminants(sub_section, contaminant_percent, is_maxquant):
    """
    Draw potential contaminants plot per group.

    Args:
        sub_section: MultiQC sub-section for adding plots
        contaminant_percent: Dictionary with contaminant percentage data
        is_maxquant: Boolean indicating if data is from MaxQuant
    """
    draw_config = {
        "id": "potential_contaminants_per_group",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Potential Contaminants Per File",
        "tt_decimals": 2,
        "ylab": "Percent [%]",
    }

    bar_html = bargraph.plot(data=contaminant_percent, pconfig=draw_config)
    bar_html = remove_subtitle(bar_html)

    if is_maxquant:
        description_text = "Potential contaminants per group from proteinGroups.txt."
        help_text = """
            External protein contamination should be controlled for, therefore MaxQuant ships with a 
            comprehensive, yet customizable protein contamination database, which is searched by MaxQuant 
            by default. 
            
            A contamination plot derived from the proteinGroups (PG) table, showing the 
            fraction of total protein intensity attributable to contaminants. 
            
            Note that this plot is based on experimental groups, and therefore may not correspond 1:1 to Raw files.
            """
    else:
        description_text = "Potential contaminants per file from mzTab."
        help_text = """
            A contamination plot derived from the Peptide section of the mzTab file, showing the 
            fraction of total intensity attributable to contaminants.
            """

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=2,
        description=description_text,
        helptext=help_text,
    )


def draw_top_n_contaminants(sub_section, top_contaminants_data):
    """
    Draw top N contaminants plot per raw file.

    Args:
        sub_section: MultiQC sub-section for adding plots
        top_contaminants_data: Dictionary containing plot data and categories
    """
    draw_config = {
        "id": "top_contaminants_per_raw_file",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Top5 Contaminants Per Raw File",
        "tt_decimals": 2,
        "ylab": "Identified [%]",
    }

    bar_html = bargraph.plot(
        top_contaminants_data["plot_data"],
        cats=top_contaminants_data["cats"],
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description="The five most abundant external protein contaminants by Raw file",
        helptext="""
            pmultiqc will explicitly show the five most abundant external protein contaminants 
            (as detected via MaxQuant's contaminants FASTA file) by Raw file, and summarize the 
            remaining contaminants as 'other'. This allows to track down which proteins exactly 
            contaminate your sample. Low contamination is obviously better.
            
            If you see less than 5 contaminants, it either means there are actually less, or that 
            one (or more) of the shortened contaminant names subsume multiple of the top5 
            contaminants (since they start with the same prefix).
            """,
    )


def remove_subtitle(plot_html):
    """
    Remove subtitle from plot HTML.

    Args:
        plot_html: Plot HTML object

    Returns:
        Modified plot HTML object
    """
    for dataset in plot_html.datasets:
        if "subtitle" in dataset.dconfig:
            dataset.dconfig["subtitle"] = ""

        title_text = ""
        if dataset.layout and "title" in dataset.layout:
            title_text = dataset.layout["title"].get("text", "")
            if "<br><sup>" in title_text:
                dataset.layout["title"]["text"] = title_text.split("<br><sup>")[0]

    return plot_html
