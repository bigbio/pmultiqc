"""
Retention time-related plotting functions for PMultiQC

This module contains functions for plotting retention time analysis and oversampling.
"""

from multiqc.plots import bargraph, linegraph
from pmultiqc.modules.core.section_groups import add_sub_section
from .contaminants import remove_subtitle


def draw_ids_rt_count(sub_section, rt_count_data, report_type):
    """
    Draw IDs over retention time plot.

    Args:
        sub_section: MultiQC sub-section for adding plots
        rt_count_data: Dictionary containing retention time count data
        report_type: String indicating the type of report
    """
    draw_config = {
        "id": "IDs_over_RT",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "IDs over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "Count",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(data=rt_count_data, pconfig=draw_config)
    linegraph_html = remove_subtitle(linegraph_html)

    if report_type == "maxquant":
        description_text = "Distribution of retention time, derived from the evidence table."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion, 
            and does not include potential contaminants.
            """
    elif report_type == "dia":
        description_text = "Distribution of retention time, derived from the main report."
        help_text = """
            [DIA-NN: main report] Distribution of retention time (RT) for each run.
            """
    elif report_type == "mzIdentML":
        description_text = "Distribution of retention time, derived from the mzIdentML (or mzML)."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion.
            """
    else:
        description_text = "Distribution of retention time, derived from the mzTab."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion.
            """

    help_text += """
            <p>This plot allows to judge column occupancy over retention time.
            Ideally, the LC gradient is chosen such that the number of identifications (here, after FDR filtering) is
            uniform over time, to ensure consistent instrument duty cycles. Sharp peaks and uneven distribution of
            identifications over time indicate potential for LC gradient optimization.
            See [Moruz 2014, DOI: 10.1002/pmic.201400036](https://pubmed.ncbi.nlm.nih.gov/24700534/) for details.</p>
"""

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=1,
        description=description_text,
        helptext=help_text,
    )


def draw_oversampling(sub_section, oversampling, oversampling_plot, is_maxquant):
    """
    Draw oversampling plot.

    Args:
        sub_section: MultiQC sub-section for adding plots
        oversampling: Dictionary containing oversampling data
        oversampling_plot: Plot configuration for oversampling
        is_maxquant: Boolean indicating if data is from MaxQuant
    """
    if is_maxquant:
        draw_config = {
            "id": "oversampling_distribution",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "MS/MS Counts Per 3D-peak",
            "tt_decimals": 2,
            "ylab": "MS/MS Counts Per 3D-peak [%]",
        }

        bar_html = bargraph.plot(
            data=oversampling["plot_data"], cats=oversampling["cats"], pconfig=draw_config
        )

    else:
        draw_config = {
            "id": "oversampling_distribution",
            "cpswitch": True,
            "cpswitch_c_active": False,
            "title": "MS/MS Counts Per 3D-peak",
            "ylab": "MS/MS Counts Per 3D-peak [%]",
            "tt_decimals": 0,
        }

        bar_html = bargraph.plot(data=oversampling, cats=oversampling_plot, pconfig=draw_config)

    bar_html = remove_subtitle(bar_html)

    helptext = """
                For high complexity samples, oversampling of individual 3D-peaks automatically leads to 
                undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. 
                Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic 
                exclusion windows for data independent acquisitions.
                """
    if is_maxquant:
        helptext += "<p>If DIA-Data: this metric is skipped.</p>"

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=7,
        description="""
            An oversampled 3D-peak is defined as a peak whose peptide ion 
            (same sequence and same charge state) was identified by at least two distinct MS2 spectra 
            in the same Raw file.
            """,
        helptext=helptext,
    )
