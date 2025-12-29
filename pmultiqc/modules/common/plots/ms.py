from multiqc.plots import bargraph, linegraph, table

from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.plots.general import plot_html_check


def draw_ms_information(
        sub_sections,
        ms1_tic=None,
        ms1_bpc=None,
        ms1_peaks=None,
        general_stats_data=None
):
    if ms1_tic:
        ms1_tic_config = {
            "id": "ms1_tic",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "Total Ion Chromatograms",
            "ylab": "Ion Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
            "save_data_file": False,
        }
        ms1_tic_html = linegraph.plot(ms1_tic, ms1_tic_config)
        ms1_tic_html = plot_html_check(ms1_tic_html)
        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_tic_html,
            order=1,
            description="MS1 quality control information extracted from the spectrum files.",
            helptext="""
                This plot displays Total Ion Chromatograms (TICs) derived from MS1 scans across all analyzed samples. 
                The x-axis represents retention time, and the y-axis shows the total ion intensity at each time point. 
                Each colored trace corresponds to a different sample. The TIC provides a global view of the ion signal
                throughout the LC-MS/MS run, reflecting when compounds elute from the chromatography column. 
            """,
        )

    if ms1_bpc:
        ms1_bpc_config = {
            "id": "ms1_bpc",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "MS1 Base Peak Chromatograms",
            "ylab": "Ion Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
            "save_data_file": False,
        }
        ms1_bpc_html = linegraph.plot(ms1_bpc, ms1_bpc_config)
        ms1_bpc_html = plot_html_check(ms1_bpc_html)
        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_bpc_html,
            order=2,
            description="MS1 base peak chromatograms extracted from the spectrum files.",
            helptext="""
                The Base Peak Chromatogram (BPC) displays the intensity of the most abundant ion at each retention time point.
            """,
        )

    if ms1_peaks:
        ms1_peaks_config = {
            "id": "ms1_peaks",
            "tt_label": "<b>{point.x} Peak Count:</b> {point.y}",
            "title": "MS1 Peaks",
            "ylab": "Peak Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
            "save_data_file": False,
        }
        ms1_peaks_html = linegraph.plot(ms1_peaks, ms1_peaks_config)
        ms1_peaks_html = plot_html_check(ms1_peaks_html)
        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_peaks_html,
            order=3,
            description="MS1 Peaks from the spectrum files",
            helptext="""
                This plot shows the number of peaks detected in MS1 scans over the course of each sample run.
            """,
        )

    if general_stats_data:
        tconfig = {
            "id": "ms_general_stats",
            "title": "General stats for MS1 information",
            "save_file": False,
            "raw_data_fn": "ms_general_stats_table",
            "no_violin": True,
            "save_data_file": False,
        }
        headers = {
            "Sample": {
                "title": "Sample [Spectra File]",
                "description": "",
                "scale": False,
            },
            "AcquisitionDateTime": {
                "title": "Acquisition Date Time",
            },
            "log10(TotalCurrent)": {
                "title": "log10(Total Current)",
                "format": "{:,.4f}",
            },
            "log10(ScanCurrent)": {
                "title": "log10(Scan Current)",
                "format": "{:,.4f}",
            },
        }

        table_html = table.plot(general_stats_data, headers=headers, pconfig=tconfig)
        add_sub_section(
            sub_section=sub_sections,
            plot=table_html,
            order=4,
            description="General stats for MS1 information extracted from the spectrum files.",
            helptext="""
                This table presents general statistics for MS1 information extracted from mass spectrometry data files.
            """,
        )


def draw_peak_intensity_distribution(
        sub_sections,
        peak_distribution_plot,
        ms_info,
):
    pconfig = {
        "id": "peak_intensity_distribution",
        "title": "Peak Intensity Distribution",
        "cpswitch": False,
        "stacking": "group",
        "logswitch": True,
        "logswitch_active": True,
        "logswitch_label": "Log10 Scale",
        "ylab": "Count",
        "tt_decimals": 0,
        "save_data_file": False,
    }
    cats = peak_distribution_plot.dict["cats"]
    bar_html = bargraph.plot(ms_info["peak_distribution"], cats, pconfig)
    bar_html = plot_html_check(bar_html)
    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=2,
        description="""
            Histogram of ion intensity vs. frequency for all MS2 spectra.
        """,
        helptext="""
            High number of low intensity noise peaks expected; disproportionate high signal peaks may indicate issues.
        """,
    )


def draw_precursor_charge_distribution(sub_sections, charge_plot=None, ms_info=None):
    if not charge_plot or not ms_info or "charge_distribution" not in ms_info:
        return

    pconfig = {
        "id": "distribution_of_precursor_charges",
        "title": "Distribution of Precursor Charges",
        "cpswitch": True,
        "tt_decimals": 0,
        "ylab": "Count",
        "save_data_file": False,
    }
    cats = charge_plot.dict["cats"]
    bar_html = bargraph.plot(ms_info["charge_distribution"], cats, pconfig)
    bar_html = plot_html_check(bar_html)
    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=5,
        description="Bar chart of precursor ion charge distribution.",
        helptext="Use to identify potential ionization problems or unexpected distributions.",
    )


def draw_peaks_per_ms2(sub_sections, peaks_ms2_plot, ms_info):
    pconfig = {
        "id": "peaks_per_ms2",
        "cpswitch": False,
        "title": "Number of Peaks per MS/MS spectrum",
        "stacking": "group",
        "logswitch": True,
        "logswitch_active": True,
        "logswitch_label": "Log10 Scale",
        "ylab": "Count",
        "tt_decimals": 0,
        "save_data_file": False,
    }
    cats = peaks_ms2_plot.dict["cats"]
    bar_html = bargraph.plot(ms_info["peaks_per_ms2"], cats, pconfig)
    bar_html = plot_html_check(bar_html)
    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=1,
        description="Histogram of number of peaks per MS/MS spectrum.",
        helptext="Too few peaks may indicate poor fragmentation; many peaks could indicate noisy spectra.",
    )