"""
MS analysis plotting functions for PMultiQC

This module contains functions for plotting MS1 and MS2 analysis results.
"""

from multiqc.plots import bargraph, linegraph, table
from pmultiqc.modules.core.section_groups import add_sub_section
from .contaminants import remove_subtitle


def draw_peaks_per_ms2(sub_sections, peaks_ms2_plot, ms_info):
    """
    Draw peaks per MS2 plot.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        peaks_ms2_plot: Plot data object containing peaks per MS2 statistics
        ms_info: Dictionary containing MS information
    """
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
    }

    cats = peaks_ms2_plot.dict["cats"]
    bar_html = bargraph.plot(ms_info["peaks_per_ms2"], cats, pconfig)
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=1,
        description="""
            This chart represents a histogram containing the number of peaks per MS/MS spectrum 
            in a given experiment.
            """,
        helptext="""
            This chart assumes centroid data. Too few peaks can identify poor fragmentation 
                or a detector fault, as opposed to a large number of peaks representing very noisy spectra. 
                This chart is extensively dependent on the pre-processing steps performed to the spectra 
                (centroiding, deconvolution, peak picking approach, etc).
                """,
        )


def draw_peak_intensity_distribution(sub_sections, peak_distribution_plot, ms_info):
    """
    Draw peak intensity distribution plot.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        peak_distribution_plot: Plot data object containing peak distribution statistics
        ms_info: Dictionary containing MS information
    """
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
    }
    
    cats = peak_distribution_plot.dict["cats"]
    bar_html = bargraph.plot(ms_info["peak_distribution"], cats, pconfig)
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=2,
        description="""
            This is a histogram representing the ion intensity vs.
            the frequency for all MS2 spectra in a whole given experiment.
            It is possible to filter the information for all, identified and unidentified spectra.
            This plot can give a general estimation of the noise level of the spectra.
            """,
        helptext="""
            Generally, one should expect to have a high number of low intensity noise peaks with a low number 
            of high intensity signal peaks. 
            A disproportionate number of high signal peaks may indicate heavy spectrum pre-filtering or 
            potential experimental problems. In the case of data reuse this plot can be useful in 
            identifying the requirement for pre-processing of the spectra prior to any downstream analysis. 
            The quality of the identifications is not linked to this data as most search engines perform internal 
            spectrum pre-processing before matching the spectra. Thus, the spectra reported are not 
            necessarily pre-processed since the search engine may have applied the pre-processing step 
            internally. This pre-processing is not necessarily reported in the experimental metadata.
            """,
    )


def draw_ms_information(sub_sections, ms1_tic=None, ms1_bpc=None, ms1_peaks=None, ms1_general_stats=None):
    """
    Draw MS1 information plots including TIC, BPC, peaks, and general stats.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        ms1_tic: Dictionary containing TIC data
        ms1_bpc: Dictionary containing BPC data
        ms1_peaks: Dictionary containing MS1 peaks data
        ms1_general_stats: Dictionary containing general MS1 statistics
    """
    if ms1_tic:
        ms1_tic_config = {
            "id": "ms1_tic",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "Total Ion Chromatograms",
            "ylab": "Ion Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
        }
        ms1_tic_html = linegraph.plot(ms1_tic, ms1_tic_config)
        ms1_tic_html = remove_subtitle(ms1_tic_html)

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
                Key aspects to assess include:
                
                * Overall intensity pattern: A consistent baseline and similar peak profiles across samples indicate good reproducibility.
                * Major peak alignment: Prominent peaks appearing at similar retention times suggest stable chromatographic performance.
                * Signal-to-noise ratio: High peaks relative to baseline noise reflect better sensitivity.
                * Chromatographic resolution: Sharp, well-separated peaks indicate effective separation.
                * Signal drift: A gradual decline in signal intensity across the run may point to source contamination or chromatography issues.
                
                Deviations such as shifted retention times, missing peaks, or inconsistent intensities may signal problems 
                in sample preparation, LC conditions, or mass spectrometer performance that require further investigation.
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
        }

        ms1_bpc_html = linegraph.plot(ms1_bpc, ms1_bpc_config)
        ms1_bpc_html = remove_subtitle(ms1_bpc_html)

        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_bpc_html,
            order=2,
            description="MS1 base peak chromatograms extracted from the spectrum files.",
            helptext="""
                The Base Peak Chromatogram (BPC) displays the intensity of the most abundant ion at each retention
                time point across your LC-MS run. Unlike the Total Ion Chromatogram (TIC) which shows the summed
                intensity of all ions, the BPC highlights the strongest signals, providing better visualization
                of compounds with high abundance while reducing baseline noise. This makes it particularly useful
                for identifying major components in complex samples, monitoring dominant species, and providing
                clearer peak visualization when signal-to-noise ratio is a concern. Comparing BPC patterns across
                samples allows you to evaluate consistency in the detection of high-abundance compounds
                and can reveal significant variations in sample composition or instrument performance.
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
        }

        ms1_peaks_html = linegraph.plot(ms1_peaks, ms1_peaks_config)
        ms1_peaks_html = remove_subtitle(ms1_peaks_html)

        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_peaks_html,
            order=3,
            description="MS1 Peaks from the spectrum files",
            helptext="""
                This plot shows the number of peaks detected in MS1 scans over the course of each sample run.
                The x-axis represents retention time (in minutes), while the y-axis displays the number of
                distinct ion signals (peaks) identified in each MS1 scan. The MS1 peak count reflects spectral
                complexity and provides insight into instrument performance during the LC-MS analysis.
                Key aspects to consider include:
                
                * Overall pattern: Peak counts typically increase during the elution of complex mixtures and
                decrease during column washing or re-equilibration phases.
                * Peak density: Higher counts suggest more complex spectra, potentially indicating a greater
                number of compounds present at that time point.
                * Peak Consistency across samples: Similar profiles among replicates or related samples
                indicate good analytical reproducibility.
                * Sudden drops: Abrupt decreases in peak count may point to transient ionization issues,
                spray instability, or chromatographic disruptions.
                * Baseline values: The minimum peak count observed reflects the level of background noise
                or instrument sensitivity in the absence of eluting compounds.
                
                Monitoring MS1 peak counts complements total ion chromatogram (TIC) and base peak chromatogram
                (BPC) data, offering an additional layer of quality control related to signal complexity,
                instrument stability, and sample composition.
                """,
        )

    if ms1_general_stats:
        tconfig = {
            "id": "ms_general_stats",
            "title": "General stats for MS1 information",
            "only_defined_headers": True,
            "col1_header": "File",
        }
        headers = {
            "File": {
                "title": "File",
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
        table_html = table.plot(ms1_general_stats, headers=headers, pconfig=tconfig)

        add_sub_section(
            sub_section=sub_sections,
            plot=table_html,
            order=4,
            description="General stats for MS1 information extracted from the spectrum files.",
            helptext="""
                This table presents general statistics for MS1 information extracted from mass spectrometry data files."
                It displays MS runs with their acquisition dates and times.
                For each file, the table shows two key metrics: TotalCurrent (the sum of all MS1 ion intensities throughout the run) and ScanCurrent (the sum of MS2 ion intensities).
                These values provide a quick overview of the total ion signals detected during both survey scans (MS1) and fragmentation scans (MS2),
                allowing for comparison of overall signal intensity across samples. Consistent TotalCurrent and ScanCurrent values across similar samples typically
                indicate good reproducibility in the mass spectrometry analysis, while significant variations may suggest issues with sample preparation,
                instrument performance, or ionization efficiency. The blue shading helps visualize the relative intensity differences between samples.
                """,
        )
