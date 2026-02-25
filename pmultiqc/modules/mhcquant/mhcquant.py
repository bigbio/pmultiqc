"""MultiQC pmultiqc plugin module for nf-core/mhcquant pipeline."""

import os

from multiqc.plots import linegraph, box, table, bargraph

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.core.section_groups import add_sub_section, add_group_modules
from pmultiqc.modules.common.logging import get_logger

log = get_logger("pmultiqc.modules.mhcquant.mhcquant")

# Files to look for relative to the directory containing multiqc_general_stats.txt
DATA_FILES = {
    "general_stats": "multiqc_general_stats.txt",
    "chromatogram": "multiqc_chromatogram.txt",
    "histogram_mz": "multiqc_histogram_mz.txt",
    "histogram_rt": "multiqc_histogram_rt.txt",
    "histogram_im": "multiqc_histogram_im.txt",
    "histogram_scores": "multiqc_histogram_scores.txt",
    "length_dist": "multiqc_length_dist.txt",
    "mass_error": "multiqc_mass_error.txt",
    "scores_xcorr": "multiqc_scores_xcorr.txt",
    "peptide_intensity": "multiqc_peptide_intensity.txt",
    "percolator_plot": "percolator_plot.txt",
}


def _parse_histogram_file(filepath):
    """Parse a TSV histogram file with header row into linegraph data.

    Format: Sample<tab>bin1<tab>bin2<tab>...
            name<tab>val1<tab>val2<tab>...

    Returns: {sample_name: {x_float: y_float, ...}}, or None on failure.
    """
    if not os.path.isfile(filepath):
        return None
    with open(filepath, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    bins = []
    for b in header[1:]:
        b = b.strip()
        if not b:
            continue
        try:
            bins.append(float(b))
        except ValueError:
            log.warning("Skipping non-numeric histogram bin label %r in %s", b, filepath)
    if not bins:
        return None
    data = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        sample = parts[0]
        row = {}
        for i, b in enumerate(bins):
            val_str = parts[i + 1].strip() if i + 1 < len(parts) else ""
            if val_str:
                try:
                    row[b] = float(val_str)
                except ValueError:
                    pass
        if row:
            data[sample] = row
    return data if data else None


def _parse_headerless_box_file(filepath):
    """Parse a TSV file with no header into boxplot data.

    Format: sample_name<tab>val1<tab>val2<tab>...

    Returns: {sample_name: [float, ...]}, or None on failure.
    """
    if not os.path.isfile(filepath):
        return None
    with open(filepath, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    if not lines:
        return None
    data = {}
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        sample = parts[0]
        values = []
        for v in parts[1:]:
            v = v.strip()
            if not v:
                continue
            try:
                values.append(float(v))
            except ValueError:
                log.warning("Skipping non-numeric value %r in %s", v, filepath)
                return None
        if values:
            data[sample] = values
    return data if data else None


def _parse_general_stats(filepath):
    """Parse general stats TSV into table data.

    Format: Sample<tab>col1<tab>col2<tab>...
            sample<tab>val1<tab>val2<tab>...

    Returns: (data_dict, headers_dict) for table.plot(), or (None, None) on failure.
    """
    if not os.path.isfile(filepath):
        return None, None
    with open(filepath, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    if len(lines) < 2:
        return None, None
    header = lines[0].split("\t")
    col_names = [c.strip() for c in header[1:] if c.strip()]
    if not col_names:
        return None, None
    data = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        sample = parts[0]
        row = {}
        for i, col in enumerate(col_names):
            val = parts[i + 1].strip() if i + 1 < len(parts) else ""
            if not val:
                continue
            try:
                row[col] = float(val) if "." in val else int(val)
            except (ValueError, TypeError):
                row[col] = val
        data[sample] = row
    headers = {col: {"title": col} for col in col_names}
    return data, headers


def _parse_percolator_median_weights(filepath):
    """Parse percolator_plot.txt into grouped bargraph data.

    Format: Sample<tab>group1<tab>group2<tab>...
            feature1<tab>val1<tab>val2<tab>...

    Returns: (data_dict, group_names) for bargraph.plot(), or (None, None) on failure.
    where data_dict = {feature: {group: median_weight}}
    """
    if not os.path.isfile(filepath):
        return None, None
    with open(filepath, encoding="utf-8") as fh:
        lines = [line.rstrip("\n") for line in fh if line.strip()]
    if len(lines) < 2:
        return None, None
    header = lines[0].split("\t")
    groups = [g.strip() for g in header[1:] if g.strip()]
    if not groups:
        return None, None

    data = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        feature = parts[0].strip()
        row = {}
        for i, group in enumerate(groups):
            v = parts[i + 1].strip() if i + 1 < len(parts) else ""
            if v and v.lower() != "nan":
                try:
                    row[group] = float(v)
                except (ValueError, TypeError):
                    log.warning("Skipping non-numeric percolator weight %r in %s", v, filepath)
        if row:
            data[feature] = row

    if not data:
        return None, None

    # Sort features from most negative to most positive weight
    sorted_data = dict(
        sorted(data.items(), key=lambda item: sum(item[1].values()))
    )

    return sorted_data, groups


class MhcquantModule(BasePMultiqcModule):

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):
        super().__init__(find_log_files_func, sub_sections, heatmap_colors)
        self.data = {}
        self.data_root = None

    def get_data(self):
        log.info("Starting mhcquant data recognition and processing...")

        # Search analysis directories for multiqc_general_stats.txt.
        # We walk the directory tree directly because MultiQC's find_log_files()
        # ignores multiqc_data/ by default and wildcard patterns (*.txt) can
        # consume our anchor file before our pattern is checked.
        from multiqc import config as mqc_config

        for analysis_dir in mqc_config.analysis_dir:
            for root, dirs, fnames in os.walk(analysis_dir):
                if "multiqc_general_stats.txt" in fnames:
                    self.data_root = root
                    log.info(f"Found mhcquant data directory: {self.data_root}")
                    break
            if self.data_root is not None:
                break

        if self.data_root is None:
            log.error("mhcquant multiqc_general_stats.txt not found.")
            return False

        # Load all data files
        for key, filename in DATA_FILES.items():
            filepath = os.path.join(self.data_root, filename)
            if not os.path.isfile(filepath):
                log.info(f"Optional file not found: {filename}")
                continue

            if key == "general_stats":
                table_data, table_headers = _parse_general_stats(filepath)
                if table_data:
                    self.data["general_stats"] = table_data
                    self.data["general_stats_headers"] = table_headers
            elif key in ("mass_error", "scores_xcorr", "peptide_intensity"):
                parsed = _parse_headerless_box_file(filepath)
                if parsed:
                    self.data[key] = parsed
            elif key == "percolator_plot":
                bar_data, groups = _parse_percolator_median_weights(filepath)
                if bar_data:
                    self.data["percolator_median"] = bar_data
                    self.data["percolator_groups"] = groups
            else:
                # Histogram-style files (chromatogram, histogram_mz, etc.)
                parsed = _parse_histogram_file(filepath)
                if parsed:
                    self.data[key] = parsed

        log.info(f"Loaded {len(self.data)} data sections from mhcquant results.")
        return "general_stats" in self.data

    def draw_plots(self):
        log.info("Generating mhcquant plots...")

        # Match mhcquant's global config: boxplot_boxpoints: false
        from multiqc import config as mqc_config
        mqc_config.boxplot_boxpoints = False

        # General Statistics table -> summary section
        if "general_stats" in self.data:
            pconfig = {
                "id": "general_stats",
                "title": "General Statistics",
                "save_file": False,
                "no_violin": True,
                "save_data_file": False,
            }
            table_html = table.plot(
                self.data["general_stats"],
                self.data["general_stats_headers"],
                pconfig,
            )
            add_sub_section(
                sub_section=self.sub_sections["summary"],
                plot=table_html,
                order=1,
                description="General statistics from the mhcquant pipeline.",
            )

        # MS1 Chromatogram -> ms1 section
        if "chromatogram" in self.data:
            pconfig = {
                "id": "mhcquant_chromatogram",
                "title": "MS1 Chromatogram",
                "xlab": "Retention time [min]",
                "ylab": "Intensity",
                "xmin": 0,
                "logswitch": True,
                "logswitch_active": True,
                "tt_decimals": 0,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["chromatogram"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["ms1"],
                plot=line_html,
                order=1,
                description=(
                    "An MS1 chromatogram is a plot of the intensity of precursor ions (MS1) "
                    "detected over time, typically with retention time (RT) on the x-axis and "
                    "ion intensity on the y-axis. It represents the total ion current (TIC) from "
                    "the MS1 scan, which is useful for monitoring the overall elution profile of "
                    "compounds during a chromatographic separation."
                ),
            )

        # Mass-to-Charge -> ms1 section
        if "histogram_mz" in self.data:
            pconfig = {
                "id": "histogram_mz",
                "title": "Mass-to-Charge Histogram",
                "xlab": "m/z",
                "ylab": "Frequency",
                "tt_decimals": 0,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["histogram_mz"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["ms1"],
                plot=line_html,
                order=2,
                description=(
                    "This plot displays the distribution of precursor ion mass-to-charge (m/z) "
                    "values detected in the dataset. It helps identify the typical m/z ranges "
                    "and assess the overall coverage and instrument performance during acquisition."
                ),
            )

        # Retention time -> rt_qc section
        if "histogram_rt" in self.data:
            pconfig = {
                "id": "histogram_rt",
                "title": "Retention time Histogram",
                "xlab": "Retention time",
                "ylab": "Frequency",
                "tt_decimals": 0,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["histogram_rt"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["rt_qc"],
                plot=line_html,
                order=1,
                description=(
                    "This plot shows the distribution of peptide retention times across the "
                    "chromatographic gradient. It provides insight into peptide elution profiles "
                    "and helps assess gradient performance and sample complexity."
                ),
            )

        # Ion mobility distribution -> ms1 section (optional, timsTOF only)
        if "histogram_im" in self.data:
            pconfig = {
                "id": "mhcquant_histogram_im",
                "title": "Ion Mobility Distribution",
                "xlab": "Ion Mobility (1/K0)",
                "ylab": "Frequency",
                "tt_decimals": 0,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["histogram_im"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["ms1"],
                plot=line_html,
                order=3,
                description="Distribution of ion mobility values (timsTOF data only).",
            )

        # Q-Value -> search_engine section
        if "histogram_scores" in self.data:
            pconfig = {
                "id": "histogram_scores",
                "title": "Percolator q-value Histogram",
                "xlab": "q-value",
                "ylab": "Frequency",
                "tt_decimals": 0,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["histogram_scores"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=line_html,
                order=1,
                description=(
                    "A q-value histogram in Percolator visualizes q-value distributions for "
                    "PSMs, helping set the FDR threshold (typically 1%). It ensures only "
                    "high-confidence identifications are retained, with a sharp drop-off "
                    "beyond the threshold."
                ),
            )

        # Peptide lengths -> identification section
        if "length_dist" in self.data:
            pconfig = {
                "id": "length_dist",
                "title": "Distribution of peptide lengths",
                "xlab": "Length",
                "ylab": "Frequency",
                "tt_decimals": 3,
                "save_data_file": False,
            }
            line_html = linegraph.plot(data=self.data["length_dist"], pconfig=pconfig)
            add_sub_section(
                sub_section=self.sub_sections["identification"],
                plot=line_html,
                order=1,
                description=(
                    "This plot shows the distribution of peptide sequence lengths to assess "
                    "variability and typical length ranges."
                ),
            )

        # Fragment mass error -> mass_error section (boxplot)
        if "mass_error" in self.data:
            pconfig = {
                "id": "mhcquant_mass_error",
                "title": "Fragment mass error",
                "xlab": "Mass deviation [Da/ppm]",
                "tt_decimals": 5,
                "boxpoints": False,
                "save_data_file": False,
            }
            box_html = box.plot(
                list_of_data_by_sample=self.data["mass_error"], pconfig=pconfig
            )
            add_sub_section(
                sub_section=self.sub_sections["mass_error"],
                plot=box_html,
                order=1,
                description=(
                    "The fragment mass error describes the differences between observed and "
                    "theoretical fragment ion masses. It helps assessing the accuracy and "
                    "calibration of the mass spectrometer by visualizing the error range and "
                    "potential deviations. Depending on your configurations this error is "
                    "plotted in Dalton (Da) or parts-per-million (ppm)."
                ),
            )

        # Xcorr-Value -> search_engine section (boxplot)
        if "scores_xcorr" in self.data:
            pconfig = {
                "id": "scores_xcorr",
                "title": "Comet Xcorr Distribution",
                "xlab": "Xcorr",
                "ylab": "Frequency",
                "tt_decimals": 3,
                "boxpoints": False,
                "save_data_file": False,
            }
            box_html = box.plot(
                list_of_data_by_sample=self.data["scores_xcorr"], pconfig=pconfig
            )
            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=box_html,
                order=2,
                description=(
                    "The Comet Xcorr boxplot visualizes the cross-correlation (Xcorr) score "
                    "distribution for PSMs, indicating match quality. Higher Xcorr values "
                    "signify better peptide-spectrum matches. This helps set a threshold to "
                    "filter out low-confidence identifications."
                ),
            )

        # Peptide Intensity -> quantification section (boxplot, optional)
        if "peptide_intensity" in self.data:
            pconfig = {
                "id": "peptide_intensity",
                "title": "Peptide Intensity Distribution",
                "xlab": "log2(Intensity)",
                "tt_decimals": 2,
                "boxpoints": False,
                "save_data_file": False,
            }
            box_html = box.plot(
                list_of_data_by_sample=self.data["peptide_intensity"], pconfig=pconfig
            )
            add_sub_section(
                sub_section=self.sub_sections["quantification"],
                plot=box_html,
                order=1,
                description=(
                    "Comparison of log2-transformed peptide intensities across samples. "
                    "Helps assess data consistency and detect potential outliers or batch effects."
                ),
            )

        # Percolator Median Feature Weights -> search_engine section (bargraph)
        if "percolator_median" in self.data:
            cats = {g: {"name": g} for g in self.data["percolator_groups"]}
            pconfig = {
                "id": "percolator_median_weights",
                "title": "Median Feature Weights",
                "ylab": "Weight",
                "tt_decimals": 4,
                "save_data_file": False,
            }
            bar_html = bargraph.plot(
                data=self.data["percolator_median"],
                cats=cats,
                pconfig=pconfig,
            )
            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=bar_html,
                order=3,
                description=(
                    "For each feature, the associated weight is the median value calculated "
                    "over all input samples."
                ),
                helptext=(
                    "The bar plot illustrates the median weights for each feature, which have "
                    "been assigned to categories based on their source (e.g. psm_file, ms2pip, "
                    "deeplc, im2deep)."
                ),
            )

        # Assemble section groups and finalize report layout
        section_group_dict = {
            "experiment_sub_section": self.sub_sections["experiment"],
            "summary_sub_section": self.sub_sections["summary"],
            "identification_sub_section": self.sub_sections["identification"],
            "search_engine_sub_section": self.sub_sections["search_engine"],
            "contaminants_sub_section": self.sub_sections["contaminants"],
            "quantification_sub_section": self.sub_sections["quantification"],
            "ms1_sub_section": self.sub_sections["ms1"],
            "ms2_sub_section": self.sub_sections["ms2"],
            "mass_error_sub_section": self.sub_sections["mass_error"],
            "rt_qc_sub_section": self.sub_sections["rt_qc"],
        }

        add_group_modules(section_group_dict, "")
        log.info("mhcquant plot generation completed.")
