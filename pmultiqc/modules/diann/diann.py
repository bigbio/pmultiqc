""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
import os

import pandas as pd

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.common import ms_io
from pmultiqc.modules.common.common_utils import (
    parse_sdrf,
    get_ms_path,
    parse_mzml,
    aggregate_general_stats
)
from pmultiqc.modules.common.dia_utils import (
    parse_diann_report,
    parse_diann_version,
    draw_diann_metadata_table
)
from pmultiqc.modules.common.plots.general import draw_exp_design
from pmultiqc.modules.common.plots.id import (
    draw_summary_protein_ident_table,
    draw_identi_num,
    draw_num_pep_per_protein,
    draw_identification,
    draw_long_trends,
    draw_peptide_length_distribution
)
from pmultiqc.modules.common.plots.ms import (
    draw_peak_intensity_distribution,
    draw_peaks_per_ms2,
    draw_ms_information
)
from pmultiqc.modules.core.section_groups import add_group_modules
from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.diann.diann")


class DiannModule(BasePMultiqcModule):

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        super().__init__(find_log_files_func, sub_sections, heatmap_colors)


    def get_data(self):

        log.info("Starting data recognition and processing...")

        self.enable_sdrf = False
        for f in self.find_log_files("pmultiqc/sdrf", filecontents=False):
            self.sdrf = os.path.join(f["root"], f["fn"])

            parse_sdrf(self.sdrf)

            # experimental_design.tsv is the default output name
            # experimental_design.tsv will be in the folder where pmultiqc is executed.
            self.exp_design = "experimental_design.tsv"
            self.enable_sdrf = True

        self.sample_df = pd.DataFrame()
        self.file_df = pd.DataFrame()
        self.is_bruker = False
        self.is_multi_conditions = False
        self.ms_with_psm = list()

        # draw the experimental design
        if self.enable_sdrf:
            (
                self.sample_df,
                self.file_df,
                self.exp_design_runs,
                self.is_bruker,
                self.is_multi_conditions
            ) = draw_exp_design(
                self.sub_sections["experiment"],
                self.exp_design
            )

        (
            self.ms_info_path,
            self.read_ms_info,
            self.ms_paths
        ) = get_ms_path(self.find_log_files)

        # DIA-NN report file path
        diann_report_path = None
        for file_type in ["pmultiqc/diann_report_tsv", "pmultiqc/diann_report_parquet"]:
            for f in self.find_log_files(file_type, filecontents=False):
                diann_report_path = os.path.join(f["root"], f["fn"])
            if diann_report_path:
                break

        if diann_report_path:
            self.diann_report_path = diann_report_path
            self.enable_dia = True
        else:
            log.error("DIANN report not found. Please check your data!")
            return False

        # DIA-NN log file for version extraction
        self.diann_version = None
        for key in ["pmultiqc/diann_log_txt", "pmultiqc/diann_log"]:
            for f in self.find_log_files(key, filecontents=False):
                log_path = os.path.join(f["root"], f["fn"])
                self.diann_version = parse_diann_version(log_path)

                if self.diann_version:
                    log.info(f"DIA-NN version detected: {self.diann_version}")
                    break

            if self.diann_version:
                break

        (
            self.mzml_table,
            self.mzml_peaks_ms2_plot,
            self.mzml_peak_distribution_plot,
            self.ms_info,
            _,  # total_ms2_spectra
            _,  # mzml_ms_df
            _,  # heatmap_charge_score
            _,  # mzml_charge_plot
            self.ms1_tic,
            self.ms1_bpc,
            self.ms1_peaks,
            self.ms1_general_stats,
            self.current_sum_by_run,
            self.long_trends
        ) = parse_mzml(
            is_bruker=self.is_bruker,
            read_ms_info=self.read_ms_info,
            ms_info_path=self.ms_info_path,
            ms_with_psm=self.ms_with_psm,
            enable_dia=self.enable_dia,
            ms_paths=self.ms_paths
        )

        log.info("Data recognition and processing completed.")
        
        return True

    def draw_plots(self):

        self.total_peptide_count = 0
        self.total_protein_quantified = 0
        self.cal_num_table_data = {}
        self.quantms_modified = {}

        # Draw DIA-NN metadata table (version info) in the experiment section
        if self.diann_version:
            draw_diann_metadata_table(self.sub_sections["experiment"], self.diann_version)

        general_stats_data = aggregate_general_stats(
            ms1_general_stats=self.ms1_general_stats,
            current_sum_by_run=self.current_sum_by_run,
            sdrf_file_df=self.file_df
        )

        draw_ms_information(
            self.sub_sections["ms1"],
            self.ms1_tic,
            self.ms1_bpc,
            self.ms1_peaks,
            general_stats_data
        )

        (
            self.total_protein_quantified,
            self.total_peptide_count,
            self.pep_plot,
            self.peptide_search_score,
            self.ms_with_psm,
            self.cal_num_table_data,
            self.quantms_modified,
            self.ms_without_psm,
            self.peptide_length
        ) = parse_diann_report(
            sub_sections=self.sub_sections,
            diann_report_path=self.diann_report_path,
            heatmap_color_list=self.heatmap_color_list,
            sample_df=self.sample_df,
            file_df=self.file_df,
            ms_with_psm=self.ms_with_psm,
            quantms_modified=self.quantms_modified,
            ms_paths=self.ms_paths
        )

        draw_summary_protein_ident_table(
            sub_sections=self.sub_sections["summary"],
            use_two_columns=self.enable_dia,
            total_peptide_count=self.total_peptide_count,
            total_protein_quantified=self.total_protein_quantified
        )

        draw_identi_num(
            sub_sections=self.sub_sections["summary"],
            enable_sdrf=self.enable_sdrf,
            is_multi_conditions=self.is_multi_conditions,
            sample_df=self.sample_df,
            file_df=self.file_df,
            cal_num_table_data=self.cal_num_table_data
        )

        draw_num_pep_per_protein(
            self.sub_sections["identification"],
            self.pep_plot
        )

        if len(self.ms_info_path) > 0 and not self.is_bruker:

            draw_peaks_per_ms2(
                self.sub_sections["ms2"],
                self.mzml_peaks_ms2_plot,
                self.ms_info
            )

            draw_peak_intensity_distribution(
                self.sub_sections["ms2"],
                self.mzml_peak_distribution_plot,
                self.ms_info
            )

        draw_identification(
            self.sub_sections["identification"],
            cal_num_table_data=self.cal_num_table_data,
            quantms_modified=self.quantms_modified
        )

        if self.long_trends:
            draw_long_trends(
                sub_sections=self.sub_sections,
                long_trends_data=self.long_trends
            )

        if self.peptide_length:
            draw_peptide_length_distribution(
                sub_section=self.sub_sections["identification"],
                plot_data=self.peptide_length
            )

        if self.enable_sdrf:
            ms_io.del_openms_convert_tsv()

        self.section_group_dict = {
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

        add_group_modules(self.section_group_dict, "")
        try:
            from pmultiqc.export.mzqc_exporter import MzQcExporter
            from multiqc.utils import config
            
            output_directory = getattr(config, "output_dir", "./")
            
            # Helper to safely extract raw data from MultiQC Plot objects (like Histogram)
            def get_raw_plot_data(plot_obj):
                if plot_obj is None:
                    return None
                # If it's a MultiQC Plot object, grab its underlying dataset
                if hasattr(plot_obj, "data"):
                    return plot_obj.data
                if hasattr(plot_obj, "plot_data"):
                    return plot_obj.plot_data
                return plot_obj

            def clean_modifications(mods_data):
                if not mods_data:
                    return None
                
                # If it's a MultiQC plot config dict, it often has 'datasets' or 'samples'
                if isinstance(mods_data, dict):
                    # Try to extract actual numeric data instead of UI config
                    if "datasets" in mods_data:
                        # Standard MultiQC plot structure
                        cleaned = {}
                        for series in mods_data["datasets"]:
                            series_name = series.get("name", "Unknown")
                            series_data = series.get("data", [])
                            cleaned[series_name] = series_data
                        return cleaned
                    
                    if "data" in mods_data:
                        return mods_data["data"]
                        
                return mods_data

            diann_payload = {
                "peptide_id_count": self.total_peptide_count,
                "protein_group_count": self.total_protein_quantified,
                "psm_count": self.ms_with_psm,
                "peptide_length_data": get_raw_plot_data(self.peptide_length),
                "peptide_per_protein_data": get_raw_plot_data(self.pep_plot),
                "modifications": clean_modifications(self.quantms_modified),  # 👈 Cleaned!
                "identification_summary": self.cal_num_table_data
            }
            
            exporter = MzQcExporter(
                pipeline_name="DIA-NN",
                raw_data=diann_payload, 
                output_dir=output_directory
            )
            
            mzqc_metrics = exporter._parse_diann(diann_payload)
            self.log.info(f"mzQC: Successfully extracted {len(mzqc_metrics)} metrics.")
            
            saved_file_path = exporter.export_to_file(mzqc_metrics, filename="diann_qc.mzQC")
            self.log.info(f"mzQC: Generated output saved directly to: {saved_file_path}")
            
        except Exception as e:
            self.log.warning(f"mzQC: Metric extraction or export failed: {e}")