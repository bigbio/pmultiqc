""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
import os

import pandas as pd

from pmultiqc.modules.common import ms_io
from pmultiqc.modules.common.common_utils import (
    parse_sdrf,
    get_ms_path,
    parse_mzml
)
from pmultiqc.modules.common.dia_utils import parse_diann_report
from pmultiqc.modules.common.plots.general import draw_exp_design
from pmultiqc.modules.common.plots.id import draw_quantms_identification
from pmultiqc.modules.common.plots.id import (
    draw_summary_protein_ident_table,
    draw_quantms_identi_num,
    draw_num_pep_per_protein
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


class DiannModule:

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections
        self.heatmap_color_list = heatmap_colors

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
            self.ms1_general_stats
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
        self.cal_num_table_data = dict()
        self.quantms_modified = dict()

        draw_ms_information(
            self.sub_sections["ms1"],
            self.ms1_tic,
            self.ms1_bpc,
            self.ms1_peaks,
            self.ms1_general_stats
        )

        (
            self.total_protein_quantified,
            self.total_peptide_count,
            self.pep_plot,
            self.peptide_search_score,
            self.ms_with_psm,
            self.cal_num_table_data,
            self.quantms_modified,
            self.ms_without_psm
        ) = parse_diann_report(
            sub_sections=self.sub_sections,
            diann_report_path=self.diann_report_path,
            heatmap_color_list=self.heatmap_color_list,
            sample_df=self.sample_df,
            file_df=self.file_df,
            ms_with_psm=self.ms_with_psm,
            cal_num_table_data=self.cal_num_table_data,
            quantms_modified=self.quantms_modified,
            ms_paths=self.ms_paths
        )

        draw_summary_protein_ident_table(
            sub_sections=self.sub_sections["summary"],
            enable_dia=self.enable_dia,
            total_peptide_count=self.total_peptide_count,
            total_protein_quantified=self.total_protein_quantified
        )

        draw_quantms_identi_num(
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

        draw_quantms_identification(
            self.sub_sections["identification"],
            cal_num_table_data=self.cal_num_table_data,
            mzml_table=self.mzml_table,
            quantms_modified=self.quantms_modified
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