from __future__ import absolute_import

import copy
import itertools
import json
import os
import re
import sqlite3
from collections import OrderedDict
from datetime import datetime
from functools import reduce
from operator import itemgetter

import numpy as np
import pandas as pd
from multiqc import config
from multiqc.plots import (
    table,
    bargraph,
    linegraph,
    box
)
from pyopenms import AASequence
from sdrf_pipelines.openms.openms import UnimodDatabase

from typing import Dict, List
from multiqc.plots.table_object import InputRow
from multiqc.types import SampleGroup, SampleName

from . import sparklines
from pmultiqc.modules.common.dia_utils import parse_diann_report
from pmultiqc.modules.common.common_utils import (
    parse_sdrf,
    get_ms_path,
    evidence_rt_count,
    evidence_calibrated_mass_error,
    parse_mzml,
    cal_num_table_at_sample,
    aggregate_msms_identified_rate,
    summarize_modifications,
    group_charge,
    aggregate_general_stats,
    sum_matching_dict_values
)
from pmultiqc.modules.common import ms_io
from pmultiqc.modules.common.ms import idxml as ms_idxml
from pmultiqc.modules.common.plots.ms import (
    draw_ms_information,
    draw_peak_intensity_distribution,
    draw_precursor_charge_distribution,
    draw_peaks_per_ms2,
)
from pmultiqc.modules.common.plots.id import (
    draw_potential_contaminants,
    draw_top_n_contaminants,
    draw_ids_rt_count,
    draw_delta_mass_da_ppm,
    draw_quantms_identification,
    draw_oversampling,
    draw_num_pep_per_protein,
    draw_charge_state,
    draw_summary_protein_ident_table,
    draw_quantms_identi_num
)
from pmultiqc.modules.common.plots.general import (
    draw_heatmap,
    plot_html_check,
    draw_exp_design,
    plot_data_check
)
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.stats import qual_uniform
from pmultiqc.modules.core.section_groups import (
    add_group_modules,
    add_sub_section
)
from pmultiqc.modules.common.logging import get_logger

log = get_logger("pmultiqc.modules.quantms")


class QuantMSModule:

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections
        self.heatmap_color_list = heatmap_colors

        self.exp_design_runs = None
        self.mzml_peak_distribution_plot = None
        self.mzml_charge_plot = None
        self.mzml_peaks_ms2_plot = None
        self.mzml_peak_distribution_plot_1 = None
        self.mzml_charge_plot_1 = None
        self.mzml_peaks_ms2_plot_1 = None
        self.msgf_label = None
        self.comet_label = None
        self.sage_label = None
        self.prot_search_score = None
        self.pep_plot = None
        self.peptide_search_score = None
        self.oversampling_plot = None
        self.psm_table_html = None
        self.protein_quantification_table_html = None
        self.pep_plot = None
        self.peptide_search_score = None
        self.ms_without_psm = None
        self.ms_with_psm = list()
        self.total_protein_identified = 0
        self.cal_num_table_data = dict()
        self.oversampling = dict()
        self.identified_spectrum = dict()
        self.total_ms2_spectra_identified = 0
        self.total_peptide_count = 0
        self.total_ms2_spectra = 0
        self.enable_dia = False
        self.is_bruker = False
        self.read_ms_info = False
        self.heatmap_charge_score = dict()
        self.missed_clevages_heatmap_score = dict()
        self.id_rt_score = dict()
        self.heatmap_over_sampling_score = dict()
        self.heatmap_pep_missing_score = dict()
        self.missed_cleavages_var_score = dict()
        self.pep_table_exists = False
        self.ms_info = dict()
        self.ms_info["charge_distribution"] = dict()
        self.ms_info["peaks_per_ms2"] = dict()
        self.ms_info["peak_distribution"] = dict()
        self.quantms_missed_cleavages = dict()
        self.quantms_modified = dict()
        self.identified_msms_spectra = dict()
        self.mztab_charge_state = dict()
        self.quantms_ids_over_rt = dict()
        self.quantms_pep_intensity = dict()
        self.quantms_contaminant_percent = dict()
        self.quantms_top_contaminant_percent = dict()
        self.quantms_mass_error = dict()

        self.enable_exp = False
        self.enable_sdrf = False
        self.msstats_input_valid = False

        self.psm_table = dict()
        self.mzml_peptide_map = dict()
        self.peptide_map_by_sample = dict()
        self.pep_quant_table = dict()
        self.mzml_table = OrderedDict()
        self.search_engine = OrderedDict()
        self.xcorr_hist_range = {"start": 0, "end": 5, "step": 0.1}
        self.hyper_hist_range = {"start": 0, "end": 5, "step": 0.1}
        self.spec_evalue_hist_range = {"start": 0, "end": 20, "step": 0.4}
        self.pep_hist_range = {
            "start": 0,
            "end": 1,
            "low_thresh": 0.3,
            "high_thresh": 0.7,
            "low_step": 0.03,
            "high_step": 0.08,
        }
        self.total_protein_quantified = 0
        self.out_csv_data = dict()
        self.ml_spec_ident_final = dict()
        self.heatmap_con_score = dict()
        self.heatmap_pep_intensity = {}
        self.ms1_tic = dict()
        self.ms1_bpc = dict()
        self.ms1_peaks = dict()
        self.ms1_general_stats = dict()
        self.is_multi_conditions = False
        self.sample_df = pd.DataFrame()
        self.file_df = pd.DataFrame()

    def get_data(self):

        log.info("Starting data recognition and processing...")

        if config.output_dir:
            if os.path.exists(config.output_dir):
                self.con = sqlite3.connect(os.path.join(config.output_dir, "quantms.db"))
            else:
                os.makedirs(config.output_dir)
                self.con = sqlite3.connect(os.path.join(config.output_dir, "quantms.db"))
        else:
            self.con = sqlite3.connect("./quantms.db")

        self.cur = self.con.cursor()
        self.cur.execute("drop table if exists PROTQUANT")
        self.con.commit()

        self.cur.execute("drop table if exists PEPQUANT")
        self.con.commit()

        self.cur.execute("drop table if exists PSM")
        self.con.commit()

        # TODO what if multiple are found??
        for f in self.find_log_files("pmultiqc/exp_design", filecontents=False):
            self.exp_design = os.path.join(f["root"], f["fn"])
            self.enable_exp = True

        if not self.enable_exp:

            for f in self.find_log_files("pmultiqc/sdrf", filecontents=False):
                self.sdrf = os.path.join(f["root"], f["fn"])

                parse_sdrf(
                    self.sdrf,
                    config.kwargs["keep_raw"],
                    config.kwargs["condition"],
                )

                # experimental_design.tsv is the default output name
                # experimental_design.tsv will be in the folder where pmultiqc is executed.
                self.exp_design = "experimental_design.tsv"
                self.enable_sdrf = True

        # draw the experimental design
        if self.enable_exp or self.enable_sdrf:
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

        # Please note that this section covers only the DIA part of quantms. For DIANN, refer to diann.py.
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

        if not self.enable_dia:
            for f in self.find_log_files("pmultiqc/mztab", filecontents=False):
                self.out_mztab_path = os.path.join(f["root"], f["fn"])
                self.parse_out_mztab()

        (
            self.mzml_table,
            self.mzml_peaks_ms2_plot,
            self.mzml_peak_distribution_plot,
            self.ms_info,
            self.total_ms2_spectra,
            _,  # mzml_ms_df
            self.heatmap_charge_score,
            self.mzml_charge_plot,
            self.ms1_tic,
            self.ms1_bpc,
            self.ms1_peaks,
            self.ms1_general_stats,
            self.current_sum_by_run
        ) = parse_mzml(
            is_bruker=self.is_bruker,
            read_ms_info=self.read_ms_info,
            ms_info_path=self.ms_info_path,
            ms_with_psm=self.ms_with_psm,
            identified_spectrum=self.identified_spectrum,
            enable_dia=self.enable_dia,
            ms_paths=self.ms_paths
        )

        self.idx_paths = []
        for idx_file in self.find_log_files("pmultiqc/idXML", filecontents=False):
            self.idx_paths.append(os.path.join(idx_file["root"], idx_file["fn"]))

        for msstats_input in self.find_log_files("pmultiqc/msstats", filecontents=False):
            self.msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])
            self.msstats_input_valid = True

        log.info("Data recognition and processing completed.")

        return True

    def draw_plots(self):

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

        # quantms: DIA
        if self.enable_dia:

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
                quantms_modified=self.quantms_modified,
                ms_paths=self.ms_paths,
                msstats_input_valid=self.msstats_input_valid
            )

            draw_summary_protein_ident_table(
                sub_sections=self.sub_sections["summary"],
                enable_dia=self.enable_dia,
                total_peptide_count=self.total_peptide_count,
                total_protein_quantified=self.total_protein_quantified
            )

            draw_quantms_identi_num(
                sub_sections=self.sub_sections["summary"],
                enable_exp=self.enable_exp,
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
        # quantms: LFQ or TMT
        else:

            if not config.kwargs["ignored_idxml"]:
                self.parse_idxml(self.mzml_table)
            self.cal_heat_map_score()

            heatmap_data, heatmap_xnames, heatmap_ynames = self.calculate_heatmap()
            draw_heatmap(
                self.sub_sections["summary"],
                self.heatmap_color_list,
                heatmap_data,
                heatmap_xnames,
                heatmap_ynames,
                False,
            )

            draw_summary_protein_ident_table(
                sub_sections=self.sub_sections["summary"],
                enable_dia=self.enable_dia,
                total_peptide_count=self.total_peptide_count,
                total_protein_quantified=self.total_protein_quantified,
                total_ms2_spectra_identified=self.total_ms2_spectra_identified,
                total_ms2_spectra=self.total_ms2_spectra,
                total_protein_identified=self.total_protein_identified
            )

            draw_quantms_identi_num(
                self.sub_sections["summary"],
                self.enable_exp,
                self.enable_sdrf,
                self.is_multi_conditions,
                self.sample_df,
                self.file_df,
                self.cal_num_table_data
            )

            draw_num_pep_per_protein(
                self.sub_sections["identification"],
                self.pep_plot
            )

            spectrum_tracking_data, spectrum_tracking_headers = aggregate_spectrum_tracking(
                mzml_table=self.mzml_table,
                peptide_map_by_sample=self.peptide_map_by_sample,
                sdrf_file_df=self.file_df
            )

            draw_mzml_ms(
                sub_section=self.sub_sections["ms2"],
                spectrum_tracking=spectrum_tracking_data,
                header_cols=spectrum_tracking_headers
            )

            if not config.kwargs["ignored_idxml"]:
                self.draw_search_engine()

            draw_precursor_charge_distribution(
                self.sub_sections["ms2"],
                charge_plot=self.mzml_charge_plot,
                ms_info=self.ms_info
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

            draw_oversampling(
                self.sub_sections["ms2"],
                self.oversampling,
                self.oversampling_plot.dict["cats"],
                False,
            )
            self.draw_delta_mass()

        msms_identified_rate = None
        if self.mzml_table and self.identified_msms_spectra:
            msms_identified_rate = aggregate_msms_identified_rate(
                self.mzml_table,
                self.identified_msms_spectra,
                self.file_df
            )

        draw_quantms_identification(
            self.sub_sections["identification"],
            cal_num_table_data=self.cal_num_table_data,
            quantms_missed_cleavages=self.quantms_missed_cleavages,
            quantms_modified=self.quantms_modified,
            msms_identified_rate=msms_identified_rate,
        )

        self.draw_quantms_contaminants()
        self.draw_quantms_quantification()
        self.draw_quantms_msms_section()
        self.draw_quantms_time_section()

        if self.msstats_input_valid:
            self.parse_msstats_input()

        if config.kwargs["quantification_method"] == "spectral_counting":
            # Add a report section with psm table plot from mzTab for spectral counting
            add_sub_section(
                sub_section=self.sub_sections["identification"],
                plot=self.psm_table_html,
                order=2,
                description="This plot shows the information of peptide spectrum matches",
                helptext="""
                        This table shows the information of peptide spectrum matches from mzTab PSM section.
                        """,
            )

        if self.enable_sdrf:
            ms_io.del_openms_convert_tsv()

        # TODO draw protein quantification from mzTab in the future with Protein and peptide tables from mzTab
        # currently only draw protein tabel for spectral counting
        if (
                not self.msstats_input_valid
                and config.kwargs["quantification_method"] == "spectral_counting"
        ):
            log.warning("MSstats input file not found!")
            add_sub_section(
                sub_section=self.sub_sections["quantification"],
                plot=self.protein_quantification_table_html,
                order=2,
                description="""
                    This plot shows the quantification information of proteins in the final result (mainly the mzTab file).
                    """,
                helptext="""
                    The quantification information (Spectral Counting) of proteins is obtained from the mzTab file. 
                    The table shows the quantitative level and distribution of proteins in different study variables and run.

                    * Peptides_Number: The number of peptides for each protein.
                    * Average Spectral Counting: Average spectral counting of each protein across all conditions with NA=0 or NA ignored.
                    * Spectral Counting in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Average spectral counting of replicates.

                    Click `Show replicates` to switch to bar plots of counting in each replicate.
                    """,
            )

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

        self.css = {
            "assets/css/quantms.css": os.path.join(
                os.path.dirname(__file__), "assets", "css", "quantms.css"
            )
        }
        self.js = {
            "assets/js/quantms.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "quantms.js"
            ),
            "assets/js/highcharts.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "highcharts.js"
            ),
            "assets/js/axios.min.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "axios.min.js"
            ),
            "assets/js/sql-optimized.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "sql-optimized.js"
            ),
        }

    def calculate_heatmap(self):

        heat_map_score = []
        if self.pep_table_exists:
            xnames = [
                "Contaminants",
                "Peptide Intensity",
                "Charge",
                "Missed Cleavages",
                "Missed Cleavages Var",
                "ID rate over RT",
                "MS2 OverSampling",
                "Pep Missing Values",
            ]
            ynames = []
            for k, v in self.heatmap_con_score.items():
                if k in self.ms_with_psm:
                    ynames.append(k)
                    heat_map_score.append(
                        [
                            v,
                            self.heatmap_pep_intensity[k],
                            self.heatmap_charge_score[k],
                            self.missed_clevages_heatmap_score[k],
                            self.missed_cleavages_var_score[k],
                            self.id_rt_score[k],
                            self.heatmap_over_sampling_score[k],
                            self.heatmap_pep_missing_score[k],
                        ]
                    )
        else:
            xnames = [
                "Charge",
                "Missed Cleavages",
                "Missed Cleavages Var",
                "ID rate over RT",
                "MS2 OverSampling",
                "Pep Missing Values",
            ]
            ynames = []
            for k, _ in self.heatmap_charge_score.items():
                if k in self.ms_with_psm:
                    ynames.append(k)
                    heat_map_score.append(
                        [
                            self.heatmap_charge_score[k],
                            self.missed_clevages_heatmap_score[k],
                            self.missed_cleavages_var_score[k],
                            self.id_rt_score[k],
                            self.heatmap_over_sampling_score[k],
                            self.heatmap_pep_missing_score[k],
                        ]
                    )
        return heat_map_score, xnames, ynames

    def draw_delta_mass(self):

        delta_mass = self.delta_mass
        delta_mass_percent = {
            "target": {
                k: v / sum(delta_mass["target"].values()) for k, v in delta_mass["target"].items()
            }
        }

        if delta_mass["decoy"]:

            delta_mass_percent["decoy"] = {
                k: v / sum(delta_mass["decoy"].values()) for k, v in delta_mass["decoy"].items()
            }

            x_values = list(delta_mass["target"].keys()) + list(delta_mass["decoy"].keys())

            range_threshold = 10
            if max(abs(x) for x in x_values) > range_threshold:
                range_abs = range_threshold
            else:
                range_abs = 1
            range_step = (max(x_values) - min(x_values)) * 0.05

            if max(abs(x) for x in x_values) > range_abs:

                delta_mass_range = dict()
                delta_mass_range["target"] = {
                    k: v for k, v in delta_mass["target"].items() if abs(k) <= range_abs
                }
                delta_mass_range["decoy"] = {
                    k: v for k, v in delta_mass["decoy"].items() if abs(k) <= range_abs
                }

                delta_mass_percent_range = dict()
                delta_mass_percent_range["target"] = {
                    k: v for k, v in delta_mass_percent["target"].items() if abs(k) <= range_abs
                }
                delta_mass_percent_range["decoy"] = {
                    k: v for k, v in delta_mass_percent["decoy"].items() if abs(k) <= range_abs
                }

                x_values_adj = list(delta_mass_range["target"].keys()) + list(
                    delta_mass_range["decoy"].keys()
                )
                range_step_adj = (max(x_values_adj) - min(x_values_adj)) * 0.05

                data_label = [
                    {
                        "name": f"Count (range: -{range_abs} to {range_abs})",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                        "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                        "ymin": 0,
                    },
                    {
                        "name": f"Relative Frequency (range: -{range_abs} to {range_abs})",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                        "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                        "ymin": 0,
                    },
                    {
                        "name": "Count (All Data)",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                    {
                        "name": "Relative Frequency (All Data)",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                ]
                pconfig = {
                    "id": "delta_mass",
                    "colors": {"target": "#b2df8a", "decoy": "#DC143C"},
                    "title": "Delta m/z",
                    "xlab": "Experimental m/z - Theoretical m/z",
                    "data_labels": data_label,
                    "style": "lines+markers",
                    "showlegend": True,
                    "save_data_file": False,
                }
                line_html = linegraph.plot(
                    [delta_mass_range, delta_mass_percent_range, delta_mass, delta_mass_percent],
                    pconfig,
                )

            else:
                data_label = [
                    {
                        "name": "Count (All Data)",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                    {
                        "name": "Relative Frequency (All Data)",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta elative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                ]
                pconfig = {
                    "id": "delta_mass",
                    "colors": {"target": "#b2df8a", "decoy": "#DC143C"},
                    "title": "Delta m/z",
                    "xlab": "Experimental m/z - Theoretical m/z",
                    "data_labels": data_label,
                    "style": "lines+markers",
                    "showlegend": True,
                    "save_data_file": False,
                }
                line_html = linegraph.plot([delta_mass, delta_mass_percent], pconfig)
        # no decoy
        else:
            delta_mass = {k: v for k, v in delta_mass.items() if k not in ["decoy"]}

            x_values = list(delta_mass["target"].keys())

            range_threshold = 10
            if max(abs(x) for x in x_values) > range_threshold:
                range_abs = range_threshold
            else:
                range_abs = 1
            range_step = (max(x_values) - min(x_values)) * 0.05

            if max(abs(x) for x in x_values) > range_abs:

                delta_mass_range = {
                    "target": {
                        k: v for k, v in delta_mass["target"].items() if abs(k) <= range_abs
                    }
                }

                delta_mass_percent_range = {
                    "target": {
                        k: v
                        for k, v in delta_mass_percent["target"].items()
                        if abs(k) <= range_abs
                    }
                }

                x_values_adj = list(delta_mass_range["target"].keys())
                range_step_adj = (max(x_values_adj) - min(x_values_adj)) * 0.05

                data_label = [
                    {
                        "name": f"Count (range: -{range_abs} to {range_abs})",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                        "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                        "ymin": 0,
                    },
                    {
                        "name": f"Relative Frequency (range: -{range_abs} to {range_abs})",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                        "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                        "ymin": 0,
                    },
                    {
                        "name": "Count (All Data)",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                    {
                        "name": "Relative Frequency (All Data)",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                ]
                pconfig = {
                    "id": "delta_mass_da",
                    "colors": {"target": "#b2df8a"},
                    "title": "Delta Mass [Da]",
                    "xlab": "Experimental m/z - Theoretical m/z",
                    "data_labels": data_label,
                    "style": "lines+markers",
                    "save_data_file": False,
                }
                line_html = linegraph.plot(
                    [delta_mass_range, delta_mass_percent_range, delta_mass, delta_mass_percent],
                    pconfig,
                )

            else:
                data_label = [
                    {
                        "name": "Count (All Data)",
                        "ylab": "Count",
                        "tt_label": "{point.x} Mass delta counts: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                    {
                        "name": "Relative Frequency (All Data)",
                        "ylab": "Relative Frequency",
                        "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                        "xmax": max(abs(x) for x in x_values) + range_step,
                        "xmin": -(max(abs(x) for x in x_values) + range_step),
                        "ymin": 0,
                    },
                ]
                pconfig = {
                    "id": "delta_mass_da",
                    "colors": {"target": "#b2df8a"},
                    "title": "Delta Mass [Da]",
                    "xlab": "Experimental m/z - Theoretical m/z",
                    "data_labels": data_label,
                    "style": "lines+markers",
                    "save_data_file": False,
                }
                line_html = linegraph.plot([delta_mass, delta_mass_percent], pconfig)

        line_html = plot_html_check(line_html)

        add_sub_section(
            sub_section=self.sub_sections["mass_error"],
            plot=line_html,
            order=1,
            description="""
                This chart represents the distribution of the relative frequency of experimental 
                precursor ion mass (m/z) - theoretical precursor ion mass (m/z).
                """,
            helptext="""
                Mass deltas close to zero reflect more accurate identifications and also 
                that the reporting of the amino acid modifications and charges have been done accurately. 
                This plot can highlight systematic bias if not centered on zero. 
                Other distributions can reflect modifications not being reported properly. 
                Also it is easy to see the different between the target and the decoys identifications.
                """,
        )

    def draw_search_engine(self):

        # Create scores summary plot
        [msgf_labels, comet_labels, sage_labels] = self.search_engine["data_label"]["score_label"]

        spec_e_pconfig = {
            "id": "spectral_e_values",  # ID used for the table
            "cpswitch": True,
            "title": "Summary of Spectral E-values",
            "xlab": "MSGF -lg(SpecEvalue) ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "",
            "tt_decimals": 0,
            "data_labels": msgf_labels,
            "save_data_file": False,
        }

        xcorr_pconfig = {
            "id": "cross_correlation_scores",  # ID used for the table
            "cpswitch": True,
            "title": "Summary of cross-correlation scores",
            "xlab": "Comet xcorr ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "",
            "tt_decimals": 0,
            "data_labels": comet_labels,
            "save_data_file": False,
        }

        hyper_pconfig = {
            "id": "summary_of_hyperscore",  # ID used for the table
            "cpswitch": True,
            "title": "Summary of Hyperscore",
            "xlab": "Sage hyperscore ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "",
            "tt_decimals": 0,
            "data_labels": sage_labels,
            "save_data_file": False,
        }

        bar_cats = OrderedDict()
        bar_cats["target"] = {"name": "target"}
        bar_cats["decoy"] = {"name": "decoy"}
        bar_cats["target+decoy"] = {"name": "target+decoy"}

        spec_e_cats = [bar_cats] * len(self.search_engine["SpecE"])
        xcorr_cats = [bar_cats] * len(self.search_engine["xcorr"])
        hyper_cats = [bar_cats] * len(self.search_engine["hyper"])
        pep_cats = [bar_cats] * len(self.search_engine["PEPs"])

        xcorr_bar_html = (
            bargraph.plot(list(self.search_engine["xcorr"].values()), xcorr_cats, xcorr_pconfig)
            if self.comet_label
            else ""
        )
        spec_e_bar_html = (
            bargraph.plot(list(self.search_engine["SpecE"].values()), spec_e_cats, spec_e_pconfig)
            if self.msgf_label
            else ""
        )
        hyper_bar_html = (
            bargraph.plot(list(self.search_engine["hyper"].values()), hyper_cats, hyper_pconfig)
            if self.sage_label
            else ""
        )

        if spec_e_bar_html != "":
            spec_e_bar_html = plot_html_check(spec_e_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=spec_e_bar_html,
                order=1,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. SpecEvalue: Spectral E-values, 
                        the search score of MSGF. The value used for plotting is -lg(SpecEvalue).
                        """,
            )

        if xcorr_bar_html != "":
            xcorr_bar_html = plot_html_check(xcorr_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=xcorr_bar_html,
                order=2,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. xcorr: cross-correlation scores, 
                        the search score of Comet. The value used for plotting is xcorr.
                        """,
            )

        if hyper_bar_html != "":
            hyper_bar_html = plot_html_check(hyper_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=hyper_bar_html,
                order=3,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. hyperscore: Hyperscore, the search 
                        score of Sage. The value used for plotting is hyperscore.
                        """,
            )

        # Create PEPs summary plot
        pep_pconfig = {
            "id": "search_engine_pep",  # ID used for the table
            "cpswitch": True,
            "title": "Summary of Search Engine PEP",
            "xlab": "PEP ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "",
            "tt_decimals": 0,
            "data_labels": self.search_engine["data_label"]["peps_label"],
            "save_data_file": False,
        }

        pep_bar_html = bargraph.plot(
            list(self.search_engine["PEPs"].values()), pep_cats, pep_pconfig
        )

        pep_bar_html = plot_html_check(pep_bar_html)

        add_sub_section(
            sub_section=self.sub_sections["search_engine"],
            plot=pep_bar_html,
            order=4,
            description="",
            helptext="This statistic is extracted from idXML files.",
        )

        # Create identified number plot
        if len(self.search_engine["data_label"]["consensus_label"]) != 0:
            consensus_pconfig = {
                "id": "consensus_summary",
                "cpswitch": True,
                "title": "Consensus Across Search Engines",
                "stacking": "normal",
                "height": 256,
                "tt_suffix": "",
                "tt_decimals": 0,
                "save_data_file": False,
            }

            consensus_bar_html = bargraph.plot(
                self.search_engine["consensus_support"],
                bar_cats,
                consensus_pconfig,
            )
            consensus_bar_html = plot_html_check(consensus_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=consensus_bar_html,
                order=5,
                description="",
                helptext="""
                    Consensus support is a measure of agreement between search engines. 
                    Every peptide sequence in the analysis has been identified by at least one search run. 
                    The consensus support defines which fraction (between 0 and 1) of the remaining search 
                    runs "supported" a peptide identification that was kept. 
                    The meaning of "support" differs slightly between algorithms: For best, worst, average 
                    and rank, each search run supports peptides that it has also identified among its top 
                    considered_hits candidates. So the "consensus support" simply gives the fraction of 
                    additional search engines that have identified a peptide. (For example, if there are 
                    three search runs, peptides identified by two of them will have a "support" of 0.5.) 
                    For the similarity-based algorithms PEPMatrix and PEPIons, the "support" for a peptide 
                    is the average similarity of the most-similar peptide from each (other) search run.
                    """,
            )

        # TODO
        # else:
        #     self.add_section(
        #         name="Summary of consensus PSMs",
        #         anchor="summary_of_consensus_PSMs",
        #         description="""#### Summary of consensus PSMs
        #         No Consensus PSMs data because of single search engine!
        #         """,
        #     )

    def cal_heat_map_score(self):
        timestamp = datetime.now().strftime("%H:%M:%S")
        log.info(f"{timestamp}: Calculating Heatmap Scores...")

        psm = self.mztab_data.spectrum_match_table
        meta_data = dict(self.mztab_data.metadata)
        if self.pep_table_exists:
            pep_table = self.mztab_data.peptide_table

            with pd.option_context("future.no_silent_downcasting", True):
                pep_table = pep_table.fillna(np.nan).infer_objects(copy=False).copy()

            study_variables = list(
                filter(
                    lambda x: re.match(r"peptide_abundance_study_variable.*?", x) is not None,
                    pep_table.columns.tolist(),
                )
            )

            pep_df_need_cols = ["accession", "opt_global_cv_MS:1002217_decoy_peptide", "spectra_ref"] + study_variables
            pep_table = pep_table[pep_df_need_cols].copy()

            pep_table.loc[:, "stand_spectra_ref"] = pep_table.apply(
                lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
                axis=1,
            )

            pep_table["average_intensity"] = pep_table[study_variables].mean(axis=1, skipna=True)

            # Contaminants
            if len(pep_table[pep_table["accession"].str.contains(config.kwargs["contaminant_affix"])]) > 0:

                self.quantms_contaminant_percent = self.cal_quantms_contaminant_percent(
                    pep_table[["average_intensity", "stand_spectra_ref", "accession"]].copy()
                )

                self.quantms_top_contaminant_percent = self.top_n_contaminant_percent(
                    pep_table[["average_intensity", "stand_spectra_ref", "accession"]].copy(), 5
                )
            else:
                log.warning(f"No contaminants found matching affix '{config.kwargs['contaminant_affix']}'")

            pep_intensity_by_run = dict()
            for name, group in pep_table.groupby("stand_spectra_ref"):

                contaminant_sum = (
                    group[group["accession"].str.contains(config.kwargs["contaminant_affix"])][
                        study_variables
                    ]
                    .sum(axis=0)
                    .sum()
                )
                all_sum = group[study_variables].sum(axis=0).sum()
                if all_sum == 0:
                    self.heatmap_con_score[name] = 0.0
                else:
                    self.heatmap_con_score[name] = 1.0 - (contaminant_sum / all_sum)

                if config.kwargs["remove_decoy"]:

                    pep_median = np.nanmedian(
                        group[(group["opt_global_cv_MS:1002217_decoy_peptide"] == 0)][
                            study_variables
                        ].to_numpy()
                    )

                    pep_intensity_by_run[name] = stat_pep_intensity(
                        group[
                            group["opt_global_cv_MS:1002217_decoy_peptide"] == 0
                        ]["average_intensity"]
                    )

                else:
                    pep_median = np.nanmedian(group[study_variables].to_numpy())

                    pep_intensity_by_run[name] = stat_pep_intensity(
                        group["average_intensity"]
                    )

                self.heatmap_pep_intensity[name] = np.minimum(
                    1.0, pep_median / (2 ** 23)
                )  # Threshold

            pep_table = pep_table.merge(
                right=self.file_df[["Sample", "Run"]].drop_duplicates(),
                left_on="stand_spectra_ref",
                right_on="Run"
            )

            pep_intensity_by_sample = dict()
            pep_table["Sample"] = pep_table["Sample"].astype(int)
            for name, group in pep_table.groupby("Sample", sort=True):

                if config.kwargs["remove_decoy"]:
                    pep_intensity_by_sample[f"Sample {str(name)}"] = stat_pep_intensity(
                        group[
                            group["opt_global_cv_MS:1002217_decoy_peptide"] == 0
                        ]["average_intensity"]
                    )

                else:
                    pep_intensity_by_sample[f"Sample {str(name)}"] = stat_pep_intensity(
                        group["average_intensity"]
                    )

            self.quantms_pep_intensity = [pep_intensity_by_run, pep_intensity_by_sample]

        #  HeatMapMissedCleavages
        global_peps = set(psm["opt_global_cv_MS:1000889_peptidoform_sequence"])
        global_peps_count = len(global_peps)
        if (
            config.kwargs["remove_decoy"]
            and "opt_global_cv_MS:1002217_decoy_peptide" in psm.columns
        ):
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] == 0].copy()

        psm_need_cols = ["spectra_ref", "opt_global_cv_MS:1000889_peptidoform_sequence", "sequence", "retention_time"]
        psm = psm[psm_need_cols].copy()

        psm.loc[:, "stand_spectra_ref"] = psm.apply(
            lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
            axis=1,
        )

        enzyme_list = [i for i in meta_data.values() if str(i).startswith("enzyme:")]
        enzyme = enzyme_list[0].split(":")[1] if len(enzyme_list) == 1 else "Trypsin"
        psm.loc[:, "missed_cleavages"] = psm.apply(
            lambda x: self.cal_miss_cleavages(x["sequence"], enzyme), axis=1
        )

        missed_cleavages_by_run = dict()

        # Calculate the ID RT Score
        for name, group in psm.groupby("stand_spectra_ref"):
            sc = group["missed_cleavages"].value_counts()

            missed_cleavages_by_run[name] = sc.to_dict()

            mis_0 = sc[0] if 0 in sc else 0
            self.missed_clevages_heatmap_score[name] = mis_0 / sc[:].sum()
            self.id_rt_score[name] = qual_uniform(group["retention_time"])

            #  For HeatMapOverSamplingScore
            self.heatmap_over_sampling_score[name] = self.oversampling[name]["1"] / np.sum(
                list(self.oversampling[name].values())
            )

            # For HeatMapPepMissingScore
            id_fraction = (
                    len(
                        set(group["opt_global_cv_MS:1000889_peptidoform_sequence"]).intersection(
                            global_peps
                        )
                    )
                    / global_peps_count
            )
            self.heatmap_pep_missing_score[name] = np.minimum(1.0, id_fraction)

        psm = psm.merge(
            right=self.file_df[["Sample", "Run"]].drop_duplicates(),
            left_on="stand_spectra_ref",
            right_on="Run"
        )

        psm["Sample"] = psm["Sample"].astype(int)
        missed_cleavages_by_sample = dict()
        for name, group in psm.groupby("Sample", sort=True):
            sc = group["missed_cleavages"].value_counts()
            missed_cleavages_by_sample[f"Sample {str(name)}"] = sc.to_dict()

        self.quantms_missed_cleavages = {
            "sdrf_samples": missed_cleavages_by_sample,
            "ms_runs": missed_cleavages_by_run,
        }

        median = np.median(list(self.missed_clevages_heatmap_score.values()))
        self.missed_cleavages_var_score = dict(
            zip(
                self.missed_clevages_heatmap_score.keys(),
                list(
                    map(
                        lambda v: 1 - np.abs(v - median),
                        self.missed_clevages_heatmap_score.values(),
                    )
                ),
            )
        )
        timestamp = datetime.now().strftime("%H:%M:%S")
        log.info(f"{timestamp}: Done calculating Heatmap Scores.")

    # if missed.cleavages is not given, it is assumed that Trypsin was used for digestion
    @staticmethod
    def cal_miss_cleavages(sequence, enzyme):
        if enzyme == "Trypsin/P":
            miss_cleavages = len(sequence[:-1]) - len(
                sequence[:-1].replace("K", "").replace("R", "").replace("P", "")
            )
        elif enzyme == "Arg-C":
            miss_cleavages = len(sequence[:-1]) - len(sequence[:-1].replace("R", ""))
        elif enzyme == "Asp-N":
            miss_cleavages = len(sequence[:-1]) - len(
                sequence[:-1].replace("B", "").replace("D", "")
            )
        elif enzyme == "Chymotrypsin":
            miss_cleavages = len(sequence[:-1]) - len(
                sequence[:-1].replace("F", "").replace("W", "").replace("Y", "").replace("L", "")
            )
        elif enzyme == "Lys-C":
            miss_cleavages = len(sequence[:-1]) - len(sequence[:-1].replace("K", ""))
        else:
            miss_cleavages = len(sequence[:-1]) - len(
                sequence[:-1].replace("K", "").replace("R", "")
            )
        return miss_cleavages

    def parse_idxml(self, mzml_table):
        # Instantiate the reader directly and parse
        reader = ms_idxml.IdXMLReader(
            file_paths=self.idx_paths,
            mzml_table=mzml_table,
            xcorr_hist_range=self.xcorr_hist_range,
            hyper_hist_range=self.hyper_hist_range,
            spec_evalue_hist_range=self.spec_evalue_hist_range,
            pep_hist_range=self.pep_hist_range,
            ml_spec_ident_final=self.ml_spec_ident_final,
            mzml_peptide_map=self.mzml_peptide_map,
            remove_decoy=config.kwargs["remove_decoy"],
        )
        reader.parse()
        self.search_engine = reader.search_engine
        self.msgf_label = reader.msgf_label
        self.comet_label = reader.comet_label
        self.sage_label = reader.sage_label

        # mass spectrum files sorted based on experimental file
        for spectrum_name in self.exp_design_runs:
            self.mzml_table[spectrum_name] = mzml_table[spectrum_name]

    def parse_out_mztab(self):

        from pmultiqc.modules.common.ms.mztab import MzTabReader

        mztab_reader = MzTabReader(file_path=self.out_mztab_path)
        mztab_reader.parse()
        self.mztab_data = mztab_reader.mztab_data
        psm = mztab_reader.psm
        pep_table = mztab_reader.pep_table
        self.delta_mass = mztab_reader.delta_mass
        prot = mztab_reader.prot
        prot_abundance_cols = mztab_reader.prot_abundance_cols
        meta_data = mztab_reader.meta_data
        self.ms_with_psm = mztab_reader.ms_with_psm
        self.total_protein_identified = mztab_reader.total_protein_identified
        self.total_protein_quantified = mztab_reader.total_protein_quantified

        self.pep_plot = Histogram("Number of peptides per proteins", plot_category="frequency")

        # There probably are no shared peptides in the final quant results. We still do it to be safe.
        # There are duplicates peptide-protein mapping in peptide table due to different feature (charge and RT)
        if config.kwargs["quantification_method"] == "spectral_counting":
            counts_per_acc = (
                psm.drop_duplicates("sequence")["accession"]
                .str.split(",")
                .explode()
                .value_counts()
            )
        else:
            self.pep_table_exists = True
            # TODO the following assumes that we always only look
            peptide_score = pep_table[
                ["opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]"]
            ]
            self.peptide_search_score = (
                peptide_score.groupby("opt_global_cv_MS:1000889_peptidoform_sequence")
                .agg("min")["best_search_engine_score[1]"]
                .to_dict()
            )
            counts_per_acc = (
                pep_table.drop_duplicates("sequence")["accession"]
                .str.split(",")
                .explode()
                .value_counts()
            )

        counts_per_acc.apply(self.pep_plot.add_value)
        # for c in counts_per_acc:
        #    self.pep_plot.addValue(c)
        categories = OrderedDict()
        categories["Frequency"] = {
            "name": "Frequency",
            "description": "Number of identified peptides per protein.",
        }
        self.pep_plot.to_dict(percentage=True, cats=categories)

        ml_spec_ident_final = {}

        # Modifications Name
        modifi_pattern = re.compile(r"UNIMOD:\d+")
        unimod_data = UnimodDatabase()

        def get_unimod_modification(modifis):
            if isinstance(modifis, str):
                modifi_matches = modifi_pattern.findall(modifis)
                mod_list = list()
                for mod in set(modifi_matches):
                    mod_list.append(unimod_data.get_by_accession(mod.upper()).get_name())
                return ",".join(set(mod_list))
            return "Unmodified"

        psm["Modifications"] = psm["modifications"].apply(get_unimod_modification)

        mod_plot_by_run = dict()
        modified_cats = list()

        data_per_run = dict()
        num_table_at_run = dict()

        if config.kwargs["remove_decoy"]:
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] == 0].copy()

        for m, group in psm.groupby("filename"):
            # m = os.path.basename(m)

            # Modifications
            mod_plot_dict, modified_cat = summarize_modifications(
                group[["sequence", "charge", "Modifications"]].drop_duplicates()
            )
            mod_plot_by_run[m] = mod_plot_dict
            modified_cats.extend(modified_cat)

            # Identified MS2 Spectra Raw File:
            self.identified_msms_spectra[m] = {"Identified": len(set(group["spectra_ref"]))}

            # Each loop resets the instance.
            self.oversampling_plot = Histogram(
                "MS/MS counts per 3D-peak", plot_category="frequency", breaks=[1, 2, 3]
            )

            with pd.option_context("future.no_silent_downcasting", True):
                group = group.fillna(pd.NA).infer_objects(copy=False)

            for i, j in group.groupby(["sequence", "charge", "modifications"]):
                self.oversampling_plot.add_value(len(j["spectra_ref"].unique()))

            self.oversampling_plot.to_dict()
            self.oversampling[m] = self.oversampling_plot.dict["data"]

            proteins = set(group["accession"])
            peptides = set(group["opt_global_cv_MS:1000889_peptidoform_sequence"])
            unique_peptides = set(
                group[group["unique"] == 1]["opt_global_cv_MS:1000889_peptidoform_sequence"]
            )

            self.identified_spectrum[m] = list(
                map(lambda x: x.split(":")[1], group["spectra_ref"])
            )
            self.mzml_peptide_map[m] = list(
                set(
                    group[group["opt_global_q-value"] <= 0.01]["sequence"].tolist()
                )
            )

            if None in proteins:
                proteins.remove(None)

            # TODO this is not really the number of proteins but the number of protein groups
            num_table_at_run[m] = {"protein_num": len(proteins)}
            num_table_at_run[m]["peptide_num"] = len(peptides)
            num_table_at_run[m]["unique_peptide_num"] = len(unique_peptides)

            modified_pep = list(
                filter(lambda x: re.match(r".*?\(.*\).*?", x) is not None, peptides)
            )
            num_table_at_run[m]["modified_peptide_num"] = len(modified_pep)

            if not self.file_df.empty:
                data_per_run[m] = {
                    "proteins": proteins,
                    "peptides": peptides,
                    "unique_peptides": unique_peptides,
                    "modified_peps": modified_pep
                }

            ml_spec_ident_final[m] = len(set(self.identified_spectrum[m]))

        num_table_at_sample = cal_num_table_at_sample(self.file_df, data_per_run)

        self.cal_num_table_data = {
            "sdrf_samples": num_table_at_sample,
            "ms_runs": num_table_at_run
        }

        # Pipeline Spectrum Tracking (Peptides quantified by Sample)
        self.peptide_map_by_sample = get_peptide_map_by_sample(
            peptide_map_by_run=self.mzml_peptide_map,
            sdrf_file_df=self.file_df
        )

        # Modifications
        mod_plot_by_sample = sample_level_modifications(
            df=psm,
            sdrf_file_df=self.file_df
        )

        self.quantms_modified["plot_data"] = [mod_plot_by_run, mod_plot_by_sample]
        self.quantms_modified["cats"] = list(
            sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x))
        )

        # Charge-state of Per File
        self.mztab_charge_state = cal_charge_state(
            psm=psm,
            sdrf_file_df=self.file_df
        )

        # IDs over RT
        quantms_rt_file_df = psm[["filename", "retention_time"]].copy()
        quantms_rt_file_df["retention time"] = quantms_rt_file_df["retention_time"] / 60
        quantms_rt_file_df.rename(columns={"filename": "raw file"}, inplace=True)
        self.quantms_ids_over_rt = evidence_rt_count(quantms_rt_file_df)

        # Delta Mass [ppm]
        mass_error = psm[["filename", "calc_mass_to_charge", "exp_mass_to_charge"]].copy()
        mass_error["mass error [ppm]"] = (
                                                 (mass_error["exp_mass_to_charge"] - mass_error["calc_mass_to_charge"])
                                                 / mass_error["calc_mass_to_charge"]
                                         ) * 1e6
        mass_error.rename(columns={"filename": "raw file"}, inplace=True)
        self.quantms_mass_error = evidence_calibrated_mass_error(mass_error)

        # TODO mzMLs without PSM: experimental design information is displayed, and all quantitative information is 0
        self.ms_without_psm = set([file_prefix(i) for i in self.ms_paths]) - set(self.ms_with_psm)
        for i in self.ms_without_psm:
            self.cal_num_table_data["ms_runs"][i] = {
                "protein_num": 0,
                "peptide_num": 0,
                "unique_peptide_num": 0,
                "modified_peptide_num": 0,
            }

        target_bin_data = {}
        decoy_bin_data = {}
        # TODO This is NOT a relative difference!
        psm["relative_diff"] = psm["exp_mass_to_charge"] - psm["calc_mass_to_charge"]
        try:
            decoy_bin = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] == 1][
                "relative_diff"
            ].value_counts(sort=False, bins=1000)
            for index in decoy_bin.index:
                decoy_bin_data[float(index.mid)] = int(decoy_bin[index])
            self.delta_mass["decoy"] = decoy_bin_data
        except (KeyError, IndexError, ValueError) as e:
            timestamp = datetime.now().strftime("%H:%M:%S")
            log.info(f"{timestamp}: No decoy peptides found -> only showing target peptides. Error: {e}")

        target_bin = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1][
            "relative_diff"
        ].value_counts(sort=False, bins=1000)
        for index in target_bin.index:
            target_bin_data[float(index.mid)] = int(target_bin[index])

        self.delta_mass["target"] = target_bin_data

        # extract delta mass
        self.ml_spec_ident_final = ml_spec_ident_final
        if config.kwargs["remove_decoy"]:
            self.total_ms2_spectra_identified = len(
                set(psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1]["spectra_ref"])
            )
            self.total_peptide_count = len(
                set(psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1]["sequence"])
            )
        else:
            self.total_ms2_spectra_identified = len(set(psm["spectra_ref"]))
            self.total_peptide_count = len(set(psm["sequence"]))

        # draw PSMs table for spectral counting
        if config.kwargs["quantification_method"] == "spectral_counting" and not config.kwargs.get(
                "disable_table", True
        ):
            mztab_data_psm_full = psm[
                [
                    "opt_global_cv_MS:1000889_peptidoform_sequence",
                    "accession",
                    "search_engine_score[1]",
                    "stand_spectra_ref",
                ]
            ].copy()
            mztab_data_psm_full.rename(
                columns={
                    "opt_global_cv_MS:1000889_peptidoform_sequence": "Sequence",
                    "accession": "Accession",
                    "search_engine_score[1]": "Search_Engine_Score",
                    "stand_spectra_ref": "Spectra_Ref",
                },
                inplace=True,
            )
            mztab_data_psm_full[["Sequence", "Modification"]] = mztab_data_psm_full.apply(
                lambda x: find_modification(x["Sequence"]), axis=1, result_type="expand"
            )
            max_search_score = mztab_data_psm_full["Search_Engine_Score"].max()
            mztab_data_psm_full = mztab_data_psm_full.to_dict("index")
            headers = OrderedDict()
            headers["Sequence"] = {"name": "Sequence", "description": "Peptide Sequence"}
            headers["Modification"] = {
                "name": "Modification",
                "description": "Modification in Peptide Sequence",
            }
            headers["Accession"] = {"name": "Accession", "description": "Protein Name"}
            headers["Search_Engine_Score"] = {
                "name": "Search Engine Score",
                "format": "{:,.5e}",
                "max": max_search_score,
                "scale": False,
            }

            # upload PSMs table to sqlite database
            self.cur.execute(
                "CREATE TABLE PSM(PSM_ID INT(200), Sequence VARCHAR(200), Modification VARCHAR(100), Accession VARCHAR(100), Search_Engine_Score FLOAT(4,5), Spectra_Ref VARCHAR(100))"
            )
            self.con.commit()
            sql_col = "PSM_ID,Sequence,Modification,Accession,Search_Engine_Score,Spectra_Ref"
            sql_t = "(" + ",".join(["?"] * 6) + ")"

            # PSM_ID is index
            all_term = [
                "Sequence",
                "Modification",
                "Accession",
                "Search_Engine_Score",
                "Spectra_Ref",
            ]
            self.cur.executemany(
                "INSERT INTO PSM (" + sql_col + ") VALUES " + sql_t,
                [(k, *itemgetter(*all_term)(v)) for k, v in mztab_data_psm_full.items()],
            )
            self.con.commit()

            pconfig = {
                "id": "peptide spectrum matches",  # ID used for the table
                "table_title": "information of peptide spectrum matches",   # Title of the table. Used in the column config modal
                "sortRows": False,  # Whether to sort rows alphabetically
                "only_defined_headers": False,  # Only show columns that are defined in the headers config
                "col1_header": "PSM_ID",
                "format": "{:,.0f}",
                "no_violin": True,
                "save_data_file": False,
            }

            mztab_data_psm_init = dict(itertools.islice(mztab_data_psm_full.items(), 50))
            table_html = table.plot(mztab_data_psm_init, headers, pconfig)
            pattern = re.compile(r'<small id="peptide_spectrum_matches_numrows_text"')
            match = re.search(pattern, table_html)
            if match is None:
                log.warning("Could not find expected pattern in table HTML, using default insertion point")
                index = len(table_html)
            else:
                index = match.span()[0]
            options = "".join(f"<option>{key}</option>" for key in ["Sequence", "Modification", "Accession", "Spectra_Ref"])
            t_html = (
                    table_html[:index]
                    + '<input type="text" placeholder="search..." class="searchInput" '
                      'onkeyup="searchPsmFunction()" id="psm_search">'
                      f'<select name="psm_search_col" id="psm_search_col">{options}</select>'
            )
            table_html = (
                    t_html + "</select>" + '<button type="button" class="btn btn-default '
                                           'btn-sm" id="psm_reset" onclick="psmFirst()">Reset</button>' + table_html[
                        index:]
            )
            table_html = (
                    table_html
                    + """<div class="page_control"><span id="psmFirst">First Page</span><span
            id="psmPre"> Previous Page</span><span id="psmNext">Next Page </span><span id="psmLast">Last 
            Page</span><span id="psmPageNum"></span>Page/Total <span id="psmTotalPage"></span>Pages <input 
            type="number" name="" id="psm_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
            onkeydown="psm_page_jump()" min="1"/> </div> """
            )

            self.psm_table_html = table_html

        # TODO implement the second option no msstats and feature intensity: draw protein quantification from mzTab
        # in the future with Protein and peptide tables from mzTab.
        # Draw protein table with spectral counting from mzTab file
        if (
                not self.msstats_input_valid
                and config.kwargs["quantification_method"] == "spectral_counting"
                and not config.kwargs.get("disable_table", True)
        ):
            mztab_data_dict_prot_full = dict()
            conditions = self.sample_df.drop_duplicates(subset="MSstats_Condition")[
                "MSstats_Condition"
            ].tolist()

            def get_spectrum_count_across_rep(condition_count_dict: dict):
                spc = []
                res = copy.deepcopy(condition_count_dict)
                for c, val in condition_count_dict.items():
                    samples_spc = dict()
                    # Average spectrum counting with NA=0 ignored with replicates
                    for sn, count_value in val.items():
                        if len(np.nonzero(count_value)[0]) == 0:
                            samples_spc[sn] = 0.0
                        else:
                            samples_spc[sn] = sum(count_value) / len(np.nonzero(count_value)[0])
                    if len(np.nonzero(list(samples_spc.values()))[0]) == 0:
                        res[c] = 0
                    else:
                        res[c] = round(
                            sum(list(samples_spc.values()))
                            / len(np.nonzero(list(samples_spc.values()))[0])
                        )
                    samples_spc = dict(filter(lambda x: x[1] != 0.0, samples_spc.items()))
                    res[c + "_distribution"] = str(samples_spc).replace("'", '"')
                    spc.append(res[c])

                # Integer for average spectrum counting with NA=0 ignored across condition
                res["Average Spectrum Counting"] = round(sum(spc) / len(np.nonzero(spc)[0]))
                return res

            for index, row in prot.iterrows():
                mztab_data_dict_prot_full[index] = {}
                for abundance_col in prot_abundance_cols:
                    # map abundance assay to factor value
                    file_name = os.path.basename(
                        meta_data[
                            meta_data[
                                abundance_col.replace("protein_abundance_", "") + "-ms_run_ref"
                                ].split(",")[0]
                            + "-location"
                            ]
                    )
                    sample_name = str(
                        self.file_df[self.file_df["Run"] == os.path.splitext(file_name)[0]][
                            "Sample"
                        ].values[0]
                    )
                    condition = str(
                        self.sample_df[self.sample_df["Sample"] == sample_name][
                            "MSstats_Condition"
                        ].values[0]
                    )

                    # Consider technical replicates and biological replicates
                    if condition in mztab_data_dict_prot_full[index]:
                        if sample_name in mztab_data_dict_prot_full[index][condition]:
                            mztab_data_dict_prot_full[index][condition][sample_name].append(
                                row[abundance_col]
                            )
                        else:
                            mztab_data_dict_prot_full[index][condition] = {
                                sample_name: [row[abundance_col]]
                            }
                    else:
                        mztab_data_dict_prot_full[index][condition] = {
                            sample_name: [row[abundance_col]]
                        }

                mztab_data_dict_prot_full[index] = get_spectrum_count_across_rep(
                    mztab_data_dict_prot_full[index]
                )
                mztab_data_dict_prot_full[index]["Peptides_Number"] = int(counts_per_acc[index])

            timestamp = datetime.now().strftime("%H:%M:%S")
            log.info(f"{timestamp}: Done aggregating mzTab file {self.out_mztab_path}...")

            headers = OrderedDict()
            headers["Peptides_Number"] = {
                "name": "Number of Peptides",
                "description": "Number of peptides per proteins",
                "format": "{:,.0f}",
            }
            headers["Average Spectrum Counting"] = {
                "name": "Average Spectrum Counting",
                "description": "Average Spectrum Counting across all conditions",
                "format": "{:,.0f}",
            }

            # upload protein table to sqlite database
            self.cur.execute(
                'CREATE TABLE PROTQUANT(ProteinName VARCHAR(100), Peptides_Number INT(100), "Average Spectrum Counting" VARCHAR)'
            )
            self.con.commit()
            sql_col = 'ProteinName,Peptides_Number,"Average Spectrum Counting"'
            sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 3)) + ")"

            for s in conditions:
                self.cur.execute('ALTER TABLE PROTQUANT ADD "' + str(s) + '" VARCHAR')
                self.con.commit()
                sql_col += ', "' + str(s) + '"'
                headers[str(s)] = {"name": s}

            for s in list(map(lambda x: str(x) + "_distribution", conditions)):
                self.cur.execute('ALTER TABLE PROTQUANT ADD "' + s + '" VARCHAR(100)')
                self.con.commit()
                sql_col += ', "' + s + '"'
                headers[str(s)] = {"name": s}

            # ProteinName is index
            all_term = (
                    ["Peptides_Number", "Average Spectrum Counting"]
                    + list(map(str, conditions))
                    + list(map(lambda x: str(x) + "_distribution", conditions))
            )
            self.cur.executemany(
                "INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
                [(k, *itemgetter(*all_term)(v)) for k, v in mztab_data_dict_prot_full.items()],
            )
            self.con.commit()

            pconfig = {
                "id": "quantification_of_protein",  # ID used for the table
                "title": "quantification information of protein",
                "anchor": "",   # Title of the table. Used in the column config modal
                "save_file": False,  # Whether to save the table data to a file
                "raw_data_fn": "multiqc_quantification_of_protein_table",  # File basename to use for raw data file
                "sort_rows": False,  # Whether to sort rows alphabetically
                "only_defined_headers": False,  # Only show columns that are defined in the headers config
                "col1_header": "ProteinName",
                "no_violin": True,
                "save_data_file": False,
            }

            max_prot_intensity = 0
            mztab_data_dict_prot_init = dict(
                itertools.islice(mztab_data_dict_prot_full.items(), 50)
            )

            table_html = sparklines.plot(
                mztab_data_dict_prot_init, headers, pconfig=pconfig, max_value=max_prot_intensity
            )
            pattern = re.compile(r'<small id="quantification_of_protein_numrows_text"')
            match = re.search(pattern, table_html)
            if match is None:
                log.warning("Could not find expected pattern in protein table HTML, using default insertion point")
                index = len(table_html)
            else:
                index = match.span()[0]
            options = "".join(f"<option>{key}</option>" for key in ["ProteinName"])
            t_html = (
                    table_html[:index]
                    + '<input type="text" placeholder="search..." class="searchInput" '
                      'onkeyup="searchProtFunction()" id="prot_search">'
                      f'<select name="prot_search_col" id="prot_search_col">{options}</select>'
            )
            table_html = (
                    t_html + "</select>" + '<button type="button" class="btn btn-default '
                                           'btn-sm" id="prot_reset" onclick="protFirst()">Reset</button>' + table_html[
                        index:]
            )
            table_html = (
                    table_html
                    + """<div class="page_control"><span id="protFirst">First Page</span><span
            id="protPre"> Previous Page</span><span id="protNext">Next Page </span><span id="protLast">Last 
            Page</span><span id="protPageNum"></span>Page/Total <span id="protTotalPage"></span>Pages <input 
            type="number" name="" id="prot_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
            onkeydown="prot_page_jump()" min="1"/> </div> """
            )

            self.protein_quantification_table_html = table_html

    def parse_msstats_input(self):
        timestamp = datetime.now().strftime("%H:%M:%S")
        log.info(f"{timestamp}: Parsing MSstats input file {self.msstats_input_path}...")
        msstats_data = pd.read_csv(self.msstats_input_path)
        # TODO we probably shouldn't even write out 0-intensity values to MSstats csv
        msstats_data = msstats_data[msstats_data["Intensity"] != 0]
        msstats_data.loc[:, "BestSearchScore"] = 1 - msstats_data.loc[:, "PeptideSequence"].map(
            self.peptide_search_score
        )
        msstats_data[["PeptideSequence", "Modification"]] = msstats_data.apply(
            lambda x: find_modification(x["PeptideSequence"]), axis=1, result_type="expand"
        )

        reps_per_condition = (
            self.sample_df.groupby("MSstats_Condition")["MSstats_BioReplicate"].agg(list).to_dict()
        )
        conditions = list(self.sample_df["MSstats_Condition"].unique())
        conditions_str = [str(c) for c in conditions]
        conditions_dists = [str(c) + "_distribution" for c in conditions]

        # TODO maybe aggregating in dicts is not the fastest. We also need to parse them again for proteins later.
        #  Maybe we should just introduce new pandas columns for every bioreplicate.
        def fill_dict(g):
            d = dict.fromkeys(reps_per_condition[str(g.name)], 0)
            d.update(zip(g["BioReplicate"].astype(str), np.log10(g["Intensity"])))
            return json.dumps(d)

        def get_inty_across_bio_reps_as_str(g):
            gdict = dict.fromkeys(conditions_str, 0.0)
            gdict.update(dict.fromkeys(conditions_dists, "{}"))
            gdict["Average Intensity"] = np.log10(g["Intensity"].mean())
            gdict["BestSearchScore"] = g["BestSearchScore"].min()
            # TODO How to determine technical replicates? Should be same BioReplicate but different Fraction_Group (but fraction group is not annotated)
            grouped = g.groupby(["Condition", "BioReplicate"], as_index=False)["Intensity"].mean()
            cond_grp = grouped.groupby("Condition", group_keys=False)[
                ["BioReplicate", "Intensity"]
            ].apply(fill_dict)

            cond_grp.index = [str(c) + "_distribution" for c in cond_grp.index]
            gdict.update(cond_grp.to_dict())
            mean = g.groupby(["Condition"])["Intensity"].mean()
            cond_grp_mean = np.log10(mean)
            cond_grp_mean.index = cond_grp_mean.index.map(str)
            gdict.update(cond_grp_mean.to_dict())
            return pd.Series(gdict)

        msstats_data_pep_agg = msstats_data.groupby(
            ["PeptideSequence", "ProteinName", "Modification"]
        )[["Intensity", "BestSearchScore", "Condition", "BioReplicate"]].apply(
            get_inty_across_bio_reps_as_str
        )

        # TODO Can we guarantee that the score was always PEP? I don't think so!
        msstats_data_pep_agg.reset_index(inplace=True)
        msstats_data_pep_agg.index = msstats_data_pep_agg.index + 1
        msstats_data_dict_pep_full = msstats_data_pep_agg.to_dict("index")

        self.cur.execute(
            'CREATE TABLE PEPQUANT(PeptideID INT(100) PRIMARY KEY, PeptideSequence VARCHAR(100), Modification VARCHAR(100), ProteinName VARCHAR(100), BestSearchScore FLOAT(4,3), "Average Intensity" FLOAT(4,3))'
        )
        self.con.commit()
        sql_col = 'PeptideID,PeptideSequence,Modification,ProteinName,BestSearchScore, "Average Intensity"'
        sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 6)) + ")"

        headers = {
            "ProteinName": {
                "title": "Protein Name",
                "description": "Name/Identifier(s) of the protein (group)",
                "minrange": "200",
            },
            "PeptideSequence": {"title": "Peptide Sequence"},
            "BestSearchScore": {"title": "Best Search Score", "format": "{:,.4f}"},
            "Average Intensity": {
                "title": "Average Intensity",
                "description": "Average intensity across all conditions",
                "format": "{:,.4f}",
            },
        }

        for s in conditions:
            self.cur.execute('ALTER TABLE PEPQUANT ADD "' + str(s) + '" VARCHAR')
            self.con.commit()
            sql_col += ', "' + str(s) + '"'
            headers[str(s)] = {"title": s, "format": "{:,.4f}"}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            self.cur.execute('ALTER TABLE PEPQUANT ADD "' + s + '" VARCHAR(100)')
            self.con.commit()
            sql_col += ', "' + s + '"'

        # PeptideID is index
        all_term = (
                [
                    "PeptideSequence",
                    "Modification",
                    "ProteinName",
                    "BestSearchScore",
                    "Average Intensity",
                ]
                + list(map(str, conditions))
                + list(map(lambda x: str(x) + "_distribution", conditions))
        )
        self.cur.executemany(
            "INSERT INTO PEPQUANT (" + sql_col + ") VALUES " + sql_t,
            [(k, *itemgetter(*all_term)(v)) for k, v in msstats_data_dict_pep_full.items()],
        )
        self.con.commit()

        draw_config = {
            "namespace": "",
            "id": "peptides_quantification_table",
            "title": "Peptides Quantification Table",
            "sort_rows": False,
            "only_defined_headers": True,
            "col1_header": "PeptideID",
            "no_violin": True,
            "save_data_file": False,
        }

        # only use the first 50 lines for the table
        max_pep_intensity = 50
        table_html = table.plot(
            dict(itertools.islice(msstats_data_dict_pep_full.items(), max_pep_intensity)),
            headers=headers,
            pconfig=draw_config,
        )

        add_sub_section(
            sub_section=self.sub_sections["quantification"],
            plot=table_html,
            order=1,
            description="""
                This plot shows the quantification information of peptides in the final result (mainly the mzTab file).
                """,
            helptext="""
                The quantification information of peptides is obtained from the MSstats input file. 
                The table shows the quantitative level and distribution of peptides in different study variables, 
                run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
                the average intensity for all the fractions, runs and peptiforms.

                * BestSearchScore: It is equal to 1 - min(Q.Value) for DIA datasets. Then it is equal to 
                1 - min(best_search_engine_score[1]), which is from best_search_engine_score[1] column in mzTab 
                peptide table for DDA datasets.
                * Average Intensity: Average intensity of each peptide sequence across all conditions with NA=0 or NA ignored.
                * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of fractions, 
                and then mean intensity in technical replicates/biological replicates separately.
                """,
        )

        # Helper functions for pandas
        def json_to_dict(s):
            if isinstance(s, str):
                return json.loads(s)
            else:
                return {}

        def reducer(accumulator, element):
            for key, value in json_to_dict(element).items():
                accumulator[key] = np.log10(pow(10, accumulator.get(key, 0)) + pow(10, value))
            return accumulator

        def my_dict_sum(series):
            return json.dumps(reduce(reducer, series, {}))

        def total_intensity(series):
            total = 0.0
            for intensity in series:
                total += pow(10, intensity)
            return np.log10(total)

        def unique_count(series):
            return len(series.unique())

        agg_funs = dict.fromkeys(conditions_dists, my_dict_sum)
        agg_funs.update(dict.fromkeys(conditions_str, total_intensity))
        agg_funs["PeptideSequence"] = unique_count
        agg_funs["Average Intensity"] = total_intensity
        msstats_data_prot = msstats_data_pep_agg.groupby("ProteinName").agg(
            agg_funs
        )  # .reset_index()
        del msstats_data_pep_agg
        msstats_data_prot.rename(columns={"PeptideSequence": "Peptides_Number"}, inplace=True)
        msstats_data_prot.reset_index(inplace=True)
        msstats_data_prot.index = msstats_data_prot.index + 1
        msstats_data_dict_prot_full = msstats_data_prot.to_dict("index")

        msstats_data_dict_prot_init = dict(
            itertools.islice(msstats_data_dict_prot_full.items(), 50)
        )

        headers = {
            "ProteinName": {
                "title": "Protein Name",
                "description": "Name/Identifier(s) of the protein (group)",
            },
            "Peptides_Number": {
                "title": "Number of Peptides",
                "description": "Number of peptides per proteins",
                "format": "{:,.0f}",
            },
            "Average Intensity": {
                "title": "Average Intensity",
                "description": "Average intensity across all conditions",
                "format": "{:,.4f}",
            },
        }

        # upload protein table to sqlite database
        self.cur.execute(
            'CREATE TABLE PROTQUANT(ProteinID INT(100), ProteinName VARCHAR(100), Peptides_Number INT(100), "Average Intensity" VARCHAR)'
        )
        self.con.commit()
        sql_col = 'ProteinID,ProteinName,Peptides_Number,"Average Intensity"'
        sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 4)) + ")"

        for s in conditions:
            self.cur.execute('ALTER TABLE PROTQUANT ADD "' + str(s) + '" VARCHAR')
            self.con.commit()
            sql_col += ', "' + str(s) + '"'
            headers[str(s)] = {"title": s, "format": "{:,.4f}"}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            self.cur.execute('ALTER TABLE PROTQUANT ADD "' + s + '" VARCHAR(100)')
            self.con.commit()
            sql_col += ', "' + s + '"'

        # ProteinID is index
        all_term = (
                ["ProteinName", "Peptides_Number", "Average Intensity"]
                + list(map(str, conditions))
                + list(map(lambda x: str(x) + "_distribution", conditions))
        )
        self.cur.executemany(
            "INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
            [(k, *itemgetter(*all_term)(v)) for k, v in msstats_data_dict_prot_full.items()],
        )
        self.con.commit()

        draw_config = {
            "namespace": "",
            "id": "protein_quant_result",
            "title": "Protein Quantification Table",
            "save_file": False,
            "sort_rows": False,
            "only_defined_headers": True,
            "col1_header": "ProteinID",
            "no_violin": True,
            "save_data_file": False,
        }
        table_html = table.plot(msstats_data_dict_prot_init, headers=headers, pconfig=draw_config)
        add_sub_section(
            sub_section=self.sub_sections["quantification"],
            plot=table_html,
            order=2,
            description="""
                This plot shows the quantification information of proteins in the final result (mainly the mzTab file).
                """,
            helptext="""
                The quantification information of proteins is obtained from the msstats input file. 
                The table shows the quantitative level and distribution of proteins in different study variables and run.

                * Peptides_Number: The number of peptides for each protein.
                * Average Intensity: Average intensity of each protein across all conditions with NA=0 or NA ignored.
                * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
                """,
        )

    def cal_quantms_contaminant_percent(self, pep_df):

        pep_df["is_contaminant"] = pep_df["accession"].str.contains(config.kwargs["contaminant_affix"], na=False)
        group_stats = pep_df.groupby("stand_spectra_ref").agg(
            total_intensity=("average_intensity", "sum"),
            cont_intensity=("average_intensity", lambda x: x[pep_df.loc[x.index, "is_contaminant"]].sum()),
        )

        group_stats["contaminant_percent"] = (
                group_stats["cont_intensity"] / group_stats["total_intensity"] * 100
        )

        result_dict = dict()
        for k, v in dict(zip(group_stats.index, group_stats["contaminant_percent"])).items():
            result_dict[k] = {"Potential Contaminants": v}

        return result_dict

    def top_n_contaminant_percent(self, pep_df, top_n):

        not_cont_tag = "NOT_CONTAM"
        is_contaminant = pep_df["accession"].str.contains(config.kwargs["contaminant_affix"], na=False)
        pep_df.loc[:, "cont_accession"] = np.where(is_contaminant, pep_df["accession"], not_cont_tag)

        pep_contaminant_df = pep_df[pep_df["cont_accession"] != not_cont_tag].copy()
        contaminant_df = (
            pep_contaminant_df.groupby("cont_accession", as_index=False)["average_intensity"]
            .sum()
            .sort_values(by="average_intensity", ascending=False)
        )

        top_contaminants = list(contaminant_df.head(top_n).cont_accession)

        plot_dict = dict()
        plot_cats = list()

        for file_name, group in pep_df.groupby("stand_spectra_ref"):
            contaminant_rows = group[group["cont_accession"] != not_cont_tag].copy()
            contaminant_rows.loc[
                ~contaminant_rows["cont_accession"].isin(top_contaminants), "cont_accession"
            ] = "Other"

            cont_df = (
                contaminant_rows.groupby("cont_accession", as_index=False)["average_intensity"]
                .sum()
                .sort_values(by="average_intensity", ascending=False)
                .reset_index(drop=True)
            )
            cont_df["contaminant_percent"] = (
                                                     cont_df["average_intensity"] / group["average_intensity"].sum()
                                             ) * 100

            plot_dict[file_name] = dict(
                zip(cont_df["cont_accession"], cont_df["contaminant_percent"])
            )
            plot_cats.extend(cont_df["cont_accession"].tolist())

        plot_dict = {k: v for k, v in plot_dict.items() if v}

        if not plot_dict:
            return None

        plot_cats = list(set(plot_cats))
        if "Other" in plot_cats:
            plot_cats = [x for x in plot_cats if x != "Other"] + ["Other"]

        result_dict = dict()
        result_dict["plot_data"] = plot_dict
        result_dict["cats"] = plot_cats

        return result_dict

    def draw_quantms_contaminants(self):

        # 1.Potential Contaminants per Group
        if self.quantms_contaminant_percent:
            draw_potential_contaminants(
                self.sub_sections["contaminants"], self.quantms_contaminant_percent, False
            )

        # 2.Top5 Contaminants per Raw file
        if self.quantms_top_contaminant_percent:
            draw_top_n_contaminants(
                self.sub_sections["contaminants"], self.quantms_top_contaminant_percent
            )

    def draw_quantms_quantification(self):

        # 1.Peptide Intensity Distribution
        if self.quantms_pep_intensity:
            draw_config = {
                "id": "peptide_intensity_distribution_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "Peptide Intensity Distribution",
                "tt_decimals": 2,
                "xlab": "log2(Intensity)",
                "data_labels": ["by Run", "by Sample"],
                "sort_samples": False,
                "save_data_file": False,
            }
            box_html = box.plot(self.quantms_pep_intensity, pconfig=draw_config)

            # box_html.flat
            box_html = plot_data_check(
                plot_data=self.quantms_pep_intensity,
                plot_html=box_html,
                log_text="pmultiqc.modules.quantms.quantms",
                function_name="draw_quantms_quantification"
            )
            box_html = plot_html_check(box_html)

            add_sub_section(
                sub_section=self.sub_sections["quantification"],
                plot=box_html,
                order=5,
                description="Peptide intensity per file from mzTab.",
                helptext="""
                    Calculate the average of peptide_abundance_study_variable[1-n] values for each peptide from the 
                    peptide table in the mzTab file, and then apply a log2 transformation.
                    """,
            )

    def draw_quantms_msms_section(self):

        # 1.Charge-state of Per File
        if self.mztab_charge_state:
            draw_charge_state(self.sub_sections["ms2"], self.mztab_charge_state, "")

    def draw_quantms_time_section(self):

        # 1.IDs over RT
        if self.quantms_ids_over_rt:
            draw_ids_rt_count(
                self.sub_sections["rt_qc"], self.quantms_ids_over_rt, ""
            )

        # 2.Delta Mass [ppm]
        if self.quantms_mass_error:
            draw_delta_mass_da_ppm(
                self.sub_sections["mass_error"], self.quantms_mass_error, "quantms_ppm"
            )

    def __del__(self):
        """Cleanup method to close SQLite connection."""
        if hasattr(self, 'con') and self.con:
            try:
                self.con.close()
            except Exception:
                pass  # Ignore errors during cleanup


def draw_mzml_ms(sub_section, spectrum_tracking, header_cols):

    pconfig = {
        "id": "pipeline_spectrum_tracking",  # ID used for the table
        "title": "Pipeline Spectrum Tracking",  # Title of the table. Used in the column config modal
        "save_file": False,  # Whether to save the table data to a file
        "raw_data_fn": "multiqc_spectrum_tracking_table",  # File basename to use for raw data file
        "save_data_file": False,
    }

    headers = OrderedDict()
    headers["MS1_Num"] = {
        "title": "#MS1 Spectra",
        "description": "Number of MS1 spectra",
    }
    headers["MS2_Num"] = {
        "title": "#MS2 Spectra",
        "description": "Number of MS2 spectra",
    }

    if "MSGF" in header_cols:
        headers["MSGF"] = {
            "description": "Number of spectra identified by MSGF search engine",
        }
    if "Comet" in header_cols:
        headers["Comet"] = {
            "description": "Number of spectra identified by Comet search engine",
        }
    if "Sage" in header_cols:
        headers["Sage"] = {
            "description": "Number of spectra identified by Sage search engine",
        }
    headers["num_quant_psms"] = {
        "title": "#PSMs from quant. peptides",
        "description": "Number of reliable PSMs from peptides IDs used in quantification",
    }
    headers["num_quant_peps"] = {
        "title": "#Peptides quantified",
        "description": "Number of quantified peptides that passed final protein and peptide FDR thresholds.",
    }
    table_html = table.plot(spectrum_tracking, headers, pconfig)

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=4,
        description="This plot shows the tracking of the number of spectra along the quantms pipeline",
        helptext="""
            This table shows the changes in the number of spectra corresponding to each input file 
            during the pipeline operation. And the number of peptides finally identified and quantified is obtained from 
            the PSM table in the mzTab file. You can also remove decoys with the `remove_decoy` parameter.:

            * MS1_Num: The number of MS1 spectra extracted from mzMLs
            * MS2_Num: The number of MS2 spectra extracted from mzMLs
            * MSGF: The Number of spectra identified by MSGF search engine
            * Comet: The Number of spectra identified by Comet search engine
            * Sage: The Number of spectra identified by Sage search engine
            * PSMs from quant. peptides: extracted from PSM table in mzTab file
            * Peptides quantified: extracted from PSM table in mzTab file
            """,
    )


def find_modification(peptide):
    peptide = str(peptide)
    pattern = re.compile(r"\((.*?)\)")
    original_mods = pattern.findall(peptide)
    peptide = re.sub(r"\(.*?\)", ".", peptide)
    position = [i.start() for i in re.finditer(r"\.", peptide)]
    for j in range(1, len(position)):
        position[j] -= j

    for k, mod in enumerate(original_mods):
        original_mods[k] = str(position[k]) + "-" + mod

    original_mods = ",".join(str(i) for i in original_mods) if len(original_mods) > 0 else "nan"

    return AASequence.fromString(peptide).toUnmodifiedString(), original_mods


def cal_charge_state(psm, sdrf_file_df):

    charge_state_df = psm[["filename", "charge"]].copy()

    charge_state_df = charge_state_df.merge(
        right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
        left_on="filename",
        right_on="Run"
    )

    charge_state_df["Sample"] = charge_state_df["Sample"].astype(int)

    charge_state_by_run = group_charge(charge_state_df, "filename", "charge")
    charge_state_by_sample = group_charge(charge_state_df, "Sample", "charge")

    return {
        "plot_data": [
            charge_state_by_run.to_dict(orient="index"),
            charge_state_by_sample.to_dict(orient="index"),
        ],
        "cats": sorted(
            set(charge_state_by_run.columns)
            | set(charge_state_by_sample.columns),
            key=int
        )
    }


def sample_level_modifications(df, sdrf_file_df):

    psm = df[
        ["filename", "sequence", "charge", "Modifications"]
    ].copy()

    psm = psm.merge(
        right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
        left_on="filename",
        right_on="Run"
    )

    psm["Sample"] = psm["Sample"].astype(int)

    mod_plot = dict()
    for sample, group in psm.groupby("Sample", sort=True):

        mod_plot_dict, _ = summarize_modifications(
            group[["sequence", "charge", "Modifications"]].drop_duplicates()
        )
        mod_plot[f"Sample {str(sample)}"] = mod_plot_dict

    return mod_plot


def get_peptide_map_by_sample(peptide_map_by_run, sdrf_file_df):

    if sdrf_file_df.empty:
        return {}
    else:
        file_df = sdrf_file_df.copy()
        file_df["Sample"] = file_df["Sample"].astype(int)

        peptide_map_by_sample = dict()
        for sample, group in file_df.groupby("Sample", sort=True):

            peptide_by_sample = list()
            for run in group["Run"].unique():

                run_data = peptide_map_by_run.get(run)
                if not run_data:
                    continue

                peptide_by_sample.extend(run_data)

            peptide_map_by_sample[str(sample)] = len(set(peptide_by_sample))

        return peptide_map_by_sample


def aggregate_spectrum_tracking(
    mzml_table,
    peptide_map_by_sample,
    sdrf_file_df
):

    header_cols = [
        "MS1_Num", "MS2_Num", "MSGF", "Comet", "Sage", "num_quant_psms", "num_quant_peps"
    ]

    for i in header_cols:
        if any([i in v for k, v in mzml_table.items()]):
            pass
        else:
            header_cols.remove(i)

    if sdrf_file_df.empty:

        rows_by_group = dict()
        for sample, value in mzml_table.items():
            for h in header_cols:
                rows_by_group[sample] = {
                    h: value.get(h, "-")
                }
    else:

        rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

        for sample in sorted(
            sdrf_file_df["Sample"].drop_duplicates().tolist(),
            key=lambda x: (str(x).isdigit(), int(x) if str(x).isdigit() else str(x).lower()),
        ):

            file_df_sample = sdrf_file_df[sdrf_file_df["Sample"] == sample].copy()

            sample_data_temp = dict()

            for h in header_cols:

                if h == "num_quant_peps":
                    sample_data_temp[h] = peptide_map_by_sample.get(sample, "-")
                    continue

                sample_data_temp[h] = sum_matching_dict_values(
                    sum_by_run=mzml_table,
                    value_col=h,
                    file_df_by_sample=file_df_sample
                )

            row_data: List[InputRow] = []
            row_data.append(
                InputRow(
                    sample=SampleName(f"Sample {str(sample)}"),
                    data=sample_data_temp,
                )
            )
            for _, row in file_df_sample.iterrows():

                run_data_temp = mzml_table.get(row["Run"], {})

                run_data = dict()
                for h in header_cols:
                    run_data[h] = run_data_temp.get(h, "-")

                row_data.append(
                    InputRow(
                        sample=SampleName(row["Run"]),
                        data=run_data,
                    )
                )
            group_name: SampleGroup = SampleGroup(sample)
            rows_by_group[group_name] = row_data

    return rows_by_group, header_cols


def stat_pep_intensity(intensities: pd.Series):

    stat_result = np.where(
        (intensities.notna()) & (intensities > 0),
        np.log2(intensities).astype(float),
        1.0
    )

    return stat_result.tolist()
