#!/usr/bin/env python

""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
import logging
from collections import OrderedDict
import itertools
from datetime import datetime
from operator import itemgetter
from multiqc import config

from sdrf_pipelines.openms.openms import OpenMS, UnimodDatabase
from multiqc.plots import table, bargraph, linegraph, box
from multiqc.plots.table_object import InputRow
from multiqc.types import SampleGroup, SampleName
from typing import Dict, List

import pandas as pd
from functools import reduce
import re
from pyteomics import mztab, mzid, mgf
from pyopenms import OpenMSBuildInfo, AASequence
import os
import sqlite3
import numpy as np
import copy
import json

from . import sparklines
from .ms_functions import get_ms_qc_info
from ..common import ms_io, common_plots
from ..common.histogram import Histogram
from ..common.calc_utils import qualUniform
from ..common.file_utils import file_prefix
from ..core.section_groups import add_group_modules, add_sub_section
from ..maxquant.maxquant_utils import (
    mod_group_percentage,
    evidence_rt_count,
    evidence_calibrated_mass_error
)
from .quantms_plots import (
    draw_dia_intensitys,
    draw_dia_ms2s,
    draw_dia_time_mass
)


# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

log.info("pyopenms has: " + str(OpenMSBuildInfo().getOpenMPMaxNumThreads()) + " threads.")

class QuantMSModule:

    def __init__(
            self,
            find_log_files_func,
            sub_sections,
            heatmap_colors
        ):

        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections
        self.heatmap_color_list = heatmap_colors


        self.ms_with_psm = list()
        self.total_protein_identified = 0
        self.cal_num_table_data = dict()
        self.oversampling = dict()
        self.identified_spectrum = dict()
        self.delta_mass = dict()
        self.Total_ms2_Spectral_Identified = 0
        self.Total_Peptide_Count = 0
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
        

        self.ms_paths = []
        for mzml_current_file in self.find_log_files("pmultiqc/mzML", filecontents=False):
            self.ms_paths.append(os.path.join(mzml_current_file["root"], mzml_current_file["fn"]))

        self.ms_info_path = []

        # Check if this is a display for mzid identification results
        if config.kwargs.get("mzid_plugin", False):

            self.mzid_peptide_map = dict()
            self.ms_without_psm = dict()

            self.mgf_paths = []
            for mgf_file in self.find_log_files("pmultiqc/mgf", filecontents=False):
                self.mgf_paths.append(os.path.join(mgf_file["root"], mgf_file["fn"]))
            self.mgf_paths.sort()

            self.mzid_paths = []
            for mzid_file in self.find_log_files("pmultiqc/mzid", filecontents=False):
                self.mzid_paths.append(os.path.join(mzid_file["root"], mzid_file["fn"]))
            self.mzid_paths.sort()

            mzid_psm = self.parse_out_mzid()

            if self.mgf_paths:
                self.parse_out_mgf()
            elif self.ms_paths:
                _ = self.parse_mzml()

            self.mzid_cal_heat_map_score(mzid_psm)

            heatmap_data, heatmap_xnames, heatmap_ynames = self.calculate_heatmap()
            common_plots.draw_heatmap(
                self.sub_sections["summary"],
                self.heatmap_color_list,
                heatmap_data,
                heatmap_xnames,
                heatmap_ynames,
                False
            )

            self.draw_summary_protein_ident_table()
            self.draw_mzid_identi_num()
            self.draw_num_pep_per_protein()
            self.draw_precursor_charge_distribution()
            self.draw_peaks_per_ms2()
            self.draw_peak_intensity_distribution()
            common_plots.draw_oversampling(
                self.sub_sections["ms2"],
                self.oversampling,
                self.oversampling_plot.dict["cats"],
                False
            )

        else:

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

            self.enable_exp = False
            self.enable_sdrf = False
            self.msstats_input_valid = False
            # TODO what if multiple are found??
            for f in self.find_log_files("pmultiqc/exp_design", filecontents=False):
                self.exp_design = os.path.join(f["root"], f["fn"])
                self.enable_exp = True

            if not self.enable_exp:

                for f in self.find_log_files("pmultiqc/sdrf", filecontents=False):
                    self.sdrf = os.path.join(f["root"], f["fn"])
                    OpenMS().openms_convert(
                        self.sdrf,
                        config.kwargs["raw"],
                        False,
                        True,
                        False,
                        config.kwargs["condition"],
                    )
                    # experimental_design.tsv is the default output name
                    self.exp_design = os.path.join(f["root"], f["experimental_design.tsv"])
                    self.enable_sdrf = True

            self.psm_table = dict()
            self.mzml_peptide_map = dict()
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
                "end": 1,
                "low_step": 0.03,
                "high_step": 0.08
                }
            self.total_protein_quantified = 0
            self.out_csv_data = dict()
            self.mL_spec_ident_final = dict()
            self.heatmap_con_score = dict()
            self.heatmap_pep_intensity = {}
            self.ms1_tic = dict()
            self.ms1_bpc = dict()
            self.ms1_peaks = dict()
            self.ms1_general_stats = dict()
            self.is_multi_conditions = False

            # draw the experimental design
            if self.enable_exp or self.enable_sdrf:
                self.draw_exp_design()

            for ms_info in self.find_log_files("pmultiqc/ms_info", filecontents=False):
                self.ms_info_path.append(os.path.join(ms_info["root"], ms_info["fn"]))
            self.ms_info_path.sort()
            
            if len(self.ms_info_path) > 0:
                self.read_ms_info = True
                self.ms_paths = [
                    file_prefix(i).replace("_ms_info", ".mzML") for i in self.ms_info_path
                ]

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

            mt = self.parse_mzml()
            self.idx_paths = []
            for idx_file in self.find_log_files("pmultiqc/idXML", filecontents=False):
                self.idx_paths.append(os.path.join(idx_file["root"], idx_file["fn"]))

            self.draw_ms_information()
            if self.enable_dia:
                self.parse_diann_report()
                self.draw_summary_protein_ident_table()
                self.draw_quantms_identi_num()
                self.draw_num_pep_per_protein()
                if len(self.ms_info_path) > 0 and not self.is_bruker:
                    self.draw_precursor_charge_distribution()
                    self.draw_peaks_per_ms2()
                    self.draw_peak_intensity_distribution()
            else:
                if not config.kwargs["ignored_idxml"]:
                    self.parse_idxml(mt)
                self.cal_heat_map_score()
                
                heatmap_data, heatmap_xnames, heatmap_ynames = self.calculate_heatmap()
                common_plots.draw_heatmap(
                    self.sub_sections["summary"],
                    self.heatmap_color_list,
                    heatmap_data,
                    heatmap_xnames,
                    heatmap_ynames,
                    False
                )

                self.draw_summary_protein_ident_table()
                self.draw_quantms_identi_num()
                self.draw_num_pep_per_protein()
                if not config.kwargs["ignored_idxml"]:
                    self.draw_mzml_ms()
                    self.draw_search_engine()
                self.draw_precursor_charge_distribution()
                self.draw_peaks_per_ms2()
                self.draw_peak_intensity_distribution()
                common_plots.draw_oversampling(
                    self.sub_sections["ms2"],
                    self.oversampling,
                    self.oversampling_plot.dict["cats"],
                    False
                )
                self.draw_delta_mass()

            self.draw_quantms_identification(mt)
            self.draw_quantms_contaminants()
            self.draw_quantms_quantification()
            self.draw_quantms_msms_section()
            self.draw_quantms_time_section()

            for msstats_input in self.find_log_files("pmultiqc/msstats", filecontents=False):
                self.msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])
                self.msstats_input_valid = True
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
                            """
                )

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
                        """
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
                "time_mass_sub_section": self.sub_sections["time_mass"],
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
            for k, v in self.heatmap_charge_score.items():
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

    def draw_exp_design(self):
        # Currently this only supports the OpenMS two-table format (default in quantms pipeline)
        # One table format would actually be even easier. You can just use pandas.read_tsv
        self.sample_df, self.file_df = read_openms_design(self.exp_design)

        self.exp_design_runs = np.unique(self.file_df["Run"].tolist())

        if self.file_df["Spectra_Filepath"][0].endswith((".d", ".d.tar")):
            self.is_bruker = True

        # Create table plot
        pattern = r'^(\w+=[^=;]+)(;\w+=[^=;]+)*$'
        self.is_multi_conditions = all(self.sample_df["MSstats_Condition"].apply(lambda x: bool(re.match(pattern, str(x)))))

        rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

        if self.is_multi_conditions:
            for sample in sorted(self.sample_df["Sample"].tolist(), key=lambda x: int(x)):
                file_df_sample = self.file_df[self.file_df["Sample"] == sample].copy()
                sample_df_slice = self.sample_df[self.sample_df["Sample"] == sample].copy()
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
                        sample=SampleName(sample),
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
            headers = {}
            headers["Sample"] = {
                "title": "Sample [Spectra File]",
                "description": "",
                "scale": False,
            }
            for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                headers["MSstats_Condition_" + str(k)] ={
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
            for sample in sorted(self.sample_df["Sample"].tolist(), key=lambda x: int(x)):
                file_df_sample = self.file_df[self.file_df["Sample"] == sample].copy()
                sample_df_slice = self.sample_df[self.sample_df["Sample"] == sample].copy()
                row_data: List[InputRow] = []
                row_data.append(
                    InputRow(
                        sample=SampleName(sample),
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
            sub_section=self.sub_sections["experiment"],
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

    def draw_summary_protein_ident_table(self):
        headers = OrderedDict()
        if self.enable_dia:
            summary_table = {
                self.Total_Peptide_Count: {"#Proteins Quantified": self.total_protein_quantified}
            }
            col_header = "#Peptides Quantified"
        else:
            summary_table = {
                self.total_ms2_spectra: {
                    "#Identified MS2 Spectra": self.Total_ms2_Spectral_Identified
                }
            }
            coverage = self.Total_ms2_Spectral_Identified / self.total_ms2_spectra * 100
            summary_table[self.total_ms2_spectra]["%Identified MS2 Spectra"] = coverage
            summary_table[self.total_ms2_spectra][
                "#Peptides Identified"
            ] = self.Total_Peptide_Count
            summary_table[self.total_ms2_spectra][
                "#Proteins Identified"
            ] = self.total_protein_identified

            if not config.kwargs.get("mzid_plugin", False):
                summary_table[self.total_ms2_spectra][
                    "#Proteins Quantified"
                ] = self.total_protein_quantified

            headers["#Identified MS2 Spectra"] = {
                "description": "Total number of MS/MS spectra identified",
            }
            headers["%Identified MS2 Spectra"] = {
                "description": "Percentage of Identified MS/MS Spectra",
                "format": "{:,.2f}",
                "suffix": "%",
            }
            col_header = "#MS2 Spectra"

        # Create table plot
        pconfig = {
            "id": "identification_summary_table",  # ID used for the table
            "title": "Summary Table",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_summary_table_table",  # File basename to use for raw data file
            "sort_rows": False,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": col_header,
            # 'format': '{:,.0f}',  # The header used for the first column
            "scale": "Set1",
        }

        table_html = table.plot(summary_table, headers, pconfig)

        if config.kwargs.get("mzid_plugin", False):
            description_str = "This plot shows the summary statistics of the submitted data."
            # TODO: add description here @Yasset
            helptext_str = """
                This plot shows the summary statistics of the submitted data.
                """
        else:
            description_str = "This table shows the quantms pipeline summary statistics."
            # TODO: add description here @Yasset
            helptext_str = """
                This table shows the quantms pipeline summary statistics.
                """

        add_sub_section(
            sub_section=self.sub_sections["summary"],
            plot=table_html,
            order=1,
            description=description_str,
            helptext=helptext_str
        )

    def draw_ms_information(self):

        if self.ms1_tic:
            ms1_tic_config = {
                "id": "ms1_tic",
                "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
                "title": "Total Ion Chromatograms",
                "ylab": "Ion Count",
                "xlab": "Retention Time (min)",
                "ymin": 0,
                "showlegend": True,
            }
            ms1_tic_html = linegraph.plot(self.ms1_tic, ms1_tic_config)

            ms1_tic_html = common_plots.remove_subtitle(ms1_tic_html)

            add_sub_section(
                sub_section=self.sub_sections["ms1"],
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
                    """
            )

        if self.ms1_bpc:

            ms1_bpc_config = {
                "id": "ms1_bpc",
                "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
                "title": "MS1 Base Peak Chromatograms",
                "ylab": "Ion Count",
                "xlab": "Retention Time (min)",
                "ymin": 0,
                "showlegend": True,
            }

            ms1_bpc_html = linegraph.plot(self.ms1_bpc, ms1_bpc_config)
            ms1_bpc_html = common_plots.remove_subtitle(ms1_bpc_html)

            add_sub_section(
                sub_section=self.sub_sections["ms1"],
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
                    """
            )

        if self.ms1_peaks:

            ms1_peaks_config = {
                "id": "ms1_peaks",
                "tt_label": "<b>{point.x} Peak Count:</b> {point.y}",
                "title": "MS1 Peaks",
                "ylab": "Peak Count",
                "xlab": "Retention Time (min)",
                "ymin": 0,
                "showlegend": True,
            }

            ms1_peaks_html = linegraph.plot(self.ms1_peaks, ms1_peaks_config)
            ms1_peaks_html = common_plots.remove_subtitle(ms1_peaks_html)

            add_sub_section(
                sub_section=self.sub_sections["ms1"],
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
                    number of compounds present at that time point."
                    * Peak Consistency across samples: Similar profiles among replicates or related samples
                    indicate good analytical reproducibility.
                    * Sudden drops: Abrupt decreases in peak count may point to transient ionization issues,
                    spray instability, or chromatographic disruptions.
                    * Baseline values: The minimum peak count observed reflects the level of background noise
                    or instrument sensitivity in the absence of eluting compounds.
                    
                    Monitoring MS1 peak counts complements total ion chromatogram (TIC) and base peak chromatogram
                    (BPC) data, offering an additional layer of quality control related to signal complexity,
                    instrument stability, and sample composition.
                    """
            )

        if self.ms1_general_stats:
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
            table_html = table.plot(self.ms1_general_stats, headers=headers, pconfig=tconfig)

            add_sub_section(
                sub_section=self.sub_sections["ms1"],
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
                    """
            )

    def draw_quantms_identi_num(self):
        # Create table data
        rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

        if self.enable_exp:
                
            if self.is_multi_conditions:
                for sample in sorted(self.sample_df["Sample"].tolist(), key=lambda x: int(x)):
                    file_df_sample = self.file_df[self.file_df["Sample"] == sample].copy()
                    sample_df_slice = self.sample_df[self.sample_df["Sample"] == sample].copy()
                    row_data: List[InputRow] = []

                    sample_data = {}
                    for k, v in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                        sample_data["MSstats_Condition_" + str(k)] = v
                    
                    sample_data["MSstats_BioReplicate"] = sample_df_slice["MSstats_BioReplicate"].iloc[0]
                    sample_data["Fraction"] = ""
                    sample_data["Peptide_Num"] = ""
                    sample_data["Unique_Peptide_Num"] = ""
                    sample_data["Modified_Peptide_Num"] = ""
                    sample_data["Protein_Num"] = ""

                    row_data.append(
                        InputRow(
                            sample=SampleName(sample),
                            data=sample_data,
                        )
                    )

                    for _, row in file_df_sample.iterrows():

                        sample_data = {}
                        for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                            sample_data["MSstats_Condition_" + str(k)] = ""

                        sample_data["Fraction"] = row["Fraction"]
                        sample_data["Peptide_Num"] = self.cal_num_table_data[row["Run"]]["peptide_num"]
                        sample_data["Unique_Peptide_Num"] = self.cal_num_table_data[row["Run"]]["unique_peptide_num"]
                        sample_data["Modified_Peptide_Num"] = self.cal_num_table_data[row["Run"]]["modified_peptide_num"]
                        sample_data["Protein_Num"] = self.cal_num_table_data[row["Run"]]["protein_num"]

                        row_data.append(
                            InputRow(
                                sample=SampleName(row["Run"]),
                                data=sample_data,
                            )
                        ) 
                    group_name: SampleGroup = SampleGroup(sample)
                    rows_by_group[group_name] = row_data
                    
                headers = {}
                for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                    headers["MSstats_Condition_" + str(k)] ={
                        "title": "MSstats Condition: " + str(k),
                        "description": "",
                        "scale": False,
                    }
                headers["Fraction"] = {
                    "title": "Fraction",
                    "description": "Fraction Identifier",
                    "scale": False,
                }
                headers["Peptide_Num"] = {
                    "title": "#Peptide IDs",
                    "description": "The number of identified PSMs in the pipeline",
                }
                headers["Unique_Peptide_Num"] = {
                    "title": "#Unambiguous Peptide IDs",
                    "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
                }
                headers["Modified_Peptide_Num"] = {
                    "title": "#Modified Peptide IDs",
                    "description": "Number of modified identified peptides in the pipeline",
                }
                headers["Protein_Num"] = {
                    "title": "#Protein (group) IDs",
                    "description": "The number of identified protein(group)s in the pipeline",
                }
            else:
                for sample in sorted(self.sample_df["Sample"].tolist(), key=lambda x: int(x)):
                    file_df_sample = self.file_df[self.file_df["Sample"] == sample].copy()
                    sample_df_slice = self.sample_df[self.sample_df["Sample"] == sample].copy()
                    row_data: List[InputRow] = []
                    row_data.append(
                        InputRow(
                            sample=SampleName(sample),
                            data={
                                "MSstats_Condition": sample_df_slice["MSstats_Condition"].iloc[0],
                                "Fraction": "",
                                "Peptide_Num": "",
                                "Unique_Peptide_Num": "",
                                "Modified_Peptide_Num": "",
                                "Protein_Num": "",
                            },
                        )
                    )
                    for _, row in file_df_sample.iterrows():
                        row_data.append(
                            InputRow(
                                sample=SampleName(row["Run"]),
                                data={
                                    "MSstats_Condition": "",
                                    "Fraction": row["Fraction"],
                                    "Peptide_Num": self.cal_num_table_data[row["Run"]]["peptide_num"],
                                    "Unique_Peptide_Num": self.cal_num_table_data[row["Run"]]["unique_peptide_num"],
                                    "Modified_Peptide_Num": self.cal_num_table_data[row["Run"]]["modified_peptide_num"],
                                    "Protein_Num": self.cal_num_table_data[row["Run"]]["protein_num"],
                                    },
                                )
                            )
                    group_name: SampleGroup = SampleGroup(sample)
                    rows_by_group[group_name] = row_data

                headers = {
                    "MSstats_Condition": {
                        "title": "MSstats_Condition",
                        "description": "MSstats Condition",
                        "scale": False,
                    },
                    "Fraction": {
                        "title": "Fraction",
                        "description": "Fraction Identifier",
                        "scale": False,
                    },
                    "Peptide_Num": {
                        "title": "#Peptide IDs",
                        "description": "The number of identified PSMs in the pipeline",
                    },
                    "Unique_Peptide_Num": {
                        "title": "#Unambiguous Peptide IDs",
                        "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
                    },
                    "Modified_Peptide_Num": {
                        "title": "#Modified Peptide IDs",
                        "description": "Number of modified identified peptides in the pipeline",
                    },
                    "Protein_Num": {
                        "title": "#Protein (group) IDs",
                        "description": "The number of identified protein(group)s in the pipeline",
                    },
                }

        else:
            rows_by_group = dict()
            for sample, value in self.cal_num_table_data.items():
                rows_by_group[sample] = {
                    "Peptide_Num": value["peptide_num"],
                    "Unique_Peptide_Num": value["unique_peptide_num"],
                    "Modified_Peptide_Num": value["modified_peptide_num"],
                    "Protein_Num": value["protein_num"],
                }
            
            headers = {
                "Peptide_Num": {
                    "title": "#Peptide IDs",
                    "description": "The number of identified PSMs in the pipeline",
                },
                "Unique_Peptide_Num": {
                    "title": "#Unambiguous Peptide IDs",
                    "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
                },
                "Modified_Peptide_Num": {
                    "title": "#Modified Peptide IDs",
                    "description": "Number of modified identified peptides in the pipeline",
                },
                "Protein_Num": {
                    "title": "#Protein (group) IDs",
                    "description": "The number of identified protein(group)s in the pipeline",
                },
            }

        # Create table plot
        pconfig = {
            "id": "pipeline_result_statistics",
            "title": "Pipeline Result Statistics",
            "save_file": False,
            "raw_data_fn": "multiqc_pipeline_result_statistics",
            "no_violin": True,
        }
        table_html = table.plot(rows_by_group, headers, pconfig)
        add_sub_section(
            sub_section=self.sub_sections["summary"],
            plot=table_html,
            order=3,
            description="This plot shows the quantms pipeline final result.",
            helptext="""
                Including Sample Name, Possible Study Variables, identified the number of peptide in the pipeline,
                and identified the number of modified peptide in the pipeline, eg. All data in this table are obtained 
                from the out_msstats file. You can also remove the decoy with the `remove_decoy` parameter.
                """
        )

    def draw_mzid_identi_num(self):
        pconfig = {
            "id": "result statistics",  # ID used for the table
            "title": "Pipeline Result Statistics",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_result statistics_table",  # File basename to use for raw data file
            "sort_rows": True,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            "no_violin": True,
            # 'format': '{:,.0f}'  # The header used for the first column
        }
        headers = {
            "peptide_num": {
                "title": "#Peptide IDs",
                "description": "The number of identified PSMs in the pipeline",
            },
            "unique_peptide_num": {
                "title": "#Unambiguous Peptide IDs",
                "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
            },
            "modified_peptide_num": {
                "title": "#Modified Peptide IDs",
                "description": "Number of modified identified peptides in the pipeline",
            },
            "protein_num": {
                "title": "#Protein (group) IDs",
                "description": "The number of identified protein(group)s in the pipeline",
            },
        }
        table_html = table.plot(self.cal_num_table_data, headers, pconfig)
        add_sub_section(
            sub_section=self.sub_sections["summary"],
            plot=table_html,
            order=3,
            description="This plot shows the submitted results",
            helptext="""
                This plot shows the submitted results.
                Including the number of identified peptides and the number of identified modified peptides in the submitted results. 
                You can also remove the decoy with the `remove_decoy` parameter.
                """
        )

    # draw number of peptides per proteins
    def draw_num_pep_per_protein(self):

        if any([len(i) >= 100 for i in self.pep_plot.dict["data"].values()]):
            data_labels = ["Frequency", "Percentage"]
        else:
            data_labels = [
                    {
                        "name": "Frequency",
                        "ylab": "Frequency",
                        "tt_suffix": "",
                        "tt_decimals": 0,
                    },
                    {
                        "name": "Percentage",
                        "ylab": "Percentage [%]",
                        "tt_suffix": "%",
                        "tt_decimals": 2,
                    },
                ]

        pconfig = {
            "id": "number_of_peptides_per_proteins",
            "cpswitch": False,
            "title": "Number of Peptides identified Per Protein",
            "data_labels": data_labels,
        }

        bar_html = bargraph.plot(
            [self.pep_plot.dict["data"]["frequency"], self.pep_plot.dict["data"]["percentage"]],
            ["Frequency", "Percentage"],
            pconfig,
        )
        bar_html = common_plots.remove_subtitle(bar_html)

        if config.kwargs.get("mzid_plugin", False):
            description_str = (
                "This plot shows the number of peptides per protein in the submitted data"
            )
            helptext_str = """
                        Proteins supported by more peptide identifications can constitute more confident results.
                    """
        else:
            description_str = "This plot shows the number of peptides per protein in quantms pipeline final result"
            helptext_str = """
                        This statistic is extracted from the out_msstats file. Proteins supported by more peptide 
                        identifications can constitute more confident results.
                    """
        
        add_sub_section(
            sub_section=self.sub_sections["identification"],
            plot=bar_html,
            order=1,
            description=description_str,
            helptext=helptext_str
        )

    def draw_mzml_ms(self):

        pconfig = {
            "id": "pipeline_spectrum_tracking",  # ID used for the table
            "title": "Pipeline Spectrum Tracking",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_spectrum_tracking_table",  # File basename to use for raw data file
            "sort_rows": False,  # Whether to sort rows alphabetically
            "only_defined_headers": True,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            # 'format': '{:,.0f}'  # The header used for the first column
        }

        headers = OrderedDict()
        headers["MS1_Num"] = {
            "title": "#MS1 Spectra",
            "description": "Number of MS1 spectra",
            "color": "#ffffff",
        }
        headers["MS2_Num"] = {
            "title": "#MS2 Spectra",
            "description": "Number of MS2 spectra",
            "color": "#ffffff",
        }

        if any(["MSGF" in v for k, v in self.mzml_table.items()]):
            headers["MSGF"] = {
                "description": "Number of spectra identified by MSGF search engine",
                "color": "#ffffff",
            }
        if any(["Comet" in v for k, v in self.mzml_table.items()]):
            headers["Comet"] = {
                "description": "Number of spectra identified by Comet search engine",
                "color": "#ffffff",
            }
        if any(["Sage" in v for k, v in self.mzml_table.items()]):
            headers["Sage"] = {
                "description": "Number of spectra identified by Sage search engine",
                "color": "#ffffff",
            }
        headers["num_quant_psms"] = {
            "title": "#PSMs from quant. peptides",
            "description": "Number of reliable PSMs from peptides IDs used in quantification",
            "color": "#ffffff",
        }
        headers["num_quant_peps"] = {
            "title": "#Peptides quantified",
            "description": "Number of quantified peptides that passed final protein and peptide FDR thresholds.",
            "color": "#ffffff",
        }
        table_html = table.plot(self.mzml_table, headers, pconfig)

        add_sub_section(
            sub_section=self.sub_sections["ms2"],
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
                """
        )

    def draw_peak_intensity_distribution(self):
        pconfig = {
            "id": "peak_intensity_distribution",
            "title": "Peak Intensity Distribution",
            "cpswitch": True,
            "ylab": "Count",
            "tt_decimals": 0,
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_peak_distribution_plot.dict["cats"]
        else:
            cats = self.mzml_peak_distribution_plot.dict["cats"]

        bar_html = bargraph.plot(self.ms_info["peak_distribution"], cats, pconfig)
        bar_html = common_plots.remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["ms2"],
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
                """
        )

    def draw_precursor_charge_distribution(self):
        pconfig = {
            "id": "distribution_of_precursor_charges",
            "title": "Distribution of Precursor Charges",
            "cpswitch": True,
            "tt_decimals": 0,
            "ylab": "Count",
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_charge_plot.dict["cats"]
        else:
            cats = self.mzml_charge_plot.dict["cats"]

        bar_html = bargraph.plot(self.ms_info["charge_distribution"], cats, pconfig)
        bar_html = common_plots.remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["ms2"],
            plot=bar_html,
            order=5,
            description="""
                This is a bar chart representing the distribution of the precursor ion charges for a given whole experiment.
                """,
            helptext="""
                This information can be used to identify potential ionization problems 
                including many 1+ charges from an ESI ionization source or an unexpected 
                distribution of charges. MALDI experiments are expected to contain almost exclusively 1+ 
                charged ions. An unexpected charge distribution may furthermore be caused by specific search 
                engine parameter settings such as limiting the search to specific ion charges.
                """
        )

    def draw_peaks_per_ms2(self):
        pconfig = {
            "id": "peaks_per_ms2",
            "cpswitch": True,
            "title": "Number of Peaks per MS/MS spectrum",
            "ylab": "Count",
            "tt_decimals": 0,
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_peaks_ms2_plot.dict["cats"]
        else:
            cats = self.mzml_peaks_ms2_plot.dict["cats"]

        bar_html = bargraph.plot(self.ms_info["peaks_per_ms2"], cats, pconfig)
        bar_html = common_plots.remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["ms2"],
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
                """
        )



    def draw_delta_mass(self):

        delta_mass = self.delta_mass
        delta_mass_percent = {
            "target": {k: v / sum(delta_mass["target"].values()) for k, v in delta_mass["target"].items()}
        }
        
        if delta_mass["decoy"]:

            delta_mass_percent["decoy"] = {k: v / sum(delta_mass["decoy"].values()) for k, v in delta_mass["decoy"].items()}

            x_values = list(delta_mass["target"].keys()) + list(delta_mass["decoy"].keys())

            range_threshold = 10
            if max(abs(x) for x in x_values) > range_threshold:
                range_abs = range_threshold
            else:
                range_abs = 1
            range_step = (max(x_values) - min(x_values)) * 0.05

            if max(abs(x) for x in x_values) > range_abs:

                delta_mass_range = dict()
                delta_mass_range["target"] = {k: v for k, v in delta_mass["target"].items() if abs(k) <= range_abs}
                delta_mass_range["decoy"] = {k: v for k, v in delta_mass["decoy"].items() if abs(k) <= range_abs}

                delta_mass_percent_range = dict()
                delta_mass_percent_range["target"] = {k: v for k, v in delta_mass_percent["target"].items() if abs(k) <= range_abs}
                delta_mass_percent_range["decoy"] = {k: v for k, v in delta_mass_percent["decoy"].items() if abs(k) <= range_abs}
                
                x_values_adj = list(delta_mass_range["target"].keys()) + list(delta_mass_range["decoy"].keys())
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
                }
                line_html = linegraph.plot(
                    [
                        delta_mass_range,
                        delta_mass_percent_range,
                        delta_mass,
                        delta_mass_percent
                    ],
                    pconfig
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
                    "target": {k: v for k, v in delta_mass["target"].items() if abs(k) <= range_abs}
                }

                delta_mass_percent_range = {
                    "target": {k: v for k, v in delta_mass_percent["target"].items() if abs(k) <= range_abs}
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
                }
                line_html = linegraph.plot(
                    [
                        delta_mass_range,
                        delta_mass_percent_range,
                        delta_mass,
                        delta_mass_percent
                    ],
                    pconfig
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
                }
                line_html = linegraph.plot([delta_mass, delta_mass_percent], pconfig)

        add_sub_section(
            sub_section=self.sub_sections["time_mass"],
            plot=line_html,
            order=3,
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
                """
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
            if self.Comet_label
            else ""
        )
        spec_e_bar_html = (
            bargraph.plot(list(self.search_engine["SpecE"].values()), spec_e_cats, spec_e_pconfig)
            if self.MSGF_label
            else ""
        )
        hyper_bar_html = (
            bargraph.plot(list(self.search_engine["hyper"].values()), hyper_cats, hyper_pconfig)
            if self.Sage_label
            else ""
        )

        if spec_e_bar_html != "":

            spec_e_bar_html = common_plots.remove_subtitle(spec_e_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=spec_e_bar_html,
                order=1,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. SpecEvalue: Spectral E-values, 
                        the search score of MSGF. The value used for plotting is -lg(SpecEvalue).
                        """
            )

        if xcorr_bar_html != "":

            xcorr_bar_html = common_plots.remove_subtitle(xcorr_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=xcorr_bar_html,
                order=2,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. xcorr: cross-correlation scores, 
                        the search score of Comet. The value used for plotting is xcorr.
                        """
            )

        if hyper_bar_html != "":

            hyper_bar_html = common_plots.remove_subtitle(hyper_bar_html)

            add_sub_section(
                sub_section=self.sub_sections["search_engine"],
                plot=hyper_bar_html,
                order=3,
                description="",
                helptext="""
                        This statistic is extracted from idXML files. hyperscore: Hyperscore, the search 
                        score of Sage. The value used for plotting is hyperscore.
                        """
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
        }

        pep_bar_html = bargraph.plot(
            list(self.search_engine["PEPs"].values()), pep_cats, pep_pconfig
        )

        pep_bar_html = common_plots.remove_subtitle(pep_bar_html)

        add_sub_section(
            sub_section=self.sub_sections["search_engine"],
            plot=pep_bar_html,
            order=4,
            description="",
            helptext="This statistic is extracted from idXML files."
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
            }

            consensus_bar_html = bargraph.plot(
                self.search_engine["consensus_support"],
                bar_cats,
                consensus_pconfig,
            )
            consensus_bar_html = common_plots.remove_subtitle(consensus_bar_html)

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
                    """
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

    def mzid_cal_heat_map_score(self, psm):
        log.info("{}: Calculating Heatmap Scores...".format(datetime.now().strftime("%H:%M:%S")))

        # HeatMapMissedCleavages
        global_peps = psm[["PeptideSequence", "Modifications"]].drop_duplicates()
        global_peps_count = len(global_peps)

        enzyme_list = list()
        for mzid_path in self.mzid_paths:
            try:
                enzyme_iter = mzid.MzIdentML(mzid_path).iterfind("Enzyme")
                enzyme = next(enzyme_iter).get("EnzymeName", None)
                if enzyme:
                    enzyme_name = list(enzyme.keys())[0]
                else:
                    enzyme_name = "Trypsin"
                enzyme_list.append(enzyme_name)
            except StopIteration:
                enzyme_list.append("Trypsin")

        enzyme_list = list(set(enzyme_list))
        enzyme = enzyme_list[0] if len(enzyme_list) == 1 else "Trypsin"
        psm["missed_cleavages"] = psm.apply(
            lambda x: self.cal_miss_cleavages(x["PeptideSequence"], enzyme), axis=1
        )

        # Calculate the ID RT Score
        if "retention_time" not in psm.columns and self.mgf_paths:
            # MGF
            psm = pd.merge(
                psm,
                self.mgf_rtinseconds[["spectrumID", "filename", "retention_time"]],
                on=["filename", "spectrumID"],
                how="left",
            )

        for name, group in psm.groupby("stand_spectra_ref"):
            sc = group["missed_cleavages"].value_counts()

            self.quantms_missed_cleavages[name] = sc.to_dict()

            mis_0 = sc.get(0, 0)
            self.missed_clevages_heatmap_score[name] = mis_0 / sc[:].sum()
            self.id_rt_score[name] = qualUniform(group["retention_time"])

            #  For HeatMapOverSamplingScore
            self.heatmap_over_sampling_score[name] = self.oversampling[name]["1"] / np.sum(
                list(self.oversampling[name].values())
            )

            # For HeatMapPepMissingScore
            id_fraction = (
                len(
                    pd.merge(
                        global_peps,
                        group[["PeptideSequence", "Modifications"]].drop_duplicates(),
                        on=["PeptideSequence", "Modifications"],
                    ).drop_duplicates()
                )
                / global_peps_count
            )
            self.heatmap_pep_missing_score[name] = np.minimum(1.0, id_fraction)

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
        log.info(
            "{}: Done calculating Heatmap Scores.".format(datetime.now().strftime("%H:%M:%S"))
        )

    def cal_heat_map_score(self):
        log.info("{}: Calculating Heatmap Scores...".format(datetime.now().strftime("%H:%M:%S")))
        mztab_data = mztab.MzTab(self.out_mztab_path)
        psm = mztab_data.spectrum_match_table
        meta_data = dict(mztab_data.metadata)
        if self.pep_table_exists:
            pep_table = mztab_data.peptide_table

            with pd.option_context("future.no_silent_downcasting", True):
                pep_table = pep_table.fillna(np.nan).infer_objects(copy=False).copy()

            pep_table.loc[:, "stand_spectra_ref"] = pep_table.apply(
                lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
                axis=1,
            )
            study_variables = list(
                filter(
                    lambda x: re.match(r"peptide_abundance_study_variable.*?", x) is not None,
                    pep_table.columns.tolist(),
                )
            )

            pep_table["average_intensity"] = pep_table[study_variables].mean(axis=1, skipna=True)

            # Contaminants
            if len(pep_table[pep_table["accession"].str.contains(config.kwargs["contaminant_affix"])]) > 0:

                self.quantms_contaminant_percent = self.cal_quantms_contaminant_percent(
                    pep_table[["average_intensity", "stand_spectra_ref", "accession"]].copy()
                )

                self.quantms_top_contaminant_percent = self.top_n_contaminant_percent(
                    pep_table[["average_intensity", "stand_spectra_ref", "accession"]].copy(),
                    5
                )

            for name, group in pep_table.groupby("stand_spectra_ref"):

                contaminant_sum = group[
                    group["accession"].str.contains(
                        config.kwargs["contaminant_affix"]
                    )][study_variables].sum(axis=0).sum()
                all_sum = group[study_variables].sum(axis=0).sum()
                self.heatmap_con_score[name] = 1.0 - (contaminant_sum / all_sum)

                if config.kwargs["remove_decoy"]:
                    pep_median = np.nanmedian(
                        group[(group["opt_global_cv_MS:1002217_decoy_peptide"] == 0)][
                            study_variables
                        ].to_numpy()
                    )
                    self.quantms_pep_intensity[name] = group[(group["opt_global_cv_MS:1002217_decoy_peptide"] == 0)]["average_intensity"].apply(
                        lambda x: np.log2(x) if pd.notnull(x) and x > 0 else 1
                    )
                else:
                    pep_median = np.nanmedian(group[study_variables].to_numpy())
                    self.quantms_pep_intensity[name] = group["average_intensity"].apply(
                        lambda x: np.log2(x) if pd.notnull(x) and x > 0 else 1
                    )

                self.heatmap_pep_intensity[name] = np.minimum(
                    1.0, pep_median / (2**23)
                )  # Threshold

        #  HeatMapMissedCleavages
        global_peps = set(psm["opt_global_cv_MS:1000889_peptidoform_sequence"])
        global_peps_count = len(global_peps)
        if (
            config.kwargs["remove_decoy"]
            and "opt_global_cv_MS:1002217_decoy_peptide" in psm.columns
        ):
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] == 0].copy()
        psm.loc[:, "stand_spectra_ref"] = psm.apply(
            lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
            axis=1,
        )

        enzyme_list = [i for i in meta_data.values() if str(i).startswith("enzyme:")]
        enzyme = enzyme_list[0].split(":")[1] if len(enzyme_list) == 1 else "Trypsin"
        psm.loc[:, "missed_cleavages"] = psm.apply(
            lambda x: self.cal_miss_cleavages(x["sequence"], enzyme), axis=1
        )

        # Calculate the ID RT Score
        for name, group in psm.groupby("stand_spectra_ref"):
            sc = group["missed_cleavages"].value_counts()

            self.quantms_missed_cleavages[name] = sc.to_dict()

            mis_0 = sc[0] if 0 in sc else 0
            self.missed_clevages_heatmap_score[name] = mis_0 / sc[:].sum()
            self.id_rt_score[name] = qualUniform(group["retention_time"])

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
        log.info(
            "{}: Done calculating Heatmap Scores.".format(datetime.now().strftime("%H:%M:%S"))
        )

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
            cut = "F*,W*,Y*,L*,!*P"
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

    @staticmethod
    def dis_decoy(protein_name):
        if config.kwargs["decoy_affix"] not in protein_name:
            return "TARGET"
        elif protein_name.split(";") == 1:
            return "DECOY"
        else:
            if config.kwargs["affix_type"] == "prefix":
                if list(
                    filter(
                        lambda x: lambda x: not x.startswith(config.kwargs["decoy_affix"]),
                        protein_name.split(";"),
                    )
                ):
                    return "TARGET"
                return "DECOY"
            else:
                if list(
                    filter(
                        lambda x: not x.endswith(config.kwargs["decoy_affix"]),
                        protein_name.split(";"),
                    )
                ):
                    return "TARGET"
                return "DECOY"

    def parse_mzml(self):
        if self.is_bruker and self.read_ms_info:
            for file in self.ms_info_path:
                log.info(
                    "{}: Parsing ms_statistics dataframe {}...".format(
                        datetime.now().strftime("%H:%M:%S"), file
                    )
                )
                mzml_df = pd.read_csv(file, sep="\t")
                (
                    self.ms1_tic[os.path.basename(file).replace("_ms_info.tsv", "")],
                    self.ms1_bpc[os.path.basename(file).replace("_ms_info.tsv", "")],
                    self.ms1_peaks[os.path.basename(file).replace("_ms_info.tsv", "")],
                    self.ms1_general_stats[os.path.basename(file).replace("_ms_info.tsv", "")],
                ) = get_ms_qc_info(mzml_df)

                log.info(
                    "{}: Done aggregating ms_statistics dataframe {}...".format(
                        datetime.now().strftime("%H:%M:%S"), file
                    )
                )
            return

        self.mzml_peak_distribution_plot = Histogram(
            "Peak Intensity",
            plot_category="range",
            breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000],
        )

        self.mzml_charge_plot = Histogram("Precursor Charge", plot_category="frequency")

        self.mzml_peaks_ms2_plot = Histogram(
            "#Peaks per MS/MS spectrum",
            plot_category="range",
            breaks=[i for i in range(0, 1001, 100)],
        )

        # New instances are used for dictionary construction.
        self.mzml_peak_distribution_plot_1 = copy.deepcopy(self.mzml_peak_distribution_plot)
        self.mzml_charge_plot_1 = copy.deepcopy(self.mzml_charge_plot)
        self.mzml_peaks_ms2_plot_1 = copy.deepcopy(self.mzml_peaks_ms2_plot)

        self.ms_without_psm = []

        mzml_table = {}
        heatmap_charge = {}

        # Use the refactored functions from ms_io.py
        if self.read_ms_info:
            result = ms_io.read_ms_info(
                self.ms_info_path,
                self.ms_with_psm,
                self.identified_spectrum,
                self.mzml_charge_plot,
                self.mzml_peak_distribution_plot,
                self.mzml_peaks_ms2_plot,
                self.mzml_charge_plot_1,
                self.mzml_peak_distribution_plot_1,
                self.mzml_peaks_ms2_plot_1,
                self.ms_without_psm,
                self.enable_dia,
            )
            (
                mzml_table,
                heatmap_charge,
                self.total_ms2_spectra,
                self.ms1_tic,
                self.ms1_bpc,
                self.ms1_peaks,
                self.ms1_general_stats,
            ) = result
        else:
            result = ms_io.read_mzmls(
                self.ms_paths,
                self.ms_with_psm,
                self.identified_spectrum,
                self.mzml_charge_plot,
                self.mzml_peak_distribution_plot,
                self.mzml_peaks_ms2_plot,
                self.mzml_charge_plot_1,
                self.mzml_peak_distribution_plot_1,
                self.mzml_peaks_ms2_plot_1,
                self.ms_without_psm,
                self.enable_dia,
            )
            mzml_table, heatmap_charge, self.total_ms2_spectra = result

        for i in self.ms_without_psm:
            log.warning("No PSM found in '{}'!".format(i))

        self.mzml_peaks_ms2_plot.to_dict()
        self.mzml_peak_distribution_plot.to_dict()
        # Construct compound dictionaries to apply to drawing functions.
        if self.enable_dia:
            self.mzml_charge_plot.to_dict()

            self.ms_info["charge_distribution"] = {
                "Whole Experiment": self.mzml_charge_plot.dict["data"]
            }
            self.ms_info["peaks_per_ms2"] = {
                "Whole Experiment": self.mzml_peaks_ms2_plot.dict["data"]
            }
            self.ms_info["peak_distribution"] = {
                "Whole Experiment": self.mzml_peak_distribution_plot.dict["data"]
            }
        else:
            self.mzml_peaks_ms2_plot_1.to_dict()
            self.mzml_peak_distribution_plot_1.to_dict()
            self.mzml_charge_plot.to_dict()
            self.mzml_charge_plot_1.to_dict()

            self.mzml_charge_plot.dict["cats"].update(self.mzml_charge_plot_1.dict["cats"])
            charge_cats_keys = [int(i) for i in self.mzml_charge_plot.dict["cats"]]
            charge_cats_keys.sort()
            self.mzml_charge_plot.dict["cats"] = OrderedDict(
                {str(i): self.mzml_charge_plot.dict["cats"][str(i)] for i in charge_cats_keys}
            )

            self.ms_info["charge_distribution"] = {
                "identified_spectra": self.mzml_charge_plot.dict["data"],
                "unidentified_spectra": self.mzml_charge_plot_1.dict["data"],
            }
            self.ms_info["peaks_per_ms2"] = {
                "identified_spectra": self.mzml_peaks_ms2_plot.dict["data"],
                "unidentified_spectra": self.mzml_peaks_ms2_plot_1.dict["data"],
            }
            self.ms_info["peak_distribution"] = {
                "identified_spectra": self.mzml_peak_distribution_plot.dict["data"],
                "unidentified_spectra": self.mzml_peak_distribution_plot_1.dict["data"],
            }

        median = np.median(list(heatmap_charge.values()))
        self.heatmap_charge_score = dict(
            zip(
                heatmap_charge.keys(),
                list(map(lambda v: 1 - np.abs(v - median), heatmap_charge.values())),
            )
        )

        return mzml_table

    def parse_idxml(self, mzml_table):
        # Use the refactored function from ms_io.py
        result = ms_io.parse_idxml(
            self.idx_paths,
            mzml_table,
            self.xcorr_hist_range,
            self.hyper_hist_range,
            self.spec_evalue_hist_range,
            self.pep_hist_range,
            self.mL_spec_ident_final,
            self.mzml_peptide_map,
            config.kwargs["remove_decoy"],
        )

        self.search_engine, self.MSGF_label, self.Comet_label, self.Sage_label = result

        # mass spectrum files sorted based on experimental file
        for spectrum_name in self.exp_design_runs:
            self.mzml_table[spectrum_name] = mzml_table[spectrum_name]

    def parse_out_mztab(self):
        log.info(
            "{}: Parsing mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mztab_path
            )
        )
        mztab_data = mztab.MzTab(self.out_mztab_path)
        log.info(
            "{}: Done parsing mzTab file {}.".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mztab_path
            )
        )
        log.info(
            "{}: Aggregating mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mztab_path
            )
        )
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)

        self.delta_mass["target"] = dict()
        self.delta_mass["decoy"] = dict()

        # PSM table data
        psm = mztab_data.spectrum_match_table
        if len(psm) == 0:
            raise ValueError("The PSM section of mzTab is missing, please check your mzTab!")

        # Generate "opt_global_cv_MS: 1002217_DECOY_peptide" column if this column is not contained in the PSM subtable
        if "opt_global_cv_MS:1002217_decoy_peptide" not in psm.columns.values:
            psm["opt_global_cv_MS:1002217_decoy_peptide"] = psm.apply(
                lambda x: 1 if self.dis_decoy(x["accession"]) == "DECOY" else 0, axis=1
            )
        # map to spectrum file name in experimental design file
        psm["stand_spectra_ref"] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(":")[0] + "-location"])
            + ":"
            + x.spectra_ref.split(":")[1],
            axis=1,
        )
        psm["filename"] = psm.apply(
            lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
            axis=1,
        )
        self.ms_with_psm = psm["filename"].unique().tolist()

        prot = mztab_data.protein_table
        self.prot_search_score = dict()

        prot_abundance_cols = list(
            filter(
                lambda x: re.match(r"protein_abundance_assay.*?", x) is not None,
                prot.columns.tolist(),
            )
        )
        opt_cols = list(filter(lambda x: x.startswith("opt_"), prot.columns.tolist()))
        score_cols = list(
            filter(lambda x: x.startswith("best_search_engine_score"), prot.columns.tolist())
        )
        # TODO in theory we do not need accession since the index is the accession
        fixed_cols = [
            "accession",
            "description",
            "taxid",
            "species",
            "database",
            "database_version",
            "search_engine",
            "ambiguity_members",
            "modifications",
            "protein_coverage",
        ]

        prot = prot[fixed_cols + score_cols + prot_abundance_cols + opt_cols]

        # We only need the overall protein (group) scores and abundances. Currently we do not care about details of single proteins (length, description,...)
        prot = prot[prot["opt_global_result_type"] != "protein_details"].copy()

        if config.kwargs["remove_decoy"]:
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1].copy()
            # TODO do we really want to remove groups that contain a single decoy? I would say ALL members need to be decoy.
            prot = prot[~prot["accession"].str.contains(config.kwargs["decoy_affix"])]

        prot.dropna(subset=["ambiguity_members"], inplace=True)

        prot["protein_group"] = prot["ambiguity_members"].apply(lambda x: x.replace(",", ";"))

        self.total_protein_identified = len(prot.index)

        prot.dropna(how="all", subset=prot_abundance_cols, inplace=True)
        self.total_protein_quantified = len(prot.index)

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
            del peptide_score
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
        
        psm["Modifications"] = psm["modifications"].apply(lambda x: get_unimod_modification(x))

        mod_plot_dict = dict()
        modified_cats = list()

        for m, group in psm.groupby("filename"):
            m = os.path.basename(m)

            if config.kwargs["remove_decoy"]:
                group = group[group["opt_global_cv_MS:1002217_decoy_peptide"] == 0]

            # Modifications
            mod_group_processed = mod_group_percentage(group[["sequence", "charge", "Modifications"]].drop_duplicates())
            mod_plot_dict[m] = dict(
                zip(mod_group_processed["modifications"], mod_group_processed["percentage"])
            )
            modified_cats.extend(mod_group_processed["modifications"])

            # Identified MS2 Spectra Raw File:
            self.identified_msms_spectra[m] = {
                "Identified": len(set(group["spectra_ref"]))
            }

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
            self.mzml_peptide_map[m] = list(set(group["sequence"].tolist()))

            if None in proteins:
                proteins.remove(None)

            ## TODO this is not really the number of proteins but the number of protein groups
            self.cal_num_table_data[m] = {"protein_num": len(proteins)}
            self.cal_num_table_data[m]["peptide_num"] = len(peptides)
            self.cal_num_table_data[m]["unique_peptide_num"] = len(unique_peptides)

            modified_pep = list(
                filter(lambda x: re.match(r".*?\(.*\).*?", x) is not None, peptides)
            )
            self.cal_num_table_data[m]["modified_peptide_num"] = len(modified_pep)

            ml_spec_ident_final[m] = len(set(self.identified_spectrum[m]))

        # Modifications
        self.quantms_modified["plot_data"] = mod_plot_dict
        self.quantms_modified["cats"] = list(sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x)))

        # Charge-state of Per File
        charge_state_df = psm.groupby(["filename", "charge"]).size().unstack(fill_value=0)
        charge_state_df.rename(columns=lambda x: str(x), inplace=True)
        self.mztab_charge_state = {
            "plot_data": charge_state_df.to_dict(orient="index"),
            "cats": charge_state_df.columns.tolist()
        }

        # IDs over RT
        quantms_rt_file_df = psm[["filename", "retention_time"]].copy()
        quantms_rt_file_df["retention time"] = quantms_rt_file_df["retention_time"] / 60
        quantms_rt_file_df.rename(columns={"filename": "raw file"}, inplace=True)
        self.quantms_ids_over_rt = evidence_rt_count(quantms_rt_file_df)

        # Delta Mass [ppm]
        mass_error = psm[["filename", "calc_mass_to_charge", "exp_mass_to_charge"]].copy()
        mass_error["mass error [ppm]"] = (
            (mass_error["exp_mass_to_charge"] - mass_error["calc_mass_to_charge"]) / mass_error["calc_mass_to_charge"]
        ) * 1e6
        mass_error.rename(columns={"filename": "raw file"}, inplace=True)
        self.quantms_mass_error = evidence_calibrated_mass_error(mass_error)

        # TODO mzMLs without PSM: experimental design information is displayed, and all quantitative information is 0
        self.ms_without_psm = set([file_prefix(i) for i in self.ms_paths]) - set(
            self.ms_with_psm
        )
        for i in self.ms_without_psm:
            self.cal_num_table_data[i] = {
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
        except Exception as e:
            log.info(
                "{}: No decoy peptides found -> only showing target peptides.".format(
                    datetime.now().strftime("%H:%M:%S")
                )
            )

        target_bin = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1][
            "relative_diff"
        ].value_counts(sort=False, bins=1000)
        for index in target_bin.index:
            target_bin_data[float(index.mid)] = int(target_bin[index])

        self.delta_mass["target"] = target_bin_data

        # extract delta mass
        self.mL_spec_ident_final = ml_spec_ident_final
        if config.kwargs["remove_decoy"]:
            self.Total_ms2_Spectral_Identified = len(
                set(psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1]["spectra_ref"])
            )
            self.Total_Peptide_Count = len(
                set(psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1]["sequence"])
            )
        else:
            self.Total_ms2_Spectral_Identified = len(set(psm["spectra_ref"]))
            self.Total_Peptide_Count = len(set(psm["sequence"]))

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
                "table_title": "information of peptide spectrum matches",
                # Title of the table. Used in the column config modal
                "save_file": False,  # Whether to save the table data to a file
                "sortRows": False,  # Whether to sort rows alphabetically
                "only_defined_headers": False,  # Only show columns that are defined in the headers config
                "col1_header": "PSM_ID",
                "format": "{:,.0f}",
                "no_violin": True,
            }

            mztab_data_psm_init = dict(itertools.islice(mztab_data_psm_full.items(), 50))
            table_html = table.plot(mztab_data_psm_init, headers, pconfig)
            pattern = re.compile(r'<small id="peptide_spectrum_matches_numrows_text"')
            index = re.search(pattern, table_html).span()[0]
            t_html = (
                table_html[:index]
                + '<input type="text" placeholder="search..." class="searchInput" '
                'onkeyup="searchPsmFunction()" id="psm_search">'
                '<select name="psm_search_col" id="psm_search_col">'
            )
            for key in ["Sequence", "Modification", "Accession", "Spectra_Ref"]:
                t_html += "<option>" + key + "</option>"
            table_html = (
                t_html + "</select>" + '<button type="button" class="btn btn-default '
                'btn-sm" id="psm_reset" onclick="psmFirst()">Reset</button>' + table_html[index:]
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

            log.info(
                "{}: Done aggregating mzTab file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.out_mztab_path
                )
            )

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
                "anchor": "",
                # Title of the table. Used in the column config modal
                "save_file": False,  # Whether to save the table data to a file
                "raw_data_fn": "multiqc_quantification_of_protein_table",  # File basename to use for raw data file
                "sort_rows": False,  # Whether to sort rows alphabetically
                "only_defined_headers": False,  # Only show columns that are defined in the headers config
                "col1_header": "ProteinName",
                "no_violin": True,
            }

            max_prot_intensity = 0
            mztab_data_dict_prot_init = dict(
                itertools.islice(mztab_data_dict_prot_full.items(), 50)
            )

            table_html = sparklines.plot(
                mztab_data_dict_prot_init, headers, pconfig=pconfig, max_value=max_prot_intensity
            )
            pattern = re.compile(r'<small id="quantification_of_protein_numrows_text"')
            index = re.search(pattern, table_html).span()[0]
            t_html = (
                table_html[:index]
                + '<input type="text" placeholder="search..." class="searchInput" '
                'onkeyup="searchProtFunction()" id="prot_search">'
                '<select name="prot_search_col" id="prot_search_col">'
            )
            for key in ["ProteinName"]:
                t_html += "<option>" + key + "</option>"
            table_html = (
                t_html + "</select>" + '<button type="button" class="btn btn-default '
                'btn-sm" id="prot_reset" onclick="protFirst()">Reset</button>' + table_html[index:]
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

    def parse_out_mgf(self):
        def get_spectrum_id(spectrum_title, index):
            if "scan=" in spectrum_title:
                spectrum_id = spectrum_title
            else:
                spectrum_id = "index=" + str(index)
            return spectrum_id

        self.mgf_peak_distribution_plot = Histogram(
            "Peak Intensity",
            plot_category="range",
            breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000],
        )
        self.mgf_charge_plot = Histogram("Precursor Charge", plot_category="frequency")
        self.mgf_peaks_ms2_plot = Histogram(
            "#Peaks per MS/MS spectrum",
            plot_category="range",
            breaks=[i for i in range(0, 1001, 100)],
        )

        self.mgf_peak_distribution_plot_1 = copy.deepcopy(self.mgf_peak_distribution_plot)
        self.mgf_charge_plot_1 = copy.deepcopy(self.mgf_charge_plot)
        self.mgf_peaks_ms2_plot_1 = copy.deepcopy(self.mgf_peaks_ms2_plot)

        heatmap_charge = dict()
        mgf_rtinseconds = {"spectrumID": [], "title": [], "filename": [], "retention_time": []}

        for m in self.mgf_paths:
            log.info("{}: Parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
            mgf_data = mgf.MGF(m)
            log.info(
                "{}: Done parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )
            m = file_prefix(m)
            log.info(
                "{}: Aggregating MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )

            charge_2 = 0

            for i, spectrum in enumerate(mgf_data):
                charge_state = int(spectrum.get("params", {}).get("charge", [])[0])
                if charge_state == 2:
                    charge_2 += 1

                peak_per_ms2 = len(spectrum["m/z array"])
                base_peak_intensity = (
                    max(spectrum["intensity array"])
                    if len(spectrum["intensity array"]) > 0
                    else None
                )

                raw_title = spectrum.get("params", {}).get("title", [])
                mgf_rtinseconds["title"].append(raw_title)
                title = get_spectrum_id(raw_title, i)
                mgf_rtinseconds["spectrumID"].append(title)
                mgf_rtinseconds["filename"].append(m)

                rtinseconds = float(spectrum.get("params", {}).get("rtinseconds", None))
                mgf_rtinseconds["retention_time"].append(rtinseconds)

                if m in self.ms_with_psm:
                    if title in self.identified_spectrum[m]:
                        self.mgf_charge_plot.add_value(charge_state)
                        self.mgf_peak_distribution_plot.add_value(base_peak_intensity)
                        self.mgf_peaks_ms2_plot.add_value(peak_per_ms2)
                    else:
                        self.mgf_charge_plot_1.add_value(charge_state)
                        self.mgf_peak_distribution_plot_1.add_value(base_peak_intensity)
                        self.mgf_peaks_ms2_plot_1.add_value(peak_per_ms2)
                else:
                    if m not in self.ms_without_psm:
                        self.ms_without_psm.append(m)
            ms2_number = i + 1

            heatmap_charge[m] = charge_2 / ms2_number
            self.total_ms2_spectra = self.total_ms2_spectra + ms2_number
            log.info(
                "{}: Done aggregating MGF file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )

        self.mgf_rtinseconds = pd.DataFrame(mgf_rtinseconds)

        for i in self.ms_without_psm:
            log.warning("No PSM found in '{}'!".format(i))

        self.mgf_peaks_ms2_plot.to_dict()
        self.mgf_peak_distribution_plot.to_dict()
        self.mgf_peaks_ms2_plot_1.to_dict()
        self.mgf_peak_distribution_plot_1.to_dict()
        self.mgf_charge_plot.to_dict()
        self.mgf_charge_plot_1.to_dict()

        self.mgf_charge_plot.dict["cats"].update(self.mgf_charge_plot_1.dict["cats"])
        charge_cats_keys = [int(i) for i in self.mgf_charge_plot.dict["cats"]]
        charge_cats_keys.sort()
        self.mgf_charge_plot.dict["cats"] = OrderedDict(
            {str(i): self.mgf_charge_plot.dict["cats"][str(i)] for i in charge_cats_keys}
        )

        self.ms_info["charge_distribution"] = {
            "identified_spectra": self.mgf_charge_plot.dict["data"],
            "unidentified_spectra": self.mgf_charge_plot_1.dict["data"],
        }
        self.ms_info["peaks_per_ms2"] = {
            "identified_spectra": self.mgf_peaks_ms2_plot.dict["data"],
            "unidentified_spectra": self.mgf_peaks_ms2_plot_1.dict["data"],
        }
        self.ms_info["peak_distribution"] = {
            "identified_spectra": self.mgf_peak_distribution_plot.dict["data"],
            "unidentified_spectra": self.mgf_peak_distribution_plot_1.dict["data"],
        }

        median = np.median(list(heatmap_charge.values()))
        self.heatmap_charge_score = dict(
            zip(
                heatmap_charge.keys(),
                list(map(lambda v: 1 - np.abs(v - median), heatmap_charge.values())),
            )
        )

    def parse_out_mzid(self):
        def parse_location(location):
            if "\\" in location:
                location = location.replace("\\", "/")
            return os.path.basename(location)

        def process_modification(modification):
            if not isinstance(modification, list):
                modifications = None
            else:
                modifi_list = list()
                for i in modification:
                    if i.get("name", None) is not None:
                        modifi_list.append(str(i.get("location")) + "-" + i.get("name", None))
                    elif i.get("cross-link receiver", None) is not None:
                        modifi_list.append(str(i.get("location")) + "-CrossLinkReceiver")
                modifications = ";".join(modifi_list)
            return modifications

        def read_mzid(m):
            log.info(
                "{}: Parsing MzIdentML file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )
            mzid_data = mzid.MzIdentML(m)
            if len(mzid_data) == 0:
                raise ValueError("Please check your MzIdentML", m)
            log.info(
                "{}: Done parsing MzIdentML file {}.".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )
            m = file_prefix(m)
            log.info(
                "{}: Aggregating MzIdentML file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )

            mzid_dict = list()
            for mzid_tmp in mzid_data:
                mzid_tmp_part = {
                    k: v for k, v in mzid_tmp.items() if k not in ["SpectrumIdentificationItem"]
                }
                for spectrum_item in mzid_tmp.get("SpectrumIdentificationItem", []):
                    spectrum_item_part = {
                        k: v
                        for k, v in spectrum_item.items()
                        if k not in ["PeptideEvidenceRef", "PeptideSequence"]
                    }
                    if (
                        spectrum_item_part.get("rank") == 1
                        and spectrum_item_part.get("peptide passes threshold", "true") == "true"
                        and spectrum_item.get("passThreshold", False)
                    ):
                        for peptide_ref in spectrum_item.get("PeptideEvidenceRef", []):
                            if "name" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_name"] = peptide_ref.pop("name")
                            if "location" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_location"] = peptide_ref.pop(
                                    "location"
                                )
                            if "FileFormat" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_FileFormat"] = peptide_ref.pop(
                                    "FileFormat"
                                )
                            mzid_dict.append(
                                {**mzid_tmp_part, **spectrum_item_part, **peptide_ref}
                            )

            mzid_t = pd.DataFrame(mzid_dict)
            log.info(
                "{}: Done aggregating MzIdentML file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )

            return mzid_t

        mzid_table = pd.DataFrame()
        for mzid_path in self.mzid_paths:
            mzid_table = pd.concat([mzid_table, read_mzid(mzid_path)], axis=0, ignore_index=True)

        search_engines = ["SEQUEST:xcorr", "Mascot:score", "PEAKS:peptideScore", "xi:score"]
        mzid_table.rename(
            columns=lambda x: "search_engine_score" if x in search_engines else x, inplace=True
        )

        if "retention time" in mzid_table.columns:
            mzid_table.rename(columns={"retention time": "retention_time"}, inplace=True)

        if "location" in mzid_table.columns:
            mzid_table["location"] = mzid_table["location"].apply(parse_location)

        if "PeptideEvidenceRef_location" in mzid_table.columns:
            mzid_table["PeptideEvidenceRef_location"] = mzid_table[
                "PeptideEvidenceRef_location"
            ].apply(parse_location)

        mzid_table["stand_spectra_ref"] = mzid_table.apply(
            lambda x: file_prefix(x.location), axis=1
        )
        mzid_table["Modifications"] = mzid_table["Modification"].apply(
            lambda x: process_modification(x)
        )

        if "cross-link spectrum identification item" in mzid_table.columns:
            mzid_table["CrossLinked_Peptide"] = ~mzid_table[
                "cross-link spectrum identification item"
            ].isna()
        else:
            mzid_table["CrossLinked_Peptide"] = False

        mzid_table["pep_to_prot_unique"] = mzid_table.groupby(["spectrumID", "PeptideSequence"])[
            "accession"
        ].transform(lambda x: len(set(x)) == 1)
        mzid_table["accession_group"] = mzid_table.groupby(["spectrumID", "PeptideSequence"])[
            "accession"
        ].transform(lambda x: ";".join(x.unique()))

        if "isDecoy" not in mzid_table.columns:
            mzid_table["isDecoy"] = False

        mzid_table["filename"] = mzid_table.apply(lambda x: file_prefix(x.location), axis=1)
        self.ms_with_psm = mzid_table["filename"].unique().tolist()

        # TODO remove_decoy
        if config.kwargs["remove_decoy"]:
            mzid_table = mzid_table[~mzid_table["isDecoy"]]

        self.total_protein_identified = mzid_table["accession_group"].nunique()

        self.pep_plot = Histogram("Number of peptides per proteins", plot_category="frequency")

        counts_per_acc = (
            mzid_table[["PeptideSequence", "accession"]]
            .drop_duplicates()["accession"]
            .value_counts()
        )
        counts_per_acc.apply(self.pep_plot.add_value)

        categories = OrderedDict()
        categories["Frequency"] = {
            "name": "Frequency",
            "description": "Number of identified peptides per protein.",
        }
        self.pep_plot.to_dict(percentage=True, cats=categories)

        if "retention_time" in mzid_table.columns:
            psm = (
                mzid_table[
                    [
                        "spectrumID",
                        "PeptideSequence",
                        "chargeState",
                        "Modifications",
                        "accession_group",
                        "experimentalMassToCharge",
                        "calculatedMassToCharge",
                        "isDecoy",
                        "pep_to_prot_unique",
                        "search_engine_score",
                        "stand_spectra_ref",
                        "location",
                        "filename",
                        "retention_time",
                    ]
                ]
                .drop_duplicates()
                .reset_index()
            )
        else:
            psm = (
                mzid_table[
                    [
                        "spectrumID",
                        "PeptideSequence",
                        "chargeState",
                        "Modifications",
                        "accession_group",
                        "experimentalMassToCharge",
                        "calculatedMassToCharge",
                        "isDecoy",
                        "pep_to_prot_unique",
                        "search_engine_score",
                        "stand_spectra_ref",
                        "location",
                        "filename",
                    ]
                ]
                .drop_duplicates()
                .reset_index()
            )

        for m, group in psm.groupby("filename"):
            self.oversampling_plot = Histogram(
                "MS/MS counts per 3D-peak", plot_category="frequency", breaks=[1, 2, 3]
            )
            for _, j in group.groupby(["PeptideSequence", "chargeState", "Modifications"]):
                self.oversampling_plot.add_value(j["spectrumID"].nunique())

            self.oversampling_plot.to_dict()
            self.oversampling[m] = self.oversampling_plot.dict["data"]

            proteins = set(group["accession_group"])
            peptides = group[["PeptideSequence", "Modifications"]].drop_duplicates()

            unique_peptides = group[group["pep_to_prot_unique"]][
                ["PeptideSequence", "Modifications"]
            ].drop_duplicates()

            self.identified_spectrum[m] = group["spectrumID"].drop_duplicates().tolist()
            self.mzid_peptide_map[m] = list(set(group["PeptideSequence"].tolist()))

            if None in proteins:
                proteins.remove(None)

            self.cal_num_table_data[m] = {"protein_num": len(proteins)}
            self.cal_num_table_data[m]["peptide_num"] = len(peptides)
            self.cal_num_table_data[m]["unique_peptide_num"] = len(unique_peptides)

            modified_pep = peptides.dropna(subset=["Modifications"])
            self.cal_num_table_data[m]["modified_peptide_num"] = len(modified_pep)

        if self.mgf_paths:
            self.ms_without_psm = set([file_prefix(i) for i in self.mgf_paths]) - set(
                self.ms_with_psm
            )
        elif self.ms_paths:
            self.ms_without_psm = set([file_prefix(i) for i in self.ms_paths]) - set(
                self.ms_with_psm
            )
        for i in self.ms_without_psm:
            self.cal_num_table_data[file_prefix(i)] = {
                "protein_num": 0,
                "peptide_num": 0,
                "unique_peptide_num": 0,
                "modified_peptide_num": 0,
            }

        self.Total_ms2_Spectral_Identified = psm["spectrumID"].nunique()
        self.Total_Peptide_Count = psm["PeptideSequence"].nunique()

        return psm

    def parse_diann_report(self):

        log.info("Parsing {}...".format(self.diann_report_path))

        # parse DIA-NN report data
        if os.path.splitext(self.diann_report_path)[1] == ".tsv":
            report_data = pd.read_csv(self.diann_report_path, header=0, sep="\t")
        else:
            report_data = pd.read_parquet(self.diann_report_path)

        # Draw: Standard Deviation of Intensity
        if "Precursor.Normalised" in report_data.columns:
            draw_dia_intensitys(self.sub_sections["quantification"], report_data)
        
        draw_dia_ms2s(self.sub_sections["ms2"], report_data)
        draw_dia_time_mass(self.sub_sections["time_mass"], report_data)

        pattern = re.compile(r"\(.*?\)")
        report_data["sequence"] = report_data.apply(
            lambda x: re.sub(pattern, "", x["Modified.Sequence"]), axis=1
        )
        self.total_protein_quantified = len(set(report_data["Protein.Names"]))
        self.Total_Peptide_Count = len(set(report_data["sequence"]))
        protein_pep_map = report_data.groupby("Protein.Group").sequence.apply(list).to_dict()

        self.pep_plot = Histogram("number of peptides per proteins", plot_category="frequency")

        for _, peps in protein_pep_map.items():
            number = len(set(peps))
            self.pep_plot.add_value(number)

        self.peptide_search_score = dict()
        pattern = re.compile(r"\((.*?)\)")
        unimod_data = UnimodDatabase()
        for peptide, group in report_data.groupby("Modified.Sequence"):
            origianl_mods = re.findall(pattern, peptide)
            for mod in set(origianl_mods):
                name = unimod_data.get_by_accession(mod.upper()).get_name()
                peptide = peptide.replace(mod, name)
            if peptide.startswith("("):
                peptide = peptide + "."

            self.peptide_search_score[peptide] = np.min(group["Q.Value"])

        categorys = OrderedDict()
        categorys["Frequency"] = {
            "name": "Frequency",
            "description": "number of peptides per proteins",
        }
        self.pep_plot.to_dict(percentage=True, cats=categorys)

        # Modifications Name
        mod_pattern = re.compile(r"\((.*?)\)")
        def find_diann_modified(peptide):
            if isinstance(peptide, str):
                mods = mod_pattern.findall(peptide)
                if mods:
                    mod_type = [
                        unimod_data.get_by_accession(mod.upper()).get_name() for mod in set(mods)
                    ]
                    return ",".join(mod_type)
                else:
                    return "Unmodified"
            return None

        report_data["Modifications"] = report_data["Modified.Sequence"].apply(lambda x: find_diann_modified(x))

        mod_plot_dict = dict()
        modified_cats = list()

        for run_file, group in report_data.groupby("Run"):
            self.ms_with_psm.append(run_file)

            # Modifications
            mod_group_processed = mod_group_percentage(group.drop_duplicates())
            mod_plot_dict[run_file] = dict(
                zip(mod_group_processed["modifications"], mod_group_processed["percentage"])
            )
            modified_cats.extend(mod_group_processed["modifications"])

            self.cal_num_table_data[run_file] = {"protein_num": len(set(group["Protein.Ids"]))}
            self.cal_num_table_data[run_file]["peptide_num"] = len(set(group["sequence"]))
            peptides = set(group["Modified.Sequence"])
            modified_pep = list(
                filter(lambda x: re.match(r".*?\(.*?\).*?", x) is not None, peptides)
            )
            group_peptides = group.groupby("sequence")["Protein.Group"].apply(list).to_dict()
            unique_peptides = [
                pep for pep, prots in group_peptides.items() if len(set(prots)) == 1
            ]
            self.cal_num_table_data[run_file]["unique_peptide_num"] = len(unique_peptides)
            self.cal_num_table_data[run_file]["modified_peptide_num"] = len(modified_pep)

        self.quantms_modified["plot_data"] = mod_plot_dict
        self.quantms_modified["cats"] = list(sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x)))

        self.ms_without_psm = set([file_prefix(i) for i in self.ms_paths]) - set(
            self.ms_with_psm
        )
        for i in self.ms_without_psm:
            log.warning("No PSM found in '{}'!".format(i))

        for i in self.ms_without_psm:
            self.cal_num_table_data[i] = {
                "protein_num": 0,
                "peptide_num": 0,
                "unique_peptide_num": 0,
                "modified_peptide_num": 0,
            }

    def parse_msstats_input(self):
        log.info(
            "{}: Parsing MSstats input file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.msstats_input_path
            )
        )
        msstats_data = pd.read_csv(self.msstats_input_path)
        ## TODO we probably shouldn't even write out 0-intensity values to MSstats csv
        msstats_data = msstats_data[-(msstats_data["Intensity"] == 0)]
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
        cond_and_dist_cols = conditions_str + conditions_dists

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
            ## TODO How to determine technical replicates? Should be same BioReplicate but different Fraction_Group (but fraction group is not annotated)
            grouped = (
                g.groupby(["Condition", "BioReplicate"], as_index=False)["Intensity"]
                .mean()
            )
            cond_grp = (
                grouped
                .groupby("Condition", group_keys=False)[["BioReplicate", "Intensity"]]
                .apply(fill_dict)
            )

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

        del msstats_data
        ## TODO Can we guarantee that the score was always PEP? I don't think so!
        msstats_data_pep_agg.reset_index(inplace=True)
        msstats_data_pep_agg.index = msstats_data_pep_agg.index + 1
        msstats_data_dict_pep_full = msstats_data_pep_agg.to_dict("index")
        msstats_data_dict_pep_init = dict(itertools.islice(msstats_data_dict_pep_full.items(), 50))

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
            "save_file": False,
            "sort_rows": False,
            "only_defined_headers": True,
            "col1_header": "PeptideID",
            "no_violin": True,
        }

        # only use the first 50 lines for the table
        max_pep_intensity = 50
        table_html = table.plot(
            dict(itertools.islice(msstats_data_dict_pep_init.items(), max_pep_intensity)),
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
                """
        )

        # Helper functions for pandas
        def json_to_dict(s):
            if type(s) is str:
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

        max_prot_intensity = 0
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
                """
        )


    def cal_quantms_contaminant_percent(self, pep_df):

        group_stats = pep_df.groupby("stand_spectra_ref").agg(
            total_intensity=("average_intensity", "sum"),
            cont_intensity=("average_intensity", lambda x: x[pep_df["accession"].str.contains("CONT")].sum())
        )

        group_stats["contaminant_percent"] = (
            group_stats["cont_intensity"] / group_stats["total_intensity"] * 100
        )

        result_dict = dict()
        for k, v in dict(zip(group_stats.index, group_stats["contaminant_percent"])).items():
            result_dict[k] = {"Potential Contaminants": v}

        return result_dict

    def top_n_contaminant_percent(self, pep_df, top_n):

        not_cont_tag = "NOT_CONT"
        pep_df["cont_accession"] = pep_df["accession"].apply(
            lambda x: x.replace("CONTAMINANT_", "") if x.startswith("CONTAMINANT_") else not_cont_tag
        )

        pep_contaminant_df = pep_df[pep_df["cont_accession"] != not_cont_tag].copy()
        contaminant_df = (
            pep_contaminant_df
            .groupby("cont_accession", as_index=False)["average_intensity"]
            .sum()
            .sort_values(by="average_intensity", ascending=False)
        )

        top_contaminants = list(contaminant_df.head(top_n).cont_accession)

        plot_dict = dict()
        plot_cats = list()
        
        for file_name, group in pep_df.groupby("stand_spectra_ref"):

            contaminant_rows = group[group["cont_accession"] != not_cont_tag].copy()
            contaminant_rows.loc[~contaminant_rows["cont_accession"].isin(top_contaminants), "cont_accession"] = (
                "Other"
            )

            cont_df = (
                contaminant_rows
                .groupby("cont_accession", as_index=False)["average_intensity"]
                .sum()
                .sort_values(by="average_intensity", ascending=False)
                .reset_index(drop=True)
            )
            cont_df["contaminant_percent"] = (cont_df["average_intensity"] / group["average_intensity"].sum()) * 100

            plot_dict[file_name] = dict(zip(cont_df["cont_accession"], cont_df["contaminant_percent"]))
            plot_cats.extend(cont_df["cont_accession"].tolist())

        plot_cats = list(set(plot_cats))
        if "Other" in plot_cats:
            plot_cats_ordered = [x for x in plot_cats if x != "Other"] + [x for x in plot_cats if x == "Other"]

        result_dict = dict()
        result_dict["plot_data"] = plot_dict
        result_dict["cats"] = plot_cats_ordered

        return result_dict


    def draw_quantms_identification(self, mzml_table):

        # 1.ProteinGroups Count
        draw_config = {
            "id": "protein_group_count",
            "cpswitch": False,
            "title": "ProteinGroups Count",
            "tt_decimals": 0,
            "ylab": "Count",
        }

        if self.cal_num_table_data:
            protein_count = {sample: {"Count": info["protein_num"]} for sample, info in self.cal_num_table_data.items()}
            peptide_count = {sample: {"Count": info["peptide_num"]} for sample, info in self.cal_num_table_data.items()}
        else:
            return

        bar_html = bargraph.plot(
            protein_count,
            pconfig=draw_config,
        )
        bar_html = common_plots.remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["identification"],
            plot=bar_html,
            order=3,
            description="Number of protein groups per raw file.",
            helptext="""
                Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
                """
        )

        # 2.Peptide ID Count
        draw_config = {
            "id": "peptide_id_count",
            "cpswitch": False,
            "title": "Peptide ID Count",
            "tt_decimals": 0,
            "ylab": "Count",
        }
        bar_html = bargraph.plot(
            peptide_count,
            pconfig=draw_config,
        )
        bar_html = common_plots.remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["identification"],
            plot=bar_html,
            order=4,
            description="""
                Number of unique (i.e. not counted twice) peptide sequences including modifications per Raw file.
                """,
            helptext="""
                Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
                """
        )

        # 3.Missed Cleavages Per Raw File
        if self.quantms_missed_cleavages:
            mc_group = {}
            for sample, counts in self.quantms_missed_cleavages.items():
                mc_group[sample] = {"0": 0, "1": 0, ">=2": 0}
                for mc, count in counts.items():
                    if mc == 0:
                        mc_group[sample]["0"] += count
                    elif mc == 1:
                        mc_group[sample]["1"] += count
                    else:
                        mc_group[sample][">=2"] += count

            mc_group_ratio = {}
            for sample, counts in mc_group.items():
                total = sum(counts.values())
                if total == 0:
                    mc_group_ratio[sample] = {"0": 0, "1": 0, ">=2": 0}
                    continue
                mc_group_ratio[sample] = {
                    group: count / total * 100
                    for group, count in counts.items()
                }

            mc_data = {
                "plot_data": mc_group_ratio,
                "cats": ["0", "1", ">=2"]
            }
            common_plots.draw_msms_missed_cleavages(
                self.sub_sections["identification"],
                mc_data,
                False
            )
    
        # 4.Modifications Per Raw File
        if self.quantms_modified:
            common_plots.draw_modifications(
                self.sub_sections["identification"],
                self.quantms_modified
            )

        # 5.MS/MS Identified per Raw File
        if self.identified_msms_spectra and mzml_table:

            msms_identified_rate = dict()

            for m in self.identified_msms_spectra.keys():
                identified_ms2 = self.identified_msms_spectra[m].get("Identified", 0)
                all_ms2 = mzml_table.get(m, {}).get("MS2_Num", 0)

                if all_ms2 > 0:
                    msms_identified_rate[m] = {
                        "Identified Rate": (identified_ms2 / all_ms2) * 100
                    }

            common_plots.draw_ms_ms_identified(
                self.sub_sections["identification"],
                msms_identified_rate
            )
    
    def draw_quantms_contaminants(self):

        # 1.Potential Contaminants per Group
        if self.quantms_contaminant_percent:
            common_plots.draw_potential_contaminants(
                self.sub_sections["contaminants"],
                self.quantms_contaminant_percent,
                False
            )

        # 2.Top5 Contaminants per Raw file
        if self.quantms_top_contaminant_percent:
            common_plots.draw_top_n_contaminants(
                self.sub_sections["contaminants"],
                self.quantms_top_contaminant_percent
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
            }
            box_html = box.plot(self.quantms_pep_intensity, pconfig=draw_config)
            box_html = common_plots.remove_subtitle(box_html)

            add_sub_section(
                sub_section=self.sub_sections["quantification"],
                plot=box_html,
                order=5,
                description="Peptide intensity per file from mzTab.",
                helptext="""
                    Calculate the average of peptide_abundance_study_variable[1-n] values for each peptide from the 
                    peptide table in the mzTab file, and then apply a log2 transformation.
                    """
            )

    def draw_quantms_msms_section(self):

        # 1.Charge-state of Per File
        if self.mztab_charge_state:
            common_plots.draw_charge_state(
                self.sub_sections["ms2"],
                self.mztab_charge_state,
                False
            )

    def draw_quantms_time_section(self):

        # 1.IDs over RT
        if self.quantms_ids_over_rt:
            common_plots.draw_ids_rt_count(
                self.sub_sections["time_mass"],
                self.quantms_ids_over_rt,
                ""
            )

        # 2.Delta Mass [ppm]
        if self.quantms_mass_error:
            common_plots.draw_delta_mass_da_ppm(
                self.sub_sections["time_mass"],
                self.quantms_mass_error,
                "quantms_ppm"
            )

def read_openms_design(desfile):
    with open(desfile, "r") as f:
        data = f.readlines()
        s_row = False
        f_table = []
        s_table = []
        for row in data:
            if row == "\n":
                continue
            if "MSstats_Condition" in row:
                s_row = True
                s_header = row.replace("\n", "").split("\t")
            elif s_row:
                s_table.append(row.replace("\n", "").split("\t"))
            elif "Spectra_Filepath" in row:
                f_header = row.replace("\n", "").split("\t")
            else:
                f_table.append(row.replace("\n", "").split("\t"))

        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table["Run"] = f_table.apply(
            lambda x: file_prefix(x["Spectra_Filepath"]), axis=1
        )
        s_data_frame = pd.DataFrame(s_table, columns=s_header)

    return s_data_frame, f_table

def find_modification(peptide):
    peptide = str(peptide)
    pattern = re.compile(r"\((.*?)\)")
    original_mods = pattern.findall(peptide)
    peptide = re.sub(r"\(.*?\)", ".", peptide)
    position = [i.start() for i in re.finditer(r"\.", peptide)]
    for j in range(1, len(position)):
        position[j] -= j

    for k in range(0, len(original_mods)):
        original_mods[k] = str(position[k]) + "-" + original_mods[k]

    original_mods = ",".join(str(i) for i in original_mods) if len(original_mods) > 0 else "nan"

    return AASequence.fromString(peptide).toUnmodifiedString(), original_mods

def condition_split(conditions):
    items = conditions.split(';')
    key_value_pairs = [item.split('=') for item in items if '=' in item]

    result_dict = {k.strip(): v.strip() for k, v in key_value_pairs}
    return result_dict
