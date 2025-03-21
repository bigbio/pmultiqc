#!/usr/bin/env python

""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
from collections import OrderedDict
import itertools
from datetime import datetime
from operator import itemgetter
import logging
from multiqc import config, report

from multiqc import BaseMultiqcModule
from sdrf_pipelines.openms.openms import OpenMS, UnimodDatabase
from multiqc.plots import table, bargraph, linegraph, heatmap, box, scatter
from multiqc.utils.mqc_colour import mqc_colour_scale

import pandas as pd
from functools import reduce
import re
from pyteomics import mztab, mzid, mgf
from pyopenms import IdXMLFile, MzMLFile, MSExperiment, OpenMSBuildInfo, AASequence
import os
import sqlite3
import numpy as np
import math
import copy
import json

from .histogram import Histogram
from . import sparklines
from .ms_functions import get_ms_qc_info
from . import maxquant

# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

if config.output_dir:
    if os.path.exists(config.output_dir):
        con = sqlite3.connect(os.path.join(config.output_dir, "quantms.db"))
    else:
        os.makedirs(config.output_dir)
        con = sqlite3.connect(os.path.join(config.output_dir, "quantms.db"))
else:
    con = sqlite3.connect("./quantms.db")

cur = con.cursor()
cur.execute("drop table if exists PROTQUANT")
con.commit()

cur.execute("drop table if exists PEPQUANT")
con.commit()

cur.execute("drop table if exists PSM")
con.commit()

log.info("pyopenms has: " + str(OpenMSBuildInfo().getOpenMPMaxNumThreads()) + " threads.")


class QuantMSModule(BaseMultiqcModule):
    @staticmethod
    def file_prefix(path):
        try:
            return os.path.splitext(os.path.basename(path))[0]
        except:
            raise SystemExit("Illegal file path: {path}")

    def __init__(self):

        # Initialise the parent module Class object
        super(QuantMSModule, self).__init__(
            name="pmultiqc",
            target="pmultiqc",
            anchor="pmultiqc",
            href="https://github.com/bigbio/pmultiqc",
            info=""" is a MultiQC module to show the pipeline performance of mass spectrometry based quantification 
            pipelines such as <a href='https://nf-co.re/quantms'>nf-core/quantms</a>, <a href='https://www.maxquant.org/'>MaxQuant</a>""",
        )

        # Halt execution if we've disabled the plugin
        if config.kwargs.get("disable_plugin", True):
            return None

        # Parse MaxQuant results
        if config.kwargs.get("parse_maxquant", False):

            self.maxquant_paths = self.maxquant_file_path()

            # parameters.txt
            if "parameters" in self.maxquant_paths.keys():
                log.info(
                    "{}: Parsing parameters file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
                    )
                )
                get_parameter_dicts = maxquant.get_parameters(
                    file_path=self.maxquant_paths["parameters"]
                )
                log.info(
                    "{}: Completed the processing of the parameters file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
                    )
                )
            else:
                get_parameter_dicts = {"parameters_tb_dict": None}

            # proteinGroups.txt
            if "proteinGroups" in self.maxquant_paths.keys():
                log.info(
                    "{}: Parsing proteinGroups file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["proteinGroups"]
                    )
                )
                get_protegroups_dicts = maxquant.get_protegroups(
                    file_path=self.maxquant_paths["proteinGroups"]
                )
                log.info(
                    "{}: Completed the processing of the proteinGroups file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["proteinGroups"]
                    )
                )
            else:
                get_protegroups_dicts = {
                    "pg_contaminant": None,
                    "pg_intensity_distri": None,
                    "pg_lfq_intensity_distri": None,
                    "raw_intensity_pca": None,
                    "lfq_intensity_pca": None,
                }

            # summary.txt
            if "summary" in self.maxquant_paths.keys():
                log.info(
                    "{}: Parsing summary file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
                    )
                )
                ms_ms_identified = maxquant.get_summary(file_path=self.maxquant_paths["summary"])
                log.info(
                    "{}: Completed the processing of the summary file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
                    )
                )
            else:
                ms_ms_identified = None

            # evidence.txt
            if "evidence" in self.maxquant_paths.keys():
                log.info(
                    "{}: Parsing evidence file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
                    )
                )
                get_evidence_dicts = maxquant.get_evidence(
                    file_path=self.maxquant_paths["evidence"]
                )
                log.info(
                    "{}: Completed the processing of the evidence file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
                    )
                )
            else:
                get_evidence_dicts = {
                    "top_contaminants": None,
                    "peptide_intensity": None,
                    "charge_counts": None,
                    "modified_percentage": None,
                    "rt_counts": None,
                    "evidence_df": None,
                    "peak_rt": None,
                    "oversampling": None,
                    "uncalibrated_mass_error": None,
                    "calibrated_mass_error": None,
                    "peptide_id_count": None,
                    "protein_group_count": None,
                }

            # msms.txt
            if "msms" in self.maxquant_paths.keys():
                log.info(
                    "{}: Parsing msms file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
                    )
                )
                get_msms_dicts = maxquant.get_msms(
                    file_path=self.maxquant_paths["msms"],
                    evidence_df=get_evidence_dicts["evidence_df"],
                )
                log.info(
                    "{}: Completed the processing of the msms file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
                    )
                )
            else:
                get_msms_dicts = {"missed_cleavages": None}

            # msScans.txt or msmsScans.txt
            if "msmsScans" in self.maxquant_paths.keys():
                msms_scans_file = "msmsScans"
            elif "msScans" in self.maxquant_paths.keys():
                msms_scans_file = "msScans"
            else:
                msms_scans_file = None
            if msms_scans_file:
                log.info(
                    "{}: Parsing msScans file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths[msms_scans_file]
                    )
                )
                get_msms_scans_dicts = maxquant.get_msms_scans(
                    file_path=self.maxquant_paths[msms_scans_file]
                )
                log.info(
                    "{}: Completed the processing of the msScans file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), self.maxquant_paths[msms_scans_file]
                    )
                )
            else:
                get_msms_scans_dicts = {
                    "ion_injec_time_RT": None,
                    "top_n": None,
                    "top_overRT": None,
                }

            # Parameters
            if get_parameter_dicts["parameters_tb_dict"]:
                self.draw_parameters(get_parameter_dicts["parameters_tb_dict"])

            # Intensity
            # self.draw_intensity_his(get_protegroups_dicts['pg_intensity_distri']['histogram'], fig_type='intensity')
            if get_protegroups_dicts["pg_intensity_distri"]:
                self.draw_intensity_box(
                    get_protegroups_dicts["pg_intensity_distri"]["box"], fig_type="intensity"
                )
            # self.draw_intensity_his(get_protegroups_dicts['pg_lfq_intensity_distri']['histogram'], fig_type='lfq_intensity')
            if get_protegroups_dicts["pg_lfq_intensity_distri"]:
                self.draw_intensity_box(
                    get_protegroups_dicts["pg_lfq_intensity_distri"]["box"],
                    fig_type="lfq_intensity",
                )
            if get_protegroups_dicts["raw_intensity_pca"]:
                self.draw_pg_pca(
                    get_protegroups_dicts["raw_intensity_pca"], fig_type="raw_intensity"
                )
            if get_protegroups_dicts["lfq_intensity_pca"]:
                self.draw_pg_pca(
                    get_protegroups_dicts["lfq_intensity_pca"], fig_type="lfq_intensity"
                )
            if ms_ms_identified:
                self.draw_summary(ms_ms_identified)
            # self.draw_intensity_his(get_evidence_dicts['peptide_intensity']['histogram'], fig_type='peptide_intensity')
            if get_evidence_dicts["peptide_intensity"]:
                self.draw_intensity_box(
                    get_evidence_dicts["peptide_intensity"]["box"], fig_type="peptide_intensity"
                )

            # Contaminants
            if get_protegroups_dicts["pg_contaminant"]:
                self.draw_pg_contaminants(get_protegroups_dicts["pg_contaminant"])
            if get_evidence_dicts["top_contaminants"]:
                self.draw_evidence_top_contaminants(get_evidence_dicts["top_contaminants"])

            if get_evidence_dicts["charge_counts"]:
                self.draw_evidence_charge(get_evidence_dicts["charge_counts"])
            if get_evidence_dicts["modified_percentage"]:
                self.draw_evidence_modified(get_evidence_dicts["modified_percentage"])

            if get_evidence_dicts["peptide_id_count"]:
                self.draw_evidence_peptide_id_count(get_evidence_dicts["peptide_id_count"])
            if get_evidence_dicts["protein_group_count"]:
                self.draw_evidence_protein_group_count(get_evidence_dicts["protein_group_count"])

            if get_evidence_dicts["oversampling"]:
                self.draw_evidence_oversampling(get_evidence_dicts["oversampling"])

            if get_msms_dicts["missed_cleavages"]:
                self.draw_msms_missed_cleavages(get_msms_dicts["missed_cleavages"])

            if get_evidence_dicts["rt_counts"]:
                self.draw_evidence_rt_count(get_evidence_dicts["rt_counts"])
            if get_evidence_dicts["peak_rt"]:
                self.draw_evidence_peak_width_rt(get_evidence_dicts["peak_rt"])

            # Uncalibrated/Calibrated Mass Error
            if get_evidence_dicts["uncalibrated_mass_error"]:
                self.draw_mass_error_box(
                    get_evidence_dicts["uncalibrated_mass_error"], fig_type="uncalibrated"
                )
            if get_evidence_dicts["calibrated_mass_error"]:
                self.draw_mass_error_box(
                    get_evidence_dicts["calibrated_mass_error"], fig_type="calibrated"
                )

            # TopN
            if get_msms_scans_dicts["top_n"]:
                self.draw_msms_scans_top_n(get_msms_scans_dicts["top_n"])
            if get_msms_scans_dicts["top_overRT"]:
                self.draw_msms_scans_top_over_rt(get_msms_scans_dicts["top_overRT"])
            if get_msms_scans_dicts["ion_injec_time_RT"]:
                self.draw_msms_scans_ion_injec_time_rt(get_msms_scans_dicts["ion_injec_time_RT"])

            return None

        self.ms_with_psm = list()
        self.Total_Protein_Identified = 0
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
        self.MissedCleavages_heatmap_score = dict()
        self.ID_RT_score = dict()
        self.HeatmapOverSamplingScore = dict()
        self.HeatmapPepMissingScore = dict()
        self.MissedCleavagesVar_score = dict()
        self.pep_table_exists = False
        self.ms_info = dict()
        self.ms_info["charge_distribution"] = dict()
        self.ms_info["peaks_per_ms2"] = dict()
        self.ms_info["peak_distribution"] = dict()

        self.ms_paths = []
        for mzML_file in self.find_log_files("quantms/mzML", filecontents=False):
            self.ms_paths.append(os.path.join(mzML_file["root"], mzML_file["fn"]))

        self.ms_info_path = []

        # Check if this is a display for mzid identification results
        if config.kwargs.get("mzid_plugin", False):

            self.mzid_peptide_map = dict()
            self.ms_without_psm = dict()

            self.mgf_paths = []
            for mgf_file in self.find_log_files("quantms/mgf", filecontents=False):
                self.mgf_paths.append(os.path.join(mgf_file["root"], mgf_file["fn"]))

            self.mzid_paths = []
            for mzid_file in self.find_log_files("quantms/mzid", filecontents=False):
                self.mzid_paths.append(os.path.join(mzid_file["root"], mzid_file["fn"]))

            mzid_psm = self.parse_out_mzid()

            if self.mgf_paths:
                self.parse_out_mgf()
            elif self.ms_paths:
                _ = self.parse_mzml()

            self.mzid_cal_heat_map_score(mzid_psm)
            self.draw_heatmap()
            self.draw_summary_protein_ident_table()
            self.draw_mzid_identi_num()
            self.draw_num_pep_per_protein()
            self.draw_precursor_charge_distribution()
            self.draw_peaks_per_ms2()
            self.draw_peak_intensity_distribution()
            self.draw_oversampling()

        else:
            self.enable_exp = False
            self.enable_sdrf = False
            self.msstats_input_valid = False
            # TODO what if multiple are found??
            for f in self.find_log_files("quantms/exp_design", filecontents=False):
                self.exp_design = os.path.join(f["root"], f["fn"])
                self.enable_exp = True

            if not self.enable_exp:
                for f in self.find_log_files("quantms/sdrf", filecontents=False):
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

            # TODO in theory the module can work without the design. We just need to remove the according sections!
            if not self.enable_sdrf and not self.enable_exp:
                raise AttributeError(
                    "Neither exp_design nor sdrf can be found! Please provide or correct your multiqc_config.yml."
                )

            self.PSM_table = dict()
            self.mzml_peptide_map = dict()
            self.pep_quant_table = dict()
            self.mzml_table = OrderedDict()
            self.search_engine = OrderedDict()
            self.XCORR_HIST_RANGE = {"start": 0, "end": 5, "step": 0.1}
            self.HYPER_HIST_RANGE = {"start": 0, "end": 5, "step": 0.1}
            self.SPECEVALUE_HIST_RANGE = {"start": 0, "end": 20, "step": 0.4}
            self.PEP_HIST_RANGE = {"start": 0, "end": 1, "step": 0.02}
            self.Total_Protein_Quantified = 0
            self.out_csv_data = dict()
            self.mL_spec_ident_final = dict()
            self.delta_mass_percent = dict()
            self.heatmap_con_score = dict()
            self.heatmap_pep_intensity = {}
            self.exp_design_table = dict()
            self.ms1_tic = dict()
            self.ms1_bpc = dict()
            self.ms1_peaks = dict()
            self.ms1_general_stats = dict()
            # parse input data
            # draw the experimental design
            self.draw_exp_design()

            for ms_info in self.find_log_files("quantms/ms_info", filecontents=False):
                self.ms_info_path.append(os.path.join(ms_info["root"], ms_info["fn"]))
                self.ms_info_path.sort()
                if len(self.ms_info_path) > 0:
                    self.read_ms_info = True
                    self.ms_paths = [
                        self.file_prefix(i).replace("_ms_info", ".mzML") for i in self.ms_info_path
                    ]

            for f in self.find_log_files("quantms/mztab", filecontents=False):
                self.out_mzTab_path = os.path.join(f["root"], f["fn"])
                self.parse_out_mzTab()

            for report in self.find_log_files("quantms/diann_report", filecontents=False):
                self.diann_report_path = os.path.join(report["root"], report["fn"])
                self.enable_dia = True

            mt = self.parse_mzml()
            self.idx_paths = []
            for idx_file in self.find_log_files("quantms/idXML", filecontents=False):
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
                self.CalHeatMapScore()
                self.draw_heatmap()
                self.draw_summary_protein_ident_table()
                self.draw_quantms_identi_num()
                self.draw_num_pep_per_protein()
                if not config.kwargs["ignored_idxml"]:
                    self.draw_mzml_ms()
                    self.draw_search_engine()
                self.draw_precursor_charge_distribution()
                self.draw_peaks_per_ms2()
                self.draw_peak_intensity_distribution()
                self.draw_oversampling()

            # TODO what if multiple are found??
            # if config.kwargs.get('disable_table', True):
            #     log.info("Skip protein/peptide table plot!")
            # else:
            for msstats_input in self.find_log_files("quantms/msstats", filecontents=False):
                self.msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])
                self.msstats_input_valid = True
                self.parse_msstats_input()

            if config.kwargs["quantification_method"] == "spectral_counting":
                # Add a report section with psm table plot from mzTab for spectral counting
                self.add_section(
                    name="Peptide Spectrum Matches",
                    anchor="peptide_spectrum_matches",
                    description="This plot shows the information of peptide spectrum matches",
                    helptext="""
                            This table shows the information of peptide spectrum matches from mzTab PSM section.
                            """,
                    plot=self.psm_table_html,
                )

            # TODO draw protein quantification from mzTab in the future with Protein and peptide tables from mzTab
            # currently only draw protein tabel for spectral counting
            if (
                    not self.msstats_input_valid
                    and config.kwargs["quantification_method"] == "spectral_counting"
            ):
                log.warning("MSstats input file not found!")
                self.add_section(
                    name="Protein Quantification Table",
                    anchor="protein_quant_result",
                    description="This plot shows the quantification information of proteins"
                                " in the final result (mainly the mzTab file).",
                    helptext="""
                            The quantification information (Spectral Counting) of proteins is obtained from the mzTab file. 
                            The table shows the quantitative level and distribution of proteins in different study variables and run.

                            * Peptides_Number: The number of peptides for each protein.
                            * Average Spectral Counting: Average spectral counting of each protein across all conditions with NA=0 or NA ignored.
                            * Spectral Counting in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Average spectral counting of replicates.

                            Click `Show replicates` to switch to bar plots of counting in each replicate.
                            """,
                    plot=self.protein_quantification_table_html,
                )

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

    def draw_heatmap(self):
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
                            self.MissedCleavages_heatmap_score[k],
                            self.MissedCleavagesVar_score[k],
                            self.ID_RT_score[k],
                            self.HeatmapOverSamplingScore[k],
                            self.HeatmapPepMissingScore[k],
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
                            self.MissedCleavages_heatmap_score[k],
                            self.MissedCleavagesVar_score[k],
                            self.ID_RT_score[k],
                            self.HeatmapOverSamplingScore[k],
                            self.HeatmapPepMissingScore[k],
                        ]
                    )

        pconfig = {
            "id": "HeatMap",
            "title": "Performance Overview",  # Plot title - should be in format "Module Name: Plot Title"
            "xlab": "metrics",  # X-axis title
            "ylab": "RawName",  # Y-axis title
            "square": False,
        }

        hm_html = heatmap.plot(heat_map_score, xnames, ynames, pconfig)
        # Add a report section with the heatmap plot
        self.add_section(
            name="HeatMap",
            anchor="quantms_heatmap",
            description="This heatmap shows a performance overview of the pipeline",
            helptext="""
                    This plot shows the pipeline performance overview. Some metrics are calculated.
                    
                    * Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                        1 = no contaminants
                    * Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                        i.e. reaching 2^21 of 2^23 gives score 0.25.
                    * Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative 
                        Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at 
                        tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), 
                        in-source fragmentation, missed cleavages and buffer irregularities can cause a shift 
                        (see Bittremieux 2017, DOI: 10.1002/mas.21544 ).
                    * Heatmap score [MC]: the fraction (0% - 100%) of fully cleaved peptides per Raw file
                    * Heatmap score [MC Var]: each Raw file is scored for its deviation from the ‘average’ digestion 
                        state of the current study.
                    * Heatmap score [ID rate over RT]: Judge column occupancy over retention time. 
                        Ideally, the LC gradient is chosen such that the number of identifications 
                        (here, after FDR filtering) is uniform over time, to ensure consistent instrument duty cycles. 
                        Sharp peaks and uneven distribution of identifications over time indicate potential for LC gradient 
                        optimization.Scored using ‘Uniform’ scoring function. i.e. constant receives good score, extreme shapes are bad.
                    * Heatmap score [MS2 Oversampling]: The percentage of non-oversampled 3D-peaks. An oversampled 
                        3D-peak is defined as a peak whose peptide ion (same sequence and same charge state) was 
                        identified by at least two distinct MS2 spectra in the same Raw file. For high complexity samples, 
                        oversampling of individual 3D-peaks automatically leads to undersampling or even omission of other 3D-peaks, 
                        reducing the number of identified peptides.
                    * Heatmap score [Pep Missing]: Linear scale of the fraction of missing peptides.
                    
                    """,
            plot=hm_html,
        )

    def draw_exp_design(self):
        # Currently this only supports the OpenMS two-table format (default in quantms pipeline)
        # One table format would actually be even easier. You can just use pandas.read_tsv
        self.sample_df, self.file_df = read_openms_design(self.exp_design)

        if self.file_df["Spectra_Filepath"][0].endswith((".d", ".d.tar")):
            self.is_bruker = True

        for file in np.unique(self.file_df["Run"].tolist()):
            file_index = self.file_df[self.file_df["Run"] == file]
            self.exp_design_table[file] = {
                "Fraction_Group": file_index["Fraction_Group"].tolist()[0]
            }
            self.exp_design_table[file]["Fraction"] = file_index["Fraction"].tolist()[0]
            self.exp_design_table[file]["Label"] = "|".join(file_index["Label"])
            sample = file_index["Sample"].tolist()
            self.exp_design_table[file]["Sample"] = "|".join(sample)
            self.exp_design_table[file]["MSstats_Condition"] = ",".join(
                [
                    row["MSstats_Condition"]
                    for _, row in self.sample_df.iterrows()
                    if row["Sample"] in sample
                ]
            )
            self.exp_design_table[file]["MSstats_BioReplicate"] = "|".join(
                [
                    row["MSstats_BioReplicate"]
                    for _, row in self.sample_df.iterrows()
                    if row["Sample"] in sample
                ]
            )

        # Create table plot
        pconfig = {
            "id": "Experimental_Design",  # ID used for the table
            "title": "Experimental Design",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_Experimental_Design_table",  # File basename to use for raw data file
            "sort_rows": False,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            "no_violin": True,
            # 'format': '{:,.0f}',  # The header used for the first column
        }
        headers = OrderedDict()
        set3_scale = mqc_colour_scale(name="Set3")
        maxnr = len(self.file_df.index)
        set3_colors = set3_scale.get_colours(name="Set3")
        colors = dict((str(i + 1), set3_colors[i % len(set3_colors)]) for i in range(maxnr))

        headers["Fraction_Group"] = {
            "description": "Fraction_Group",
            "bgcols": colors,
        }
        headers["Fraction"] = {
            "description": "Fraction identifier",
            "bgcols": colors,
        }
        headers["Label"] = {
            "description": "Label",
            "bgcols": colors,
        }
        headers["Sample"] = {
            "description": "Sample Name",
            "bgcols": colors,
        }
        headers["MSstats_Condition"] = {
            "description": "MSstats Condition",
            "bgcols": colors,
        }
        headers["MSstats_BioReplicate"] = {
            "description": "MSstats BioReplicate",
            "bgcols": colors,
        }
        table_html = table.plot(self.exp_design_table, headers, pconfig)
        # table_html = sparklines.plot(self.exp_design_table, headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
            name="Experimental Design",
            anchor="exp_design",
            description="This table shows the design of the experiment. I.e., which files and channels correspond to which sample/condition/fraction.",
            helptext="""
            You can see details about it in 
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            """,
            plot=table_html,
        )

    def draw_summary_protein_ident_table(self):
        headers = OrderedDict()
        if self.enable_dia:
            summary_table = {
                self.Total_Peptide_Count: {"#Proteins Quantified": self.Total_Protein_Quantified}
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
            ] = self.Total_Protein_Identified

            if not config.kwargs.get("mzid_plugin", False):
                summary_table[self.total_ms2_spectra][
                    "#Proteins Quantified"
                ] = self.Total_Protein_Quantified

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
            "id": "identification summary table",  # ID used for the table
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
            description_str = "This plot shows the summary statistics of the submitted data"
            helptext_str = """
                        This plot shows the summary statistics of the submitted data
                        """
        else:
            description_str = "This table shows the quantms pipeline summary statistics"
            helptext_str = """
                        This table shows the quantms pipeline summary statistics
                        """

        # Add a report section with the table
        self.add_section(
            name="Summary Table",
            anchor="quantms_summary_table",
            description=description_str,
            helptext=helptext_str,
            plot=table_html,
        )

    def draw_ms_information(self):
        ms1_tic_config = {
            "id": "ms1_tic",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "Total Ion Chromatograms",
            "ylab": "Ion Count",
            "ymin": 0,
        }

        ms1_tic_html = linegraph.plot(self.ms1_tic, ms1_tic_config)

        ms1_bpc_config = {
            "id": "ms1_bpc",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "MS1 BPC",
            "ylab": "Ion Count",
            "ymin": 0,
        }

        ms1_bpc_html = linegraph.plot(self.ms1_bpc, ms1_bpc_config)

        ms1_peaks_config = {
            "id": "ms1_peaks",
            "tt_label": "<b>{point.x} Peak Count:</b> {point.y}",
            "title": "MS1 Peaks",
            "ylab": "Peak Count",
            "ymin": 0,
        }

        ms1_peaks_html = linegraph.plot(self.ms1_peaks, ms1_peaks_config)

        tconfig = {
            "id": "ms_general_stats",
            "title": "General Stats",
            "only_defined_headers": False,
            "col1_header": "File",
        }
        table_html = table.plot(self.ms1_general_stats, pconfig=tconfig)

        self.add_section(
            name="MS1 Information",
            anchor="ms1_information",
            description="""MS1 quality control information extracted from the spectrum files.""",
            plot=ms1_tic_html,
        )

        self.add_section(
            description="""#### MS1 base peak chromatograms extracted from the spectrum files
            """,
            plot=ms1_bpc_html,
        )
        self.add_section(
            description="""#### MS1 Peaks from the spectrum files
            """,
            plot=ms1_peaks_html,
        )
        self.add_section(
            description="""#### General stats for MS1 information from .d files
            """,
            plot=table_html,
        )

    def draw_quantms_identi_num(self):
        # Create table plot
        pconfig = {
            "id": "result statistics",  # ID used for the table
            "title": "pipeline result statistics",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_result statistics_table",  # File basename to use for raw data file
            "sort_rows": True,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            "no_violin": True,
            # 'format': '{:,.0f}'  # The header used for the first column
        }
        colors = dict((str(i + 1), "#ffffff") for i in range(len(self.file_df.index)))
        headers = OrderedDict()

        ## TODO BIG! We desperately need rules for naming columns!

        headers["sample_name"] = {
            "title": "Sample Name",
            "description": "Sample identifier",
            "bgcols": colors,
        }
        headers["condition"] = {
            "title": "Condition",
            "description": "Combination of possible study variables",
            "bgcols": colors,
        }
        headers["fraction"] = {
            "title": "Fraction",
            "description": "Fraction identifier",
            "bgcols": colors,
        }
        headers["peptide_num"] = {
            "title": "#Peptide IDs",
            "description": "The number of identified PSMs in the pipeline",
            "color": "#ffffff",
        }
        headers["unique_peptide_num"] = {
            "title": "#Unambiguous Peptide IDs",
            "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database.",
            "color": "#ffffff",
        }
        headers["modified_peptide_num"] = {
            "title": "#Modified Peptide IDs",
            "description": "Number of modified identified peptides in the pipeline",
            "color": "#ffffff",
        }
        headers["protein_num"] = {
            "title": "#Protein (group) IDs",
            "description": "The number of identified protein(group)s in the pipeline",
            "color": "#ffffff",
        }
        table_html = table.plot(self.cal_num_table_data, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Pipeline Result Statistics",
            anchor="Pipeline Result Statistics",
            description="This plot shows the quantms pipeline final result",
            helptext="""
            This plot shows the quantms pipeline final result.
            Including Sample Name、Possible Study Variables、identified the number of peptide in the pipeline、
            and identified the number of modified peptide in the pipeline, eg. All data in this table are obtained 
            from the out_msstats file. You can also remove the decoy with the `remove_decoy` parameter.
            """,
            plot=table_html,
        )

    def draw_mzid_identi_num(self):
        pconfig = {
            "id": "result statistics",  # ID used for the table
            "title": "pipeline result statistics",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_result statistics_table",  # File basename to use for raw data file
            "sort_rows": True,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            "no_violin": True,
            # 'format': '{:,.0f}'  # The header used for the first column
        }

        headers = OrderedDict()
        headers["peptide_num"] = {
            "title": "#Peptide IDs",
            "description": "The number of identified PSMs in the pipeline",
            "color": "#ffffff",
        }
        headers["unique_peptide_num"] = {
            "title": "#Unambiguous Peptide IDs",
            "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database.",
            "color": "#ffffff",
        }
        headers["modified_peptide_num"] = {
            "title": "#Modified Peptide IDs",
            "description": "Number of modified identified peptides in the pipeline",
            "color": "#ffffff",
        }
        headers["protein_num"] = {
            "title": "#Protein (group) IDs",
            "description": "The number of identified protein(group)s in the pipeline",
            "color": "#ffffff",
        }

        table_html = table.plot(self.cal_num_table_data, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Submitted Result Statistics",
            anchor="Submitted Result Statistics",
            description="This plot shows the submitted results",
            helptext="""
            This plot shows the submitted results.
            Including the number of identified peptides and the number of identified modified peptides in the submitted results. 
            You can also remove the decoy with the `remove_decoy` parameter.
            """,
            plot=table_html,
        )

    # draw number of peptides per proteins
    def draw_num_pep_per_protein(self):
        # Create table plot
        if any([len(i) >= 100 for i in self.pep_plot.dict["data"].values()]):
            data_labels = ["Frequency", "Percentage"]
        else:
            data_labels = [
                {
                    "name": "Frequency",
                    "ylab": "Frequency",
                    "title": "Number of Peptides identified Per Protein",
                },
                {"name": "Percentage", "ylab": "Percentage(%)"},
            ]
        pconfig = {
            "id": "number_of_peptides_per_proteins",  # ID used for the table
            "cpswitch": False,
            "title": "Number of Peptides identified per Protein",
            "xlab": "Number of Peptides",
            "tt_suffix": "",
            "tt_decimals": 2,
            "data_labels": data_labels,
        }
        bar_html = bargraph.plot(
            [self.pep_plot.dict["data"]["frequency"], self.pep_plot.dict["data"]["percentage"]],
            ["Frequency", "Percentage"],
            pconfig,
        )

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

        # Add a report section with the line plot
        self.add_section(
            name="Number of Peptides identified Per Protein",
            anchor="number_of_peptides_per_proteins",
            description=description_str,
            helptext=helptext_str,
            plot=bar_html,
        )

    def draw_mzml_ms(self):

        pconfig = {
            "id": "mzML_tracking",  # ID used for the table
            "title": "Pipeline spectrum tracking",  # Title of the table. Used in the column config modal
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

        # Add a report section with the line plot
        self.add_section(
            name="Spectra Tracking",
            anchor="spectra_tracking",
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
            plot=table_html,
        )

    def draw_peak_intensity_distribution(self):
        pconfig = {
            "id": "peak_intensity_distribution",
            "cpswitch": False,
            "title": "Peak Intensity Distribution",
            "ylab": "Count",
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_peak_distribution_plot.dict["cats"]
        else:
            cats = self.mzml_peak_distribution_plot.dict["cats"]
        bar_html = bargraph.plot(self.ms_info["peak_distribution"], cats, pconfig)
        self.add_section(
            name="Peak Intensity Distribution",
            anchor="Peak Intensity Distribution",
            description="""This is a histogram representing the ion intensity vs.
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
            plot=bar_html,
        )

    def draw_precursor_charge_distribution(self):
        pconfig = {
            "id": "Precursor_Ion_Charge_Distribution",
            "cpswitch": False,
            "title": "Precursor Ion Charge Distribution",
            "ylab": "Count",
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_charge_plot.dict["cats"]
        else:
            cats = self.mzml_charge_plot.dict["cats"]
        bar_html = bargraph.plot(self.ms_info["charge_distribution"], cats, pconfig)
        self.add_section(
            name="Distribution of precursor charges",
            anchor="Distribution of precursor charges",
            description="""This is a bar chart representing the distribution of the precursor ion charges 
                        for a given whole experiment. 
                        """,
            helptext="""This information can be used to identify potential ionization problems 
                        including many 1+ charges from an ESI ionization source or an unexpected 
                        distribution of charges. MALDI experiments are expected to contain almost exclusively 1+ 
                        charged ions. An unexpected charge distribution may furthermore be caused by specific search 
                        engine parameter settings such as limiting the search to specific ion charges.
                    """,
            plot=bar_html,
        )

    def draw_peaks_per_ms2(self):
        pconfig = {
            "id": "peaks_per_ms2",
            "cpswitch": False,
            "title": "Number of Peaks per MS/MS spectrum",
            "ylab": "Count",
        }
        if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
            cats = self.mgf_peaks_ms2_plot.dict["cats"]
        else:
            cats = self.mzml_peaks_ms2_plot.dict["cats"]
        bar_html = bargraph.plot(self.ms_info["peaks_per_ms2"], cats, pconfig)
        self.add_section(
            name="Number of Peaks per MS/MS spectrum",
            anchor="Number of Peaks per MS/MS spectrum",
            description="""This chart represents a histogram containing the number of peaks per MS/MS spectrum 
            in a given experiment. This chart assumes centroid data. Too few peaks can identify poor fragmentation 
            or a detector fault, as opposed to a large number of peaks representing very noisy spectra. 
            This chart is extensively dependent on the pre-processing steps performed to the spectra 
            (centroiding, deconvolution, peak picking approach, etc).
                        """,
            helptext="""
                    """,
            plot=bar_html,
        )

    def draw_oversampling(self):
        # Create bar plot
        pconfig = {
            "id": "Oversampling_Distribution",
            "cpswitch": False,
            "title": "MS2 counts per 3D-peak",
            "ylab": "Count",
        }
        cats = self.oversampling_plot.dict["cats"]
        bar_html = bargraph.plot(self.oversampling, cats, pconfig)
        # Add a report section with the bar plot
        self.add_section(
            name="Oversampling Distribution",
            anchor="Oversampling (MS/MS counts per 3D-peak)",
            description="""An oversampled 3D-peak is defined as a peak whose peptide ion 
                (same sequence and same charge state) was identified by at least two distinct MS2 spectra 
                in the same Raw file. 
                                """,
            helptext="""For high complexity samples, oversampling of individual 3D-peaks automatically leads to 
                undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. 
                Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic 
                exclusion windows for data independent acquisitions.
                
                * Heatmap score [EVD: MS2 Oversampling]: The percentage of non-oversampled 3D-peaks.
                            """,
            plot=bar_html,
        )

    def draw_delta_mass(self):

        self.delta_mass_percent["target"] = dict(
            zip(
                self.delta_mass["target"].keys(),
                list(
                    map(
                        lambda v: v / len(self.delta_mass["target"]),
                        self.delta_mass["target"].values(),
                    )
                ),
            )
        )
        if "decoy" in self.delta_mass.keys():
            self.delta_mass_percent["decoy"] = dict(
                zip(
                    self.delta_mass["decoy"].keys(),
                    list(
                        map(
                            lambda v: v / len(self.delta_mass["decoy"]),
                            self.delta_mass["decoy"].values(),
                        )
                    ),
                )
            )
            lineconfig = {
                # Building the plot
                "id": "delta_mass",  # HTML ID used for plot
                "tt_label": "<b>{point.x} Delta mass error:</b> {point.y}",
                # Plot configuration
                "title": "Delta m/z",  # Plot title - should be in format "Module Name: Plot Title"
                "xlab": "Experimental m/z - Theoretical m/z",  # X axis label
                "ylab": "Relative Frequency",  # Y axis label
                "colors": {"target": "#b2df8a", "decoy": "#DC143C"},
                "xmax": max(
                    list(self.delta_mass["decoy"].keys())
                    + (list(self.delta_mass["target"].keys()))
                )
                        + 0.01,
                "xmin": min(
                    list(self.delta_mass["decoy"].keys())
                    + (list(self.delta_mass["target"].keys()))
                )
                        - 0.01,
                "data_labels": [
                    {
                        "name": "Counts",
                        "ylab": "Count",
                        "tt_label": "<b>{point.x} Mass delta counts</b>: {point.y}",
                    },
                    {
                        "name": "Relative Frequency",
                        "ylab": "Relative Frequency",
                        "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}",
                    },
                ],
            }
        else:
            lineconfig = {
                # Building the plot
                "id": "delta_mass",  # HTML ID used for plot
                "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}",
                # Plot configuration
                "title": "Delta m/z",  # Plot title - should be in format "Module Name: Plot Title"
                "xlab": "Experimental m/z - Theoretical m/z",  # X axis label
                "ylab": "Relative Frequency",  # Y axis label
                "colors": {"target": "#b2df8a"},
                "xmax": max(list(self.delta_mass["target"].keys())) + 0.01,
                "xmin": min((list(self.delta_mass["target"].keys()))) - 0.01,
                "data_labels": [
                    {
                        "name": "Counts",
                        "ylab": "Count",
                        "tt_label": "<b>{point.x} Mass delta counts</b>: {point.y}",
                    },
                    {
                        "name": "Relative Frequency",
                        "ylab": "Relative Frequency",
                        "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}",
                    },
                ],
            }
        line_html = linegraph.plot([self.delta_mass, self.delta_mass_percent], lineconfig)

        self.add_section(
            name="Delta Mass",
            anchor="delta_mass",
            description="""This chart represents the distribution of the relative frequency of experimental 
                            precursor ion mass (m/z) - theoretical precursor ion mass (m/z). 
                            """,
            helptext="""
                    Mass deltas close to zero reflect more accurate identifications and also 
                    that the reporting of the amino acid modifications and charges have been done accurately. 
                    This plot can highlight systematic bias if not centered on zero. 
                    Other distributions can reflect modifications not being reported properly. 
                    Also it is easy to see the different between the target and the decoys identifications.
                    """,
            plot=line_html,
        )

    def draw_search_engine(self):
        # Add a report section with multiple histograms
        self.add_section(
            name="Summary of Search Engine Scores",
            anchor="search_engine_summary",
            description="These plots contain search scores and PEPs counts for different search engines in different files, "
                        "and they also contain a summary of the consensus PSMs if two or more search engines are used",
            helptext="""
                        This statistic is extracted from idXML files. 
                    """,
        )
        # Create scores summary plot
        [MSGF_labels, Comet_labels, Sage_labels] = self.search_engine["data_label"]["score_label"]

        spec_e_pconfig = {
            "id": "search_scores_summary",  # ID used for the table
            "cpswitch": False,
            "title": "Summary of Spectral E-values",
            "xlab": "MSGF -lg(SpecEvalue) ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "%",
            "tt_decimals": 0,
            "data_labels": MSGF_labels,
        }

        xcorr_pconfig = {
            "id": "search_scores_summary",  # ID used for the table
            "cpswitch": False,
            "title": "Summary of cross-correlation scores",
            "xlab": "Comet xcorr ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "%",
            "tt_decimals": 0,
            "data_labels": Comet_labels,
        }

        hyper_pconfig = {
            "id": "search_scores_summary",  # ID used for the table
            "cpswitch": False,
            "title": "Summary of Hyperscore",
            "xlab": "Sage hyperscore ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "%",
            "tt_decimals": 0,
            "data_labels": Sage_labels,
        }

        bar_cats = OrderedDict()
        bar_cats["target"] = {"name": "target", "color": "#2b908f"}
        bar_cats["decoy"] = {"name": "decoy", "color": "#90ed7d"}
        bar_cats["target+decoy"] = {"name": "target+decoy", "color": "#434348"}

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
            self.add_section(
                description="""#### SpecEvalue Description
                SpecEvalue : Spectral E-values, the search score of MSGF. The value used for plotting is -lg(SpecEvalue).
                """,
                plot=spec_e_bar_html,
            )

        if xcorr_bar_html != "":
            self.add_section(
                description="""#### xcorr description
                xcorr : cross-correlation scores, the search score of Comet. The value used for plotting is xcorr.
                """,
                plot=xcorr_bar_html,
            )

        if hyper_bar_html != "":
            self.add_section(
                description="""#### hyperscore description
                hyperscore : Hyperscore, the search score of Sage. The value used for plotting is hyperscore.
                """,
                plot=hyper_bar_html,
            )

        # Create PEPs summary plot
        pep_pconfig = {
            "id": "search_engine_PEP",  # ID used for the table
            "cpswitch": False,
            "title": "Summary of Search Engine PEP",
            "xlab": "PEP ranges",
            "stacking": "normal",
            "height": 550,
            "tt_suffix": "%",
            "tt_decimals": 0,
            "data_labels": self.search_engine["data_label"]["PEPs_label"],
        }

        pep_bar_html = bargraph.plot(
            list(self.search_engine["PEPs"].values()), pep_cats, pep_pconfig
        )

        self.add_section(
            description="""#### Summary of Posterior Error Probabilities (PEP)""",
            plot=pep_bar_html,
        )
        # Create identified number plot
        if len(self.search_engine["data_label"]["consensus_label"]) != 0:
            consensus_support_cats = [bar_cats] * len(self.search_engine["consensus_support"])
            consensus_pconfig = {
                "id": "consensus_summary",  # ID used for the table
                "cpswitch": False,
                "title": "Consensus Across Search Engines",
                "stacking": "normal",
                "height": 256,
                "tt_suffix": "%",
                "tt_decimals": 0,
                "data_labels": self.search_engine["data_label"]["consensus_label"],
            }

            consensus_bar_html = bargraph.plot(
                list(self.search_engine["consensus_support"].values()),
                consensus_support_cats,
                consensus_pconfig,
            )

            self.add_section(
                description="""#### Summary of consensus support for PSMs
                Consensus support is a measure of agreement between search engines. Every peptide sequence in the analysis has been
                identified by at least one search run. The consensus support defines which fraction (between 0 and 1) of the remaining
                search runs "supported" a peptide identification that was kept. The meaning of "support" differs slightly between
                algorithms: For best, worst, average and rank, each search run supports peptides that it has also identified among its
                top considered_hits candidates. So the "consensus support" simply gives the fraction of additional search engines that
                have identified a peptide. (For example, if there are three search runs, peptides identified by two of them will have a
                "support" of 0.5.) For the similarity-based algorithms PEPMatrix and PEPIons, the "support" for a peptide is the average
                similarity of the most-similar peptide from each (other) search run.
                """,
                plot=consensus_bar_html,
            )
        else:
            self.add_section(
                description="""#### Summary of consensus PSMs
                No Consensus PSMs data because of single search engine!
                """
            )

    def mzid_cal_heat_map_score(self, psm):
        log.info("Calculating Heatmap Scores...")

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
            lambda x: self.cal_MissCleavages(x["PeptideSequence"], enzyme), axis=1
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
            mis_0 = sc.get(0, 0)
            self.MissedCleavages_heatmap_score[name] = mis_0 / sc[:].sum()

            x = group["retention_time"] / np.sum(group["retention_time"])
            n = len(group["retention_time"])
            y = np.sum(x) / n
            worst = ((1 - y) ** 0.5) * 1 / n + (y ** 0.5) * (n - 1) / n
            sc = np.sum(np.abs(x - y) ** 0.5) / n
            if worst == 0:
                self.ID_RT_score[name] = 1.0
            else:
                self.ID_RT_score[name] = float((worst - sc) / worst)

            #  For HeatMapOverSamplingScore
            self.HeatmapOverSamplingScore[name] = self.oversampling[name]["1"] / np.sum(
                list(self.oversampling[name].values())
            )

            # For HeatMapPepMissingScore
            idFraction = (
                    len(
                        pd.merge(
                            global_peps,
                            group[["PeptideSequence", "Modifications"]].drop_duplicates(),
                            on=["PeptideSequence", "Modifications"],
                        ).drop_duplicates()
                    )
                    / global_peps_count
            )
            self.HeatmapPepMissingScore[name] = np.minimum(1.0, idFraction)

        median = np.median(list(self.MissedCleavages_heatmap_score.values()))
        self.MissedCleavagesVar_score = dict(
            zip(
                self.MissedCleavages_heatmap_score.keys(),
                list(
                    map(
                        lambda v: 1 - np.abs(v - median),
                        self.MissedCleavages_heatmap_score.values(),
                    )
                ),
            )
        )
        log.info("Done calculating Heatmap Scores.")

    def CalHeatMapScore(self):
        log.info("Calculating Heatmap Scores...")
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        psm = mztab_data.spectrum_match_table
        meta_data = dict(mztab_data.metadata)
        if self.pep_table_exists:
            pep_table = mztab_data.peptide_table
            pep_table = pep_table.fillna(np.nan)
            pep_table.loc[:, "stand_spectra_ref"] = pep_table.apply(
                lambda x: self.file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
                axis=1,
            )
            study_variables = list(
                filter(
                    lambda x: re.match(r"peptide_abundance_study_variable.*?", x) is not None,
                    pep_table.columns.tolist(),
                )
            )

            for name, group in pep_table.groupby("stand_spectra_ref"):

                self.heatmap_con_score[name] = 1.0 - np.sum(
                    np.sum(
                        group[group["accession"].str.contains(config.kwargs["contaminant_affix"])][
                            study_variables
                        ]
                    )
                ) / np.sum(np.sum(group[study_variables]))

                if config.kwargs["remove_decoy"]:
                    pep_median = np.nanmedian(
                        group[(group["opt_global_cv_MS:1002217_decoy_peptide"] == 0)][
                            study_variables
                        ].to_numpy()
                    )
                else:
                    pep_median = np.nanmedian(group[study_variables].to_numpy())
                self.heatmap_pep_intensity[name] = np.minimum(
                    1.0, pep_median / (2 ** 23)
                )  # Threshold

        #  HeatMapMissedCleavages
        global_peps = set(psm["opt_global_cv_MS:1000889_peptidoform_sequence"])
        global_peps_count = len(global_peps)
        if (
                config.kwargs["remove_decoy"]
                and "opt_global_cv_MS:1002217_decoy_peptide" in psm.columns
        ):
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] == 0]
        psm.loc[:, "stand_spectra_ref"] = psm.apply(
            lambda x: self.file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
            axis=1,
        )

        enzyme_list = [i for i in meta_data.values() if str(i).startswith("enzyme:")]
        enzyme = enzyme_list[0].split(":")[1] if len(enzyme_list) == 1 else "Trypsin"
        psm.loc[:, "missed_cleavages"] = psm.apply(
            lambda x: self.cal_MissCleavages(x["sequence"], enzyme), axis=1
        )

        # Calculate the ID RT Score
        for name, group in psm.groupby("stand_spectra_ref"):
            sc = group["missed_cleavages"].value_counts()
            mis_0 = sc[0] if 0 in sc else 0
            self.MissedCleavages_heatmap_score[name] = mis_0 / sc[:].sum()

            x = group["retention_time"] / np.sum(group["retention_time"])
            n = len(group["retention_time"])
            y = np.sum(x) / n
            worst = ((1 - y) ** 0.5) * 1 / n + (y ** 0.5) * (n - 1) / n
            sc = np.sum(np.abs(x - y) ** 0.5) / n
            if worst == 0:
                self.ID_RT_score[name] = 1.0
            else:
                self.ID_RT_score[name] = float((worst - sc) / worst)

            #  For HeatMapOverSamplingScore
            self.HeatmapOverSamplingScore[name] = self.oversampling[name]["1"] / np.sum(
                list(self.oversampling[name].values())
            )

            # For HeatMapPepMissingScore
            idFraction = (
                    len(
                        set(group["opt_global_cv_MS:1000889_peptidoform_sequence"]).intersection(
                            global_peps
                        )
                    )
                    / global_peps_count
            )
            self.HeatmapPepMissingScore[name] = np.minimum(1.0, idFraction)

        median = np.median(list(self.MissedCleavages_heatmap_score.values()))
        self.MissedCleavagesVar_score = dict(
            zip(
                self.MissedCleavages_heatmap_score.keys(),
                list(
                    map(
                        lambda v: 1 - np.abs(v - median),
                        self.MissedCleavages_heatmap_score.values(),
                    )
                ),
            )
        )
        log.info("Done calculating Heatmap Scores.")

    # if missed.cleavages is not given, it is assumed that Trypsin was used for digestion
    @staticmethod
    def cal_MissCleavages(sequence, enzyme):
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
    def dis_decoy(ProteinName):
        if config.kwargs["decoy_affix"] not in ProteinName:
            return "TARGET"
        elif ProteinName.split(";") == 1:
            return "DECOY"
        else:
            if config.kwargs["affix_type"] == "prefix":
                if list(
                        filter(
                            lambda x: lambda x: not x.startswith(config.kwargs["decoy_affix"]),
                            ProteinName.split(";"),
                        )
                ):
                    return "TARGET"
                return "DECOY"
            else:
                if list(
                        filter(
                            lambda x: not x.endswith(config.kwargs["decoy_affix"]),
                            ProteinName.split(";"),
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

        mzml_table = {}
        heatmap_charge = {}
        self.ms_without_psm = []

        def add_ms_values(info_df, ms_name):
            charge_state = int(info_df["Charge"]) if info_df["Charge"] is not None else None
            base_peak_intensity = (
                float(info_df["Base_Peak_Intensity"])
                if info_df["Base_Peak_Intensity"] is not None
                else None
            )
            peak_per_ms2 = int(info_df["MS_peaks"]) if info_df["MS_peaks"] is not None else None

            if self.enable_dia:
                self.mzml_charge_plot.add_value(charge_state)
                self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                return

            if ms_name in self.ms_with_psm:
                if info_df["SpectrumID"] in self.identified_spectrum[ms_name]:
                    self.mzml_charge_plot.add_value(charge_state)
                    self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                    self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                else:
                    self.mzml_charge_plot_1.add_value(charge_state)
                    self.mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                    self.mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
            else:
                if ms_name not in self.ms_without_psm:
                    self.ms_without_psm.append(ms_name)

        def read_mzmls():
            exp = MSExperiment()
            for m in self.ms_paths:
                ms1_number = 0
                ms2_number = 0
                log.info(
                    "{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
                )
                MzMLFile().load(m, exp)
                log.info(
                    "{}: Done parsing mzML file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), m
                    )
                )
                m = self.file_prefix(m)
                log.info(
                    "{}: Aggregating mzML file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), m
                    )
                )
                charge_2 = 0
                for i in exp:
                    if i.getMSLevel() == 1:
                        ms1_number += 1
                    elif i.getMSLevel() == 2:
                        ms2_number += 1
                        charge_state = i.getPrecursors()[0].getCharge()
                        peaks_tuple = i.get_peaks()
                        peak_per_ms2 = len(peaks_tuple[0])
                        if i.getMetaValue("base peak intensity"):
                            base_peak_intensity = i.getMetaValue("base peak intensity")
                        else:
                            base_peak_intensity = (
                                max(peaks_tuple[1]) if len(peaks_tuple[1]) > 0 else None
                            )

                        if charge_state == 2:
                            charge_2 += 1

                        if self.enable_dia:
                            self.mzml_charge_plot.add_value(charge_state)
                            self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                            self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                            continue

                        if m in self.ms_with_psm:
                            if i.getNativeID() in self.identified_spectrum[m]:
                                self.mzml_charge_plot.add_value(charge_state)
                                self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                                self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                            else:
                                self.mzml_charge_plot_1.add_value(charge_state)
                                self.mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                                self.mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
                        else:
                            if m not in self.ms_without_psm:
                                self.ms_without_psm.append(m)

                heatmap_charge[m] = charge_2 / ms2_number
                self.total_ms2_spectra = self.total_ms2_spectra + ms2_number
                mzml_table[m] = {"MS1_Num": ms1_number}
                mzml_table[m]["MS2_Num"] = ms2_number
                log.info(
                    "{}: Done aggregating mzML file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), m
                    )
                )

        def read_ms_info():
            for file in self.ms_info_path:
                log.info(
                    "{}: Parsing ms_statistics dataframe {}...".format(
                        datetime.now().strftime("%H:%M:%S"), file
                    )
                )
                mzml_df = pd.read_parquet(file)
                m = self.file_prefix(file).replace("_ms_info", "")
                if m not in mzml_table:
                    mzml_table[m] = dict.fromkeys(["MS1_Num", "MS2_Num", "Charge_2"], 0)
                charge_group = mzml_df.groupby("Charge").size()
                ms_level_group = mzml_df.groupby("MSLevel").size()
                charge_2 = charge_group[2] if 2 in charge_group else 0
                ms1_number = int(ms_level_group[1]) if 1 in ms_level_group else 0
                ms2_number = int(ms_level_group[2]) if 2 in ms_level_group else 0
                self.total_ms2_spectra = self.total_ms2_spectra + ms2_number
                mzml_table[m].update({"MS1_Num": mzml_table[m]["MS1_Num"] + ms1_number})
                mzml_table[m].update({"MS2_Num": mzml_table[m]["MS2_Num"] + ms2_number})
                mzml_table[m].update({"Charge_2": mzml_table[m]["Charge_2"] + charge_2})

                (
                    self.ms1_tic[os.path.basename(file).replace("_ms_info.parquet", "")],
                    self.ms1_bpc[os.path.basename(file).replace("_ms_info.parquet", "")],
                    self.ms1_peaks[os.path.basename(file).replace("_ms_info.parquet", "")],
                    self.ms1_general_stats[os.path.basename(file).replace("_ms_info.parquet", "")],
                ) = get_ms_qc_info(mzml_df)

                group = mzml_df[mzml_df["MSLevel"] == 2]
                del mzml_df
                group.apply(add_ms_values, args=(m,), axis=1)

                for m in mzml_table.keys():
                    heatmap_charge[m] = mzml_table[m]["Charge_2"] / mzml_table[m]["MS2_Num"]
                log.info(
                    "{}: Done aggregating ms_statistics dataframe {}...".format(
                        datetime.now().strftime("%H:%M:%S"), file
                    )
                )

        read_ms_info() if self.read_ms_info else read_mzmls()
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
        consensus_paths = []
        for raw_id in self.idx_paths:
            if "consensus" in os.path.split(raw_id)[1]:
                consensus_paths.append(raw_id)
                self.idx_paths.remove(raw_id)

        self.MSGF_label, self.Comet_label, self.Sage_label = False, False, False
        self.search_engine = {
            "SpecE": OrderedDict(),
            "xcorr": OrderedDict(),
            "hyper": OrderedDict(),
            "PEPs": OrderedDict(),
            "consensus_support": OrderedDict(),
            "data_label": OrderedDict(),
        }
        SpecE_label, xcorr_label, hyper_label, PEPs_label, consensus_label = [], [], [], [], []

        for raw_id in self.idx_paths:

            # log.info("Parsing search result file {}...".format(raw_id))
            log.info(
                "{}: Parsing search result file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), raw_id
                )
            )

            protein_ids = []
            peptide_ids = []
            IdXMLFile().load(raw_id, protein_ids, peptide_ids)
            raw_id = self.file_prefix(raw_id)
            ## TODO I would use the QC functionality of pyopenms. Should be much faster.
            if config.kwargs["remove_decoy"]:
                identified_num = len(
                    set(
                        [
                            i.getMetaValue("spectrum_reference")
                            for i in peptide_ids
                            if i.getHits()[0].getMetaValue("target_decoy") == "target"
                        ]
                    )
                )
            else:
                identified_num = len(peptide_ids)

            ## TODO make clear if this is before or after first step of filtering
            ms_name = self.file_prefix(
                protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8")
            )
            search_engine = protein_ids[0].getSearchEngine()

            self.search_engine["SpecE"][raw_id] = OrderedDict()
            self.search_engine["xcorr"][raw_id] = OrderedDict()
            self.search_engine["hyper"][raw_id] = OrderedDict()
            self.search_engine["PEPs"][raw_id] = OrderedDict()

            xcorr_breaks = list(
                np.arange(
                    self.XCORR_HIST_RANGE["start"],
                    self.XCORR_HIST_RANGE["end"] + self.XCORR_HIST_RANGE["step"],
                    self.XCORR_HIST_RANGE["step"],
                ).round(2)
            )

            hyper_breaks = list(
                np.arange(
                    self.HYPER_HIST_RANGE["start"],
                    self.HYPER_HIST_RANGE["end"] + self.HYPER_HIST_RANGE["step"],
                    self.HYPER_HIST_RANGE["step"],
                ).round(2)
            )

            SpecE_breaks = list(
                np.arange(
                    self.SPECEVALUE_HIST_RANGE["start"],
                    self.SPECEVALUE_HIST_RANGE["end"] + self.SPECEVALUE_HIST_RANGE["step"],
                    self.SPECEVALUE_HIST_RANGE["step"],
                ).round(2)
            )
            SpecE_breaks.append(float("inf"))
            SpecE_breaks.sort()

            PEP_breaks = list(
                np.arange(
                    self.PEP_HIST_RANGE["start"],
                    self.PEP_HIST_RANGE["end"] + self.PEP_HIST_RANGE["step"],
                    self.PEP_HIST_RANGE["step"],
                ).round(2)
            )

            bar_stacks = ["target", "decoy", "target+decoy"]
            Xcorr = Histogram(
                "Comet cross-correlation score",
                plot_category="range",
                stacks=bar_stacks,
                breaks=xcorr_breaks,
            )
            Hyper = Histogram(
                "Sage hyperscore", plot_category="range", stacks=bar_stacks, breaks=hyper_breaks
            )
            SpecE = Histogram(
                "MSGF spectral E-value",
                plot_category="range",
                stacks=bar_stacks,
                breaks=SpecE_breaks,
            )
            PEP = Histogram(
                "Posterior error probability",
                plot_category="range",
                stacks=bar_stacks,
                breaks=PEP_breaks,
            )
            Consensus_support = Histogram(
                "Consensus PSM number", plot_category="frequency", stacks=bar_stacks
            )

            if search_engine == "MSGF+" or "msgf" in raw_id:
                mzml_table[ms_name]["MSGF"] = identified_num
                self.MSGF_label = True
                SpecE_label.append({"name": raw_id, "ylab": "Counts"})
                PEPs_label.append({"name": raw_id, "ylab": "Counts"})
                for peptide_id in peptide_ids:
                    for hit in peptide_id.getHits():
                        spec_e = (
                            hit.getMetaValue("SpecEvalue-score")
                            if hit.getMetaValue("SpecEvalue-score")
                            else hit.getMetaValue("MS:1002052")
                        )
                        logSpecE = -math.log(spec_e, 10)
                        pep = (
                            hit.getMetaValue("MS:1001493")
                            if hit.getMetaValue("MS:1001493")
                            else hit.getScore()
                        )
                        SpecE.add_value(logSpecE, stack=hit.getMetaValue("target_decoy"))
                        PEP.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                SpecE.to_dict()
                PEP.to_dict()
                self.search_engine["SpecE"][raw_id] = SpecE.dict["data"]
                self.search_engine["PEPs"][raw_id] = PEP.dict["data"]

            elif search_engine == "Comet" or "comet" in raw_id:
                self.Comet_label = True
                mzml_table[ms_name]["Comet"] = identified_num
                xcorr_label.append({"name": raw_id, "ylab": "Counts"})
                PEPs_label.append({"name": raw_id, "ylab": "Counts"})
                for peptide_id in peptide_ids:
                    for hit in peptide_id.getHits():
                        xcorr = hit.getMetaValue("MS:1002252")
                        pep = (
                            hit.getMetaValue("MS:1001493")
                            if hit.getMetaValue("MS:1001493")
                            else hit.getScore()
                        )
                        Xcorr.add_value(xcorr, stack=hit.getMetaValue("target_decoy"))
                        PEP.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                Xcorr.to_dict()
                PEP.to_dict()
                self.search_engine["xcorr"][raw_id] = Xcorr.dict["data"]
                self.search_engine["PEPs"][raw_id] = PEP.dict["data"]

            elif search_engine == "Sage" or "sage" in raw_id:
                self.Sage_label = True
                mzml_table[ms_name]["Sage"] = identified_num
                hyper_label.append({"name": raw_id, "ylab": "Counts"})
                PEPs_label.append({"name": raw_id, "ylab": "Counts"})
                for peptide_id in peptide_ids:
                    for hit in peptide_id.getHits():
                        hyper = hit.getMetaValue("hyperscore")
                        pep = (
                            hit.getMetaValue("MS:1001493")
                            if hit.getMetaValue("MS:1001493")
                            else hit.getScore()
                        )
                        Hyper.add_value(hyper, stack=hit.getMetaValue("target_decoy"))
                        PEP.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                Hyper.to_dict()
                PEP.to_dict()
                self.search_engine["hyper"][raw_id] = Hyper.dict["data"]
                self.search_engine["PEPs"][raw_id] = PEP.dict["data"]

            else:
                mzml_table[ms_name][search_engine] = identified_num

            mzml_table[ms_name]["num_quant_psms"] = (
                self.mL_spec_ident_final[ms_name]
                if ms_name in self.mL_spec_ident_final.keys()
                else 0
            )
            mzml_table[ms_name]["num_quant_peps"] = (
                len(self.mzml_peptide_map[ms_name])
                if ms_name in self.mL_spec_ident_final.keys()
                else 0
            )

        for raw_id in consensus_paths:
            log.info("Parsing consensus file {}...".format(raw_id))
            protein_ids = []
            peptide_ids = []
            IdXMLFile().load(raw_id, protein_ids, peptide_ids)
            raw_id = self.file_prefix(raw_id)
            consensus_label.append({"name": raw_id, "ylab": "Counts"})
            self.search_engine["consensus_support"][raw_id] = OrderedDict()

            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    support = hit.getMetaValue("consensus_support")
                    Consensus_support.add_value(support, stack=hit.getMetaValue("target_decoy"))
            Consensus_support.to_dict()

            for i in Consensus_support.dict["data"].keys():
                self.search_engine["consensus_support"][raw_id][i] = Consensus_support.dict[
                    "data"
                ][i]

        self.search_engine["data_label"] = {
            "score_label": [SpecE_label, xcorr_label, hyper_label],
            "PEPs_label": PEPs_label,
            "consensus_label": consensus_label,
        }

        # mass spectrum files sorted based on experimental file
        for spectrum_name in self.exp_design_table.keys():
            self.mzml_table[spectrum_name] = mzml_table[spectrum_name]

    def parse_out_mzTab(self):
        log.info(
            "{}: Parsing mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mzTab_path
            )
        )
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        log.info(
            "{}: Done parsing mzTab file {}.".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mzTab_path
            )
        )
        log.info(
            "{}: Aggregating mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.out_mzTab_path
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
            lambda x: self.file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
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
        prot = prot[prot["opt_global_result_type"] != "protein_details"]

        if config.kwargs["remove_decoy"]:
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1]
            # TODO do we really want to remove groups that contain a single decoy? I would say ALL members need to be decoy.
            prot = prot[~prot["accession"].str.contains(config.kwargs["decoy_affix"])]

        prot["protein_group"] = prot["ambiguity_members"].apply(lambda x: x.replace(",", ";"))

        self.Total_Protein_Identified = len(prot.index)

        prot.dropna(how="all", subset=prot_abundance_cols, inplace=True)
        self.Total_Protein_Quantified = len(prot.index)

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

        mL_spec_ident_final = {}

        for m, group in psm.groupby("filename"):
            m = os.path.basename(m)
            self.cal_num_table_data[m] = {"sample_name": self.exp_design_table[m]["Sample"]}
            self.cal_num_table_data[m]["condition"] = self.exp_design_table[m]["MSstats_Condition"]
            self.cal_num_table_data[m]["fraction"] = self.exp_design_table[m]["Fraction"]

            if config.kwargs["remove_decoy"]:
                group = group[group["opt_global_cv_MS:1002217_decoy_peptide"] == 0]
            # Each loop resets the instance.
            self.oversampling_plot = Histogram(
                "MS/MS counts per 3D-peak", plot_category="frequency", breaks=[1, 2, 3]
            )

            group.fillna(pd.NA, inplace=True)
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
            self.cal_num_table_data[m]["protein_num"] = len(proteins)
            self.cal_num_table_data[m]["peptide_num"] = len(peptides)
            self.cal_num_table_data[m]["unique_peptide_num"] = len(unique_peptides)

            modified_pep = list(
                filter(lambda x: re.match(r".*?\(.*\).*?", x) is not None, peptides)
            )
            self.cal_num_table_data[m]["modified_peptide_num"] = len(modified_pep)

            mL_spec_ident_final[m] = len(set(self.identified_spectrum[m]))

        # TODO mzMLs without PSM: experimental design information is displayed, and all quantitative information is 0
        self.ms_without_psm = set([self.file_prefix(i) for i in self.ms_paths]) - set(
            self.ms_with_psm
        )
        for i in self.ms_without_psm:
            self.cal_num_table_data[i] = {
                "sample_name": self.exp_design_table[i]["Sample"],
                "condition": self.exp_design_table[i]["MSstats_Condition"],
                "fraction": self.exp_design_table[i]["Fraction"],
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
        self.mL_spec_ident_final = mL_spec_ident_final
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
            cur.execute(
                "CREATE TABLE PSM(PSM_ID INT(200), Sequence VARCHAR(200), Modification VARCHAR(100), Accession VARCHAR(100), Search_Engine_Score FLOAT(4,5), Spectra_Ref VARCHAR(100))"
            )
            con.commit()
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
            cur.executemany(
                "INSERT INTO PSM (" + sql_col + ") VALUES " + sql_t,
                [(k, *itemgetter(*all_term)(v)) for k, v in mztab_data_psm_full.items()],
            )
            con.commit()

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

            def getSpectrumCountAcrossRep(condition_count_dict: dict):
                Spc = []
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
                    Spc.append(res[c])

                # Integer for average spectrum counting with NA=0 ignored across condition
                res["Average Spectrum Counting"] = round(sum(Spc) / len(np.nonzero(Spc)[0]))
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

                mztab_data_dict_prot_full[index] = getSpectrumCountAcrossRep(
                    mztab_data_dict_prot_full[index]
                )
                mztab_data_dict_prot_full[index]["Peptides_Number"] = int(counts_per_acc[index])

            log.info(
                "{}: Done aggregating mzTab file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.out_mzTab_path
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
            cur.execute(
                'CREATE TABLE PROTQUANT(ProteinName VARCHAR(100), Peptides_Number INT(100), "Average Spectrum Counting" VARCHAR)'
            )
            con.commit()
            sql_col = 'ProteinName,Peptides_Number,"Average Spectrum Counting"'
            sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 3)) + ")"

            for s in conditions:
                cur.execute('ALTER TABLE PROTQUANT ADD "' + str(s) + '" VARCHAR')
                con.commit()
                sql_col += ', "' + str(s) + '"'
                headers[str(s)] = {"name": s}

            for s in list(map(lambda x: str(x) + "_distribution", conditions)):
                cur.execute('ALTER TABLE PROTQUANT ADD "' + s + '" VARCHAR(100)')
                con.commit()
                sql_col += ', "' + s + '"'
                headers[str(s)] = {"name": s}

            # ProteinName is index
            all_term = (
                    ["Peptides_Number", "Average Spectrum Counting"]
                    + list(map(str, conditions))
                    + list(map(lambda x: str(x) + "_distribution", conditions))
            )
            cur.executemany(
                "INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
                [(k, *itemgetter(*all_term)(v)) for k, v in mztab_data_dict_prot_full.items()],
            )
            con.commit()

            pconfig = {
                "id": "quantification_of_protein",  # ID used for the table
                "title": "quantification information of protein",
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

    def parse_out_mgf(self):
        def get_spectrumID(SpectrumTitle, index):
            if "scan=" in SpectrumTitle:
                spectrumID = SpectrumTitle
            else:
                spectrumID = "index=" + str(index)
            return spectrumID

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
            m = self.file_prefix(m)
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
                title = get_spectrumID(raw_title, i)
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

        def Process_Modification(Modification):
            if not isinstance(Modification, list):
                modifications = None
            else:
                modifi_list = list()
                for i in Modification:
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
            m = self.file_prefix(m)
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
            lambda x: self.file_prefix(x.location), axis=1
        )
        mzid_table["Modifications"] = mzid_table["Modification"].apply(
            lambda x: Process_Modification(x)
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

        mzid_table["filename"] = mzid_table.apply(lambda x: self.file_prefix(x.location), axis=1)
        self.ms_with_psm = mzid_table["filename"].unique().tolist()

        # TODO remove_decoy
        if config.kwargs["remove_decoy"]:
            mzid_table = mzid_table[~mzid_table["isDecoy"]]

        self.Total_Protein_Identified = mzid_table["accession_group"].nunique()

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

            # unique_peptides = group[group["pep_to_prot_unique"] == True][["PeptideSequence", "Modifications"]].drop_duplicates()
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
            self.ms_without_psm = set([self.file_prefix(i) for i in self.mgf_paths]) - set(
                self.ms_with_psm
            )
        elif self.ms_paths:
            self.ms_without_psm = set([self.file_prefix(i) for i in self.ms_paths]) - set(
                self.ms_with_psm
            )
        for i in self.ms_without_psm:
            self.cal_num_table_data[self.file_prefix(i)] = {
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
        pattern = re.compile(r"\(.*?\)")
        report_data = pd.read_csv(self.diann_report_path, header=0, sep="\t")
        report_data["sequence"] = report_data.apply(
            lambda x: re.sub(pattern, "", x["Modified.Sequence"]), axis=1
        )
        self.Total_Protein_Quantified = len(set(report_data["Protein.Names"]))
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

        for run_file, group in report_data.groupby("File.Name"):
            run_file = self.file_prefix(run_file)
            self.ms_with_psm.append(run_file)
            self.cal_num_table_data[run_file] = {
                "sample_name": self.exp_design_table[run_file]["Sample"]
            }
            self.cal_num_table_data[run_file]["condition"] = self.exp_design_table[run_file][
                "MSstats_Condition"
            ]
            self.cal_num_table_data[run_file]["fraction"] = self.exp_design_table[run_file][
                "Fraction"
            ]
            self.cal_num_table_data[run_file]["protein_num"] = len(set(group["Protein.Ids"]))
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

        self.ms_without_psm = set([self.file_prefix(i) for i in self.ms_paths]) - set(
            self.ms_with_psm
        )
        for i in self.ms_without_psm:
            log.warning("No PSM found in '{}'!".format(i))

        for i in self.ms_without_psm:
            self.cal_num_table_data[i] = {
                "sample_name": self.exp_design_table[i]["Sample"],
                "condition": self.exp_design_table[i]["MSstats_Condition"],
                "fraction": self.exp_design_table[i]["Fraction"],
                "protein_num": 0,
                "peptide_num": 0,
                "unique_peptide_num": 0,
                "modified_peptide_num": 0,
            }

    def parse_msstats_input(self):
        log.info("Parsing MSstats input file " + self.msstats_input_path)
        msstats_data = pd.read_csv(self.msstats_input_path)
        ## TODO we probably shouldn't even write out 0-intensity values to MSstats csv
        msstats_data = msstats_data[-(msstats_data["Intensity"] == 0)]
        msstats_data.loc[:, "BestSearchScore"] = 1 - msstats_data.loc[:, "PeptideSequence"].map(
            self.peptide_search_score
        )
        msstats_data[["PeptideSequence", "Modification"]] = msstats_data.apply(
            lambda x: find_modification(x["PeptideSequence"]), axis=1, result_type="expand"
        )

        max_pep_intensity = 0.0

        repsPerCondition = (
            self.sample_df.groupby("MSstats_Condition")["MSstats_BioReplicate"].agg(list).to_dict()
        )
        conditions = list(self.sample_df["MSstats_Condition"].unique())
        conditions_str = [str(c) for c in conditions]
        conditions_dists = [str(c) + "_distribution" for c in conditions]
        cond_and_dist_cols = conditions_str + conditions_dists

        # TODO maybe aggregating in dicts is not the fastest. We also need to parse them again for proteins later.
        #  Maybe we should just introduce new pandas columns for every bioreplicate.
        def fillDict(g):
            d = dict.fromkeys(repsPerCondition[str(g.name)], 0)
            d.update(zip(g["BioReplicate"].astype(str), np.log10(g["Intensity"])))
            return json.dumps(d)

        def getIntyAcrossBioRepsAsStr(g):
            gdict = dict.fromkeys(conditions_str, 0.0)
            gdict.update(dict.fromkeys(conditions_dists, "{}"))
            gdict["Average Intensity"] = np.log10(g["Intensity"].mean())
            gdict["BestSearchScore"] = g["BestSearchScore"].min()
            ## TODO How to determine technical replicates? Should be same BioReplicate but different Fraction_Group (but fraction group is not annotated)
            condGrp = (
                g.groupby(["Condition", "BioReplicate"])["Intensity"]
                .mean()
                .reset_index()
                .groupby("Condition")
                .apply(fillDict)
            )
            condGrp.index = [str(c) + "_distribution" for c in condGrp.index]
            gdict.update(condGrp.to_dict())
            mean = g.groupby(["Condition"])["Intensity"].mean()
            condGrpMean = np.log10(mean)
            condGrpMean.index = condGrpMean.index.map(str)
            gdict.update(condGrpMean.to_dict())
            return pd.Series(gdict)

        msstats_data_pep_agg = msstats_data.groupby(
            ["PeptideSequence", "ProteinName", "Modification"]
        ).apply(
            getIntyAcrossBioRepsAsStr
        )  # .unstack()
        del msstats_data
        ## TODO Can we guarantee that the score was always PEP? I don't think so!
        msstats_data_pep_agg.reset_index(inplace=True)
        msstats_data_pep_agg.index = msstats_data_pep_agg.index + 1
        msstats_data_dict_pep_full = msstats_data_pep_agg.to_dict("index")
        msstats_data_dict_pep_init = dict(itertools.islice(msstats_data_dict_pep_full.items(), 50))

        cur.execute(
            'CREATE TABLE PEPQUANT(PeptideID INT(100) PRIMARY KEY, PeptideSequence VARCHAR(100), Modification VARCHAR(100), ProteinName VARCHAR(100), BestSearchScore FLOAT(4,3), "Average Intensity" FLOAT(4,3))'
        )
        con.commit()
        sql_col = 'PeptideID,PeptideSequence,Modification,ProteinName,BestSearchScore, "Average Intensity"'
        sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 6)) + ")"

        headers = OrderedDict()
        headers = {  # 'PeptideID': {'name': 'PeptideID'}, # this is the index
            "PeptideSequence": {"rid": "PeptideSequence"},
            "Modification": {"rid": "Modification"},
            "ProteinName": {"rid": "ProteinName"},
            "BestSearchScore": {"rid": "BestSearchScore", "format": "{:,.5e}"},
            "Average Intensity": {"rid": "Average Intensity", "format": "{:,.3f}"},
        }

        for s in conditions:
            cur.execute('ALTER TABLE PEPQUANT ADD "' + str(s) + '" VARCHAR')
            con.commit()
            sql_col += ', "' + str(s) + '"'
            headers[str(s)] = {"rid": s, "format": "{:,.5f}"}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute('ALTER TABLE PEPQUANT ADD "' + s + '" VARCHAR(100)')
            con.commit()
            sql_col += ', "' + s + '"'
            # we need a thousandsSep_format otherwise commas will be replaced
            headers[str(s)] = {"rid": s}

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
        cur.executemany(
            "INSERT INTO PEPQUANT (" + sql_col + ") VALUES " + sql_t,
            [(k, *itemgetter(*all_term)(v)) for k, v in msstats_data_dict_pep_full.items()],
        )
        con.commit()

        pconfig = {
            "id": "quantification_of_peptides",  # ID used for the table
            "title": "quantification information of peptides",
            # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_quantification_of_peptides_table",  # File basename to use for raw data file
            "sort_rows": False,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "PeptideID",
            "no_violin": True,
        }

        # only use the first 50 lines for the table
        table_html = sparklines.plot(
            data=msstats_data_dict_pep_init,
            headers=headers,
            pconfig=pconfig,
            max_value=max_pep_intensity,
        )
        pattern = re.compile(r'<small id="quantification_of_peptides_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = (
                table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" '
                                     'onkeyup="searchQuantFunction()" id="quant_search">'
                                     '<select name="quant_search_col" id="quant_search_col">'
        )
        for key in ["ProteinName", "PeptideSequence", "Modification", "PeptideID"]:
            t_html += "<option>" + key + "</option>"
        table_html = (
                t_html + "</select>" + '<button type="button" class="btn btn-default '
                                       'btn-sm" id="quant_reset" onclick="quantFirst()">Reset</button>' + table_html[
                                                                                                          index:]
        )
        table_html = (
                table_html
                + """<div class="page_control"><span id="quantFirst">First Page</span><span 
        id="quantPre"> Previous Page</span><span id="quantNext">Next Page </span><span id="quantLast">Last 
        Page</span><span id="quantPageNum"></span>Page/Total <span id="quantTotalPage"></span>Pages <input 
        type="number" name="" id="pep_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
        onkeydown="quant_page_jump()" min="1"/> </div> """
        )

        # Add a report section with the line plot
        self.add_section(
            name="Peptides Quantification Table",
            anchor="peptides_quant_result",
            description="This plot shows the quantification information of peptides"
                        " in the final result (mainly the mzTab file).",
            helptext="""
                    The quantification information of peptides is obtained from the MSstats input file. 
                    The table shows the quantitative level and distribution of peptides in different study variables, run and peptiforms. The distribution show all the intensity values in a bar plot above and below the average intensity for all the fractions, runs and peptiforms.

                    * BestSearchScore: It is equal to 1 - min(Q.Value) for DIA datasets. Then it is equal to 1 - min(best_search_engine_score[1]), which is from best_search_engine_score[1] column in mzTab peptide table for DDA datasets.
                    * Average Intensity: Average intensity of each peptide sequence across all conditions with NA=0 or NA ignored.
                    * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of fractions, and then mean intensity in technical replicates/biological replicates separately.
                    Click `Show replicates` to switch to bar plots for every replicate.

                    """,
            plot=table_html,
        )

        report.plot_by_id["peptides_quant_result"] = table_html

        # Helper functions for pandas
        def jsonToDict(s):
            if type(s) is str:
                return json.loads(s)
            else:
                return {}

        def reducer(accumulator, element):
            for key, value in jsonToDict(element).items():
                accumulator[key] = np.log10(pow(10, accumulator.get(key, 0)) + pow(10, value))
            return accumulator

        def myDictSum(series):
            return json.dumps(reduce(reducer, series, {}))

        def totalIntensity(series):
            total = 0.0
            for intensity in series:
                total += pow(10, intensity)
            return np.log10(total)

        def uniqueCount(series):
            return len(series.unique())

        max_prot_intensity = 0
        agg_funs = dict.fromkeys(conditions_dists, myDictSum)
        agg_funs.update(dict.fromkeys(conditions_str, totalIntensity))
        agg_funs["PeptideSequence"] = uniqueCount
        agg_funs["Average Intensity"] = totalIntensity
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

        headers = OrderedDict()
        # This will be the rowheader/index
        # headers['ProteinID'] = {
        #    'name': 'Protein ID',
        #    'description': 'Count(s) of the protein (group)',
        #    'format': '{:,.0f}'
        # }
        headers["ProteinName"] = {
            "rid": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
            "format": "{:,.0f}",
        }
        headers["Peptides_Number"] = {
            "rid": "Number of Peptides",
            "description": "Number of peptides per proteins",
            "format": "{:,.0f}",
        }
        headers["Average Intensity"] = {
            "rid": "Average Intensity",
            "description": "Average intensity across all conditions",
            "format": "{:,.3f}",
        }

        # upload protein table to sqlite database
        cur.execute(
            'CREATE TABLE PROTQUANT(ProteinID INT(100), ProteinName VARCHAR(100), Peptides_Number INT(100), "Average Intensity" VARCHAR)'
        )
        con.commit()
        sql_col = 'ProteinID,ProteinName,Peptides_Number,"Average Intensity"'
        sql_t = "(" + ",".join(["?"] * (len(conditions) * 2 + 4)) + ")"

        for s in conditions:
            cur.execute('ALTER TABLE PROTQUANT ADD "' + str(s) + '" VARCHAR')
            con.commit()
            sql_col += ', "' + str(s) + '"'
            headers[str(s)] = {"rid": s}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute('ALTER TABLE PROTQUANT ADD "' + s + '" VARCHAR(100)')
            con.commit()
            sql_col += ', "' + s + '"'
            headers[str(s)] = {"rid": s}

        # ProteinID is index
        all_term = (
                ["ProteinName", "Peptides_Number", "Average Intensity"]
                + list(map(str, conditions))
                + list(map(lambda x: str(x) + "_distribution", conditions))
        )
        cur.executemany(
            "INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
            [(k, *itemgetter(*all_term)(v)) for k, v in msstats_data_dict_prot_full.items()],
        )
        con.commit()

        pconfig = {
            "id": "quantification_of_protein",  # ID used for the table
            "title": "quantification information of protein",
            # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_quantification_of_protein_table",  # File basename to use for raw data file
            "sort_rows": False,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "ProteinID",
            # 'format': '{:,.5f}',
            "no_violin": True,
        }

        table_html = sparklines.plot(
            msstats_data_dict_prot_init, headers, pconfig=pconfig, max_value=max_prot_intensity
        )
        pattern = re.compile(r'<small id="quantification_of_protein_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = (
                table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" '
                                     'onkeyup="searchProtFunction()" id="prot_search">'
                                     '<select name="prot_search_col" id="prot_search_col">'
        )
        for key in ["ProteinName", "Peptides_Number", "ProteinID"]:
            t_html += "<option>" + key + "</option>"
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

        self.add_section(
            name="Protein Quantification Table",
            anchor="protein_quant_result",
            description="This plot shows the quantification information of proteins"
                        " in the final result (mainly the mzTab file).",
            helptext="""
                    The quantification information of proteins is obtained from the msstats input file. 
                    The table shows the quantitative level and distribution of proteins in different study variables and run.

                    * Peptides_Number: The number of peptides for each protein.
                    * Average Intensity: Average intensity of each protein across all conditions with NA=0 or NA ignored.
                    * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
                    
                    Click `Show replicates` to switch to bar plots of quantities in each replicate.
                    """,
            plot=table_html,
        )
        report.plot_by_id["protein_quant_result"] = table_html
        self.sections[-1].plot_id = "protein_quant_result"
        self.sections[-2].plot_id = "peptides_quant_result"

    # MaxQuant Files
    def maxquant_file_path(self):
        maxquant_files = []
        for maxquant_file in self.find_log_files("quantms/maxquant_result", filecontents=False):
            maxquant_files.append(maxquant_file["fn"])

        # MaxQuant file check.
        # TODO "mqpar.xml" ?

        # msmsScans.txt or msScans.txt
        # required_files = ["parameters.txt", "summary.txt", "proteinGroups.txt", "evidence.txt", "msms.txt", "msmsScans.txt"]
        # required_files = ["parameters.txt", "summary.txt", "proteinGroups.txt", "evidence.txt", "msms.txt", "msScans.txt"]
        # if not set(required_files).issubset(set(maxquant_files)):
        #     raise ValueError("Required files are missing. Please check your MaxQuant results!")

        required_files = [
            "parameters.txt",
            "summary.txt",
            "proteinGroups.txt",
            "evidence.txt",
            "msms.txt",
            "msScans.txt",
            "msmsScans.txt",
        ]

        maxquant_paths = {}
        for maxquant_file in self.find_log_files("quantms/maxquant_result", filecontents=False):
            if maxquant_file["fn"] in required_files:
                f_path = os.path.join(maxquant_file["root"], maxquant_file["fn"])
                maxquant_paths[self.file_prefix(f_path)] = f_path

        return maxquant_paths

    # MaxQuant parameters table
    def draw_parameters(self, parameter_table):
        draw_config = {
            "namespace": "",
            "id": "parameters",
            "title": "Parameters",
            "save_file": False,
            "sort_rows": False,
            "only_defined_headers": True,
            "col1_header": "No.",
            "no_violin": True,
        }
        headers = {"Parameter": {"title": "Parameter"}, "Value": {"title": "Value"}}
        table_html = table.plot(parameter_table, headers=headers, pconfig=draw_config)
        self.add_section(
            name="Parameters",
            anchor="parameters",
            # description='Parameters',
            helptext="""
            MaxQuant parameters, extracted from parameters.txt, summarizes the settings used for the MaxQuant analysis. 
            Key parameters are MaxQuant version, Re-quantify, Match-between-runs and mass search tolerances. 
            A list of protein database files is also provided, allowing to track database completeness 
            and database version information (if given in the filename).
            """,
            plot=table_html,
        )

    # MaxQuant Fig 1
    def draw_pg_contaminants(self, contaminant_percent):
        draw_config = {
            "id": "potential_contaminants_per_group",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Potential Contaminants per Group",
            "tt_decimals": 2,
            "ylab": "Percent [%]",
        }
        bar_html = bargraph.plot(contaminant_percent, pconfig=draw_config)
        self.add_section(
            name="Potential Contaminants per Group",
            anchor="potential_contaminants_per_group",
            description="Potential contaminants per group from proteinGroups.txt.",
            helptext="""
            External protein contamination should be controlled for, therefore MaxQuant ships with a 
            comprehensive, yet customizable protein contamination database, which is searched by MaxQuant 
            by default. A contamination plot derived from the proteinGroups (PG) table, showing the 
            fraction of total protein intensity attributable to contaminants. 
            
            Note that this plot is based on experimental groups, and therefore may not correspond 1:1 to Raw files.
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 2, Fig 4
    # def draw_intensity_his(self, distribution_histogram, fig_type):

    #     # 'intensity'
    #     if fig_type == 'intensity':
    #         draw_config = {
    #             'id': 'intensity_distribution_histogram',
    #             'cpswitch': False,
    #             'cpswitch_c_active': False,
    #             'title': 'Intensity Distribution',
    #             'tt_decimals': 0,
    #             'data_labels': ['Sample', 'Contaminants'],
    #             'ylab': 'Count'
    #         }
    #         bar_html = bargraph.plot(distribution_histogram, pconfig=draw_config)
    #         self.add_section(
    #             name='Intensity Distribution (Histogram)',
    #             anchor='intensity_distribution_histogram',
    #             # description='Description of ~Intensity Distribution',
    #             # helptext='',
    #             plot=bar_html
    #         )

    #     # 'LFQ intensity'
    #     elif fig_type == 'lfq_intensity':
    #         draw_config = {
    #             'id': 'lfq_intensity_distribution_histogram',
    #             'cpswitch': False,
    #             'cpswitch_c_active': False,
    #             'title': 'LFQ Intensity Distribution',
    #             'tt_decimals': 0,
    #             'data_labels': ['Sample', 'Contaminants'],
    #             'ylab': 'Count'
    #         }
    #         bar_html = bargraph.plot(distribution_histogram, pconfig=draw_config)
    #         self.add_section(
    #             name='LFQ Intensity Distribution (Histogram)',
    #             anchor='lfq_intensity_distribution_histogram',
    #             # description='Description of ~LFQ Intensity Distribution',
    #             # helptext='',
    #             plot=bar_html
    #         )

    #     # 'peptide intensity'
    #     elif fig_type == 'peptide_intensity':
    #         draw_config = {
    #             'id': 'peptide_intensity_distribution_histogram',
    #             'cpswitch': False,
    #             'cpswitch_c_active': False,
    #             'title': 'Peptide Intensity Distribution',
    #             'tt_decimals': 0,
    #             'data_labels': ['Sample', 'Contaminants'],
    #             'ylab': 'Count'
    #         }
    #         bar_html = bargraph.plot(distribution_histogram, pconfig=draw_config)
    #         self.add_section(
    #             name='Peptide Intensity Distribution (Histogram)',
    #             anchor='peptide_intensity_distribution_histogram',
    #             # description='Description of ~Peptide Intensity Distribution',
    #             # helptext='',
    #             plot=bar_html
    #         )

    def draw_intensity_box(self, distribution_box, fig_type):

        # 'intensity'
        if fig_type == "intensity":
            draw_config = {
                "id": "intensity_distribution_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "Intensity Distribution",
                "tt_decimals": 2,
                "data_labels": ["Sample", "Contaminants"],
                "xlab": "log2(Intensity)",
            }
            box_html = box.plot(distribution_box, pconfig=draw_config)
            self.add_section(
                # name='Intensity Distribution (Box Plot)',
                name="Intensity Distribution",
                anchor="intensity_distribution_box",
                # description='Intensity boxplots by experimental groups.',
                helptext="""
                Intensity boxplots by experimental groups. Groups are user-defined during MaxQuant configuration. 
                This plot displays a (customizable) threshold line for the desired mean intensity of proteins. 
                Groups which underperform here, are likely to also suffer from a worse MS/MS id rate and higher 
                contamination due to the lack of total protein loaded/detected. If possible, all groups should 
                show a high and consistent amount of total protein. 
                
                The height of the bar correlates to the number of proteins with non-zero abundance.
                """,
                plot=box_html,
            )

        # 'LFQ intensity'
        elif fig_type == "lfq_intensity":
            draw_config = {
                "id": "lfq_intensity_distribution_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "LFQ Intensity Distribution",
                "tt_decimals": 2,
                "data_labels": ["Sample", "Contaminants"],
                "xlab": "log2(Intensity)",
            }
            box_html = box.plot(distribution_box, pconfig=draw_config)
            self.add_section(
                # name='LFQ Intensity Distribution (Box Plot)',
                name="LFQ Intensity Distribution",
                anchor="lfq_intensity_distribution_box",
                description="Label-free quantification (LFQ) intensity boxplots by experimental groups.",
                helptext="""
                Label-free quantification (LFQ) intensity boxplots by experimental groups. Groups are user-defined 
                during MaxQuant configuration. This plot displays a (customizable) threshold line for the desired 
                mean of LFQ intensity of proteins. Raw files which underperform in *Raw* intensity, are likely to 
                show an *increased* mean here, since only high-abundance proteins are recovered and quantifyable by 
                MaxQuant in this Raw file. The remaining proteins are likely to receive an LFQ value of 0 (i.e. do 
                not contribute to the distribution).
                
                The height of the bar correlates to the number of proteins with non-zero abundance.
                """,
                plot=box_html,
            )

        # 'peptide intensity'
        elif fig_type == "peptide_intensity":
            draw_config = {
                "id": "peptide_intensity_distribution_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "Peptide Intensity Distribution",
                "tt_decimals": 2,
                "data_labels": ["Sample", "Contaminants"],
                "xlab": "log2(Intensity)",
            }
            box_html = box.plot(distribution_box, pconfig=draw_config)
            self.add_section(
                # name='Peptide Intensity Distribution (Box Plot)',
                name="Peptide Intensity Distribution",
                anchor="peptide_intensity_distribution_box",
                description="Peptide precursor intensity per Raw file from evidence.txt WITHOUT match-between-runs evidence.",
                helptext="""
                Peptide precursor intensity per Raw file from evidence.txt WITHOUT match-between-runs evidence. 
                Low peptide intensity usually goes hand in hand with low MS/MS identifcation rates and unfavourable 
                signal/noise ratios, which makes signal detection harder. Also instrument acquisition time increases 
                for trapping instruments. 
                Failing to reach the intensity threshold is usually due to unfavorable column conditions, inadequate 
                column loading or ionization issues. If the study is not a dilution series or pulsed SILAC experiment, 
                we would expect every condition to have about the same median log-intensity (of 2<sup>%1.1f</sup>). 
                The relative standard deviation (RSD) gives an indication about reproducibility across files and should 
                be below 5%%.
                """,
                plot=box_html,
            )

    # MaxQuant Fig 6: PCA
    def draw_pg_pca(self, pca_data, fig_type):

        # fig_type: 'raw_intensity' or 'lfq_intensity'

        if fig_type == "raw_intensity":
            fig_id = "pca_of_raw_intensity"
            fig_title = "PCA of Raw Intensity"
        if fig_type == "lfq_intensity":
            fig_id = "pca_of_lfq_intensity"
            fig_title = "PCA of LFQ Intensity"

        draw_config = {
            "id": fig_id,
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": fig_title,
            "xlab": "PC #1",
            "ylab": "PC #2",
        }
        scatter_html = scatter.plot(pca_data, pconfig=draw_config)
        self.add_section(
            name=fig_title,
            anchor=fig_id,
            description="[Excludes Contaminants] Principal components plots of experimental groups (as defined during MaxQuant configuration).",
            helptext="""
            This plot is shown only if more than one experimental group was defined. 
            If LFQ was activated in MaxQuant, an additional PCA plot for LFQ intensities is shown. 
            Similarly, if iTRAQ/TMT reporter intensities are detected. 
            Since experimental groups and Raw files do not necessarily correspond 1:1, this plot cannot use the abbreviated 
            Raw file names, but instead must rely on automatic shortening of group names.
            """,
            plot=scatter_html,
        )

    # MaxQuant Fig 6
    def draw_summary(self, msms_identified_percent):
        draw_config = {
            "id": "msms_identified_per_raw_file",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "MS/MS Identified per Raw File",
            "tt_decimals": 2,
            "ylab": "MS/MS Identified [%]",
        }
        bar_html = bargraph.plot(msms_identified_percent, pconfig=draw_config)
        self.add_section(
            name="MS/MS Identified per Raw File",
            anchor="msms_identified_per_raw_file",
            description="MS/MS identification rate per Raw file from summary.txt.",
            helptext="MS/MS identification rate per Raw file from summary.txt.",
            plot=bar_html,
        )

    # MaxQuant Fig 7
    def draw_evidence_top_contaminants(self, top_contaminants_data):
        draw_config = {
            "id": "top_contaminants_per_raw_file",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Top5 Contaminants per Raw file",
            "tt_decimals": 2,
            "ylab": "Identified [%]",
        }
        bar_html = bargraph.plot(
            top_contaminants_data["plot_data"],
            cats=top_contaminants_data["cats"],
            pconfig=draw_config,
        )
        self.add_section(
            name="Top5 Contaminants per Raw file",
            anchor="top_contaminants_per_raw_file",
            # description='',
            helptext="""
            pmultiqc will explicitly show the five most abundant external protein contaminants 
            (as detected via MaxQuant's contaminants FASTA file) by Raw file, and summarize the 
            remaining contaminants as 'other'. This allows to track down which proteins exactly 
            contaminate your sample. Low contamination is obviously better.
            
            If you see less than 5 contaminants, it either means there are actually less, or that 
            one (or more) of the shortened contaminant names subsume multiple of the top5 
            contaminants (since they start with the same prefix).
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 9
    def draw_evidence_charge(self, charge_data):
        draw_config = {
            "id": "Charge_state_of_per_Raw_file",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Charge-state of per Raw file",
            "tt_decimals": 0,
            "ylab": "Count",
        }
        bar_html = bargraph.plot(
            charge_data["plot_data"], cats=charge_data["cats"], pconfig=draw_config
        )
        self.add_section(
            name="Charge-state of per Raw file",
            anchor="Charge_state_of_per_Raw_file",
            description="The distribution of the charge-state of the precursor ion, excluding potential contaminants.",
            helptext="The distribution of the charge-state of the precursor ion, excluding potential contaminants.",
            plot=bar_html,
        )

    # MaxQuant Fig 10
    def draw_evidence_modified(self, modified_data):
        draw_config = {
            "id": "modifications_per_raw_file",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Modifications per Raw file",
            "stacking": "group",
            "tt_decimals": 3,
            "ylab": "Occurence [%]",
        }
        bar_html = bargraph.plot(
            modified_data["plot_data"], cats=modified_data["cats"], pconfig=draw_config
        )
        self.add_section(
            name="Modifications per Raw file",
            anchor="modifications_per_raw_file",
            description="Compute an occurence table of modifications (e.g. Oxidation (M)) for all peptides, including the unmodified.",
            helptext="Post-translational modifications contained within the identified peptide sequence.",
            plot=bar_html,
        )

    # MaxQuant Fig 11
    def draw_evidence_rt_count(self, rt_count_data):
        draw_config = {
            "id": "IDs_over_RT",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "IDs over RT",
            "ymin": 0,
            "tt_decimals": 3,
            "ylab": "Counts",
            "xlab": "Retention time [min]",
            "showlegend": True,
        }
        linegraph_html = linegraph.plot(rt_count_data, pconfig=draw_config)
        self.add_section(
            name="IDs over RT",
            anchor="IDs_over_RT",
            description="Distribution of retention time, derived from the evidence table.",
            helptext="The uncalibrated retention time in minutes in the elution profile of the precursor ion, and does not include potential contaminants.",
            plot=linegraph_html,
        )

    # MaxQuant Fig 12
    def draw_evidence_peak_width_rt(self, peak_rt_data):
        draw_config = {
            "id": "peak_width_over_RT",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Peak width over RT",
            "ymin": 0,
            "tt_decimals": 3,
            "ylab": "Retention length [median]",
            "xlab": "Retention time [min]",
            "showlegend": True,
        }
        linegraph_html = linegraph.plot(peak_rt_data, pconfig=draw_config)
        self.add_section(
            name="Peak width over RT",
            anchor="peak_width_over_RT",
            description="Distribution of widths of peptide elution peaks, derived from the evidence table.",
            helptext="""
            The distribution of the widths of peptide elution peaks, 
            derived from the evidence table and excluding potential contaminants, 
            is one parameter of optimal and reproducible chromatographic separation.
            """,
            plot=linegraph_html,
        )

    # MaxQuant Fig 13
    def draw_evidence_oversampling(self, oversampling_data):
        draw_config = {
            "id": "oversampling",
            "cpswitch": False,
            "cpswitch_c_active": False,
            # 'title': 'Oversampling (MS/MS counts per 3D-peak)',
            "title": "Oversampling",
            "tt_decimals": 2,
            "ylab": "MS/MS counts per 3D-peak [%]",
        }
        bar_html = bargraph.plot(
            oversampling_data["plot_data"], cats=oversampling_data["cats"], pconfig=draw_config
        )
        self.add_section(
            # name='Oversampling (MS/MS counts per 3D-peak)',
            name="Oversampling",
            anchor="oversampling",
            description="""
            An oversampled 3D-peak is defined as a peak whose peptide ion (same sequence 
            and same charge state) was identified by at least two distinct MS2 spectra in the same Raw file.
            """,
            helptext="""
            For high complexity samples, oversampling of individual 3D-peaks automatically leads to undersampling 
            or even omission of other 3D-peaks, reducing the number of identified peptides. Oversampling occurs 
            in low-complexity samples or long LC gradients, as well as undersized dynamic exclusion windows 
            for data independent acquisitions.
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 15
    def draw_mass_error_box(self, mass_error_data, fig_type):
        # fig_type: 'uncalibrated' or 'calibrated'
        # 'uncalibrated'
        if fig_type == "uncalibrated":
            draw_config = {
                "id": "uncalibrated_mass_error_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "Uncalibrated Mass Error",
                "tt_decimals": 2,
                "xlab": "Mass Error [ppm]",
            }
            box_html = box.plot(mass_error_data, pconfig=draw_config)
            self.add_section(
                name="Uncalibrated Mass Error",
                anchor="uncalibrated_mass_error_box",
                description="[Excludes Contaminants] Mass accurary before calibration.",
                helptext="""
                Mass error of the uncalibrated mass-over-charge value of the precursor ion in comparison 
                to the predicted monoisotopic mass of the identified peptide sequence. 
                """,
                plot=box_html,
            )
        # calibrated
        if fig_type == "calibrated":
            draw_config = {
                "id": "calibrated_mass_error_box",
                "cpswitch": False,
                "cpswitch_c_active": False,
                "title": "Calibrated Mass Error",
                "tt_decimals": 2,
                "xlab": "Mass Error [ppm]",
            }
            box_html = box.plot(mass_error_data, pconfig=draw_config)
            self.add_section(
                name="Calibrated Mass Error",
                anchor="Calibrated_mass_error_box",
                description="Mass accuracy after calibration (Excludes Contaminants).",
                helptext="""
                Mass error of the recalibrated mass-over-charge value of the precursor ion in comparison 
                to the predicted monoisotopic mass of the identified peptide sequence in parts per million. 
                """,
                plot=box_html,
            )

    # MaxQuant Fig 16
    def draw_evidence_peptide_id_count(self, peptide_id_count_data):
        if peptide_id_count_data["title_value"]:
            fig_title = "Peptide ID Count" + " [" + peptide_id_count_data["title_value"] + "]"
        else:
            fig_title = "Peptide ID Count"
        draw_config = {
            "id": "peptide_id_count",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": fig_title,
            "tt_decimals": 0,
            "ylab": "Count",
        }
        bar_html = bargraph.plot(
            peptide_id_count_data["plot_data"],
            cats=peptide_id_count_data["cats"],
            pconfig=draw_config,
        )
        self.add_section(
            name="Peptide ID Count",
            anchor="peptide_id_count",
            description="""
            [Excludes Contaminants] Number of unique (i.e. not counted twice) peptide sequences including modifications (after FDR) per Raw file.
            """,
            helptext="""
            If MBR was enabled, three categories ('Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)'
            are shown, so the user can judge the gain that MBR provides. Peptides in the 'Genuine + Transferred' category 
            were identified within the Raw file by MS/MS, but at the same time also transferred to this Raw file using MBR. 
            This ID transfer can be correct (e.g. in case of different charge states), or incorrect -- see MBR-related 
            metrics to tell the difference. 
            Ideally, the 'Genuine + Transferred' category should be rather small, the other two should be large.

            If MBR would be switched off, you can expect to see the number of peptides corresponding to 'Genuine (Exclusive)' + 'Genuine + Transferred'. 
            In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics),
            MBR should be switched off for the Raw files which are affected (could be a few or all). 
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 17
    def draw_evidence_protein_group_count(self, protein_group_count_data):
        if protein_group_count_data["title_value"]:
            fig_title = (
                    "ProteinGroups Count" + " [" + protein_group_count_data["title_value"] + "]"
            )
        else:
            fig_title = "ProteinGroups Count"
        draw_config = {
            "id": "protein_group_count",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": fig_title,
            "tt_decimals": 0,
            "ylab": "Count",
        }
        bar_html = bargraph.plot(
            protein_group_count_data["plot_data"],
            cats=protein_group_count_data["cats"],
            pconfig=draw_config,
        )
        self.add_section(
            name="ProteinGroups Count",
            anchor="protein_group_count",
            description="""
            [Excludes Contaminants] Number of Protein groups (after FDR) per Raw file.
            """,
            helptext="""
            If MBR was enabled, three categories ('Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)' 
            are shown, so the user can judge the gain that MBR provides. Here, 'Transferred (Exclusive)' means that this protein group 
            has peptide evidence which originates only from transferred peptide IDs. The quantification is (of course) always from the 
            local Raw file. 
            Proteins in the 'Genuine + Transferred' category have peptide evidence from within the Raw file by MS/MS, but at the same time 
            also peptide IDs transferred to this Raw file using MBR were used. It is not unusual to see the 'Genuine + Transferred' category be the 
            rather large, since a protein group usually has peptide evidence from both sources. 
            To see of MBR worked, it is better to look at the two MBR-related metrics.

            If MBR would be switched off, you can expect to see the number of protein groups corresponding to 'Genuine (Exclusive)' + 'Genuine + Transferred'. 
            In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics), 
            MBR should be switched off for the Raw files which are affected (could be a few or all).
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 18
    def draw_msms_missed_cleavages(self, missed_cleavages_data):
        draw_config = {
            "id": "missed_cleavages_per_raw_file",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Missed Cleavages per Raw file",
            "tt_decimals": 2,
            "ylab": "Missed Cleavages [%]",
        }
        bar_html = bargraph.plot(
            missed_cleavages_data["plot_data"],
            cats=missed_cleavages_data["cats"],
            pconfig=draw_config,
        )
        self.add_section(
            name="Missed Cleavages per Raw file",
            anchor="missed_cleavages_per_raw_file",
            description="[Excludes Contaminants] Missed Cleavages per Raw file.",
            helptext="""
            Under optimal digestion conditions (high enzyme grade etc.), only few missed cleavages (MC) are expected. In 
            general, increased MC counts also increase the number of peptide signals, thus cluttering the available 
            space and potentially provoking overlapping peptide signals, biasing peptide quantification.
            Thus, low MC counts should be favored. Interestingly, it has been shown recently that 
            incorporation of peptides with missed cleavages does not negatively influence protein quantification (see 
            [Chiva, C., Ortega, M., and Sabido, E. Influence of the Digestion Technique, Protease, and Missed 
            Cleavage Peptides in Protein Quantitation. J. Proteome Res. 2014, 13, 3979-86](https://doi.org/10.1021/pr500294d) ). 
            However this is true only if all samples show the same degree of digestion. High missed cleavage values 
            can indicate for example, either a) failed digestion, b) a high (post-digestion) protein contamination, or 
            c) a sample with high amounts of unspecifically degraded peptides which are not digested by trypsin. 

            If MC>=1 is high (>20%) you should increase the missed cleavages settings in MaxQuant and compare the number of peptides.
            Usually high MC correlates with bad identification rates, since many spectra cannot be matched to the forward database.

            In the rare case that 'no enzyme' was specified in MaxQuant, neither scores nor plots are shown.
            """,
            plot=bar_html,
        )

    # MaxQuant Fig 19
    def draw_msms_scans_ion_injec_time_rt(self, ion_injec_time_rt_data):
        draw_config = {
            "id": "ion_injection_time_over_rt",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Ion Injection Time over RT",
            "ymin": 0,
            "tt_decimals": 2,
            "ylab": "Ion injection time [ms]",
            "xlab": "Retention time [min]",
            "showlegend": True,
        }
        linegraph_html = linegraph.plot(ion_injec_time_rt_data, pconfig=draw_config)
        self.add_section(
            name="Ion Injection Time over RT",
            anchor="ion_injection_time_over_rt",
            # description='Ion Injection Time over RT.',
            helptext="""
            Ion injection time score - should be as low as possible to allow fast cycles. Correlated with peptide intensity. 
            Note that this threshold needs customization depending on the instrument used (e.g., ITMS vs. FTMS).
            """,
            plot=linegraph_html,
        )

    # MaxQuant Fig 20
    def draw_msms_scans_top_over_rt(self, top_over_RT_data):
        draw_config = {
            "id": "topn_over_rt",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "TopN over RT",
            "ymin": 0,
            "tt_decimals": 2,
            "ylab": "Highest N [median per RT bin]",
            "xlab": "Retention time [min]",
            "showlegend": True,
        }
        linegraph_html = linegraph.plot(top_over_RT_data, pconfig=draw_config)
        self.add_section(
            name="TopN over RT",
            anchor="topn_over_rt",
            description="TopN over retention time.",
            helptext="""
            TopN over retention time. Similar to ID over RT, this metric reflects the complexity of the sample 
            at any point in time. Ideally complexity should be made roughly equal (constant) by choosing a proper (non-linear) LC gradient. 
            See [Moruz 2014, DOI: 10.1002/pmic.201400036](https://pubmed.ncbi.nlm.nih.gov/24700534/) for details.
            """,
            plot=linegraph_html,
        )

    # MaxQuant Fig 21
    def draw_msms_scans_top_n(self, top_n_data):
        draw_config = {
            "id": "top_n",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "TopN",
            "stacking": "group",
            "tt_decimals": 0,
            "ylab": "Highest Scan Event",
        }
        bar_html = bargraph.plot(
            top_n_data["plot_data"], cats=top_n_data["cats"], pconfig=draw_config
        )
        self.add_section(
            name="TopN",
            anchor="top_n",
            description='This metric somewhat summarizes "TopN over RT"',
            helptext="""
            Reaching TopN on a regular basis indicates that all sections of the LC gradient 
            deliver a sufficient number of peptides to keep the instrument busy. This metric somewhat summarizes "TopN over RT".
            """,
            plot=bar_html,
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
            lambda x: QuantMSModule.file_prefix(x["Spectra_Filepath"]), axis=1
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