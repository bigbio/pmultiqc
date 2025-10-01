""" MultiQC mzIdentML plugin module """

from __future__ import absolute_import
import logging
from collections import OrderedDict
from datetime import datetime
from multiqc import config

from multiqc.plots import table, bargraph

import pandas as pd
import re
from pyteomics import mztab, mzid
from pyopenms import OpenMSBuildInfo
import os
import numpy as np

from pmultiqc.modules.common.plots import *
from pmultiqc.modules.common.id.mzid import MzIdentMLReader
from pmultiqc.modules.common.ms.mgf import MGFReader
from pmultiqc.modules.common.utils import parse_mzml
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.stats import qual_uniform
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.core.section_groups import add_group_modules, add_sub_section
from pmultiqc.modules.mzidentml.mzidentml_plots import draw_mzid_quant_table
from pmultiqc.modules.mzidentml.mzidentml_utils import (
    get_mzidentml_mzml_df,
    get_mzidentml_charge,
    get_mzid_rt_id,
    get_mzid_num_data,
)
from pmultiqc.modules.common.base_module import BaseModule


# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

log.info("pyopenms has: " + str(OpenMSBuildInfo().getOpenMPMaxNumThreads()) + " threads.")


class MzIdentMLModule(BaseModule):

    def __init__(
            self,
            find_log_files_func,
            sub_sections
        ):
        # Call parent constructor
        super().__init__(
            find_log_files_func=find_log_files_func,
            sub_sections=sub_sections,
            module_name="mzidentml"
        )

        # Initialize module-specific attributes
        self._initialize_mzidentml_data()

    def _initialize_module(self):
        """Initialize MzIdentML-specific attributes."""
        # Discover files manually to handle multiple files per type
        self.ms_paths = []
        for mzml_file in self.find_log_files("pmultiqc/mzML", filecontents=False):
            self.ms_paths.append(os.path.join(mzml_file["root"], mzml_file["fn"]))
        self.ms_paths.sort()

        self.mgf_paths = []
        for mgf_file in self.find_log_files("pmultiqc/mgf", filecontents=False):
            self.mgf_paths.append(os.path.join(mgf_file["root"], mgf_file["fn"]))
        self.mgf_paths.sort()

        self.mzid_paths = []
        for mzid_file in self.find_log_files("pmultiqc/mzid", filecontents=False):
            self.mzid_paths.append(os.path.join(mzid_file["root"], mzid_file["fn"]))
        self.mzid_paths.sort()

    def _initialize_mzidentml_data(self):
        """Initialize MzIdentML-specific data structures."""
        self.ms_with_psm = list()
        self.total_protein_identified = 0
        self.cal_num_table_data = dict()
        self.oversampling = dict()
        self.identified_spectrum = dict()
        self.total_ms2_spectral_identified = 0
        self.total_peptide_count = 0
        self.total_ms2_spectra = 0
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
        self.identified_msms_spectra = dict()
        self.mzid_peptide_map = dict()
        self.ms_without_psm = dict()

    def get_data(self):
        """Extract and process MzIdentML data."""
        # Parse mzIdentML files
        mzid_psm = self.parse_out_mzid()

        # Process MGF or mzML files
        if self.mgf_paths:
            self._process_mgf_data(mzid_psm)
        elif self.ms_paths:
            self._process_mzml_data(mzid_psm)

        # Generate plots
        self.draw_plots()

        return self.results

    def _process_mgf_data(self, mzid_psm):
        """Process MGF files and generate related plots."""
        mgf_reader = MGFReader()
        mgf_result = mgf_reader.read(
            file_paths=self.mgf_paths,
            ms_with_psm=self.ms_with_psm,
            identified_spectrum=self.identified_spectrum,
            ms_without_psm=self.ms_without_psm
        )

        # Extract results from MGFReader
        self.mgf_rtinseconds = mgf_result["mgf_rtinseconds"]
        self.mgf_charge_plot = mgf_result["mgf_charge_plot"]
        self.mgf_peak_distribution_plot = mgf_result["mgf_peak_distribution_plot"]
        self.mgf_peaks_ms2_plot = mgf_result["mgf_peaks_ms2_plot"]
        self.mgf_charge_plot_1 = mgf_result["mgf_charge_plot_1"]
        self.mgf_peak_distribution_plot_1 = mgf_result["mgf_peak_distribution_plot_1"]
        self.mgf_peaks_ms2_plot_1 = mgf_result["mgf_peaks_ms2_plot_1"]
        self.ms_info.update(mgf_result["ms_info"])
        self.heatmap_charge_score = mgf_result["heatmap_charge_score"]
        self.total_ms2_spectra = mgf_result["total_ms2_spectra"]

        self.mzid_cal_heat_map_score(mzid_psm)

    def _process_mzml_data(self, mzid_psm):
        """Process mzML files and generate related plots."""
        (
            self.mzml_table,
            self.mzml_peaks_ms2_plot,
            self.mzml_peak_distribution_plot,
            self.ms_info,
            self.total_ms2_spectra,
            mzml_ms_df,
            self.heatmap_charge_score,
            self.mzml_charge_plot,
            _,  # ms1_tic,
            _,  # ms1_bpc,
            _,  # ms1_peaks,
            _,  # ms1_general_stats
        ) = parse_mzml(
            ms_with_psm=self.ms_with_psm,
            identified_spectrum=self.identified_spectrum,
            ms_paths=self.ms_paths,
            enable_mzid=True
        )

        mzidentml_df = get_mzidentml_mzml_df(mzid_psm, mzml_ms_df)
        if len(mzidentml_df) > 0:
            draw_mzid_quant_table(self.sub_sections["quantification"], mzidentml_df)

            mzid_mzml_charge_state = get_mzidentml_charge(mzidentml_df)
            draw_charge_state(
                self.sub_sections["ms2"], mzid_mzml_charge_state, "mzIdentML"
            )

            mzid_ids_over_rt = get_mzid_rt_id(mzidentml_df)
            draw_ids_rt_count(
                self.sub_sections["rt_qc"], mzid_ids_over_rt, "mzIdentML"
            )

            (self.cal_num_table_data, self.identified_msms_spectra) = get_mzid_num_data(
                mzidentml_df
            )

            draw_quantms_identification(
                self.sub_sections["identification"],
                cal_num_table_data=self.cal_num_table_data,
                mzml_table=self.mzml_table,
                quantms_missed_cleavages=self.quantms_missed_cleavages,
                identified_msms_spectra=self.identified_msms_spectra,
            )

            self.mzid_cal_heat_map_score(mzidentml_df)

    def draw_plots(self):
        """Generate plots and add them to report sections."""
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
            total_peptide_count=self.total_peptide_count,
            total_ms2_spectral_identified=self.total_ms2_spectral_identified,
            total_ms2_spectra=self.total_ms2_spectra,
            total_protein_identified=self.total_protein_identified
        )

        self.draw_mzid_identi_num()

        draw_num_pep_per_protein(
            self.sub_sections["identification"],
            self.pep_plot
        )

        if self.mgf_paths:
            charge_plot = self.mgf_charge_plot
            peaks_ms2_plot = self.mgf_peaks_ms2_plot
            peak_distribution_plot = self.mgf_peak_distribution_plot
        else:
            charge_plot = self.mzml_charge_plot
            peaks_ms2_plot = self.mzml_peaks_ms2_plot
            peak_distribution_plot = self.mzml_peak_distribution_plot

        draw_precursor_charge_distribution(
            self.sub_sections["ms2"],
            charge_plot,
            self.ms_info
        )

        draw_peaks_per_ms2(
            self.sub_sections["ms2"],
            peaks_ms2_plot,
            self.ms_info
        )

        draw_peak_intensity_distribution(
            self.sub_sections["ms2"],
            peak_distribution_plot,
            self.ms_info
        )

        draw_oversampling(
            self.sub_sections["ms2"],
            self.oversampling,
            self.oversampling_plot.dict["cats"],
            False
        )

        # Create section groups
        self.section_group_dict = {
            "experiment_sub_section": self.sub_sections["experiment"],
            "summary_sub_section": self.sub_sections["summary"],
            "identification_sub_section": self.sub_sections["identification"],
            "quantification_sub_section": self.sub_sections["quantification"],
            "ms1_sub_section": self.sub_sections["ms1"],
            "ms2_sub_section": self.sub_sections["ms2"],
            "mass_error_sub_section": self.sub_sections["mass_error"],
            "rt_qc_sub_section": self.sub_sections["rt_qc"],
        }

        add_group_modules(self.section_group_dict, "")

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
                """,
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

            spec_e_bar_html = remove_subtitle(spec_e_bar_html)

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

            xcorr_bar_html = remove_subtitle(xcorr_bar_html)

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

            hyper_bar_html = remove_subtitle(hyper_bar_html)

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
        }

        pep_bar_html = bargraph.plot(
            list(self.search_engine["PEPs"].values()), pep_cats, pep_pconfig
        )

        pep_bar_html = remove_subtitle(pep_bar_html)

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
            }

            consensus_bar_html = bargraph.plot(
                self.search_engine["consensus_support"],
                bar_cats,
                consensus_pconfig,
            )
            consensus_bar_html = remove_subtitle(consensus_bar_html)

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

        for name, group in psm.groupby("filename"):
            sc = group["missed_cleavages"].value_counts()

            self.quantms_missed_cleavages[name] = sc.to_dict()

            mis_0 = sc.get(0, 0)
            self.missed_clevages_heatmap_score[name] = mis_0 / sc[:].sum()
            self.id_rt_score[name] = qual_uniform(group["retention_time"])

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
            if (
                len(
                    pep_table[
                        pep_table["accession"].str.contains(config.kwargs["contaminant_affix"])
                    ]
                )
                > 0
            ):
                pass

            for name, group in pep_table.groupby("stand_spectra_ref"):

                contaminant_sum = (
                    group[group["accession"].str.contains(config.kwargs["contaminant_affix"])][
                        study_variables
                    ]
                    .sum(axis=0)
                    .sum()
                )
                all_sum = group[study_variables].sum(axis=0).sum()
                self.heatmap_con_score[name] = 1.0 - (contaminant_sum / all_sum)

                if config.kwargs["remove_decoy"]:
                    pep_median = np.nanmedian(
                        group[(group["opt_global_cv_MS:1002217_decoy_peptide"] == 0)][
                            study_variables
                        ].to_numpy()
                    )
                pep_median = np.nanmedian(group[study_variables].to_numpy())

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


    def parse_out_mzid(self):

        mzid_reader = MzIdentMLReader()
        mzid_table = mzid_reader.read(self.mzid_paths)

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

        psm_cols = [
            "spectrumID",
            "PeptideSequence",
            "chargeState",
            "Modifications",
            "accession_group",
            "experimentalMassToCharge",
            "calculatedMassToCharge",
            "search_engine_score",
            "filename",
            "retention_time",
        ]

        if "retention_time" not in mzid_table.columns:
            psm_cols.remove("retention_time")

        psm = mzid_table[psm_cols].drop_duplicates().reset_index(drop=True)

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

            unique_peptides = group[["PeptideSequence", "Modifications"]].drop_duplicates()

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

        self.total_ms2_spectral_identified = psm["spectrumID"].nunique()
        self.total_peptide_count = psm["PeptideSequence"].nunique()

        return psm

    def cal_quantms_contaminant_percent(self, pep_df):
        pass

    def top_n_contaminant_percent(self, pep_df, top_n):
        pass

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
            protein_count = {
                sample: {"Count": info["protein_num"]}
                for sample, info in self.cal_num_table_data.items()
            }
            peptide_count = {
                sample: {"Count": info["peptide_num"]}
                for sample, info in self.cal_num_table_data.items()
            }
        else:
            return

        bar_html = bargraph.plot(
            protein_count,
            pconfig=draw_config,
        )
        bar_html = remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["identification"],
            plot=bar_html,
            order=3,
            description="Number of protein groups per raw file.",
            helptext="""
                Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
                """,
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
        bar_html = remove_subtitle(bar_html)

        add_sub_section(
            sub_section=self.sub_sections["identification"],
            plot=bar_html,
            order=4,
            description="""
                Number of unique (i.e. not counted twice) peptide sequences including modifications per Raw file.
                """,
            helptext="""
                Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
                """,
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
                    group: count / total * 100 for group, count in counts.items()
                }

            mc_data = {"plot_data": mc_group_ratio, "cats": ["0", "1", ">=2"]}
            draw_msms_missed_cleavages(
                self.sub_sections["identification"], mc_data, False
            )

        # 5.MS/MS Identified per Raw File
        if self.identified_msms_spectra and mzml_table:

            msms_identified_rate = dict()

            for m in self.identified_msms_spectra.keys():
                identified_ms2 = self.identified_msms_spectra[m].get("Identified", 0)
                all_ms2 = mzml_table.get(m, {}).get("MS2_Num", 0)

                if all_ms2 > 0:
                    msms_identified_rate[m] = {"Identified Rate": (identified_ms2 / all_ms2) * 100}

            draw_ms_ms_identified(
                self.sub_sections["identification"], msms_identified_rate
            )

    def draw_quantms_quantification(self):
        pass

    def draw_quantms_msms_section(self):
        # 1.Charge-state of Per File
        pass

    def draw_quantms_time_section(self):
        pass
