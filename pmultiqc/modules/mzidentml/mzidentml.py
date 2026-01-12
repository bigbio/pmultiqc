""" mzIdentML pmultiqc plugin module """

from __future__ import absolute_import

import copy
import os
from collections import OrderedDict
from datetime import datetime

import numpy as np
import pandas as pd
from multiqc import config
from multiqc.plots import table
from pyteomics import mzid, mgf

from pmultiqc.modules.common.mzidentml_utils import (
    get_mzidentml_mzml_df,
    get_mzidentml_charge,
    get_mzid_rt_id,
    draw_mzid_quant_table,
)

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.common.plots.id import (
    draw_charge_state,
    draw_ids_rt_count,
    draw_identification,
    draw_summary_protein_ident_table,
    draw_num_pep_per_protein,
    draw_oversampling,
    draw_long_trends
)

from pmultiqc.modules.common.plots.ms import (
    draw_precursor_charge_distribution,
    draw_peak_intensity_distribution,
    draw_peaks_per_ms2,
    draw_ms_information,
)
from pmultiqc.modules.common.plots.general import draw_heatmap
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.stats import qual_uniform
from pmultiqc.modules.core.section_groups import (
    add_group_modules,
    add_sub_section
)
from pmultiqc.modules.common.common_utils import aggregate_msms_identified_rate


class MzIdentMLModule(BasePMultiqcModule):

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        super().__init__(find_log_files_func, sub_sections, heatmap_colors)

        self.oversampling_plot = None
        self.pep_plot = None
        self.mgf_rtinseconds = None
        self.mgf_peaks_ms2_plot_1 = None
        self.mgf_charge_plot_1 = None
        self.mgf_peak_distribution_plot_1 = None
        self.mgf_peaks_ms2_plot = None
        self.mgf_charge_plot = None
        self.mgf_peak_distribution_plot = None
        self.mzml_ms_df = None
        self.mzml_peaks_ms2_plot_1 = None
        self.mzml_charge_plot_1 = None
        self.mzml_peak_distribution_plot_1 = None
        self.mzml_peaks_ms2_plot = None
        self.mzml_charge_plot = None
        self.mzml_peak_distribution_plot = None
        self.section_group_dict = None
        self.mzid_paths = []
        self.ms_without_psm = None
        self.mzid_peptide_map = None
        self.mgf_paths = None
        self.ms_paths = None
        self.ms_with_psm = []
        self.total_protein_identified = 0
        self.cal_num_table_data = {}
        self.oversampling = {}
        self.identified_spectrum = {}
        self.delta_mass = {}
        self.total_ms2_spectra_identified = 0
        self.total_peptide_count = 0
        self.total_ms2_spectra = 0
        self.heatmap_charge_score = {}
        self.missed_clevages_heatmap_score = {}
        self.id_rt_score = {}
        self.heatmap_over_sampling_score = {}
        self.heatmap_pep_missing_score = {}
        self.missed_cleavages_var_score = {}
        self.ms_info = {}
        self.ms_info["charge_distribution"] = {}
        self.ms_info["peaks_per_ms2"] = {}
        self.ms_info["peak_distribution"] = {}
        self.quantms_missed_cleavages = {}
        self.quantms_modified = {}
        self.identified_msms_spectra = {}
        self.ms1_tic = {}
        self.ms1_bpc = {}
        self.ms1_peaks = {}
        self.ms1_general_stats = {}
        self.long_trends = {}

    def get_data(self) -> bool | None:
        self.log.info("Start parsing the MzIdentML results and spectra files...")

        self.ms_paths = []
        for mzml_current_file in self.find_log_files("pmultiqc/mzML", filecontents=False):
            self.ms_paths.append(os.path.join(mzml_current_file["root"], mzml_current_file["fn"]))

        self.mgf_paths = []
        for mgf_file in self.find_log_files("pmultiqc/mgf", filecontents=False):
            self.mgf_paths.append(os.path.join(mgf_file["root"], mgf_file["fn"]))
        self.mgf_paths.sort()

        self.mzid_peptide_map = dict()
        self.ms_without_psm = dict()

        for mzid_file in self.find_log_files("pmultiqc/mzid", filecontents=False):
            self.mzid_paths.append(os.path.join(mzid_file["root"], mzid_file["fn"]))
        self.mzid_paths.sort()

        mzid_psm = self.parse_out_mzid()

        if self.mgf_paths:
            self.parse_out_mgf()
            self.mzid_cal_heat_map_score(mzid_psm)

        elif self.ms_paths:
            mt = self.parse_mzml()

            mzidentml_df = get_mzidentml_mzml_df(mzid_psm, self.mzml_ms_df)
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

                msms_identified_rate = aggregate_msms_identified_rate(
                    mzml_table=mt,
                    identified_msms_spectra=self.identified_msms_spectra,
                    sdrf_file_df=None
                )

                self.mzid_cal_heat_map_score(mzidentml_df)

                draw_identification(
                    self.sub_sections["identification"],
                    cal_num_table_data=self.cal_num_table_data,
                    quantms_missed_cleavages=self.quantms_missed_cleavages,
                    quantms_modified=self.quantms_modified,
                    msms_identified_rate=msms_identified_rate
                )

        return True

    def draw_plots(self) -> None:
        self.log.info("Start plotting the MzIdentML results...")

        heatmap_data, heatmap_xnames, heatmap_ynames = self.calculate_heatmap()
        draw_heatmap(
            self.sub_sections["summary"],
            self.heatmap_color_list,
            heatmap_data,
            heatmap_xnames,
            heatmap_ynames,
            "",
        )

        draw_ms_information(
            self.sub_sections["ms1"],
            self.ms1_tic,
            self.ms1_bpc,
            self.ms1_peaks,
            self.ms1_general_stats
        )

        draw_summary_protein_ident_table(
            sub_sections=self.sub_sections["summary"],
            total_peptide_count=self.total_peptide_count,
            total_ms2_spectra_identified=self.total_ms2_spectra_identified,
            total_ms2_spectra=self.total_ms2_spectra,
            total_protein_identified=self.total_protein_identified,
            enable_mzid=True
        )

        self.draw_mzid_identi_num()

        draw_num_pep_per_protein(
            self.sub_sections["identification"],
            self.pep_plot,
            True
        )

        if self.mgf_paths:
            draw_precursor_charge_distribution(
                self.sub_sections["ms2"],
                charge_plot=self.mgf_charge_plot,
                ms_info=self.ms_info
            )

            draw_peaks_per_ms2(
                self.sub_sections["ms2"],
                self.mgf_peaks_ms2_plot,
                self.ms_info
            )

            draw_peak_intensity_distribution(
                self.sub_sections["ms2"],
                self.mgf_peak_distribution_plot,
                self.ms_info
            )
        else:
            draw_precursor_charge_distribution(
                self.sub_sections["ms2"],
                charge_plot=self.mzml_charge_plot,
                ms_info=self.ms_info
            )

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
            self.sub_sections["ms2"], self.oversampling, self.oversampling_plot.dict["cats"], False
        )

        if self.long_trends:
            draw_long_trends(
                sub_sections=self.sub_sections,
                long_trends_data=self.long_trends
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

    def calculate_heatmap(self):

        heat_map_score = []
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

    def draw_mzid_identi_num(self):
        pconfig = {
            "id": "result_statistics",  # ID used for the table
            "title": "Pipeline Result Statistics",  # Title of the table. Used in the column config modal
            "save_file": False,  # Whether to save the table data to a file
            "raw_data_fn": "multiqc_result_statistics_table",  # File basename to use for raw data file
            "sort_rows": True,  # Whether to sort rows alphabetically
            "only_defined_headers": False,  # Only show columns that are defined in the headers config
            "col1_header": "Spectra File",
            "no_violin": True,
            "save_data_file": False,
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
        table_html = table.plot(self.cal_num_table_data["ms_runs"], headers, pconfig)
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

    def mzid_cal_heat_map_score(self, psm):
        self.log.info("{}: Calculating Heatmap Scores...".format(datetime.now().strftime("%H:%M:%S")))

        # HeatMapMissedCleavages
        global_peps = psm[["PeptideSequence", "Modifications"]].drop_duplicates()
        global_peps_count = len(global_peps)

        enzyme_list = list()
        for mzid_path in self.mzid_paths:
            try:
                enzyme_iter = mzid.MzIdentML(mzid_path).iterfind("Enzyme")
                enzyme = next(enzyme_iter).get("EnzymeName", None)
                if enzyme:
                    enzyme_name = next(iter(enzyme.keys()))
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

        missed_cleavages_by_run = dict()

        for name, group in psm.groupby("filename"):
            sc = group["missed_cleavages"].value_counts()

            missed_cleavages_by_run[name] = sc.to_dict()

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

        self.quantms_missed_cleavages = {
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
                strict=True
            )
        )
        self.log.info(
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

    def parse_mzml(self):

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

        # Use the refactored class from ms/mzml.py
        from pmultiqc.modules.common.ms.mzml import MzMLReader

        mzml_reader = MzMLReader(
            file_paths=self.ms_paths,
            ms_with_psm=self.ms_with_psm,
            identified_spectrum=self.identified_spectrum,
            mzml_charge_plot=self.mzml_charge_plot,
            mzml_peak_distribution_plot=self.mzml_peak_distribution_plot,
            mzml_peaks_ms2_plot=self.mzml_peaks_ms2_plot,
            mzml_charge_plot_1=self.mzml_charge_plot_1,
            mzml_peak_distribution_plot_1=self.mzml_peak_distribution_plot_1,
            mzml_peaks_ms2_plot_1=self.mzml_peaks_ms2_plot_1,
            ms_without_psm=self.ms_without_psm,
            enable_dia=False,
            enable_mzid=True
        )

        mzml_reader.parse()

        mzml_table = mzml_reader.mzml_table
        heatmap_charge = mzml_reader.heatmap_charge
        self.total_ms2_spectra = mzml_reader.total_ms2_spectra
        self.mzml_ms_df = mzml_reader.mzml_ms_df
        self.ms1_tic = mzml_reader.ms1_tic
        self.ms1_bpc = mzml_reader.ms1_bpc
        self.ms1_peaks = mzml_reader.ms1_peaks
        self.ms1_general_stats = mzml_reader.ms1_general_stats
        self.long_trends = mzml_reader.long_trends

        for i in self.ms_without_psm:
            self.log.warning("No PSM found in '{}'!".format(i))

        self.mzml_peaks_ms2_plot.to_dict()
        self.mzml_peak_distribution_plot.to_dict()
        # Construct compound dictionaries to apply to drawing functions.

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

    def parse_out_mgf(self) -> None:
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
            self.log.info("{}: Parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
            mgf_data = mgf.MGF(m)
            self.log.info(
                "{}: Done parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )
            m = file_prefix(m)
            self.log.info(
                "{}: Aggregating MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )

            charge_2 = 0
            ms2_number = 0
            for i, spectrum in enumerate(mgf_data):

                charge_state = 0
                charge_list = spectrum.get("params", {}).get("charge", [])
                if charge_list:
                    charge_state = int(charge_list[0])
                    ms2_number += 1

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

            heatmap_charge[m] = charge_2 / ms2_number if ms2_number > 0 else 0
            self.total_ms2_spectra = self.total_ms2_spectra + ms2_number
            self.log.info(
                "{}: Done aggregating MGF file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )

        self.mgf_rtinseconds = pd.DataFrame(mgf_rtinseconds)

        for i in self.ms_without_psm:
            self.log.warning("No PSM found in '{}'!".format(i))

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

        # Use the class-based MzID reader
        from pmultiqc.modules.common.ms.mzid import MzidReader

        mzid_reader = MzidReader(self.mzid_paths)
        mzid_reader.parse()
        mzid_table = mzid_reader.filtered_mzid_df

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

        num_table_at_run = dict()
        identified_msms_spectra = dict()

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

            num_table_at_run[m] = {"protein_num": len(proteins)}
            num_table_at_run[m]["peptide_num"] = len(peptides)
            num_table_at_run[m]["unique_peptide_num"] = len(unique_peptides)

            modified_pep = peptides.dropna(subset=["Modifications"])
            num_table_at_run[m]["modified_peptide_num"] = len(modified_pep)

            identified_msms_spectra[m] = {"Identified": len(set(group["spectrumID"]))}

        self.identified_msms_spectra = identified_msms_spectra

        if self.mgf_paths:
            self.ms_without_psm = set([file_prefix(i) for i in self.mgf_paths]) - set(
                self.ms_with_psm
            )
        elif self.ms_paths:
            self.ms_without_psm = set([file_prefix(i) for i in self.ms_paths]) - set(
                self.ms_with_psm
            )
        for i in self.ms_without_psm:
            num_table_at_run[file_prefix(i)] = {
                "protein_num": 0,
                "peptide_num": 0,
                "unique_peptide_num": 0,
                "modified_peptide_num": 0,
            }

        self.cal_num_table_data = {
            "ms_runs": num_table_at_run
        }

        self.total_ms2_spectra_identified = psm["spectrumID"].nunique()
        self.total_peptide_count = psm["PeptideSequence"].nunique()

        return psm
