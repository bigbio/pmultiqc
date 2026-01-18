import os
from datetime import datetime

from pmultiqc.modules.maxquant import (
    maxquant_utils,
    maxquant_io,
    maxquant_plots
)
from pmultiqc.modules.common.plots import id as id_plots
from pmultiqc.modules.common.plots.general import (
    draw_heatmap,
    draw_search_engine_scores
)
from pmultiqc.modules.core.section_groups import add_group_modules
from pmultiqc.modules.base import BasePMultiqcModule


class MaxQuantModule(BasePMultiqcModule):

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        super().__init__(find_log_files_func, sub_sections, heatmap_colors)
        self.section_group_dict = None

        self.maxquant_paths = None
        self.mq_results = None

    def get_data(self) -> bool | None:
        """Process MaxQuant data files and populate results."""
        self.maxquant_paths = maxquant_io.maxquant_file_path(self.find_log_files)

        # Process different file types
        get_parameter_dicts = self._process_parameters_file()
        self._process_sdrf_file()
        get_protegroups_dicts = self._process_protein_groups_file()
        ms_ms_identified = self._process_summary_file()
        get_evidence_dicts = self._process_evidence_file()
        get_msms_dicts = self._process_msms_file(get_evidence_dicts)
        get_msms_scans_dicts = self._process_msms_scans_file()
        maxquant_heatmap = self._calculate_heatmap(get_evidence_dicts, get_msms_dicts)

        # Store results
        self.mq_results = {
            "get_parameter_dicts": get_parameter_dicts,
            "get_protegroups_dicts": get_protegroups_dicts,
            "ms_ms_identified": ms_ms_identified,
            "get_evidence_dicts": get_evidence_dicts,
            "get_msms_dicts": get_msms_dicts,
            "get_msms_scans_dicts": get_msms_scans_dicts,
            "maxquant_heatmap": maxquant_heatmap,
        }

        return bool(self.mq_results)

    def _process_sdrf_file(self):
        """Process SDRF file if present."""

        if "sdrf" not in self.maxquant_paths.keys():
            self.log.info("{}: No SDRF file found.".format(
                    datetime.now().strftime("%H:%M:%S")
                )
            )
            return  # no SDRF file found in txt folder

        try:
            self.log.info(
                "{}: Parsing SDRF file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["sdrf"]
                )
            )
            sdrf_data = maxquant_utils.read_sdrf(self.maxquant_paths["sdrf"])
            if sdrf_data is not None:
                maxquant_plots.draw_exp_design(sdrf_data, self.sub_sections["experiment"])
        except Exception as e:
            self.log.warning(f"Error occurred while draw_exp_design: {e}")

    def _process_parameters_file(self):
        """Process parameters.txt file if present."""
        get_parameter_dicts = {"parameters_tb_dict": None}

        if "parameters" not in self.maxquant_paths.keys():
            return get_parameter_dicts

        self.log.info(
            "{}: Parsing parameters file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
            )
        )

        if os.path.getsize(self.maxquant_paths["parameters"]) == 0:
            self.log.warning("parameters.txt is empty. Please check.")
            return get_parameter_dicts

        try:
            get_parameter_dicts = maxquant_utils.get_parameters(
                file_path=self.maxquant_paths["parameters"]
            )
            self.log.info(
                "{}: Completed the processing of the parameters file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
                )
            )
        except Exception as e:
            self.log.warning(f"Error occurred while reading parameters.txt: {e}")

        return get_parameter_dicts

    def _process_protein_groups_file(self):
        """Process proteinGroups.txt file if present."""
        get_protegroups_dicts = {
            "pg_contaminant": None,
            "pg_intensity_distri": None,
            "pg_lfq_intensity_distri": None,
            "raw_intensity_pca": None,
            "lfq_intensity_pca": None,
            "protein_summary": None,
            "num_pep_per_protein_dict": None,
        }

        if "proteinGroups" not in self.maxquant_paths.keys():
            return get_protegroups_dicts

        self.log.info(
            "{}: Parsing proteinGroups file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["proteinGroups"]
            )
        )

        if os.path.getsize(self.maxquant_paths["proteinGroups"]) == 0:
            self.log.warning("proteinGroups.txt is empty. Please check.")
            return get_protegroups_dicts

        try:
            get_protegroups_dicts = maxquant_utils.get_protegroups(
                file_path=self.maxquant_paths["proteinGroups"]
            )
            self.log.info(
                "{}: Completed the processing of the proteinGroups file {}...".format(
                    datetime.now().strftime("%H:%M:%S"),
                    self.maxquant_paths["proteinGroups"],
                )
            )
        except Exception as e:
            self.log.warning(f"Error occurred while reading proteinGroups.txt: {e}")

        return get_protegroups_dicts

    def _process_summary_file(self):
        """Process summary.txt file if present."""
        if "summary" not in self.maxquant_paths.keys():
            return None

        self.log.info(
            "{}: Parsing summary file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
            )
        )

        if os.path.getsize(self.maxquant_paths["summary"]) == 0:
            self.log.warning("summary.txt is empty. Please check.")
            return None

        try:
            ms_ms_identified = maxquant_utils.get_summary(
                file_path=self.maxquant_paths["summary"]
            )
            self.log.info(
                "{}: Completed the processing of the summary file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
                )
            )
            return ms_ms_identified
        except Exception as e:
            self.log.warning(f"Error occurred while reading summary.txt: {e}")
            return None

    def _process_evidence_file(self):
        """Process evidence.txt file if present."""
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
            "summary_identified_msms_count": None,
            "summary_identified_peptides": None,
            "maxquant_delta_mass_da": None,
            "peptides_quant_table": None,
            "protein_quant_table": None,
        }

        if "evidence" not in self.maxquant_paths.keys():
            return get_evidence_dicts

        self.log.info(
            "{}: Parsing evidence file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
            )
        )

        if os.path.getsize(self.maxquant_paths["evidence"]) == 0:
            self.log.warning("evidence.txt is empty. Please check.")
            return get_evidence_dicts

        try:
            get_evidence_dicts = maxquant_utils.get_evidence(
                file_path=self.maxquant_paths["evidence"]
            )
            self.log.info(
                "{}: Completed the processing of the evidence file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
                )
            )
        except Exception as e:
            self.log.warning(f"Error occurred while reading evidence.txt: {e}")

        return get_evidence_dicts

    def _process_msms_file(self, get_evidence_dicts):
        """Process msms.txt file if present."""
        get_msms_dicts = {
            "missed_cleavages": None,
            "search_engine_scores": None,
        }

        if "msms" not in self.maxquant_paths.keys():
            return get_msms_dicts

        self.log.info(
            "{}: Parsing msms file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
            )
        )

        if os.path.getsize(self.maxquant_paths["msms"]) == 0:
            self.log.warning("msms.txt is empty. Please check.")
            return get_msms_dicts

        try:
            get_msms_dicts = maxquant_utils.get_msms(
                file_path=self.maxquant_paths["msms"],
                evidence_df=get_evidence_dicts["evidence_df"],
            )
            self.log.info(
                "{}: Completed the processing of the msms file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
                )
            )
        except Exception as e:
            self.log.warning(f"Error occurred while reading msms.txt: {e}")

        return get_msms_dicts

    def _process_msms_scans_file(self):
        """Process msScans.txt or msmsScans.txt file if present."""
        get_msms_scans_dicts = {
            "ion_injec_time_rt": None,
            "top_n": None,
            "top_over_rt": None,
            "summary_msms_spectra": None,
        }

        # Determine which file to use
        if "msmsScans" in self.maxquant_paths.keys():
            msms_scans_file = "msmsScans"
        elif "msScans" in self.maxquant_paths.keys():
            msms_scans_file = "msScans"
        else:
            return get_msms_scans_dicts

        self.log.info(
            "{}: Parsing msScans file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.maxquant_paths[msms_scans_file]
            )
        )

        if os.path.getsize(self.maxquant_paths[msms_scans_file]) == 0:
            self.log.warning(f"{msms_scans_file}.txt is empty. Please check.")
            return get_msms_scans_dicts

        try:
            get_msms_scans_dicts = maxquant_utils.get_msms_scans(
                file_path=self.maxquant_paths[msms_scans_file]
            )
            self.log.info(
                "{}: Completed the processing of the msScans file {}...".format(
                    datetime.now().strftime("%H:%M:%S"),
                    self.maxquant_paths[msms_scans_file],
                )
            )
        except Exception as e:
            self.log.warning(f"Error occurred while reading {msms_scans_file}.txt: {e}")

        return get_msms_scans_dicts

    def _calculate_heatmap(self, get_evidence_dicts, get_msms_dicts):
        """Calculate heatmap if required data is available."""
        if not (
            get_evidence_dicts.get("evidence_df") is not None
            and get_evidence_dicts.get("oversampling") is not None
            and get_msms_dicts.get("missed_cleavages") is not None
        ):
            return None

        try:
            if (
                get_evidence_dicts["oversampling"].get("plot_data") is not None
                and get_msms_dicts["missed_cleavages"].get("plot_data") is not None
            ):
                return maxquant_utils.calculate_heatmap(
                    evidence_df=get_evidence_dicts["evidence_df"],
                    oversampling=get_evidence_dicts["oversampling"]["plot_data"],
                    msms_missed_cleavages=get_msms_dicts["missed_cleavages"]["plot_data"],
                )
        except Exception as e:
            self.log.warning(f"Error occurred while calculating heatmap: {e}")

        return None

    def draw_plots(self) -> None:
        """Draw all MaxQuant plots and set up section groups."""
        # Draw different categories of plots
        self._draw_parameter_plots()
        self._draw_quantification_plots()
        self._draw_identification_plots()
        self._draw_contaminant_plots()
        self._draw_mass_error_plots()
        self._draw_summary_plots()
        self._draw_rt_qc_plots()

        # Set up section groups
        self._setup_section_groups()

    def _draw_parameter_plots(self):
        """Draw parameter-related plots."""
        if self.mq_results["get_parameter_dicts"].get("parameters_tb_dict") is not None:
            self._safe_draw(
                maxquant_plots.draw_parameters,
                self.sub_sections["experiment"],
                self.mq_results["get_parameter_dicts"]["parameters_tb_dict"],
                error_name="draw_parameters"
            )

    def _draw_quantification_plots(self):
        """Draw quantification-related plots."""
        # HeatMap
        if self.mq_results["maxquant_heatmap"]:
            self._safe_draw(
                draw_heatmap,
                self.sub_sections["summary"],
                self.heatmap_color_list,
                self.mq_results["maxquant_heatmap"],
                "",
                "",
                "maxquant",
                error_name="draw_heatmap"
            )

        # Quantification tables
        self._safe_draw_if_exists(
            maxquant_plots.draw_peptide_table,
            self.sub_sections["quantification"],
            self.mq_results["get_evidence_dicts"].get("peptides_quant_table"),
            error_name="draw_peptide_table"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_protein_table,
            self.sub_sections["quantification"],
            self.mq_results["get_evidence_dicts"].get("protein_quant_table"),
            error_name="draw_protein_table"
        )

        # Intensity plots
        self._draw_intensity_plots()

    def _draw_intensity_plots(self):
        """Draw intensity-related plots."""
        # Protein group intensity plots
        self._safe_draw_if_exists(
            maxquant_plots.draw_intensity_box,
            self.sub_sections["quantification"],
            self.mq_results["get_protegroups_dicts"].get("pg_intensity_distri", {}),
            "intensity",
            error_name="draw_intensity_box"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_intensity_box,
            self.sub_sections["quantification"],
            self.mq_results["get_protegroups_dicts"].get("pg_lfq_intensity_distri", {}),
            "lfq_intensity",
            error_name="draw_intensity_box"
        )

        # PCA plots
        self._safe_draw_if_exists(
            maxquant_plots.draw_pg_pca,
            self.sub_sections["quantification"],
            self.mq_results["get_protegroups_dicts"].get("raw_intensity_pca"),
            "raw_intensity",
            error_name="draw_pg_pca"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_pg_pca,
            self.sub_sections["quantification"],
            self.mq_results["get_protegroups_dicts"].get("lfq_intensity_pca"),
            "lfq_intensity",
            error_name="draw_pg_pca"
        )

        # Peptide intensity
        self._safe_draw_if_exists(
            maxquant_plots.draw_intensity_box,
            self.sub_sections["quantification"],
            self.mq_results["get_evidence_dicts"].get("peptide_intensity", {}),
            "peptide_intensity",
            error_name="draw_intensity_box"
        )

    def _draw_identification_plots(self):
        """Draw identification-related plots."""
        # MS/MS identified
        if self.mq_results["ms_ms_identified"]:
            self._safe_draw(
                id_plots.draw_ms_ms_identified,
                self.sub_sections["identification"],
                self.mq_results["ms_ms_identified"],
                error_name="draw_ms_ms_identified"
            )

        # Charge state
        self._safe_draw_if_exists(
            id_plots.draw_charge_state,
            self.sub_sections["ms2"],
            self.mq_results["get_evidence_dicts"].get("charge_counts"),
            "MaxQuant",
            error_name="draw_charge_state"
        )

        # Modifications
        self._safe_draw_if_exists(
            id_plots.draw_modifications,
            self.sub_sections["identification"],
            self.mq_results["get_evidence_dicts"].get("modified_percentage"),
            error_name="draw_modifications"
        )

        # Peptide and protein counts
        self._safe_draw_if_exists(
            maxquant_plots.draw_evidence_peptide_id_count,
            self.sub_sections["identification"],
            self.mq_results["get_evidence_dicts"].get("peptide_id_count"),
            "maxquant",
            error_name="draw_evidence_peptide_id_count"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_evidence_protein_group_count,
            self.sub_sections["identification"],
            self.mq_results["get_evidence_dicts"].get("protein_group_count"),
            error_name="draw_evidence_protein_group_count"
        )

        # Oversampling
        self._safe_draw_if_exists(
            id_plots.draw_oversampling,
            self.sub_sections["ms2"],
            self.mq_results["get_evidence_dicts"].get("oversampling"),
            "",
            "maxquant",
            error_name="draw_oversampling"
        )

        # Missed cleavages
        self._safe_draw_if_exists(
            id_plots.draw_msms_missed_cleavages,
            self.sub_sections["identification"],
            self.mq_results["get_msms_dicts"].get("missed_cleavages"),
            True,
            error_name="draw_msms_missed_cleavages"
        )

        # Number of peptides per protein
        self._safe_draw_if_exists(
            maxquant_plots.draw_maxquant_num_pep_pro,
            self.sub_sections["identification"],
            self.mq_results["get_protegroups_dicts"].get("num_pep_per_protein_dict"),
            error_name="draw_maxquant_num_pep_pro"
        )

        # Search engine scores
        self._safe_draw_if_exists(
            draw_search_engine_scores,
            self.sub_sections["search_engine"],
            self.mq_results["get_msms_dicts"].get("search_engine_scores"),
            "maxquant",
            error_name="draw_search_engine_scores"
        )

    def _draw_contaminant_plots(self):
        """Draw contaminant-related plots."""
        self._safe_draw_if_exists(
            id_plots.draw_potential_contaminants,
            self.sub_sections["contaminants"],
            self.mq_results["get_protegroups_dicts"].get("pg_contaminant"),
            "maxquant",
            error_name="draw_potential_contaminants"
        )

        self._safe_draw_if_exists(
            id_plots.draw_top_n_contaminants,
            self.sub_sections["contaminants"],
            self.mq_results["get_evidence_dicts"].get("top_contaminants"),
            error_name="draw_top_n_contaminants"
        )

    def _draw_mass_error_plots(self):
        """Draw mass error-related plots."""
        # Uncalibrated mass error
        self._safe_draw_if_exists(
            maxquant_plots.draw_mass_error_box,
            self.sub_sections["mass_error"],
            self.mq_results["get_evidence_dicts"].get("uncalibrated_mass_error"),
            error_name="draw_mass_error_box"
        )

        # Delta mass plots
        self._safe_draw_if_exists(
            id_plots.draw_delta_mass_da_ppm,
            self.sub_sections["mass_error"],
            self.mq_results["get_evidence_dicts"].get("maxquant_delta_mass_da"),
            "Mass Error [Da]",
            error_name="draw_delta_mass_da_ppm"
        )

        self._safe_draw_if_exists(
            id_plots.draw_delta_mass_da_ppm,
            self.sub_sections["mass_error"],
            self.mq_results["get_evidence_dicts"].get("calibrated_mass_error"),
            "Mass Error [ppm]",
            error_name="draw_delta_mass_da_ppm"
        )

    def _draw_summary_plots(self):
        """Draw summary-related plots."""
        # Summary table
        if (
            self.mq_results["get_msms_scans_dicts"].get("summary_msms_spectra") is not None
            and self.mq_results["get_evidence_dicts"].get("summary_stat") is not None
            and self.mq_results["get_protegroups_dicts"].get("protein_summary") is not None
        ):
            summary_stat = self.mq_results["get_evidence_dicts"]["summary_stat"]

            if (
                summary_stat.get("summary_identified_msms_count") is not None
                and summary_stat.get("summary_identified_peptides") is not None
            ):
                self._safe_draw(
                    maxquant_plots.draw_maxquant_summary_table,
                    self.sub_sections["summary"],
                    self.mq_results["get_msms_scans_dicts"]["summary_msms_spectra"],
                    summary_stat["summary_identified_msms_count"],
                    summary_stat["summary_identified_peptides"],
                    self.mq_results["get_protegroups_dicts"]["protein_summary"],
                    error_name="draw_maxquant_summary_table"
                )

    def _draw_rt_qc_plots(self):
        """Draw RT QC-related plots."""
        # RT counts
        self._safe_draw_if_exists(
            id_plots.draw_ids_rt_count,
            self.sub_sections["rt_qc"],
            self.mq_results["get_evidence_dicts"].get("rt_counts"),
            "maxquant",
            error_name="draw_ids_rt_count"
        )

        # Peak RT
        self._safe_draw_if_exists(
            maxquant_plots.draw_evidence_peak_width_rt,
            self.sub_sections["rt_qc"],
            self.mq_results["get_evidence_dicts"].get("peak_rt"),
            error_name="draw_evidence_peak_width_rt"
        )

        # TopN plots
        self._safe_draw_if_exists(
            maxquant_plots.draw_msms_scans_top_n,
            self.sub_sections["rt_qc"],
            self.mq_results["get_msms_scans_dicts"].get("top_n"),
            error_name="draw_msms_scans_top_n"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_msms_scans_top_over_rt,
            self.sub_sections["rt_qc"],
            self.mq_results["get_msms_scans_dicts"].get("top_over_rt"),
            error_name="draw_msms_scans_top_over_rt"
        )

        self._safe_draw_if_exists(
            maxquant_plots.draw_msms_scans_ion_injec_time_rt,
            self.sub_sections["rt_qc"],
            self.mq_results["get_msms_scans_dicts"].get("ion_injec_time_rt"),
            error_name="draw_msms_scans_ion_injec_time_rt"
        )

    def _setup_section_groups(self):
        """Set up section groups for the module."""
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

    def _safe_draw(self, draw_func, *args, error_name=None):
        """Safely execute a drawing function with error handling."""
        try:
            draw_func(*args)
        except Exception as e:
            self.log.warning(f"Error occurred while {error_name}: {e}")

    def _safe_draw_if_exists(self, draw_func, *args, error_name=None):
        """Safely execute a drawing function if data exists."""
        if args[1] is not None:  # Check if data exists (second argument is usually the data)
            self._safe_draw(draw_func, *args, error_name=error_name)