import logging
from datetime import datetime
import os

from . import maxquant_utils, maxquant_io, maxquant_plots
from ..core.section_groups import add_group_modules
from ..common import common_plots


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

class MaxQuantModule:

    def __init__(
            self,
            find_log_files_func,
            sub_sections,
            heatmap_colors
        ):

        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections

        self.maxquant_paths = None
        self.mq_results = None
        self.heatmap_color_list = heatmap_colors
    
    def get_mq_data(self):

        # parameters.txt
        get_parameter_dicts = {
            "parameters_tb_dict": None
        }

        if "parameters" in self.maxquant_paths.keys():
            log.info(
                "{}: Parsing parameters file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
                )
            )
            if os.path.getsize(self.maxquant_paths["parameters"]) == 0:
                log.warning("parameters.txt is empty. Please check.")
            else:
                try:
                    get_parameter_dicts = maxquant_utils.get_parameters(
                        file_path=self.maxquant_paths["parameters"]
                    )
                    log.info(
                        "{}: Completed the processing of the parameters file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["parameters"]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading parameters.txt: {e}")


        # proteinGroups.txt
        get_protegroups_dicts = {
            "pg_contaminant": None,
            "pg_intensity_distri": None,
            "pg_lfq_intensity_distri": None,
            "raw_intensity_pca": None,
            "lfq_intensity_pca": None,
            "protein_summary": None,
            "num_pep_per_protein_dict": None,
        }

        if "proteinGroups" in self.maxquant_paths.keys():
            log.info(
                "{}: Parsing proteinGroups file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["proteinGroups"]
                )
            )
            if os.path.getsize(self.maxquant_paths["proteinGroups"]) == 0:
                log.warning("proteinGroups.txt is empty. Please check.")
            else:
                try:
                    get_protegroups_dicts = maxquant_utils.get_protegroups(
                        file_path=self.maxquant_paths["proteinGroups"]
                    )
                    log.info(
                        "{}: Completed the processing of the proteinGroups file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["proteinGroups"]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading proteinGroups.txt: {e}")


        # summary.txt
        ms_ms_identified = None
        
        if "summary" in self.maxquant_paths.keys():
            log.info(
                "{}: Parsing summary file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
                )
            )
            if os.path.getsize(self.maxquant_paths["summary"]) == 0:
                log.warning("summary.txt is empty. Please check.")
            else:
                try:
                    ms_ms_identified = maxquant_utils.get_summary(file_path=self.maxquant_paths["summary"])
                    log.info(
                        "{}: Completed the processing of the summary file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["summary"]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading summary.txt: {e}")


        # evidence.txt
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
        
        if "evidence" in self.maxquant_paths.keys():
            log.info(
                "{}: Parsing evidence file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
                )
            )
            if os.path.getsize(self.maxquant_paths["evidence"]) == 0:
                log.warning("evidence.txt is empty. Please check.")
            else:
                try:
                    get_evidence_dicts = maxquant_utils.get_evidence(
                        file_path=self.maxquant_paths["evidence"]
                    )
                    log.info(
                        "{}: Completed the processing of the evidence file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["evidence"]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading evidence.txt: {e}")


        # msms.txt
        get_msms_dicts = {
            "missed_cleavages": None,
            "search_engine_scores": None,
        }
        if "msms" in self.maxquant_paths.keys():
            log.info(
                "{}: Parsing msms file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
                )
            )
            if os.path.getsize(self.maxquant_paths["msms"]) == 0:
                log.warning("msms.txt is empty. Please check.")
            else:
                try:
                    get_msms_dicts = maxquant_utils.get_msms(
                        file_path=self.maxquant_paths["msms"],
                        evidence_df=get_evidence_dicts["evidence_df"],
                    )
                    log.info(
                        "{}: Completed the processing of the msms file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths["msms"]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading msms.txt: {e}")


        # msScans.txt or msmsScans.txt
        get_msms_scans_dicts = {
            "ion_injec_time_rt": None,
            "top_n": None,
            "top_over_rt": None,
            "summary_msms_spectra": None,
        }

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
            if os.path.getsize(self.maxquant_paths[msms_scans_file]) == 0:
                log.warning(f"{msms_scans_file}.txt is empty. Please check.")
            else:
                try:
                    get_msms_scans_dicts = maxquant_utils.get_msms_scans(
                        file_path=self.maxquant_paths[msms_scans_file]
                    )
                    log.info(
                        "{}: Completed the processing of the msScans file {}...".format(
                            datetime.now().strftime("%H:%M:%S"), self.maxquant_paths[msms_scans_file]
                        )
                    )
                except Exception as e:
                    log.warning(f"Error occurred while reading {msms_scans_file}.txt: {e}")


        # HeatMap
        maxquant_heatmap = None
        if (
            get_evidence_dicts.get("evidence_df") is not None and
            get_evidence_dicts.get("oversampling") is not None and
            get_msms_dicts.get("missed_cleavages") is not None
        ):
            try:
                maxquant_heatmap = maxquant_utils.calculate_heatmap(
                    evidence_df=get_evidence_dicts["evidence_df"],
                    oversampling=get_evidence_dicts["oversampling"]["plot_data"],
                    msms_missed_cleavages=get_msms_dicts["missed_cleavages"]["plot_data"]
                    )
            except Exception as e:
                log.warning(f"Error occurred while calculate heatmap: {e}")  

        return {
            "get_parameter_dicts": get_parameter_dicts,
            "get_protegroups_dicts": get_protegroups_dicts,
            "ms_ms_identified": ms_ms_identified,
            "get_evidence_dicts": get_evidence_dicts,
            "get_msms_dicts": get_msms_dicts,
            "get_msms_scans_dicts": get_msms_scans_dicts,
            "maxquant_heatmap": maxquant_heatmap,
        }
    
    def get_data(self):

        self.maxquant_paths = maxquant_io.maxquant_file_path(self.find_log_files)
        self.mq_results = self.get_mq_data()

        if self.mq_results:
            return self.mq_results
        
    def draw_report_plots(self):

        # Parameters
        if self.mq_results["get_parameter_dicts"]["parameters_tb_dict"]:
            maxquant_plots.draw_parameters(
                self.sub_sections["experiment"],
                self.mq_results["get_parameter_dicts"]["parameters_tb_dict"]
            )

        # HeatMap
        if self.mq_results["maxquant_heatmap"]:
            common_plots.draw_heatmap(
                self.sub_sections["summary"],
                self.heatmap_color_list,
                self.mq_results["maxquant_heatmap"],
                "",
                "",
                True
            )

        # Peptides Quantification Table
        if self.mq_results["get_evidence_dicts"]["peptides_quant_table"]:
            maxquant_plots.draw_peptide_table(
                self.sub_sections["quantification"],
                self.mq_results["get_evidence_dicts"]["peptides_quant_table"]
            )

        # Protein Quantification Table
        if self.mq_results["get_evidence_dicts"]["protein_quant_table"]:
            maxquant_plots.draw_protein_table(
                self.sub_sections["quantification"],
                self.mq_results["get_evidence_dicts"]["protein_quant_table"]
            )

        # Intensity
        if self.mq_results["get_protegroups_dicts"]["pg_intensity_distri"]:
            maxquant_plots.draw_intensity_box(
                self.sub_sections["quantification"],
                self.mq_results["get_protegroups_dicts"]["pg_intensity_distri"]["box"],
                "intensity"
            )
        if self.mq_results["get_protegroups_dicts"]["pg_lfq_intensity_distri"]:
            maxquant_plots.draw_intensity_box(
                self.sub_sections["quantification"],
                self.mq_results["get_protegroups_dicts"]["pg_lfq_intensity_distri"]["box"],
                "lfq_intensity"
            )
        if self.mq_results["get_protegroups_dicts"]["raw_intensity_pca"]:
            maxquant_plots.draw_pg_pca(
                self.sub_sections["quantification"],
                self.mq_results["get_protegroups_dicts"]["raw_intensity_pca"],
                "raw_intensity"
            )
        if self.mq_results["get_protegroups_dicts"]["lfq_intensity_pca"]:
            maxquant_plots.draw_pg_pca(
                self.sub_sections["quantification"],
                self.mq_results["get_protegroups_dicts"]["lfq_intensity_pca"],
                "lfq_intensity"
            )
        if self.mq_results["ms_ms_identified"]:
            common_plots.draw_ms_ms_identified(
                self.sub_sections["identification"],
                self.mq_results["ms_ms_identified"]
            )
        if self.mq_results["get_evidence_dicts"]["peptide_intensity"]:
            maxquant_plots.draw_intensity_box(
                self.sub_sections["quantification"],
                self.mq_results["get_evidence_dicts"]["peptide_intensity"]["box"],
                "peptide_intensity"
            )

        # Contaminants
        if self.mq_results["get_protegroups_dicts"]["pg_contaminant"]:
            common_plots.draw_potential_contaminants(
                self.sub_sections["contaminants"],
                self.mq_results["get_protegroups_dicts"]["pg_contaminant"],
                True
            )
        if self.mq_results["get_evidence_dicts"]["top_contaminants"]:
            common_plots.draw_top_n_contaminants(
                self.sub_sections["contaminants"],
                self.mq_results["get_evidence_dicts"]["top_contaminants"]
            )

        if self.mq_results["get_evidence_dicts"]["charge_counts"]:
            common_plots.draw_charge_state(
                self.sub_sections["ms2"],
                self.mq_results["get_evidence_dicts"]["charge_counts"],
                "MaxQuant"
            )

        if self.mq_results["get_evidence_dicts"]["modified_percentage"]:
            common_plots.draw_modifications(
                self.sub_sections["identification"],
                self.mq_results["get_evidence_dicts"]["modified_percentage"]
            )

        if self.mq_results["get_evidence_dicts"]["peptide_id_count"]:
            maxquant_plots.draw_evidence_peptide_id_count(
                self.sub_sections["identification"],
                self.mq_results["get_evidence_dicts"]["peptide_id_count"]
            )

        if self.mq_results["get_evidence_dicts"]["protein_group_count"]:
            maxquant_plots.draw_evidence_protein_group_count(
                self.sub_sections["identification"],
                self.mq_results["get_evidence_dicts"]["protein_group_count"]
            )

        if self.mq_results["get_evidence_dicts"]["oversampling"]:
            common_plots.draw_oversampling(
                self.sub_sections["ms2"],
                self.mq_results["get_evidence_dicts"]["oversampling"],
                "",
                True
            )

        if self.mq_results["get_msms_dicts"]["missed_cleavages"]:
            common_plots.draw_msms_missed_cleavages(
                self.sub_sections["identification"],
                self.mq_results["get_msms_dicts"]["missed_cleavages"],
                True
            )

        if self.mq_results["get_evidence_dicts"]["rt_counts"]:
            common_plots.draw_ids_rt_count(
                self.sub_sections["time_mass"],
                self.mq_results["get_evidence_dicts"]["rt_counts"],
                "maxquant"
            )

        if self.mq_results["get_evidence_dicts"]["peak_rt"]:
            maxquant_plots.draw_evidence_peak_width_rt(
                self.sub_sections["time_mass"],
                self.mq_results["get_evidence_dicts"]["peak_rt"],
            )

        # Uncalibrated Mass Error
        if self.mq_results["get_evidence_dicts"]["uncalibrated_mass_error"]:
            maxquant_plots.draw_mass_error_box(
                self.sub_sections["time_mass"],
                self.mq_results["get_evidence_dicts"]["uncalibrated_mass_error"]
            )


        # Summary Table
        if (
            self.mq_results["get_msms_scans_dicts"]["summary_msms_spectra"]
            and
            self.mq_results["get_evidence_dicts"]["summary_stat"]["summary_identified_msms_count"]
        ):
            maxquant_plots.draw_maxquant_summary_table(
                self.sub_sections["summary"],
                self.mq_results["get_msms_scans_dicts"]["summary_msms_spectra"],
                self.mq_results["get_evidence_dicts"]["summary_stat"]["summary_identified_msms_count"],
                self.mq_results["get_evidence_dicts"]["summary_stat"]["summary_identified_peptides"],
                self.mq_results["get_protegroups_dicts"]["protein_summary"],
            )

        # Number of Peptides identified Per Protein
        if self.mq_results["get_protegroups_dicts"]["num_pep_per_protein_dict"]:
            maxquant_plots.draw_maxquant_num_pep_pro(
                self.sub_sections["identification"],
                self.mq_results["get_protegroups_dicts"]["num_pep_per_protein_dict"]
            )

        # Search Engine Scores
        if self.mq_results["get_msms_dicts"]["search_engine_scores"]:
            maxquant_plots.draw_maxquant_scores(
                self.sub_sections["search_engine"],
                self.mq_results["get_msms_dicts"]["search_engine_scores"]
            )

        # MaxQuant: Delta Mass [Da]
        if self.mq_results["get_evidence_dicts"]["maxquant_delta_mass_da"]:
            common_plots.draw_delta_mass_da_ppm(
                self.sub_sections["time_mass"],
                self.mq_results["get_evidence_dicts"]["maxquant_delta_mass_da"],
                "Mass Error [Da]"
            )

        # MaxQuant: Delta Mass [ppm]
        if self.mq_results["get_evidence_dicts"]["calibrated_mass_error"]:
            common_plots.draw_delta_mass_da_ppm(
                self.sub_sections["time_mass"],
                self.mq_results["get_evidence_dicts"]["calibrated_mass_error"],
                "Mass Error [ppm]"
            )

        # TopN
        if self.mq_results["get_msms_scans_dicts"]["top_n"]:
            maxquant_plots.draw_msms_scans_top_n(
                self.sub_sections["time_mass"],
                self.mq_results["get_msms_scans_dicts"]["top_n"]
            )
        if self.mq_results["get_msms_scans_dicts"]["top_over_rt"]:
            maxquant_plots.draw_msms_scans_top_over_rt(
                self.sub_sections["time_mass"],
                self.mq_results["get_msms_scans_dicts"]["top_over_rt"]
            )
        if self.mq_results["get_msms_scans_dicts"]["ion_injec_time_rt"]:
            maxquant_plots.draw_msms_scans_ion_injec_time_rt(
                self.sub_sections["time_mass"],
                self.mq_results["get_msms_scans_dicts"]["ion_injec_time_rt"]
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
