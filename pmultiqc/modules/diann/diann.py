import logging
import os
import pandas as pd
import re
import numpy as np
import inspect

from pmultiqc.modules.common.histogram import Histogram
from sdrf_pipelines.openms.openms import UnimodDatabase
from collections import OrderedDict
from pmultiqc.modules.common.ms_io import del_openms_convert_tsv
from pmultiqc.modules.common.common_utils import (
    get_exp_sdrf,
    get_ms_path,
    get_msstats_path,
    parse_mzml
)
from pmultiqc.modules.common.common_plots import (
    HEATMAP_COLOR_LIST,
    draw_exp_design,
    draw_summary_protein_ident_table,
    draw_quantms_identi_num,
    draw_num_pep_per_protein,
    draw_peaks_per_ms2,
    draw_peak_intensity_distribution,
    draw_ms_information,
    draw_quantms_identification
)
from pmultiqc.modules.core.section_groups import add_group_modules
from pmultiqc.modules.diann.diann_utils import (
    get_diann_path,
    draw_dia_heatmap,
    draw_dia_intensitys,
    draw_dia_ms1,
    draw_dia_ms2s,
    draw_dia_mass_error,
    draw_dia_rt_qc,
    draw_diann_quant_table,
)
from pmultiqc.modules.common.calc_utils import mod_group_percentage
from pmultiqc.modules.common.file_utils import file_prefix


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

class DiaNNModule:

    def __init__(
            self,
            find_log_files_func,
            sub_sections
        ):

        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections

        self.heatmap_color_list = HEATMAP_COLOR_LIST

        self.ms_with_psm = list()
        self.cal_num_table_data = dict()
        self.quantms_modified = dict()

        self.sample_df = pd.DataFrame()
        self.file_df = pd.DataFrame()
        self.is_bruker = False
        self.is_multi_conditions = False

        (
            self.exp_design,
            self.enable_exp,
            self.enable_sdrf
        ) = get_exp_sdrf(self.find_log_files)

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

        (
            self.diann_report_path,
            self.enable_dia
        ) = get_diann_path(self.find_log_files)

        (
            self.msstats_input_path,
            self.msstats_input_valid
        ) = get_msstats_path(self.find_log_files)
        
        (
            mzml_table,
            mzml_peaks_ms2_plot,
            mzml_peak_distribution_plot,
            ms_info,
            _,  # total_ms2_spectra
            _,  # mzml_ms_df
            _,  # heatmap_charge_score
            _,  # mzml_charge_plot
            ms1_tic,
            ms1_bpc,
            ms1_peaks,
            ms1_general_stats
        ) = parse_mzml(
            is_bruker=self.is_bruker,
            read_ms_info=self.read_ms_info,
            ms_info_path=self.ms_info_path,
            ms_with_psm=self.ms_with_psm,
            enable_dia=self.enable_dia,
            ms_paths=self.ms_paths
        )

        draw_ms_information(
            self.sub_sections["ms1"],
            ms1_tic,
            ms1_bpc,
            ms1_peaks,
            ms1_general_stats
        )

        self.parse_diann_report()

        draw_summary_protein_ident_table(
            sub_sections=self.sub_sections["summary"],
            enable_dia=self.enable_dia,
            total_peptide_count=self.total_peptide_count,
            total_protein_quantified=self.total_protein_quantified
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

        if len(self.ms_info_path) > 0 and not self.is_bruker:

            draw_peaks_per_ms2(
                self.sub_sections["ms2"],
                mzml_peaks_ms2_plot,
                ms_info
            )

            draw_peak_intensity_distribution(
                self.sub_sections["ms2"],
                mzml_peak_distribution_plot,
                ms_info
            )
  
        draw_quantms_identification(
            self.sub_sections["identification"],
            cal_num_table_data=self.cal_num_table_data,
            mzml_table=mzml_table,
            quantms_modified=self.quantms_modified
        )

        caller_frame = inspect.stack()[1].frame
        caller_instance = caller_frame.f_locals.get('self')

        if caller_instance is not None:

            from pmultiqc.modules.quantms import QuantMSModule
            if not isinstance(caller_instance, QuantMSModule):

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

        if self.enable_sdrf:
            del_openms_convert_tsv()

    def get_data_from_diann(self):
        return (
            self.peptide_search_score,
            self.sample_df,
            self.sub_sections["experiment"],
            self.sub_sections["summary"],
            self.sub_sections["identification"],
            self.sub_sections["quantification"],
            self.sub_sections["ms1"],
            self.sub_sections["ms2"],
            self.sub_sections["mass_error"],
            self.sub_sections["rt_qc"]
        )

    def parse_diann_report(self):

        log.info("Parsing {}...".format(self.diann_report_path))

        # parse DIA-NN report data
        if os.path.splitext(self.diann_report_path)[1] == ".tsv":
            report_data = pd.read_csv(
                self.diann_report_path,
                header=0,
                sep="\t",
                on_bad_lines="warn"
            )
        else:
            report_data = pd.read_parquet(self.diann_report_path)

        # "Decoy" appears only in DIA-NN 2.0 and later.
        # 0 or 1 based on whether the precursor is decoy, relevant when using --report-decoys
        if "Decoy" in report_data.columns:
            report_data = report_data[report_data["Decoy"] == 0].copy()

        # Normalisation.Factor: can be calculated as Precursor.Normalised/Precursor.Quantity
        required_cols = ["Precursor.Normalised", "Precursor.Quantity"]
        if "Normalisation.Factor" not in report_data.columns and all(
            col in report_data.columns for col in required_cols
        ):
            report_data["Normalisation.Factor"] = report_data[required_cols[0]] / report_data[required_cols[1]]

        # Draw: Standard Deviation of Intensity
        if "Precursor.Quantity" in report_data.columns:
            draw_dia_intensitys(self.sub_sections["quantification"], report_data)
            draw_dia_heatmap(self.sub_sections["summary"], report_data, self.heatmap_color_list)
        
        log.info("Draw the DIA MS1 subsection.")
        draw_dia_ms1(self.sub_sections["ms1"], report_data)

        log.info("Draw the DIA MS2 subsection.")
        draw_dia_ms2s(self.sub_sections["ms2"], report_data)

        log.info("Draw the DIA mass_error subsection.")
        draw_dia_mass_error(self.sub_sections["mass_error"], report_data)

        log.info("Draw the DIA rt_qc subsection.")
        draw_dia_rt_qc(self.sub_sections["rt_qc"], report_data)

        # Draw: Quantification Table (DIA-NN, without msstats data)
        if not self.msstats_input_valid:
            log.info("Draw the DIA quant table subsection.")
            draw_diann_quant_table(
                self.sub_sections["quantification"],
                report_data,
                self.sample_df,
                self.file_df
            )

        pattern = re.compile(r"\(.*?\)")
        report_data["sequence"] = [
            pattern.sub("", s) for s in report_data["Modified.Sequence"]
        ]

        self.total_protein_quantified = len(set(report_data["Protein.Group"]))
        self.total_peptide_count = len(set(report_data["sequence"]))

        log.info("Processing DIA pep_plot.")
        protein_pep_map = report_data.groupby("Protein.Group")["sequence"].agg(list).to_dict()
        self.pep_plot = Histogram("number of peptides per proteins", plot_category="frequency")
        for _, peps in protein_pep_map.items():
            number = len(set(peps))
            self.pep_plot.add_value(number)

        log.info("Processing DIA peptide_search_score.")
        self.peptide_search_score = dict()
        pattern = re.compile(r"\((.*?)\)")
        unimod_data = UnimodDatabase()
        for peptide, group in report_data.groupby("Modified.Sequence"):
            original_mods = re.findall(pattern, peptide)
            for mod in set(original_mods):
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
        log.info("Processing DIA Modifications.")
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

        log.info("Processing DIA mod_plot_dict.")
        mod_plot_dict = dict()
        modified_cats = list()
        for run_file, group in report_data.groupby("Run"):
            run_file = str(run_file)
            self.ms_with_psm.append(run_file)

            # Modifications
            mod_group_processed = mod_group_percentage(group.drop_duplicates())
            mod_plot_dict[run_file] = dict(
                zip(mod_group_processed["modifications"], mod_group_processed["percentage"])
            )
            modified_cats.extend(mod_group_processed["modifications"])

            self.cal_num_table_data[run_file] = {"protein_num": len(set(group["Protein.Group"]))}
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
