from __future__ import absolute_import

import os
import pandas as pd

from multiqc import config

from sdrf_pipelines.converters.openms.unimod import UnimodDatabase

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.common.common_utils import (
    parse_sdrf,
    cal_num_table_at_sample,
    summarize_modifications,
    evidence_rt_count
)
from pmultiqc.modules.common.plots.general import (
    draw_exp_design
)
from pmultiqc.modules.common.plots.id import (
    draw_identi_num,
    draw_identification,
    draw_peptide_length_distribution,
    draw_peptide_intensity,
    draw_ids_rt_count
)

from pmultiqc.modules.core.section_groups import add_group_modules
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.qpx.qpx_io import parse_qpx_parquet
from pmultiqc.modules.qpx.qpx_plot import (
    draw_summary_table,
    draw_whole_exp_charge,
    draw_qpx_ms2_charge
)
from pmultiqc.modules.qpx.qpx_utils import (
    calculate_run_stat,
    get_unimod_mod_qpx,
    get_pep_intensity,
    get_missed_cleavages
)
from pmultiqc.modules.common.ms_io import del_openms_convert_tsv


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.qpx.qpx")


class QpxModule(BasePMultiqcModule):

    def __init__(
            self,
            find_log_files_func,
            sub_sections,
            heatmap_colors
        ):
        """Initialize the QpxModule."""
        super().__init__(
            find_log_files_func,
            sub_sections,
            heatmap_colors
        )

        self.enable_exp = False
        self.enable_sdrf = False

        self.psm_df = None
        self.pg_df = None
        self.feature_df = None
        self.run_df = None
        self.psm_df_valid = False
        self.pg_df_valid = False
        self.feature_df_valid = False

        self.sample_df = pd.DataFrame()
        self.file_df = pd.DataFrame()
        self.exp_design_runs = None
        self.is_bruker = False
        self.is_multi_conditions = False

        self.cal_num_table_data = {}

    def get_data(self):

        log.info("[get_data] Starting data recognition and processing...")

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

        self.psm_df = self._read_qpx_parquets("pmultiqc/qpx_psm", "psm")
        self.pg_df = self._read_qpx_parquets("pmultiqc/qpx_pg", "pg")
        self.feature_df = self._read_qpx_parquets("pmultiqc/qpx_feature", "feature")
        self.run_df = self._read_qpx_parquets("pmultiqc/qpx_run", None)

        self.psm_df_valid = self._is_valid(self.psm_df)
        self.pg_df_valid = self._is_valid(self.pg_df)
        self.feature_df_valid = self._is_valid(self.feature_df)

        if self.psm_df_valid or self.pg_df_valid or self.feature_df_valid:
            log.info(
                f"[get_data] Data loading completed! "
                f"psm_df_valid: {self.psm_df_valid}, "
                f"pg_df_valid: {self.pg_df_valid}, "
                f"feature_df_valid: {self.feature_df_valid}."
            )
            return True
        else:
            log.warning("[get_data] No valid PSM or PG or feature data found. Files may be missing or empty.")
            return False


    def draw_plots(self):

        log.info("[draw_plots] Starting data processing and plot generation...")

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

        # Results Overview
        self._safe_draw(
            self.plot_results_overview,
            name="plot_results_overview"
        )

        # Identification Summary
        self._safe_draw(
            self.plot_id_summary,
            name="plot_id_summary"
        )

        # Quantification Analysis
        self._safe_draw(
            self.plot_quant_analysis,
            name="plot_quant_analysis"
        )

        # MS2 and Spectral Stats
        self._safe_draw(
            self.plot_ms2_stats,
            name="plot_ms2_stats"
        )

        # Mass Error Trends
        self._safe_draw(
            self.plot_mass_error,
            name="plot_mass_error"
        )

        # RT Quality Control
        self._safe_draw(
            self.plot_rt,
            name="plot_rt"
        )


        self.section_group_dict = {
            "experiment_sub_section": self.sub_sections["experiment"],
            "summary_sub_section": self.sub_sections["summary"],
            "identification_sub_section": self.sub_sections["identification"],
            # "search_engine_sub_section": self.sub_sections["search_engine"],
            # "contaminants_sub_section": self.sub_sections["contaminants"],
            "quantification_sub_section": self.sub_sections["quantification"],
            # "ms1_sub_section": self.sub_sections["ms1"],
            "ms2_sub_section": self.sub_sections["ms2"],
            "mass_error_sub_section": self.sub_sections["mass_error"],
            "rt_qc_sub_section": self.sub_sections["rt_qc"],
        }

        add_group_modules(self.section_group_dict, "")

        if self.enable_sdrf:
            del_openms_convert_tsv()

    # Results Overview
    def plot_results_overview(self):

        log.info("[Results Overview] Starting generation...")

        total_ms2_spectra_identified = 0
        total_peptide_count = 0
        total_protein_identified = 0
        total_protein_quantified = 0

        if self.psm_df_valid:

            # Summary Table & Pipeline Result Statistics
            # 'scan' is a list column (list<int32>) in psm.parquet, so .str[0] unwraps
            # the scan number; the column itself is unhashable and cannot be grouped on.
            total_ms2_spectra_identified = self.psm_df.groupby(
                ["run", self.psm_df["scan"].str[0]]
            ).ngroups
            total_peptide_count = self.psm_df["sequence"].nunique()

        if self.pg_df_valid:
            total_protein_identified = self.pg_df["anchor_protein"].nunique()

            valid_intensity = self.pg_df["intensity"].notna() & (self.pg_df["intensity"] != "")
            total_protein_quantified = self.pg_df.loc[valid_intensity, "anchor_protein"].nunique()

        # Pipeline Result Statistics
        if self.psm_df_valid and self.pg_df_valid:

            stat_at_run = dict()
            data_per_run = dict()

            pg_prots_by_run = self.pg_df.groupby("run")["anchor_protein"].apply(set).to_dict()
            for run, psm_group in self.psm_df.groupby("run"):
                run_str = str(run)
                prots = pg_prots_by_run.get(run_str, set())
                (
                    stat_at_run[run_str],
                    data_per_run[run_str]
                ) = calculate_run_stat(psm_group, prots)

            num_table_at_sample = cal_num_table_at_sample(self.file_df, data_per_run)

            self.cal_num_table_data = {
                "sdrf_samples": num_table_at_sample,
                "ms_runs": stat_at_run
            }

        if self.cal_num_table_data:
            draw_identi_num(
                sub_sections=self.sub_sections["summary"],
                enable_exp=self.enable_exp,
                enable_sdrf=self.enable_sdrf,
                is_multi_conditions=self.is_multi_conditions,
                sample_df=self.sample_df,
                file_df=self.file_df,
                cal_num_table_data=self.cal_num_table_data
            )

        draw_summary_table(
            self.sub_sections["summary"],
            total_ms2_spectra_identified,
            total_peptide_count,
            total_protein_identified,
            total_protein_quantified
        )


    # Identification Summary
    def plot_id_summary(self):

        log.info("[Identification Summary] Starting generation...")

        psm_modified = {}
        peptide_length = {}
        missed_cleavages = {}

        if self.psm_df_valid:

            psm = self.psm_df[["run", "sequence", "charge", "modifications"]].copy()

            psm["pep_length"] = psm["sequence"].str.len()

            unimod_data = UnimodDatabase()
            psm["modifications"] = psm["modifications"].apply(
                get_unimod_mod_qpx,
                args=(unimod_data,)
            )

            mod_plot_by_run = {}
            modified_cats = []

            for m, group in psm.groupby("run"):

                mod_plot_dict, modified_cat = summarize_modifications(
                    group[["sequence", "charge", "modifications"]].drop_duplicates()
                )
                mod_plot_by_run[m] = mod_plot_dict
                modified_cats.extend(modified_cat)

                peptide_length[m] = group["pep_length"].value_counts().sort_index().to_dict()

            mod_plot_by_sample = _sample_level_mods(
                df=psm,
                sdrf_file_df=self.file_df
            )

            # Modifications
            psm_modified["plot_data"] = [mod_plot_by_run, mod_plot_by_sample]
            psm_modified["cats"] = list(
                sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x))
            )

            # Missed Cleavages
            missed_cleavages = get_missed_cleavages(psm, self.run_df, self.file_df)

        draw_identification(
            sub_sections=self.sub_sections["identification"],
            cal_num_table_data=self.cal_num_table_data,
            missed_cleavages=missed_cleavages,
            modified=psm_modified,
            msms_identified_rate=None,
        )

        # Peptide Length Distribution
        if peptide_length:
           self._safe_draw(
                draw_peptide_length_distribution,
                name="draw_peptide_length_distribution",
                sub_section=self.sub_sections["identification"],
                plot_data=peptide_length
           )

    # Quantification Analysis
    def plot_quant_analysis(self):

        log.info("[Quantification Analysis] Starting generation...")

        qpx_pep_intensity = []

        if self.feature_df_valid:
            feature_tmp = self.feature_df[["run", "peptidoform", "intensities"]].copy()
            qpx_pep_intensity = get_pep_intensity(feature_tmp, self.file_df)

        if qpx_pep_intensity:
            self._safe_draw(
                draw_peptide_intensity,
                name="draw_peptide_intensity",
                sub_section=self.sub_sections["quantification"],
                plot_data=qpx_pep_intensity
            )

    # MS2 and Spectral Stats
    def plot_ms2_stats(self):
        log.info("[MS2 and Spectral Stats] Starting generation...")

        if self.feature_df_valid:
            feature_tmp = self.feature_df[["run", "charge"]].copy()

            if not feature_tmp.empty:

                self._safe_draw(
                    draw_whole_exp_charge,
                    name="draw_whole_exp_charge",
                    sub_section=self.sub_sections["ms2"],
                    df=feature_tmp
                )

                self._safe_draw(
                    draw_qpx_ms2_charge,
                    name="draw_qpx_ms2_charge",
                    sub_section=self.sub_sections["ms2"],
                    df=feature_tmp,
                    sdrf_file_df=self.file_df
                )

    # Mass Error Trends
    def plot_mass_error(self):
        log.info(
            "mass_error_ppm column is entirely empty in psm.parquet, may need separate calculation."
        )

    # RT Quality Control
    def plot_rt(self):

        log.info("[RT Quality Control] Starting generation...")

        if self.psm_df_valid:

            psm = self.psm_df[["run", "rt"]].copy()
            psm["retention time"] = psm["rt"] / 60
            psm.rename(columns={"run": "raw file"}, inplace=True)
            qpx_ids_over_rt = evidence_rt_count(psm)

            if qpx_ids_over_rt:
                self._safe_draw(
                    draw_ids_rt_count,
                    name="draw_ids_rt_count",
                    sub_section=self.sub_sections["rt_qc"],
                    rt_count_data=qpx_ids_over_rt,
                    report_type=""
                )

    def _read_qpx_parquets(self, search_pattern, qpx_type):
        """Read every parquet file matching search_pattern and concatenate them.

        A quantms.io project may split a table across several parquet files, so all
        matches are read rather than only the first one. Returns None if nothing matched.
        """
        dfs = []
        for f in self.find_log_files(search_pattern, filecontents=False):
            file_path = os.path.join(f["root"], f["fn"])
            if qpx_type is None:
                dfs.append(pd.read_parquet(file_path))
            else:
                dfs.append(parse_qpx_parquet(file_path=file_path, qpx_type=qpx_type))

        if not dfs:
            return None
        if len(dfs) == 1:
            return dfs[0]

        log.info(f"[get_data] Concatenating {len(dfs)} files for '{search_pattern}'.")
        return pd.concat(dfs, ignore_index=True)

    def _is_valid(self, df):
        return df is not None and not df.empty

    def _safe_draw(self, func, *args, name="plot", **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            log.exception(f"Failed to generate {name}")
            return None


def _sample_level_mods(df, sdrf_file_df):

    mod_plot = dict()

    if sdrf_file_df.empty:
        return mod_plot

    df_merged = df.merge(
        right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
        left_on="run",
        right_on="Run"
    )

    df_merged["Sample"] = df_merged["Sample"].astype(int)

    for sample, group in df_merged.groupby("Sample", sort=True):

        mod_plot_dict, _ = summarize_modifications(
            group[["sequence", "charge", "modifications"]].drop_duplicates()
        )
        mod_plot[f"Sample {str(sample)}"] = mod_plot_dict

    return mod_plot
