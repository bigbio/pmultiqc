import pandas as pd
import numpy as np
import re

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.fragpipe import fragpipe_io
from pmultiqc.modules.fragpipe.fragpipe_io import (
    ion_reader,
    get_ion_intensity_data,
)
from pmultiqc.modules.common.stats import (
    cal_delta_mass_dict,
    nanmedian,
    cal_hm_charge,
    qual_uniform
)
from pmultiqc.modules.common.plots.id import (
    draw_summary_protein_ident_table,
    draw_delta_mass_da_ppm,
    draw_identi_num,
    draw_identification,
    draw_ids_rt_count,
    draw_num_pep_per_protein,
    draw_peptide_intensity,
    draw_msms_missed_cleavages,
    rebuild_dict_structure,
    draw_top_n_contaminants,
    draw_potential_contaminants,
    draw_modifications
)
from pmultiqc.modules.core.section_groups import (
    add_group_modules,
    add_sub_section
)
from pmultiqc.modules.common.common_utils import (
    group_charge,
    evidence_rt_count,
    top_n_contaminant_percent,
    cal_contaminant_percent,
    mods_statistics
)
from pmultiqc.modules.common.plots.general import (
    plot_html_check,
    plot_data_check,
    stat_pep_intensity,
    search_engine_score_bins,
    draw_search_engine_scores,
    draw_heatmap
)

from collections import OrderedDict
from pmultiqc.modules.common.histogram import Histogram

from multiqc import config
from multiqc.plots import bargraph, box

from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.fragpipe.fragpipe")


NOT_CONT_TAG = "NOT_CONT"

class FragPipeModule(BasePMultiqcModule):
    """pmultiqc module for FragPipe results."""

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        super().__init__(find_log_files_func, sub_sections, heatmap_colors)

        self.delta_masses = []
        self.charge_states = []
        self.pipeline_stats = []
        self.retentions = []
        self.intensities = []
        self.missed_cleavages = []
        self.hyperscores = []
        self.contam_df = []
        self.mods = []
        self.hm_data = []

        # Ion-level intensity data from ion.tsv
        self.ion_intensity_data = None
        self.ion_sample_cols = []


    def get_data(self):

        log.info("Starting data recognition and processing...")

        self.fragpipe_files = fragpipe_io.get_fragpipe_files(self.find_log_files)

        if self.fragpipe_files is None:
            log.warning("No FragPipe files found.")
            return False

        if self.fragpipe_files.get("psm"):
            (
                self.delta_masses,
                self.charge_states,
                self.pipeline_stats,
                self.retentions,
                self.intensities,
                self.missed_cleavages,
                self.hyperscores,
                self.contam_df,
                self.mods,
                self.hm_data
            ) = self.parse_psm(
                fragpipe_files=self.fragpipe_files
            )
        else:
            log.warning("Required input not found: psm.tsv")
            return False

        # Parse ion.tsv for ion-level intensity data
        if self.fragpipe_files.get("ion"):
            self.ion_intensity_data, self.ion_sample_cols = self.parse_ion(
                fragpipe_files=self.fragpipe_files
            )

        return True

    def draw_plots(self):

        log.info("Starting to process plotting data...")

        # Delta Mass
        if self.delta_masses:
            self.draw_delta_mass(
                sub_sections=self.sub_sections["mass_error"],
                delta_masses=self.delta_masses
            )

        # Charge-state
        if self.charge_states:
            self.draw_charge_state(
                sub_section=self.sub_sections["ms2"],
                charge_states=self.charge_states
            )

        if self.pipeline_stats:
            
            # Statistics
            (
                summary_result,
                statistics_result,
                peptides_per_protein
            ) = _calculate_statistics(self.pipeline_stats)

            # Summary Table
            draw_summary_protein_ident_table(
                sub_sections=self.sub_sections["summary"],
                use_two_columns=True,
                total_peptide_count=summary_result["total_peptides"],
                total_protein_quantified=summary_result["total_proteins"]
            )

            # Pipeline Result Statistics
            draw_identi_num(
                sub_sections=self.sub_sections["summary"],
                cal_num_table_data=statistics_result
            )

            # Number of Peptides identified Per Protein
            draw_num_pep_per_protein(
                sub_sections=self.sub_sections["identification"],
                pep_plot=peptides_per_protein,
                is_fragpipe_or_mzid=True
            )

            # ProteinGroups Count & Peptide ID Count
            draw_identification(
                self.sub_sections["identification"],
                cal_num_table_data=statistics_result,
            )

        # Peptide Intensity Distribution
        if self.intensities:
            self.draw_intensity(
                sub_section=self.sub_sections["quantification"],
                intensities=self.intensities
            )

        mc_plot_data = {}

        # Missed Cleavages
        if self.missed_cleavages:
            mc_plot_data = self.draw_missed_cleavages(
                sub_section=self.sub_sections["identification"],
                missed_cleavages=self.missed_cleavages
            )

        # Summary of Hyperscore
        if self.hyperscores:
            self.draw_hyperscore(
                sub_section=self.sub_sections["search_engine"],
                hyperscores=self.hyperscores
            )

        # Contaminants
        if self.contam_df:
            self.draw_contaminants(
                sub_section=self.sub_sections["contaminants"],
                contam_df=self.contam_df
            )

        # Modifications
        if self.mods:
            self.draw_mods(
                sub_section=self.sub_sections["identification"],
                mods=self.mods
            )

        # self.hm_data
        if self.hm_data:
            self.draw_fragpipe_heatmap(
                sub_section=self.sub_sections["summary"],
                hm=self.hm_data,
                hm_color=self.heatmap_color_list,
                missed_cleavages=mc_plot_data
            )

        # IDs over RT
        if self.retentions:
            self.draw_ids_over_rt(
                sub_section=self.sub_sections["rt_qc"],
                retentions=self.retentions
            )

        # Ion-level intensity plots from ion.tsv
        if self.ion_intensity_data:
            self.draw_ion_intensity_distribution(
                sub_section=self.sub_sections["quantification"],
                intensity_data=self.ion_intensity_data
            )

        section_group_dict = {
            "summary_sub_section": self.sub_sections["summary"],
            "identification_sub_section": self.sub_sections["identification"],
            "search_engine_sub_section": self.sub_sections["search_engine"],
            "contaminants_sub_section": self.sub_sections["contaminants"],
            "quantification_sub_section": self.sub_sections["quantification"],
            "ms2_sub_section": self.sub_sections["ms2"],
            "mass_error_sub_section": self.sub_sections["mass_error"],
            "rt_qc_sub_section": self.sub_sections["rt_qc"],
        }

        add_group_modules(section_group_dict, "")

        log.info("Plotting data processing completed.")


    @staticmethod
    def parse_psm(fragpipe_files):

        delta_masses = []
        charge_states = []
        pipeline_stats = []
        retentions = []
        intensities = []
        missed_cleavages = []
        hyperscores = []
        contam_df = []
        mods = []
        hm_data = []

        for psm in fragpipe_files.get("psm", []):

            psm_df = fragpipe_io.psm_reader(psm)

            if psm_df is None or psm_df.empty:
                log.warning(f"Skipping unreadable/empty FragPipe PSM file: {psm}")
                continue

            psm_df, psm_cont_df = _mark_contaminants(
                df=psm_df,
                contam_affix=config.kwargs["contaminant_affix"]
            )
            log.info(f"Number of non-contaminant rows in {psm}: {len(psm_df)}")

            if "Delta Mass" in psm_df.columns:
                delta_masses.append(psm_df["Delta Mass"].copy())
            
            if "Charge" in psm_df.columns:
                charge_states.append(psm_df[["Run", "Charge"]].copy())

            # Pipeline Result Statistics
            stats_requires = [
                "Run", "Is Unique", "Modified Peptide", "Protein", "Peptide"
            ]
            if all(
                col in psm_df.columns
                for col in stats_requires
            ):
                pipeline_stats.append(psm_df[stats_requires].copy())

            # Peptide Intensity Distribution
            if _has_valid_column(psm_df, "Intensity"):
                log.info(f"Intensity in {psm} is available.")

                intensities.append(psm_df[["Run", "Intensity"]].copy())

            else:
                log.info(f"All intensity values in {psm} are 0 and not available")

            # Missed Cleavages
            if "Number of Missed Cleavages" in psm_df.columns:
                missed_cleavages.append(
                    psm_df[["Run", "Number of Missed Cleavages"]].copy()
                )

            # Summary of Hyperscore
            if "Hyperscore" in psm_df.columns:
                hyperscores.append(psm_df[["Run", "Hyperscore"]].copy())

            # Retention: MS2 scan's precursor retention time (in seconds)
            if "Retention" in psm_df.columns:
                retentions.append(psm_df[["Run", "Retention"]].copy())

            # Modifications
            if "Assigned Modifications" in psm_df.columns:
                mods.append(psm_df[["Run", "Assigned Modifications"]].copy())

            # HeatMap
            hm_requires = [
                "Run", "Modified Peptide", "Protein", "Peptide",
                "Intensity", "Retention", "Charge"
            ]
            if all(
                col in psm_cont_df.columns
                for col in hm_requires
            ):
                hm_data.append(psm_cont_df[hm_requires].copy())

            # Contaminants
            if _has_valid_contaminant(psm_cont_df):
                log.info(f"{psm} contains contaminants.")

                if _has_valid_column(psm_cont_df, "Intensity"):
                    log.info("Contaminant analysis...")

                    contam_df.append(
                        psm_cont_df[["Run", "Intensity", "cont_protein", "Protein"]].copy()
                    )
                else:
                    log.warning(
                        f"All intensity values in {psm} are unavailable; skipping contaminant analysis."
                    )
            else:
                log.info(f"No contaminants found in {psm}; skipping contaminant analysis.")

        return (
            delta_masses,
            charge_states,
            pipeline_stats,
            retentions,
            intensities,
            missed_cleavages,
            hyperscores,
            contam_df,
            mods,
            hm_data
        )

    # Delta Mass
    @staticmethod
    def draw_delta_mass(sub_sections, delta_masses: list):

        # Delta Mass: difference between calibrated observed peptide mass and calculated peptide mass (in Da)
        df = (
            pd.concat(delta_masses, ignore_index=True)
            .to_frame(name="Delta Mass")
        )

        df["Delta Mass"] = pd.to_numeric(df["Delta Mass"], errors="coerce")
        df = df.dropna(subset=["Delta Mass"])

        log.info(f"Number of delta mass data points: {len(df)}")

        delta_mass_da = cal_delta_mass_dict(df, "Delta Mass")

        draw_delta_mass_da_ppm(
            sub_sections, delta_mass_da, "Mass Error [Da]"
        )

        log.info("Delta mass [Da] plot generated.")


    # Charge-state
    @staticmethod
    def draw_charge_state(sub_section, charge_states: list):

        # Charge: charge state of the identified peptide

        df = pd.concat(charge_states, ignore_index=True)
        df["Charge"] = df["Charge"].astype("str")
        log.info(f"Number of charge state rows in DataFrame: {len(df)}")

        stat_data_by_run = group_charge(df, "Run", "Charge")

        bar_data = stat_data_by_run.to_dict(orient="index")

        draw_config = {
            "id": "charge_state",
            "cpswitch": True,
            "cpswitch_c_active": False,
            "title": "Charge-state",
            "tt_decimals": 0,
            "ylab": "Count",
            "save_data_file": False,
        }

        bar_html = bargraph.plot(
            data=bar_data,
            pconfig=draw_config,
        )

        bar_html = plot_html_check(bar_html)

        add_sub_section(
            sub_section=sub_section,
            plot=bar_html,
            order=6,
            description="Distribution of identified peptide charge states.",
            helptext="""
                [FragPipe: psm.tsv] Charge: charge state of the identified peptide.
                """,
        )

        log.info("Charge-state plot generated.")


    # Peptide Intensity Distribution
    @staticmethod
    def draw_intensity(sub_section, intensities: list):

        df = pd.concat(intensities, ignore_index=True)
        log.info(f"Number of intensity rows in DataFrame: {len(df)}")

        intensity_by_run = {}
        for run, group in df.groupby("Run"):
            intensity_by_run[run] = stat_pep_intensity(group["Intensity"])

        plot_data = [intensity_by_run]

        draw_peptide_intensity(
            sub_section=sub_section,
            plot_data=plot_data
        )


    # Missed Cleavages
    @staticmethod
    def draw_missed_cleavages(sub_section, missed_cleavages: list):

        df = pd.concat(missed_cleavages, ignore_index=True)
        log.info(f"Number of missed cleavages rows in DataFrame: {len(df)}")

        mc_by_run = {}
        for run, group in df.groupby("Run", sort=True):
            mc = group["Number of Missed Cleavages"].value_counts()
            mc_by_run[run] = mc.to_dict()

        re_mc_by_run = rebuild_dict_structure(mc_by_run)

        mc_plot = {
            "plot_data": re_mc_by_run,
            "cats": ["0", "1", ">=2"]
        }

        draw_msms_missed_cleavages(
            sub_section=sub_section,
            missed_cleavages_data=mc_plot,
            is_maxquant=False
        )

        return re_mc_by_run


    # Summary of Hyperscore
    @staticmethod
    def draw_hyperscore(sub_section, hyperscores: list):

        df = pd.concat(hyperscores, ignore_index=True)
        log.info(f"Number of hyperscore rows in DataFrame: {len(df)}")
        log.info(f"Maximum Hyperscore value: {df['Hyperscore'].max()}")

        plot_data = search_engine_score_bins(
            bins_start=0,
            bins_end=120,
            bins_step=3,
            df=df,
            groupby_col="Run",
            score_col="Hyperscore"
        )

        draw_search_engine_scores(
            sub_section=sub_section,
            plot_data=plot_data,
            plot_type="fragpipe"
        )


    # Contaminants
    @staticmethod
    def draw_contaminants(sub_section, contam_df: list):

        df = pd.concat(contam_df, ignore_index=True)

        num_cont = (df["cont_protein"] != NOT_CONT_TAG).sum()
        log.info(f"Number of contaminants rows in DataFrame: {num_cont}")

        contam_percent = cal_contaminant_percent(
            df=df[["Run", "Intensity", "Protein"]].copy(),
            protein_col="Protein",
            intensity_col="Intensity",
            run_col="Run",
            contam_affix=config.kwargs["contaminant_affix"]
        )

        draw_potential_contaminants(
            sub_section=sub_section,
            contaminant_percent=contam_percent,
            report_type="fragpipe"
        )

        top_contams = top_n_contaminant_percent(
            df=df,
            not_cont_tag=NOT_CONT_TAG,
            cont_tag_col="cont_protein",
            intensity_col="Intensity",
            run_col="Run",
            top_n=5,
        )

        draw_top_n_contaminants(
            sub_section=sub_section,
            top_contaminants_data=top_contams
        )


    # Modifications
    @staticmethod
    def draw_mods(sub_section, mods: list):

        df = pd.concat(mods, ignore_index=True)

        num_mods = df["Assigned Modifications"].notna().sum()
        log.info(f"Number of modifications rows in DataFrame: {num_mods}")

        if num_mods == 0:
            log.info("No modification data found; skipping modifications analysis.")
            return

        df["modifications"] = df["Assigned Modifications"].apply(_extract_modifications)

        modified_data = mods_statistics(
            df=df,
            run_col="Run"
        )
        
        draw_modifications(
            sub_section=sub_section,
            modified_data=modified_data
        )


    # HeatMap
    @staticmethod
    def draw_fragpipe_heatmap(sub_section, hm: list, hm_color: list, missed_cleavages: dict):

        if not hm:
            log.warning("No heatmap data; skipping heatmap.")
            return

        df = pd.concat(hm, ignore_index=True)

        if df.empty:
            log.warning("Heatmap DataFrame is empty; skipping heatmap.")
            return

        log.info(f"Number of HeatMap rows in DataFrame: {len(df)}")

        heatmap_cols = [
            "Contaminants",
            "Peptide Intensity",
            "Charge",
            "Missed Cleavages",
            "Missed Cleavages Var",
            "ID rate over RT",
            "Pep Missing Values",
        ]

        if not _has_valid_column(df, "Intensity"):
            heatmap_cols = [
                x
                for x in heatmap_cols
                if x not in ["Contaminants", "Peptide Intensity"]
            ]

        # Charge
        if "Charge" in heatmap_cols:
            hm_charge_all = cal_hm_charge(
                df=df,
                run_col="Run",
                charge_col="Charge"
            )

        if missed_cleavages:

            # Missed Cleavages
            mc = {
                key: value.get("0", 0) / 100
                for key, value in missed_cleavages.items()
            }

            # Missed Cleavages Var
            mc_median = np.median(list(mc.values()))
            mc_var = dict(
                zip(
                    mc.keys(),
                    list(map(lambda v: 1 - np.abs(v - mc_median), mc.values())),
                )
            )
        else:
            heatmap_cols = [
                x
                for x in heatmap_cols
                if x not in ["Missed Cleavages", "Missed Cleavages Var"]
            ]

        # 8. Pep Missing Values
        global_peps = df["Modified Peptide"].unique()
        global_peps_count = len(global_peps)

        contam_affix = config.kwargs["contaminant_affix"]

        heatmap_plot = {}
        for run, group in df.groupby("Run"):

            heatmap_plot[run] = {}

            # 1. Contaminants
            if "Contaminants" in heatmap_cols:

                is_contam = group["Protein"].str.contains(contam_affix, na=False, case=False)
                cont_intensity = group[is_contam]["Intensity"].sum()
                all_intensity = group["Intensity"].sum()

                if np.isnan(cont_intensity) or np.isnan(all_intensity) or all_intensity == 0:
                    hm_contam = 1
                else:
                    hm_contam = 1 - cont_intensity / all_intensity
                
                heatmap_plot[run]["Contaminants"] = hm_contam

            # 2. Peptide Intensity
            if "Peptide Intensity" in heatmap_cols:
                median_int = nanmedian(group["Intensity"], 0)  # if everything is NaN, use 0
                hm_intensity = np.minimum(
                    1.0, median_int / (2 ** 23)
                )  # score = 1, if intensity >= 2**23

                heatmap_plot[run]["Peptide Intensity"] = hm_intensity

            # 3. Charge
            if "Charge" in heatmap_cols:
                heatmap_plot[run]["Charge"] = hm_charge_all.get(run, 0)

            # 4. Missed Cleavages
            if "Missed Cleavages" in heatmap_cols:
                heatmap_plot[run]["Missed Cleavages"] = mc.get(run, 0)
            
            # 5. Missed Cleavages Var
            if "Missed Cleavages Var" in heatmap_cols:
                heatmap_plot[run]["Missed Cleavages Var"] = mc_var.get(run, 0)

            # 6. ID rate over RT
            if "ID rate over RT" in heatmap_cols:
                hm_rt = qual_uniform(group["Retention"])
                heatmap_plot[run]["ID rate over RT"] = hm_rt

            # 7. Pep Missing Values
            if "Pep Missing Values" in heatmap_cols:

                if global_peps_count > 0:
                    hm_pep_missing_values = np.minimum(
                        1.0,
                        len(set(global_peps) & set(group["Modified Peptide"].unique())) / global_peps_count,
                    )
                else:
                    hm_pep_missing_values = 0

                heatmap_plot[run]["Pep Missing Values"] = hm_pep_missing_values

        draw_heatmap(
            sub_sections=sub_section,
            hm_colors=hm_color,
            heatmap_data=heatmap_plot,
            heatmap_xnames="",
            heatmap_ynames="",
            report_type="fragpipe"
        )


    # IDs over RT
    @staticmethod
    def draw_ids_over_rt(sub_section, retentions: list):

        df = pd.concat(retentions, ignore_index=True)
        log.info(f"Number of retention rows in DataFrame: {len(df)}")

        df.rename(
            columns={"Run": "raw file", "Retention": "retention time"},
            inplace=True
        )

        df["retention time"] = pd.to_numeric(df["retention time"], errors="coerce") / 60.0
        df = df.dropna(subset=["retention time"])

        plot_data = evidence_rt_count(df)

        draw_ids_rt_count(
            sub_section=sub_section,
            rt_count_data=plot_data,
            report_type="fragpipe"
        )

    @staticmethod
    def parse_ion(fragpipe_files):
        """
        Parse ion.tsv files for ion-level intensity data.

        Parameters
        ----------
        fragpipe_files : dict
            Dictionary of FragPipe file paths.

        Returns
        -------
        tuple
            (ion_intensity_data, sample_cols) where:
            - ion_intensity_data: Dictionary with intensity distribution and CV data
            - sample_cols: List of sample intensity column names
        """
        ion_files = fragpipe_files.get("ion", [])

        if not ion_files:
            log.info("No ion.tsv files found.")
            return None, []

        all_intensity_data = {
            'intensity_distribution': {},
            'intensity_cv': {},
            'intensity_summary': {}
        }
        all_sample_cols = []

        for ion_file in ion_files:
            try:
                ion_df, sample_cols = ion_reader(ion_file)

                if ion_df is None or ion_df.empty:
                    log.warning(f"Skipping unreadable/empty ion.tsv file: {ion_file}")
                    continue

                log.info(f"Loaded ion.tsv with {len(ion_df)} rows and {len(sample_cols)} samples")

                intensity_data = get_ion_intensity_data(ion_df, sample_cols)

                if intensity_data:
                    # Merge intensity distributions
                    for sample, values in intensity_data.get('intensity_distribution', {}).items():
                        if sample in all_intensity_data['intensity_distribution']:
                            all_intensity_data['intensity_distribution'][sample].extend(values)
                        else:
                            all_intensity_data['intensity_distribution'][sample] = values

                    # Merge CV data
                    for sample, values in intensity_data.get('intensity_cv', {}).items():
                        if sample in all_intensity_data['intensity_cv']:
                            all_intensity_data['intensity_cv'][sample].extend(values)
                        else:
                            all_intensity_data['intensity_cv'][sample] = values

                    # Merge summary data (take first occurrence)
                    for sample, summary in intensity_data.get('intensity_summary', {}).items():
                        if sample not in all_intensity_data['intensity_summary']:
                            all_intensity_data['intensity_summary'][sample] = summary

                    all_sample_cols.extend([c for c in sample_cols if c not in all_sample_cols])

            except Exception as e:
                log.warning(f"Error parsing ion.tsv file {ion_file}: {e}")
                continue

        if not all_intensity_data['intensity_distribution']:
            log.info("No valid intensity data found in ion.tsv files.")
            return None, []

        log.info(f"Ion intensity data parsed for {len(all_sample_cols)} samples")

        return all_intensity_data, all_sample_cols

    @staticmethod
    def draw_ion_intensity_distribution(sub_section, intensity_data: dict):
        """
        Draw ion-level intensity distribution box plot.

        Parameters
        ----------
        sub_section : dict
            Section to add the plot to.
        intensity_data : dict
            Dictionary containing 'intensity_distribution' data.
        """
        distribution = intensity_data.get('intensity_distribution', {})

        if not distribution:
            log.info("No ion intensity distribution data available.")
            return

        log.info(f"Drawing ion intensity distribution for {len(distribution)} samples")

        draw_config = {
            "id": "ion_intensity_distribution_box",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Ion Intensity Distribution",
            "tt_decimals": 2,
            "xlab": "log2(Intensity)",
            "sort_samples": False,
            "save_data_file": False,
        }

        box_html = box.plot(list_of_data_by_sample=distribution, pconfig=draw_config)

        box_html = plot_data_check(
            plot_data=distribution,
            plot_html=box_html,
            log_text="pmultiqc.modules.fragpipe.fragpipe",
            function_name="draw_ion_intensity_distribution"
        )
        box_html = plot_html_check(box_html)

        add_sub_section(
            sub_section=sub_section,
            plot=box_html,
            order=6,
            description="Ion-level intensity distribution per sample from ion.tsv.",
            helptext="""
                [FragPipe: ion.tsv] This plot shows the log2-transformed ion intensity
                distribution for each sample/channel. The ion.tsv file contains precursor-level
                quantification data from IonQuant.

                For TMT experiments, each box represents a TMT channel.
                For label-free experiments, each box represents a sample/run.

                A higher median intensity and narrower distribution typically indicate
                better quantification quality. Large differences between samples may
                indicate normalization issues or batch effects.
            """,
        )

        log.info("Ion intensity distribution plot generated.")


def _calculate_statistics(pipeline_stats: list):

    df = pd.concat(pipeline_stats, ignore_index=True)
    log.info(f"Number of pipeline result statistics rows in DataFrame: {len(df)}")

    summary_data = {
        "total_proteins": len(set(df["Protein"])),
        "total_peptides": len(set(df["Peptide"]))
    }

    stats_by_run = dict()
    for run, group in df.groupby("Run"):

        unique_group = group.loc[group["Is Unique"]]

        modified_peptides = group.loc[
            group["Modified Peptide"].notna() & (group["Modified Peptide"] != ""),
            "Modified Peptide"
        ]

        stats_by_run[run] = {
            "protein_num": len(set(group["Protein"])),
            "peptide_num": len(set(group["Peptide"])),
            "unique_peptide_num": len(set(unique_group["Peptide"])),
            "modified_peptide_num": modified_peptides.nunique()
        }
    
    statistics_data = {
        "ms_runs": stats_by_run
    }

    protein_pep_map = df.groupby("Protein")["Peptide"].agg(list).to_dict()
    pep_plot = Histogram("number of peptides per proteins", plot_category="frequency")
    for _, peps in protein_pep_map.items():
        number = len(set(peps))
        pep_plot.add_value(number)

    categorys = OrderedDict()
    categorys["Frequency"] = {
        "name": "Frequency",
        "description": "number of peptides per proteins",
    }

    pep_plot.to_dict(percentage=True, cats=categorys)

    return summary_data, statistics_data, pep_plot


def _has_valid_column(df: pd.DataFrame, col: str):

    if col not in df.columns:
        return False

    values = df[col].dropna()

    return not (values == 0).all()


def _mark_contaminants(df, contam_affix="CONT"):

    is_contaminant = df["Protein"].str.contains(contam_affix, na=False, case=False)
    
    psm_cont_df = df.copy()
    psm_cont_df["cont_protein"] = np.where(is_contaminant, psm_cont_df["Protein"], NOT_CONT_TAG)
    
    psm_df = df[~is_contaminant].copy()
    
    return psm_df, psm_cont_df


def _has_valid_contaminant(df: pd.DataFrame):
    return not (df["cont_protein"] == NOT_CONT_TAG).all()


def _extract_modifications(x):

    if pd.isna(x):
        return "Unmodified"
    
    pattern = re.compile(r"N-term|\d+([A-Z])\(")

    mod_types = []
    for m in pattern.finditer(x):
        if m.group(0).startswith("N-term"):
            site = "N-term"
        else:
            site = m.group(1)

        if site not in mod_types:
            mod_types.append(site)

    if not mod_types:
        return "Unmodified"

    return ",".join(mod_types)
