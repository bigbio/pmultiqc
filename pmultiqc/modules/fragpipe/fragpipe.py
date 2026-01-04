import pandas as pd

from pmultiqc.modules.base import BasePMultiqcModule
from pmultiqc.modules.fragpipe import fragpipe_io
from pmultiqc.modules.common.stats import cal_delta_mass_dict
from pmultiqc.modules.common.plots.id import (
    draw_summary_protein_ident_table,
    draw_delta_mass_da_ppm,
    draw_identi_num,
    draw_identification,
    draw_ids_rt_count,
    draw_num_pep_per_protein
)
from pmultiqc.modules.core.section_groups import (
    add_group_modules,
    add_sub_section
)
from pmultiqc.modules.common.common_utils import (
    group_charge,
    evidence_rt_count
)
from pmultiqc.modules.common.plots.general import plot_html_check
from collections import OrderedDict
from pmultiqc.modules.common.histogram import Histogram

from multiqc.plots import bargraph

from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.fragpipe.fragpipe")


class FragPipeModule(BasePMultiqcModule):
    """pmultiqc module for FragPipe results."""

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):

        super().__init__(find_log_files_func, sub_sections, heatmap_colors)

        self.delta_masses = []
        self.charge_states = []
        self.pipeline_stats = []
        self.retentions = []


    def get_data(self):

        log.info("Starting data recognition and processing...")

        self.fragpipe_files = fragpipe_io.get_fragpipe_files(self.find_log_files)

        if self.fragpipe_files["psm"]:
            (
                self.delta_masses,
                self.charge_states,
                self.pipeline_stats,
                self.retentions
            ) = self.parse_psm(
                fragpipe_files=self.fragpipe_files
            )
        else:
            log.warning("Required input not found: psm.tsv")
            return False

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
        
        if self.retentions:
            self.draw_ids_over_rt(
                sub_section=self.sub_sections["rt_qc"],
                retentions=self.retentions
            )

        section_group_dict = {
            "summary_sub_section": self.sub_sections["summary"],
            "identification_sub_section": self.sub_sections["identification"],
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

        for psm in fragpipe_files.get("psm", []):

            psm_df = fragpipe_io.psm_reader(psm)

            if psm_df is None or psm_df.empty:
                log.warning(f"Skipping unreadable/empty FragPipe PSM file: {psm}")
                continue

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

            # Retention: MS2 scan's precursor retention time (in seconds)
            if "Retention" in psm_df.columns:
                retentions.append(psm_df[["Run", "Retention"]].copy())

        return (
            delta_masses,
            charge_states,
            pipeline_stats,
            retentions
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
