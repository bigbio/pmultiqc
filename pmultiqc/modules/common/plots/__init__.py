"""
PMultiQC Plotting Package

This package contains organized plotting functions for PMultiQC reports.
Functions are separated by logical categories for better maintainability.
"""

from .contaminants import (
    draw_potential_contaminants,
    draw_top_n_contaminants
)

from .identification import (
    draw_ms_ms_identified,
    draw_quantms_identification,
    draw_num_pep_per_protein
)

from .charge_modifications import (
    draw_charge_state,
    draw_modifications,
    draw_precursor_charge_distribution
)

from .mass_analysis import (
    draw_delta_mass_da_ppm,
    draw_msms_missed_cleavages
)

from .retention_time import (
    draw_ids_rt_count,
    draw_oversampling
)

from .tables import (
    draw_peptides_table,
    draw_protein_table,
    draw_exp_design,
    draw_summary_protein_ident_table,
    draw_quantms_identi_num
)

from .ms_analysis import (
    draw_ms_information,
    draw_peaks_per_ms2,
    draw_peak_intensity_distribution
)

from .general import (
    draw_heatmap,
    remove_subtitle,
    HEATMAP_COLOR_LIST
)

__all__ = [
    # Contaminants
    "draw_potential_contaminants",
    "draw_top_n_contaminants",

    # Identification
    "draw_ms_ms_identified",
    "draw_quantms_identification",
    "draw_num_pep_per_protein",

    # Charge & Modifications
    "draw_charge_state",
    "draw_modifications",
    "draw_precursor_charge_distribution",

    # Mass Analysis
    "draw_delta_mass_da_ppm",
    "draw_msms_missed_cleavages",

    # Retention Time
    "draw_ids_rt_count",
    "draw_oversampling",

    # Tables
    "draw_peptides_table",
    "draw_protein_table",
    "draw_exp_design",
    "draw_summary_protein_ident_table",
    "draw_quantms_identi_num",

    # MS Analysis
    "draw_ms_information",
    "draw_peaks_per_ms2",
    "draw_peak_intensity_distribution",

    # General
    "draw_heatmap",
    "remove_subtitle",
    "HEATMAP_COLOR_LIST"
]
