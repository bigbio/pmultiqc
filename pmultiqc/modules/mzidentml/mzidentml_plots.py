

# mzIdentML: Quantification Table
from pmultiqc.modules.mzidentml import mzidentml_utils
from pmultiqc.modules.quantms.quantms_plots import draw_peptides_table, draw_protein_table


def draw_mzid_quant_table(sub_section, mzid_mzml_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = mzidentml_utils.create_peptides_table(
        mzid_mzml_df
    )
    draw_peptides_table(
        sub_section,
        peptides_table,
        peptides_headers,
        "mzIdentML"
    )

    # Protein Quantification Table
    protein_table, protein_headers = mzidentml_utils.create_protein_table(
        mzid_mzml_df
    )
    draw_protein_table(
        sub_section,
        protein_table,
        protein_headers,
        "mzIdentML"
    )