from __future__ import annotations
from pathlib import Path
from datetime import datetime
from pyteomics import mztab
import pandas as pd
import os
import re

from multiqc import config
from pmultiqc.modules.common.ms.base import BaseParser
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.file_utils import file_prefix


class MzTabReader(BaseParser):
    def __init__(
        self,
        file_path: Path | str
    ) -> None:
        super().__init__([file_path])

        self.file_path = file_path

        self.pep_table: pd.DataFrame = pd.DataFrame()
        self.delta_mass: dict = dict()
        self.psm: pd.DataFrame = pd.DataFrame()
        self.prot: pd.DataFrame = pd.DataFrame()
        self.prot_abundance_cols: list = list()
        self.meta_data: dict = {}
        self.mztab_data: mztab.MzTab = None
        self.ms_with_psm: list = list()
        self.total_protein_identified: int = 0
        self.total_protein_quantified: int = 0

        self.log = get_logger("pmultiqc.modules.common.ms.mztab")

    def parse(self, **_kwargs) -> None:

        self.log.info(
            "{}: Parsing mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.file_path
            )
        )
        mztab_data = mztab.MzTab(self.file_path)
        self.log.info(
            "{}: Done parsing mzTab file {}.".format(
                datetime.now().strftime("%H:%M:%S"), self.file_path
            )
        )
        self.log.info(
            "{}: Aggregating mzTab file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.file_path
            )
        )
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)

        self.delta_mass["target"] = dict()
        self.delta_mass["decoy"] = dict()

        # PSM table data
        psm = mztab_data.spectrum_match_table
        if len(psm) == 0:
            raise ValueError("The PSM section of mzTab is missing, please check your mzTab!")

        # Generate "opt_global_cv_MS: 1002217_DECOY_peptide" column if this column is not contained in the PSM subtable
        if "opt_global_cv_MS:1002217_decoy_peptide" not in psm.columns.values:
            psm["opt_global_cv_MS:1002217_decoy_peptide"] = psm.apply(
                lambda x: 1 if self.dis_decoy(x["accession"]) == "DECOY" else 0, axis=1
            )
        # map to spectrum file name in experimental design file
        psm["stand_spectra_ref"] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(":")[0] + "-location"])
            + ":"
            + x.spectra_ref.split(":")[1],
            axis=1,
        )
        psm["filename"] = psm.apply(
            lambda x: file_prefix(meta_data[x.spectra_ref.split(":")[0] + "-location"]),
            axis=1,
        )
        self.ms_with_psm = psm["filename"].unique().tolist()

        prot = mztab_data.protein_table
        self.prot_search_score = dict()

        prot_abundance_cols = list(
            filter(
                lambda x: re.match(r"protein_abundance_assay.*?", x) is not None,
                prot.columns.tolist(),
            )
        )
        opt_cols = list(filter(lambda x: x.startswith("opt_"), prot.columns.tolist()))
        score_cols = list(
            filter(lambda x: x.startswith("best_search_engine_score"), prot.columns.tolist())
        )
        # TODO in theory we do not need accession since the index is the accession
        fixed_cols = [
            "accession",
            "description",
            "taxid",
            "species",
            "database",
            "database_version",
            "search_engine",
            "ambiguity_members",
            "modifications",
            "protein_coverage",
        ]

        keep_cols = [
            c
            for c in (fixed_cols + score_cols + prot_abundance_cols + opt_cols)
            if c in prot.columns
        ]
        prot = prot[keep_cols].copy()

        # We only need the overall protein (group) scores and abundances. Currently we do not care about details of single proteins (length, description,...)
        if "opt_global_result_type" in prot.columns:
            prot = prot[prot["opt_global_result_type"] != "protein_details"].copy()

        if config.kwargs["remove_decoy"]:
            psm = psm[psm["opt_global_cv_MS:1002217_decoy_peptide"] != 1].copy()
            # TODO do we really want to remove groups that contain a single decoy? I would say ALL members need to be decoy.
            prot = prot[~prot["accession"].str.contains(config.kwargs["decoy_affix"], regex=False)]

        prot.dropna(subset=["ambiguity_members"], inplace=True)

        prot["protein_group"] = prot["ambiguity_members"].str.replace(",", ";")

        self.total_protein_identified = len(prot.index)

        if prot_abundance_cols:
            prot.dropna(how="all", subset=prot_abundance_cols, inplace=True)

        self.total_protein_quantified = len(prot.index)

        psm["pep_length"] = psm["sequence"].str.len()

        self.mztab_data = mztab_data
        self.pep_table = pep_table
        self.psm = psm
        self.prot = prot
        self.prot_abundance_cols = prot_abundance_cols
        self.meta_data = meta_data

        return None

    @staticmethod
    def dis_decoy(protein_name):

        decoy_affix = config.kwargs["decoy_affix"]
        affix_type = config.kwargs["affix_type"]

        proteins = protein_name.split(";")

        if decoy_affix not in protein_name:
            return "TARGET"

        is_decoy = (
            lambda x: x.startswith(decoy_affix) if affix_type == "prefix" else x.endswith(decoy_affix)
        )

        if any(is_decoy(p) for p in proteins):
            return "DECOY"

        return "TARGET"
