"""
Class-based parser for idXML search results, migrated from ms_io.parse_idxml.
"""

import os
import math
from collections import OrderedDict
from datetime import datetime

import numpy as np
from pyopenms import IdXMLFile

from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.ms.base import BaseParser
from pathlib import Path

class IdXMLReader(BaseParser):

    def __init__(
        self,
        file_paths: list[str | Path],
        mzml_table,
        xcorr_hist_range,
        hyper_hist_range,
        spec_evalue_hist_range,
        pep_hist_range,
        ml_spec_ident_final,
        mzml_peptide_map,
        remove_decoy: bool = True,
    ):
        super().__init__(file_paths)
        self.mzml_table = mzml_table
        self.xcorr_hist_range = xcorr_hist_range
        self.hyper_hist_range = hyper_hist_range
        self.spec_evalue_hist_range = spec_evalue_hist_range
        self.pep_hist_range = pep_hist_range
        self.ml_spec_ident_final = ml_spec_ident_final
        self.mzml_peptide_map = mzml_peptide_map
        self.remove_decoy = remove_decoy
        self.log = get_logger("pmultiqc.modules.common.ms.idxml")
        # Outputs populated by parse()
        self.search_engine = {
            "SpecE": OrderedDict(),
            "xcorr": OrderedDict(),
            "hyper": OrderedDict(),
            "PEPs": OrderedDict(),
            "consensus_support": OrderedDict(),
            "data_label": OrderedDict(),
        }
        self.msgf_label = False
        self.comet_label = False
        self.sage_label = False

    def parse(self, **_kwargs) -> None:
        idx_paths = list(self.file_paths or [])
        mzml_table = self.mzml_table
        ml_spec_ident_final = self.ml_spec_ident_final
        mzml_peptide_map = self.mzml_peptide_map
        remove_decoy = self.remove_decoy

        idx_paths, consensus_paths = self._separate_consensus_paths(idx_paths)
        _, _, _, labels = self._init_labels()

        for raw_id in idx_paths:
            self._parse_search_file(raw_id, remove_decoy, labels, mzml_table, ml_spec_ident_final, mzml_peptide_map)

        self._parse_consensus_files(consensus_paths, labels)

        self.search_engine["data_label"] = {
            "score_label": [labels["spec_e_label"], labels["xcorr_label"], labels["hyper_label"]],
            "peps_label": labels["peps_label"],
            "consensus_label": labels["consensus_label"],
        }
        self.msgf_label = labels["msgf_label"]
        self.comet_label = labels["comet_label"]
        self.sage_label = labels["sage_label"]
        return None

    def _separate_consensus_paths(self, idx_paths: list[str | Path]) -> tuple[list[str | Path], list[str | Path]]:
        consensus_paths: list[str | Path] = []
        for raw_id in list(idx_paths):
            if "consensus" in os.path.split(raw_id)[1]:
                consensus_paths.append(raw_id)
                idx_paths.remove(raw_id)
        return idx_paths, consensus_paths

    def _init_labels(self):
        labels = {
            "spec_e_label": [],
            "xcorr_label": [],
            "hyper_label": [],
            "peps_label": [],
            "consensus_label": [],
            "msgf_label": False,
            "comet_label": False,
            "sage_label": False,
        }
        return False, False, False, labels

    def _parse_search_file(self, raw_id, remove_decoy, labels, mzml_table, ml_spec_ident_final, mzml_peptide_map):
        self.log.info(
            "{}: Parsing search result file {}...".format(
                datetime.now().strftime("%H:%M:%S"), raw_id
            )
        )
        protein_ids: list = []
        peptide_ids: list = []
        # peptide_ids = PeptideIdentificationList()     # pyopenms >= 3.5.0
        IdXMLFile().load(raw_id, protein_ids, peptide_ids)
        raw_id_name = file_prefix(raw_id)

        identified_num = self._count_identified(peptide_ids, remove_decoy)
        ms_name = file_prefix(protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8"))
        search_engine_name = protein_ids[0].getSearchEngine()

        self._ensure_engine_slots(raw_id_name)
        hist = self._build_histograms()

        if search_engine_name == "MSGF+" or "msgf" in raw_id_name:
            mzml_table[ms_name]["MSGF"] = identified_num
            labels["msgf_label"] = True
            labels["spec_e_label"].append({"name": raw_id_name, "ylab": "Counts"})
            labels["peps_label"].append({"name": raw_id_name, "ylab": "Counts"})
            self._fill_msgf_hist(peptide_ids, hist)
            self.search_engine["SpecE"][raw_id_name] = hist["spectral_e"].dict["data"]
            self.search_engine["PEPs"][raw_id_name] = hist["posterior_error"].dict["data"]

        elif search_engine_name == "Comet" or "comet" in raw_id_name:
            labels["comet_label"] = True
            mzml_table[ms_name]["Comet"] = identified_num
            labels["xcorr_label"].append({"name": raw_id_name, "ylab": "Counts"})
            labels["peps_label"].append({"name": raw_id_name, "ylab": "Counts"})
            self._fill_comet_hist(peptide_ids, hist)
            self.search_engine["xcorr"][raw_id_name] = hist["cross_corr"].dict["data"]
            self.search_engine["PEPs"][raw_id_name] = hist["posterior_error"].dict["data"]

        elif search_engine_name == "Sage" or "sage" in raw_id_name:
            labels["sage_label"] = True
            mzml_table[ms_name]["Sage"] = identified_num
            labels["hyper_label"].append({"name": raw_id_name, "ylab": "Counts"})
            labels["peps_label"].append({"name": raw_id_name, "ylab": "Counts"})
            self._fill_sage_hist(peptide_ids, hist)
            self.search_engine["hyper"][raw_id_name] = hist["hyper"].dict["data"]
            self.search_engine["PEPs"][raw_id_name] = hist["posterior_error"].dict["data"]

        else:
            mzml_table[ms_name][search_engine_name] = identified_num

        mzml_table[ms_name]["num_quant_psms"] = ml_spec_ident_final.get(ms_name, 0)
        mzml_table[ms_name]["num_quant_peps"] = len(mzml_peptide_map.get(ms_name, []))

    def _parse_consensus_files(self, consensus_paths, labels):
        bar_stacks = ["target", "decoy", "target+decoy"]
        for raw_id in consensus_paths:
            self.log.info(
                "{}: Parsing consensus file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), format(raw_id)
                )
            )
            protein_ids = []
            peptide_ids = []
            IdXMLFile().load(raw_id, protein_ids, peptide_ids)
            raw_id_name = file_prefix(raw_id)
            labels["consensus_label"].append({"name": raw_id_name, "ylab": "Counts"})
            consensus_support = Histogram(
                "Consensus PSM number", plot_category="frequency", stacks=bar_stacks
            )
            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    support = hit.getMetaValue("consensus_support")
                    consensus_support.add_value(support, stack=hit.getMetaValue("target_decoy"))
            consensus_support.to_dict()
            for i in consensus_support.dict["data"].keys():
                self.search_engine["consensus_support"][f"{raw_id_name} ({i})"] = consensus_support.dict[
                    "data"
                ][i]

    def _ensure_engine_slots(self, raw_id_name: str) -> None:
        self.search_engine["SpecE"].setdefault(raw_id_name, OrderedDict())
        self.search_engine["xcorr"].setdefault(raw_id_name, OrderedDict())
        self.search_engine["hyper"].setdefault(raw_id_name, OrderedDict())
        self.search_engine["PEPs"].setdefault(raw_id_name, OrderedDict())

    def _count_identified(self, peptide_ids, remove_decoy: bool) -> int:
        if remove_decoy:
            return len(
                set(
                    [
                        i.getMetaValue("spectrum_reference")
                        for i in peptide_ids
                        if i.getHits()[0].getMetaValue("target_decoy") == "target"
                    ]
                )
            )
        return len(peptide_ids)

    def _build_histograms(self) -> dict:
        xcorr_breaks = list(
            np.arange(
                self.xcorr_hist_range["start"],
                self.xcorr_hist_range["end"] + self.xcorr_hist_range["step"],
                self.xcorr_hist_range["step"],
            ).round(2)
        )
        hyper_breaks = list(
            np.arange(
                self.hyper_hist_range["start"],
                self.hyper_hist_range["end"] + self.hyper_hist_range["step"],
                self.hyper_hist_range["step"],
            ).round(2)
        )
        spec_e_breaks = list(
            np.arange(
                self.spec_evalue_hist_range["start"],
                self.spec_evalue_hist_range["end"] + self.spec_evalue_hist_range["step"],
                self.spec_evalue_hist_range["step"],
            ).round(2)
        )
        spec_e_breaks.append(float("inf"))
        spec_e_breaks.sort()
        pep_breaks = list(
            np.concatenate(
                [
                    np.arange(
                        self.pep_hist_range["start"],
                        self.pep_hist_range["low_thresh"],
                        self.pep_hist_range["low_step"],
                    ),
                    np.arange(
                        self.pep_hist_range["low_thresh"],
                        self.pep_hist_range["high_thresh"],
                        self.pep_hist_range["high_step"],
                    ),
                    np.arange(
                        self.pep_hist_range["high_thresh"],
                        self.pep_hist_range["end"] + 0.01,
                        self.pep_hist_range["low_step"],
                    ),
                ]
            ).round(2)
        )
        bar_stacks = ["target", "decoy", "target+decoy"]
        return {
            "cross_corr": Histogram(
                "Comet cross-correlation score", plot_category="range", stacks=bar_stacks, breaks=xcorr_breaks
            ),
            "hyper": Histogram(
                "Sage hyperscore", plot_category="range", stacks=bar_stacks, breaks=hyper_breaks
            ),
            "spectral_e": Histogram(
                "MSGF spectral E-value", plot_category="range", stacks=bar_stacks, breaks=spec_e_breaks
            ),
            "posterior_error": Histogram(
                "Posterior error probability", plot_category="range", stacks=bar_stacks, breaks=pep_breaks
            ),
        }

    def _fill_msgf_hist(self, peptide_ids, hist):
        for peptide_id in peptide_ids:
            for hit in peptide_id.getHits():
                spec_e = (
                    hit.getMetaValue("SpecEvalue-score")
                    if hit.getMetaValue("SpecEvalue-score")
                    else hit.getMetaValue("MS:1002052")
                )
                log_spec_e = -math.log(spec_e, 10)
                pep = (
                    hit.getMetaValue("MS:1001493")
                    if hit.getMetaValue("MS:1001493")
                    else hit.getScore()
                )
                hist["spectral_e"].add_value(log_spec_e, stack=hit.getMetaValue("target_decoy"))
                hist["posterior_error"].add_value(pep, stack=hit.getMetaValue("target_decoy"))
        hist["spectral_e"].to_dict()
        hist["posterior_error"].to_dict()

    def _fill_comet_hist(self, peptide_ids, hist):
        for peptide_id in peptide_ids:
            for hit in peptide_id.getHits():
                xcorr = hit.getMetaValue("MS:1002252")
                pep = (
                    hit.getMetaValue("MS:1001493")
                    if hit.getMetaValue("MS:1001493")
                    else hit.getScore()
                )
                hist["cross_corr"].add_value(xcorr, stack=hit.getMetaValue("target_decoy"))
                hist["posterior_error"].add_value(pep, stack=hit.getMetaValue("target_decoy"))
        hist["cross_corr"].to_dict()
        hist["posterior_error"].to_dict()

    def _fill_sage_hist(self, peptide_ids, hist):
        for peptide_id in peptide_ids:
            for hit in peptide_id.getHits():
                hyper_score = hit.getMetaValue("hyperscore")
                pep = (
                    hit.getMetaValue("MS:1001493")
                    if hit.getMetaValue("MS:1001493")
                    else hit.getScore()
                )
                hist["hyper"].add_value(hyper_score, stack=hit.getMetaValue("target_decoy"))
                hist["posterior_error"].add_value(pep, stack=hit.getMetaValue("target_decoy"))
        hist["hyper"].to_dict()
        hist["posterior_error"].to_dict()