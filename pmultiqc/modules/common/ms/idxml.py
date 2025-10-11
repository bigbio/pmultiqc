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
        file_paths: list[str, Path],
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

    def parse(self, **kwargs) -> None:
        # Reuse the previously implemented logic (kept below)
        idx_paths = list(self.file_paths or [])
        mzml_table = self.mzml_table
        xcorr_hist_range = self.xcorr_hist_range
        hyper_hist_range = self.hyper_hist_range
        spec_evalue_hist_range = self.spec_evalue_hist_range
        pep_hist_range = self.pep_hist_range
        ml_spec_ident_final = self.ml_spec_ident_final
        mzml_peptide_map = self.mzml_peptide_map
        remove_decoy = self.remove_decoy

        consensus_paths = []
        for raw_id in idx_paths:
            if "consensus" in os.path.split(raw_id)[1]:
                consensus_paths.append(raw_id)

        for raw_id in consensus_paths:
            if raw_id in idx_paths:
                idx_paths.remove(raw_id)

        msgf_label, comet_label, sage_label = False, False, False
        search_engine = self.search_engine
        spec_e_label, xcorr_label, hyper_label, peps_label, consensus_label = [], [], [], [], []

        for raw_id in idx_paths:
            self.log.info(
                "{}: Parsing search result file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), raw_id
                )
            )

            protein_ids = []
            peptide_ids = []
            IdXMLFile().load(raw_id, protein_ids, peptide_ids)
            raw_id_name = file_prefix(raw_id)

            if remove_decoy:
                identified_num = len(
                    set(
                        [
                            i.getMetaValue("spectrum_reference")
                            for i in peptide_ids
                            if i.getHits()[0].getMetaValue("target_decoy") == "target"
                        ]
                    )
                )
            else:
                identified_num = len(peptide_ids)

            ms_name = file_prefix(protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8"))
            search_engine_name = protein_ids[0].getSearchEngine()

            search_engine["SpecE"][raw_id_name] = OrderedDict()
            search_engine["xcorr"][raw_id_name] = OrderedDict()
            search_engine["hyper"][raw_id_name] = OrderedDict()
            search_engine["PEPs"][raw_id_name] = OrderedDict()

            xcorr_breaks = list(
                np.arange(
                    xcorr_hist_range["start"],
                    xcorr_hist_range["end"] + xcorr_hist_range["step"],
                    xcorr_hist_range["step"],
                ).round(2)
            )

            hyper_breaks = list(
                np.arange(
                    hyper_hist_range["start"],
                    hyper_hist_range["end"] + hyper_hist_range["step"],
                    hyper_hist_range["step"],
                ).round(2)
            )

            spec_e_breaks = list(
                np.arange(
                    spec_evalue_hist_range["start"],
                    spec_evalue_hist_range["end"] + spec_evalue_hist_range["step"],
                    spec_evalue_hist_range["step"],
                ).round(2)
            )
            spec_e_breaks.append(float("inf"))
            spec_e_breaks.sort()

            pep_breaks = list(
                np.concatenate(
                    [
                        np.arange(
                            pep_hist_range["start"],
                            pep_hist_range["low_thresh"],
                            pep_hist_range["low_step"],
                        ),
                        np.arange(
                            pep_hist_range["low_thresh"],
                            pep_hist_range["high_thresh"],
                            pep_hist_range["high_step"],
                        ),
                        np.arange(
                            pep_hist_range["high_thresh"],
                            pep_hist_range["end"] + 0.01,
                            pep_hist_range["low_step"],
                        ),
                    ]
                ).round(2)
            )

            bar_stacks = ["target", "decoy", "target+decoy"]
            cross_corr = Histogram(
                "Comet cross-correlation score",
                plot_category="range",
                stacks=bar_stacks,
                breaks=xcorr_breaks,
            )
            hyper = Histogram(
                "Sage hyperscore", plot_category="range", stacks=bar_stacks, breaks=hyper_breaks
            )
            spectral_e = Histogram(
                "MSGF spectral E-value",
                plot_category="range",
                stacks=bar_stacks,
                breaks=spec_e_breaks,
            )
            posterior_error = Histogram(
                "Posterior error probability",
                plot_category="range",
                stacks=bar_stacks,
                breaks=pep_breaks,
            )

            if search_engine_name == "MSGF+" or "msgf" in raw_id_name:
                mzml_table[ms_name]["MSGF"] = identified_num
                msgf_label = True
                spec_e_label.append({"name": raw_id_name, "ylab": "Counts"})
                peps_label.append({"name": raw_id_name, "ylab": "Counts"})
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
                        spectral_e.add_value(log_spec_e, stack=hit.getMetaValue("target_decoy"))
                        posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                spectral_e.to_dict()
                posterior_error.to_dict()
                search_engine["SpecE"][raw_id_name] = spectral_e.dict["data"]
                search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

            elif search_engine_name == "Comet" or "comet" in raw_id_name:
                comet_label = True
                mzml_table[ms_name]["Comet"] = identified_num
                xcorr_label.append({"name": raw_id_name, "ylab": "Counts"})
                peps_label.append({"name": raw_id_name, "ylab": "Counts"})
                for peptide_id in peptide_ids:
                    for hit in peptide_id.getHits():
                        xcorr = hit.getMetaValue("MS:1002252")
                        pep = (
                            hit.getMetaValue("MS:1001493")
                            if hit.getMetaValue("MS:1001493")
                            else hit.getScore()
                        )
                        cross_corr.add_value(xcorr, stack=hit.getMetaValue("target_decoy"))
                        posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                cross_corr.to_dict()
                posterior_error.to_dict()
                search_engine["xcorr"][raw_id_name] = cross_corr.dict["data"]
                search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

            elif search_engine_name == "Sage" or "sage" in raw_id_name:
                sage_label = True
                mzml_table[ms_name]["Sage"] = identified_num
                hyper_label.append({"name": raw_id_name, "ylab": "Counts"})
                peps_label.append({"name": raw_id_name, "ylab": "Counts"})
                for peptide_id in peptide_ids:
                    for hit in peptide_id.getHits():
                        hyper_score = hit.getMetaValue("hyperscore")
                        pep = (
                            hit.getMetaValue("MS:1001493")
                            if hit.getMetaValue("MS:1001493")
                            else hit.getScore()
                        )
                        hyper.add_value(hyper_score, stack=hit.getMetaValue("target_decoy"))
                        posterior_error.add_value(pep, stack=hit.getMetaValue("target_decoy"))

                hyper.to_dict()
                posterior_error.to_dict()
                search_engine["hyper"][raw_id_name] = hyper.dict["data"]
                search_engine["PEPs"][raw_id_name] = posterior_error.dict["data"]

            else:
                mzml_table[ms_name][search_engine_name] = identified_num

            mzml_table[ms_name]["num_quant_psms"] = (
                ml_spec_ident_final[ms_name] if ms_name in ml_spec_ident_final.keys() else 0
            )
            mzml_table[ms_name]["num_quant_peps"] = (
                len(mzml_peptide_map[ms_name]) if ms_name in ml_spec_ident_final.keys() else 0
            )

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

            consensus_label.append({"name": raw_id_name, "ylab": "Counts"})

            consensus_support = Histogram(
                "Consensus PSM number", plot_category="frequency", stacks=bar_stacks
            )

            for peptide_id in peptide_ids:
                for hit in peptide_id.getHits():
                    support = hit.getMetaValue("consensus_support")
                    consensus_support.add_value(support, stack=hit.getMetaValue("target_decoy"))
            consensus_support.to_dict()

            for i in consensus_support.dict["data"].keys():
                search_engine["consensus_support"][f"{raw_id_name} ({i})"] = consensus_support.dict[
                    "data"
                ][i]

        search_engine["data_label"] = {
            "score_label": [spec_e_label, xcorr_label, hyper_label],
            "peps_label": peps_label,
            "consensus_label": consensus_label,
        }
        # Persist flags
        self.msgf_label = msgf_label
        self.comet_label = comet_label
        self.sage_label = sage_label
        return None