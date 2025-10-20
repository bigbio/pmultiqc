from __future__ import annotations

import os
from pathlib import Path
import pandas as pd
from datetime import datetime
from pyteomics import mzid

from pmultiqc.modules.common.ms.base import BaseParser
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.file_utils import file_prefix


class MzidReader(BaseParser):
    def __init__(
            self,
            file_paths: list[str | Path],
    ) -> None:

        super().__init__(file_paths)

        # Outputs populated by parse()
        self.filtered_mzid_df: pd.DataFrame = pd.DataFrame()

        self.log = get_logger("pmultiqc.modules.common.ms.mzid")

    def parse(self, **_kwargs) -> None:

        mzid_dicts = list()
        for file_path in self.file_paths:

            self.log.info(
                "{}: Parsing MzIdentML file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file_path
                )
            )

            with mzid.MzIdentML(file_path) as mzid_data:

                if len(mzid_data) == 0:
                    raise ValueError(f"Please check your MzIdentML: {file_path}")

                self.log.info(
                    "{}: Done parsing MzIdentML file {}.".format(
                        datetime.now().strftime("%H:%M:%S"), file_path
                    )
                )
                m = file_prefix(file_path)
                self.log.info(
                    "{}: Aggregating MzIdentML file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
                )

                for mzid_tmp in mzid_data:
                    mzid_tmp_part = {
                        k: v for k, v in mzid_tmp.items() if k not in ["SpectrumIdentificationItem"]
                    }

                    for spectrum_item in mzid_tmp.get("SpectrumIdentificationItem", []):
                        spectrum_item_part = {
                            k: v
                            for k, v in spectrum_item.items()
                            if k not in ["PeptideEvidenceRef", "PeptideSequence"]
                        }

                        rank = spectrum_item_part.get("rank")
                        peptide_pass = spectrum_item_part.get("peptide passes threshold", "true") == "true"
                        pass_threshold = spectrum_item.get("passThreshold", False)

                        if rank != 1 or not peptide_pass or not pass_threshold:
                            continue

                        for peptide_ref in spectrum_item.get("PeptideEvidenceRef", []):

                            if "name" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_name"] = peptide_ref.pop("name")
                            if "location" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_location"] = peptide_ref.pop("location")
                            if "FileFormat" in peptide_ref:
                                peptide_ref["PeptideEvidenceRef_FileFormat"] = peptide_ref.pop(
                                    "FileFormat"
                                )

                            mzid_dict = {
                                **mzid_tmp_part,
                                **spectrum_item_part,
                                **peptide_ref,
                                "mzid_file_name": m,
                            }

                            need_keys = [
                                "SEQUEST:xcorr",
                                "Mascot:score",
                                "PEAKS:peptideScore",
                                "xi:score",
                                "retention time",
                                "location",
                                "Modification",
                                "spectrumID",
                                "isDecoy",
                                "accession",
                                "PeptideSequence",
                                "experimentalMassToCharge",
                                "calculatedMassToCharge",
                                "chargeState",
                                "mzid_file_name",
                                "Andromeda:score",
                                "Proteome Discoverer Delta Score",
                                "MS-GF:SpecEValue",
                            ]
                            mzid_dicts.append({k: v for k, v in mzid_dict.items() if k in need_keys})
                self.log.info(
                    "{}: Done aggregating MzIdentML file {}...".format(
                        datetime.now().strftime("%H:%M:%S"), m
                    )
                )

        mzid_df = pd.DataFrame(mzid_dicts)

        # Check columns
        check_list = [
            "spectrumID",
            "PeptideSequence",
            "chargeState",
            "accession",
            "Modification",
            "experimentalMassToCharge",
            "calculatedMassToCharge",
        ]
        missing_cols = [col for col in check_list if col not in mzid_df.columns]
        if missing_cols:
            self.log.warning(f"MzIdentML file is missing required fields: {missing_cols}")
            self.filtered_mzid_df = pd.DataFrame()
            return None

        # Filter out contaminants: start-with "cont" (case-insensitive)
        filtered_mzid_df = mzid_df[
            ~mzid_df["accession"].astype(str).str.lower().str.startswith("cont")
        ].copy()

        search_engines = [
            "SEQUEST:xcorr",
            "Mascot:score",
            "PEAKS:peptideScore",
            "xi:score",
            "Andromeda:score",
            "Proteome Discoverer Delta Score",
            "MS-GF:SpecEValue",
        ]
        filtered_mzid_df.rename(
            columns=lambda x: "search_engine_score" if x in search_engines else x, inplace=True
        )

        if "search_engine_score" not in filtered_mzid_df.columns:
            self.log.warning(
                "Please check the 'search_engine_score' field in the mzIdentML file."
            )

        if "retention time" in filtered_mzid_df.columns:
            filtered_mzid_df.rename(columns={"retention time": "retention_time"}, inplace=True)

        if "location" in filtered_mzid_df.columns:
            filtered_mzid_df["location"] = filtered_mzid_df["location"].apply(self.parse_location)

        if "Modification" in filtered_mzid_df.columns:
            filtered_mzid_df["Modifications"] = filtered_mzid_df["Modification"].apply(
                self.process_modification
            )
        else:
            filtered_mzid_df["Modifications"] = None

        filtered_mzid_df["accession_group"] = filtered_mzid_df.groupby(
            ["spectrumID", "PeptideSequence"]
        )["accession"].transform(lambda x: ";".join(x.unique()))

        if "isDecoy" not in filtered_mzid_df.columns:
            filtered_mzid_df["isDecoy"] = False

        # location: path of mzML file
        if "location" in filtered_mzid_df.columns:
            filtered_mzid_df["filename"] = filtered_mzid_df.apply(
                lambda x: file_prefix(x.location), axis=1
            )
        else:
            filtered_mzid_df["filename"] = filtered_mzid_df["mzid_file_name"]

        self.filtered_mzid_df = filtered_mzid_df

        return None

    @staticmethod
    def parse_location(location):
        if "\\" in location:
            location = location.replace("\\", "/")
        return os.path.basename(location)

    @staticmethod
    def process_modification(modification):
        if not isinstance(modification, list):
            modifications = None
        else:
            modifi_list = list()
            for i in modification:
                if i.get("name", None) is not None:
                    modifi_list.append(str(i.get("location")) + "-" + i.get("name", None))
                elif i.get("cross-link receiver", None) is not None:
                    modifi_list.append(str(i.get("location")) + "-CrossLinkReceiver")
            modifications = ";".join(modifi_list)
        return modifications