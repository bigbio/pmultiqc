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
        """Parse MzIdentML files and populate filtered_mzid_df."""
        mzid_dicts = self._extract_mzid_data()
        if not mzid_dicts:
            return

        mzid_df = pd.DataFrame(mzid_dicts)

        if not self._validate_required_columns(mzid_df):
            return

        self.filtered_mzid_df = self._process_dataframe(mzid_df)

    def _extract_mzid_data(self):
        """Extract data from MzIdentML files."""
        mzid_dicts = list()
        for file_path in self.file_paths:
            self._log_parsing_start(file_path)

            with mzid.MzIdentML(file_path) as mzid_data:
                if len(mzid_data) == 0:
                    raise ValueError(f"Please check your MzIdentML: {file_path}")

                self._log_parsing_complete(file_path)
                file_prefix_name = file_prefix(file_path)
                self._log_aggregation_start(file_prefix_name)

                for mzid_tmp in mzid_data:
                    mzid_dicts.extend(self._process_mzid_entry(mzid_tmp, file_prefix_name))

                self._log_aggregation_complete(file_prefix_name)

        return mzid_dicts

    def _log_parsing_start(self, file_path):
        """Log the start of parsing."""
        self.log.info(
            "{}: Parsing MzIdentML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), file_path
            )
        )

    def _log_parsing_complete(self, file_path):
        """Log the completion of parsing."""
        self.log.info(
            "{}: Done parsing MzIdentML file {}.".format(
                datetime.now().strftime("%H:%M:%S"), file_path
            )
        )

    def _log_aggregation_start(self, file_prefix_name):
        """Log the start of aggregation."""
        self.log.info(
            "{}: Aggregating MzIdentML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), file_prefix_name
            )
        )

    def _log_aggregation_complete(self, file_prefix_name):
        """Log the completion of aggregation."""
        self.log.info(
            "{}: Done aggregating MzIdentML file {}...".format(
                datetime.now().strftime("%H:%M:%S"), file_prefix_name
            )
        )

    def _process_mzid_entry(self, mzid_tmp, file_prefix_name):
        """Process a single MzIdentML entry."""
        mzid_tmp_part = {
            k: v for k, v in mzid_tmp.items() if k not in ["SpectrumIdentificationItem"]
        }

        mzid_dicts = []
        for spectrum_item in mzid_tmp.get("SpectrumIdentificationItem", []):
            if not self._should_process_spectrum_item(spectrum_item):
                continue

            spectrum_item_part = self._extract_spectrum_item_part(spectrum_item)

            for peptide_ref in spectrum_item.get("PeptideEvidenceRef", []):
                self._normalize_peptide_ref(peptide_ref)

                mzid_dict = {
                    **mzid_tmp_part,
                    **spectrum_item_part,
                    **peptide_ref,
                    "mzid_file_name": file_prefix_name,
                }

                mzid_dicts.append(self._filter_required_keys(mzid_dict))

        return mzid_dicts

    def _should_process_spectrum_item(self, spectrum_item):
        """Check if spectrum item should be processed."""
        spectrum_item_part = {
            k: v for k, v in spectrum_item.items()
            if k not in ["PeptideEvidenceRef", "PeptideSequence"]
        }

        rank = spectrum_item_part.get("rank")
        peptide_pass = spectrum_item_part.get("peptide passes threshold", "true") == "true"
        pass_threshold = spectrum_item.get("passThreshold", False)

        return rank == 1 and peptide_pass and pass_threshold

    def _extract_spectrum_item_part(self, spectrum_item):
        """Extract relevant parts from spectrum item."""
        return {
            k: v for k, v in spectrum_item.items()
            if k not in ["PeptideEvidenceRef", "PeptideSequence"]
        }

    def _normalize_peptide_ref(self, peptide_ref):
        """Normalize peptide reference keys."""
        if "name" in peptide_ref:
            peptide_ref["PeptideEvidenceRef_name"] = peptide_ref.pop("name")
        if "location" in peptide_ref:
            peptide_ref["PeptideEvidenceRef_location"] = peptide_ref.pop("location")
        if "FileFormat" in peptide_ref:
            peptide_ref["PeptideEvidenceRef_FileFormat"] = peptide_ref.pop("FileFormat")

    def _filter_required_keys(self, mzid_dict):
        """Filter dictionary to only include required keys."""
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
        return {k: v for k, v in mzid_dict.items() if k in need_keys}

    def _validate_required_columns(self, mzid_df):
        """Validate that required columns are present."""
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
            return False
        return True

    def _process_dataframe(self, mzid_df):
        """Process the dataframe with filtering and transformations."""
        # Filter out contaminants
        filtered_df = mzid_df[
            ~mzid_df["accession"].astype(str).str.lower().str.startswith("cont")
        ].copy()

        # Rename search engine columns
        self._rename_search_engine_columns(filtered_df)

        # Process column names and values
        self._process_column_renames(filtered_df)
        self._process_modifications(filtered_df)
        self._process_accession_groups(filtered_df)
        self._process_filename_column(filtered_df)

        return filtered_df

    def _rename_search_engine_columns(self, df):
        """Rename search engine score columns."""
        search_engines = [
            "SEQUEST:xcorr",
            "Mascot:score",
            "PEAKS:peptideScore",
            "xi:score",
            "Andromeda:score",
            "Proteome Discoverer Delta Score",
            "MS-GF:SpecEValue",
        ]
        df.rename(
            columns=lambda x: "search_engine_score" if x in search_engines else x, inplace=True
        )

        if "search_engine_score" not in df.columns:
            self.log.warning(
                "Please check the 'search_engine_score' field in the mzIdentML file."
            )

    def _process_column_renames(self, df):
        """Process column renames and location parsing."""
        if "retention time" in df.columns:
            df.rename(columns={"retention time": "retention_time"}, inplace=True)

        if "location" in df.columns:
            df["location"] = df["location"].apply(self.parse_location)

    def _process_modifications(self, df):
        """Process modification columns."""
        if "Modification" in df.columns:
            df["Modifications"] = df["Modification"].apply(self.process_modification)
        else:
            df["Modifications"] = None

    def _process_accession_groups(self, df):
        """Process accession groups."""
        df["accession_group"] = df.groupby(
            ["spectrumID", "PeptideSequence"]
        )["accession"].transform(lambda x: ";".join(x.unique()))

        if "isDecoy" not in df.columns:
            df["isDecoy"] = False

    def _process_filename_column(self, df):
        """Process filename column."""
        if "location" in df.columns:
            df["filename"] = df.apply(lambda x: file_prefix(x.location), axis=1)
        else:
            df["filename"] = df["mzid_file_name"]

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