import logging
import os

import pytest

from pmultiqc.modules.common.ms.idxml import IdXMLReader

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "resources")
IDXML_DIR = os.path.join(TEST_DATA_DIR, "lfq", "raw_ids")

HIST_RANGES = {
    "xcorr_hist_range": {"start": 0, "end": 5, "step": 0.1},
    "hyper_hist_range": {"start": 0, "end": 5, "step": 0.1},
    "spec_evalue_hist_range": {"start": 0, "end": 20, "step": 0.4},
    "pep_hist_range": {
        "start": 0, "end": 1,
        "low_thresh": 0.3, "high_thresh": 0.7,
        "low_step": 0.03, "high_step": 0.08,
    },
}


class TestQuantMS:

    def test_read_ms_info(self):
        # parquet_file = os.path.join(TEST_DATA_DIR, "quantms/a05191_ms_info.parquet")
        # quantms_module = QuantMSModule()
        # quantms_module.read_ms_info()
        logging.info("Test finished")

    @pytest.mark.parametrize("idxml_filename,expected_engine_key", [
        (
            "SF_200217_pPeptideLibrary_pool1_HCDOT_rep1_comet_feat_perc.idXML",
            "Comet",
        ),
        (
            "SF_200217_pPeptideLibrary_pool1_HCDOT_rep1_msgf_feat_perc.idXML",
            "MSGF",
        ),
    ])
    def test_idxml_reader_empty_mzml_table(self, idxml_filename, expected_engine_key):
        """Regression test: IdXMLReader must not raise KeyError when mzml_table
        starts empty (e.g. when quantms runs with mzml_statistics = false)."""
        idxml_path = os.path.join(IDXML_DIR, idxml_filename)
        mzml_table = {}

        reader = IdXMLReader(
            file_paths=[idxml_path],
            mzml_table=mzml_table,
            ml_spec_ident_final={},
            mzml_peptide_map={},
            **HIST_RANGES,
        )
        reader.parse()

        assert len(mzml_table) > 0, "mzml_table should be populated from idXML spectra_data"
        for run_name, run_data in mzml_table.items():
            assert expected_engine_key in run_data, (
                f"Expected '{expected_engine_key}' key in mzml_table['{run_name}']"
            )
            assert "num_quant_psms" in run_data
            assert "num_quant_peps" in run_data