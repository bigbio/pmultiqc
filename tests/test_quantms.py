import logging
import os

from pmultiqc.modules.quantms import QuantMSModule

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "resources")


class TestQuantMS:

    def test_read_ms_info(self):
        # parquet_file = os.path.join(TEST_DATA_DIR, "quantms/a05191_ms_info.parquet")
        # quantms_module = QuantMSModule()
        # quantms_module.read_ms_info()
        logging.info("Test finished")
