from unittest import TestCase

from pmultiqc.modules.quantms import QuantMSModule


class TestMultiqcModule(TestCase):
  path_mztab1 = 'resources/PXD005942-Sample-25-out.mzTab'
  def test_parse_out_mz_tab(self):
    print("Improve tests here")
