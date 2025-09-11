import logging
import pytest
import numpy as np
from pmultiqc.modules.common.statistics_utils import nanmedian as pm_nanmedian


class TestStatisticsUtils:
    """
    Test class for statistics module functions.

    """

    def test_nanmedian(self):

        ## no NaN's
        a = np.array([10.0, 7, 4, 3, 2, 1])
        assert pm_nanmedian(a, 0)==3.5, "nanmedian() on just numbers failed"

        ## a few NaNs (fallback not required)
        a = np.array([10.0, np.nan, 4, np.nan, 2, 1])
        assert pm_nanmedian(a, 0) == 3, "nanmedian() on mixed numbers/NaN failed"

        ## all NaNs (fallback required)
        a = np.array([np.nan, np.nan,  np.nan])
        assert pm_nanmedian(a, 0) == 0, "nanmedian() on all-NaN failed"

        ## empty array (fallback required)
        a = np.array([])
        assert pm_nanmedian(a, 0) == 0, "nanmedian() on all-NaN failed"
