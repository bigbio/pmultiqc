"""Mass-error trends for the QPX module.

quantms.io declares ``mass_error_ppm`` on both psm and feature, but it is optional and
commonly left null -- it is empty in every test project seen so far. ``calculated_mz``
and ``observed_mz`` are core fields, though, so the error is derivable:

    ppm = (observed_mz - calculated_mz) / calculated_mz * 1e6
    Da  = (observed_mz - calculated_mz) * charge

This module prefers the stored column when it carries data and falls back to deriving
it, reporting which route it took. If neither is possible the section is skipped.
"""

from __future__ import absolute_import

import numpy as np
import pandas as pd

from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.stats import cal_delta_mass_dict
from pmultiqc.modules.qpx.qpx_io import has_data


log = get_logger("pmultiqc.modules.qpx.qpx_mass_error")

# Beyond this the value is not a calibration error but a mis-assignment; excluded so a
# handful of outliers cannot flatten the histogram.
PPM_SANITY_LIMIT = 1000.0


def calculate_mass_error(psm_df):
    """Return ``(ppm_dict, da_dict, source)`` for the mass-error plots.

    Each dict is the binned structure ``draw_delta_mass_da_ppm`` expects, or None when
    that flavour cannot be produced. ``source`` is "column" or "derived" for logging.
    """
    if psm_df is None or getattr(psm_df, "empty", True):
        return None, None, None

    ppm_series, source = _ppm_series(psm_df)
    da_series = _da_series(psm_df)

    ppm_dict = _binned(ppm_series, "ppm")
    da_dict = _binned(da_series, "delta_da")

    if ppm_dict is None and da_dict is None:
        log.info(
            "[Mass Error] Skipped: neither 'mass_error_ppm' nor "
            "'calculated_mz'/'observed_mz' carry usable values."
        )
        return None, None, None

    return ppm_dict, da_dict, source


def _ppm_series(psm_df):
    """Prefer the stored ppm column; derive from m/z when it is absent or empty."""
    if has_data(psm_df, "mass_error_ppm"):
        series = pd.to_numeric(psm_df["mass_error_ppm"], errors="coerce")
        log.info("[Mass Error] Using the 'mass_error_ppm' column.")
        return series, "column"

    if has_data(psm_df, "calculated_mz", "observed_mz"):
        calculated = pd.to_numeric(psm_df["calculated_mz"], errors="coerce")
        observed = pd.to_numeric(psm_df["observed_mz"], errors="coerce")
        series = (observed - calculated) / calculated.replace(0, np.nan) * 1e6
        log.info(
            "[Mass Error] 'mass_error_ppm' is empty; deriving ppm from "
            "observed_mz/calculated_mz."
        )
        return series, "derived"

    return None, None


def _da_series(psm_df):
    """Delta mass in Da needs the charge to convert from m/z difference to mass."""
    if not has_data(psm_df, "calculated_mz", "observed_mz", "charge"):
        return None

    calculated = pd.to_numeric(psm_df["calculated_mz"], errors="coerce")
    observed = pd.to_numeric(psm_df["observed_mz"], errors="coerce")
    charge = pd.to_numeric(psm_df["charge"], errors="coerce")

    return (observed - calculated) * charge


def _binned(series, name):
    """Bin a mass-error series for plotting, or None when it has no usable spread."""
    if series is None:
        return None

    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    if name == "ppm":
        series = series[series.abs() <= PPM_SANITY_LIMIT]
    if series.empty:
        return None

    frame = pd.DataFrame({name: series})
    # cal_delta_mass_dict bins with pd.value_counts(bins=...), which needs a spread.
    if frame[name].nunique() < 2:
        return None

    return cal_delta_mass_dict(frame, name)
