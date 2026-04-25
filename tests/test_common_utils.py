import numpy as np
import pandas as pd

from pmultiqc.modules.common.common_utils import (
    ISOTOPE_STEP_DA,
    isotope_corrected_mz_delta,
)


class TestIsotopeCorrectedMzDelta:
    """Verify that the precursor mass-delta is shifted by the isotope step
    recorded per PSM in `opt_global_isotope_error`, matching how search
    engines score M+1 / M+2 precursor picks."""

    def test_no_correction_when_metadata_missing(self):
        psm = pd.DataFrame({
            "exp_mass_to_charge": [500.10, 750.00],
            "calc_mass_to_charge": [500.00, 750.00],
        })
        out = isotope_corrected_mz_delta(psm)
        np.testing.assert_allclose(out.to_numpy(), [0.10, 0.00])

    def test_zero_isotope_error_is_identity(self):
        psm = pd.DataFrame({
            "exp_mass_to_charge": [500.10, 750.00],
            "calc_mass_to_charge": [500.00, 750.00],
            "charge": [2, 3],
            "opt_global_isotope_error": [0, 0],
        })
        out = isotope_corrected_mz_delta(psm)
        np.testing.assert_allclose(out.to_numpy(), [0.10, 0.00])

    def test_isotope_step_subtracted_per_charge(self):
        # exp_mz - calc_mz is inflated by iso_err * 1.00335 / z when the
        # precursor was picked at M+N. Subtracting that recovers the true
        # delta (here, exactly zero for a perfectly-matched peptide).
        psm = pd.DataFrame({
            "exp_mass_to_charge": [
                500.00 + ISOTOPE_STEP_DA / 2,  # M+1 at z=2
                750.00 + 2 * ISOTOPE_STEP_DA / 3,  # M+2 at z=3
            ],
            "calc_mass_to_charge": [500.00, 750.00],
            "charge": [2, 3],
            "opt_global_isotope_error": [1, 2],
        })
        out = isotope_corrected_mz_delta(psm)
        np.testing.assert_allclose(out.to_numpy(), [0.0, 0.0], atol=1e-12)

    def test_handles_missing_or_zero_charge_without_crashing(self):
        psm = pd.DataFrame({
            "exp_mass_to_charge": [500.10, 750.00, 400.00],
            "calc_mass_to_charge": [500.00, 750.00, 400.00],
            "charge": [2, 0, np.nan],
            "opt_global_isotope_error": [0, 1, 1],
        })
        out = isotope_corrected_mz_delta(psm)
        # z=2 iso=0 -> unchanged; z=0 / z=NaN -> correction skipped (raw delta).
        np.testing.assert_allclose(out.to_numpy(), [0.10, 0.00, 0.00])

    def test_string_isotope_error_coerces_cleanly(self):
        psm = pd.DataFrame({
            "exp_mass_to_charge": [500.00 + ISOTOPE_STEP_DA / 2],
            "calc_mass_to_charge": [500.00],
            "charge": ["2"],
            "opt_global_isotope_error": ["1"],
        })
        out = isotope_corrected_mz_delta(psm)
        np.testing.assert_allclose(out.to_numpy(), [0.0], atol=1e-12)
