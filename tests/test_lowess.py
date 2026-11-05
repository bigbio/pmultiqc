"""Check LOWESS interpolation accuracy and finite-input handling."""

import numpy as np
import pandas as pd
import pytest
from statsmodels.nonparametric import smoothers_lowess

from pmultiqc.modules.common.dia_utils import cal_rt_irt_loess


@pytest.mark.parametrize("frac", [0.2, 0.3, 0.6])
def test_interpolated_curve_stays_close_to_exact_fit(monkeypatch, frac):
    """Keep binned curves within 0.05 RT units on seeded nonlinear inputs."""
    rng = np.random.default_rng(734)
    x = np.sort(rng.uniform(-50, 150, 2000))
    y = 30 + 0.25 * x + 3 * np.sin(x / 30) + rng.normal(0, 0.3, len(x))
    y[::100] += 10
    frame = pd.DataFrame({"Run": "a", "iRT": x, "RT": y})
    shifted = frame.assign(Run="b", iRT=x * 0.1 + 500, RT=y + 20)
    frame = pd.concat([frame, shifted], ignore_index=True).sample(frac=1, random_state=734)
    before = frame.copy(deep=True)
    approximate = cal_rt_irt_loess(frame, frac=frac)
    original_lowess = smoothers_lowess.lowess

    def exact_lowess(*args, **kwargs):
        """Run the same fitting library with interpolation disabled."""
        kwargs["delta"] = 0
        return original_lowess(*args, **kwargs)

    monkeypatch.setattr(smoothers_lowess, "lowess", exact_lowess)
    exact = cal_rt_irt_loess(frame, frac=frac)
    assert set(approximate) == {"a", "b"}
    for run, curve in approximate.items():
        assert list(curve) == list(exact[run])
        np.testing.assert_allclose(
            list(curve.values()), list(exact[run].values()), atol=0.05, rtol=0
        )
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("seed", [734, 42, 2026])
def test_separated_clusters_stay_close_to_exact_fit(monkeypatch, seed):
    """Limit interpolation error to 0.03 RT units on seeded separated clusters."""
    rng = np.random.default_rng(seed)
    x = np.concatenate([rng.uniform(0, 1, 500), rng.uniform(99, 100, 500)])
    frame = pd.DataFrame({"Run": "clustered", "iRT": x, "RT": 10 + 0.25 * x + np.sin(x)})
    original_lowess = smoothers_lowess.lowess
    fitted = []

    def record_lowess(*args, **kwargs):
        """Capture the fit before binning can average away interpolation errors."""
        result = original_lowess(*args, **kwargs)
        fitted.append(result)
        return result

    monkeypatch.setattr(smoothers_lowess, "lowess", record_lowess)
    approximate = cal_rt_irt_loess(frame)

    def exact_lowess(*args, **kwargs):
        """Capture the corresponding fit with interpolation disabled."""
        kwargs["delta"] = 0
        return record_lowess(*args, **kwargs)

    monkeypatch.setattr(smoothers_lowess, "lowess", exact_lowess)
    exact = cal_rt_irt_loess(frame)
    assert list(approximate["clustered"]) == list(exact["clustered"])
    np.testing.assert_array_equal(fitted[0][:, 0], fitted[1][:, 0])
    np.testing.assert_allclose(fitted[0][:, 1], fitted[1][:, 1], atol=0.03, rtol=0)
    np.testing.assert_allclose(
        list(approximate["clustered"].values()),
        list(exact["clustered"].values()),
        atol=0.03,
        rtol=0,
    )


def test_nonfinite_pairs_match_clean_input():
    """Discard invalid RT/iRT pairs before computing bins and delta."""
    x = np.linspace(0, 100, 200)
    clean = pd.DataFrame({"Run": "valid", "iRT": x, "RT": 2 * x + 10})
    invalid = pd.DataFrame(
        {
            "Run": ["valid"] * 5,
            "iRT": [np.nan, np.inf, -np.inf, 1000, 2000],
            "RT": [10, 10, 10, np.nan, np.inf],
        }
    )
    dirty = pd.concat([clean, invalid], ignore_index=True)
    before = dirty.copy(deep=True)
    assert cal_rt_irt_loess(dirty) == cal_rt_irt_loess(clean)
    pd.testing.assert_frame_equal(dirty, before)


@pytest.mark.parametrize("values", [[], [1.0], [1.0, 1.0], [np.nan, np.inf]])
def test_no_usable_irt_range_is_skipped(values):
    """Return no curve when there is no finite iRT range to fit."""
    frame = pd.DataFrame({"Run": ["a"] * len(values), "iRT": values, "RT": values})
    assert cal_rt_irt_loess(frame) is None


def test_degenerate_runs_do_not_hide_valid_runs(caplog):
    """Skip single-point and constant-iRT runs while retaining usable curves."""
    frame = pd.DataFrame(
        {
            "Run": ["valid"] * 20 + ["constant"] * 3 + ["single"],
            "iRT": list(range(20)) + [5, 5, 5, 8],
            "RT": list(range(20)) + [2, 3, 4, 9],
        }
    )
    frame["Run"] = pd.Categorical(
        frame["Run"], categories=["valid", "constant", "single", "unused"]
    )
    result = cal_rt_irt_loess(frame)
    assert set(result) == {"valid"}
    assert result["valid"]
    assert "constant" in caplog.text and "single" in caplog.text


def test_repeated_irt_values_produce_finite_curve():
    """Repeated precursor iRT observations retain a finite linear trend."""
    x = np.repeat(np.linspace(0, 100, 50), 4)
    frame = pd.DataFrame({"Run": "a", "iRT": x, "RT": 3 * x + 5})
    curve = cal_rt_irt_loess(frame)["a"]
    assert curve
    np.testing.assert_allclose(list(curve.values()), 3 * np.array(list(curve)) + 5, atol=1e-8)


def test_large_input_keeps_existing_guard(caplog):
    """Skip reports above the existing million-row limit before copying data."""
    frame = pd.DataFrame(index=pd.RangeIndex(1_000_001))
    assert cal_rt_irt_loess(frame) is None
    assert "Dataset too large" in caplog.text
