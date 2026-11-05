"""QPX retention-time plots tolerate unavailable measurements."""

from unittest.mock import Mock

import pandas as pd
import pytest
from numpy.testing import assert_equal
from pandas.testing import assert_frame_equal

from pmultiqc.modules.qpx.qpx import QpxModule


def _no_files(search_pattern: str, filecontents: bool = False) -> list[dict[str, str]]:
    """Provide an empty file inventory for the synthetic module."""
    del search_pattern, filecontents
    return []


def _module(frame: pd.DataFrame) -> QpxModule:
    """Create an RT-only module using the supplied identification measurements."""
    module = QpxModule(_no_files, {"rt_qc": []}, [])
    module.id_df = frame
    module.id_df_valid = True
    return module


@pytest.mark.parametrize(
    "frame",
    [
        pd.DataFrame({"run": ["run_A", "run_B"]}),
        pd.DataFrame({"run": ["run_A", "run_B"], "rt": [None, None]}),
        pd.DataFrame({"run": ["run_A", "run_B"], "rt": [float("nan"), float("nan")]}),
        pd.DataFrame(
            {"run": ["run_A", "run_B"], "rt": pd.array([pd.NA, pd.NA], dtype="Float64")}
        ),
    ],
    ids=["missing-column", "nulls", "nan", "nullable-floats"],
)
def test_no_rt_measurements_skip_the_plot(frame: pd.DataFrame, monkeypatch: pytest.MonkeyPatch) -> None:
    """Missing and all-null RT must skip drawing without a swallowed exception."""
    draw = Mock()
    monkeypatch.setattr("pmultiqc.modules.qpx.qpx.draw_ids_rt_count", draw)
    original = frame.copy(deep=True)

    # Call directly: the module's plot exception handler must not hide a failure.
    _module(frame).plot_rt()

    draw.assert_not_called()
    assert_frame_equal(frame, original)


def test_only_runs_with_reported_rt_are_plotted(monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep measured zeroes and omit runs whose RT measurements are all missing."""
    draw = Mock()
    monkeypatch.setattr("pmultiqc.modules.qpx.qpx.draw_ids_rt_count", draw)
    frame = pd.DataFrame(
        {"run": ["run_A", "run_A", "run_A", "run_B"], "rt": [0.0, 60.0, None, None]}
    )
    original = frame.copy(deep=True)

    _module(frame).plot_rt()

    draw.assert_called_once()
    counts = draw.call_args.kwargs["rt_count_data"]
    assert_equal(set(counts), {"run_A"})
    assert_equal(sum(counts["run_A"].values()), 2)
    assert_equal(counts["run_A"][0.5], 1, "zero RT is a measurement")
    assert_equal(counts["run_A"][1.5], 1, "QPX seconds become plot minutes")
    assert_frame_equal(frame, original)


def test_complete_rt_data_keep_the_existing_histograms(monkeypatch: pytest.MonkeyPatch) -> None:
    """Complete measurements retain per-run counts and seconds-to-minutes bins."""
    draw = Mock()
    monkeypatch.setattr("pmultiqc.modules.qpx.qpx.draw_ids_rt_count", draw)
    frame = pd.DataFrame(
        {"run": ["run_A", "run_A", "run_B", "run_B"], "rt": [60.0, 120.0, 180.0, 240.0]}
    )

    _module(frame).plot_rt()

    draw.assert_called_once()
    counts = draw.call_args.kwargs["rt_count_data"]
    assert_equal(set(counts), {"run_A", "run_B"})
    assert_equal(sum(counts["run_A"].values()), 2)
    assert_equal(sum(counts["run_B"].values()), 2)
    assert_equal([counts["run_A"][1.5], counts["run_A"][2.5]], [1, 1])
    assert_equal([counts["run_B"][3.5], counts["run_B"][4.5]], [1, 1])
