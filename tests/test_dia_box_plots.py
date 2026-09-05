"""DIA box plots must hand MultiQC summary statistics, not raw points, on large runs (#717).

On PXD030304 (5,798 runs) draw_dia_ms1_area serialised every raw MS1 area into the
box plot and polars panicked on a >4 GiB buffer; draw_dia_intensity_std carried
42.7 M points. Above the flat threshold both must pass {min,q1,median,q3,max,mean}.
"""

import numpy as np
import pandas as pd
import pytest

from pmultiqc.modules.common.plots import dia as dia_plots
from pmultiqc.modules.common.plots import general


@pytest.fixture
def capture_box(monkeypatch):
    seen = []
    monkeypatch.setattr(dia_plots.box, "plot", lambda list_of_data_by_sample, pconfig=None: seen.append(list_of_data_by_sample) or "html")
    monkeypatch.setattr(dia_plots, "add_sub_section", lambda **kw: None)
    monkeypatch.setattr(dia_plots, "plot_html_check", lambda h: h, raising=False)
    return seen


def _big_ms1(n_runs=20, per_run=None):
    per_run = per_run or (general.FLAT_THRESHOLD // n_runs + 1)
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "Run": np.repeat([f"run{i}" for i in range(n_runs)], per_run),
        "log_ms1_area": rng.normal(20, 2, n_runs * per_run),
    })


def test_ms1_area_uses_summary_stats_above_threshold(capture_box):
    dia_plots.draw_dia_ms1_area(None, _big_ms1())
    (data,) = capture_box
    assert data, "box.plot received no data"
    for run, stats in data.items():
        assert isinstance(stats, dict), f"{run}: raw list reached box.plot"
        assert {"min", "q1", "median", "q3", "max", "mean"} <= set(stats)


def test_ms1_area_keeps_raw_points_below_threshold(capture_box):
    dia_plots.draw_dia_ms1_area(None, _big_ms1(n_runs=2, per_run=10))
    (data,) = capture_box
    assert all(isinstance(v, list) for v in data.values()), "small reports should keep raw points"


def test_intensity_std_uses_summary_stats_above_threshold(capture_box, monkeypatch):
    n = general.FLAT_THRESHOLD + 10
    fake = [{"Sample 1": list(np.linspace(0, 1, n))}]
    monkeypatch.setattr(dia_plots, "calculate_dia_intensity_std", lambda df, sdrf: fake)
    dia_plots.draw_dia_intensity_std(None, pd.DataFrame(), pd.DataFrame())
    (data,) = capture_box
    ds = data[0] if isinstance(data, list) else data
    assert isinstance(ds["Sample 1"], dict)
    assert {"min", "q1", "median", "q3", "max", "mean"} <= set(ds["Sample 1"])
