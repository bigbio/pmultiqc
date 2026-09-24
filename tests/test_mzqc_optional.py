"""mzQC export is guarded so it remains robust when pymzqc is unavailable.

These tests pass with and without pymzqc installed. The checks that
need pymzqc itself are skipped when it is absent, and the checks for the absent
case hide the package so they also run on installations that have it.
"""

import importlib
import json
import logging
import subprocess
import sys

import pytest

import pmultiqc.export as mzqc_export


@pytest.fixture(autouse=True)
def reset_availability_cache():
    """Each test sees a fresh, uncomputed availability result."""
    mzqc_export._availability = None
    yield
    mzqc_export._availability = None


class TestIsMzqcAvailable:
    def test_true_when_pymzqc_is_installed(self, monkeypatch):
        monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
        assert mzqc_export.is_mzqc_available() is True

    def test_false_when_pymzqc_is_missing(self, monkeypatch):
        monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
        assert mzqc_export.is_mzqc_available() is False

    def test_looks_up_the_mzqc_package(self, monkeypatch):
        looked_up = []
        monkeypatch.setattr(importlib.util, "find_spec", looked_up.append)
        mzqc_export.is_mzqc_available()
        assert looked_up == ["mzqc"]

    def test_result_is_computed_once(self, monkeypatch):
        calls = []
        monkeypatch.setattr(importlib.util, "find_spec", calls.append)
        for _ in range(3):
            mzqc_export.is_mzqc_available()
        assert len(calls) == 1

    def test_missing_package_is_reported_once_with_install_hint(self, monkeypatch, caplog):
        monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)
        with caplog.at_level(logging.INFO, logger="pmultiqc.export"):
            mzqc_export.is_mzqc_available()
            mzqc_export.is_mzqc_available()

        messages = [r.getMessage() for r in caplog.records]
        assert len(messages) == 1
        assert "pip install pymzqc" in messages[0]
        # skipping is a graceful fallback when the dependency is unavailable
        assert caplog.records[0].levelno == logging.INFO

    def test_nothing_is_logged_when_pymzqc_is_installed(self, monkeypatch, caplog):
        monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
        with caplog.at_level(logging.DEBUG, logger="pmultiqc.export"):
            mzqc_export.is_mzqc_available()
        assert caplog.records == []


# Runs in a fresh interpreter: a None entry in sys.modules makes ``import mzqc``
# raise ImportError and find_spec("mzqc") return None, as on an install without
# the extra, without disturbing the modules already loaded by this test session.
_IMPORT_WITHOUT_PYMZQC = """
import importlib, sys
sys.modules["mzqc"] = None
importlib.import_module(sys.argv[1])
from pmultiqc.export import is_mzqc_available
assert is_mzqc_available() is False, "pymzqc reported as available"
leaked = sorted(m for m in sys.modules if m.startswith("mzqc."))
assert not leaked, f"mzqc submodules imported: {leaked}"
"""


class TestWithoutPymzqc:
    @pytest.mark.parametrize(
        "module_name",
        [
            "pmultiqc.export",
            "pmultiqc.modules.maxquant.maxquant",
            "pmultiqc.modules.diann.diann",
            "pmultiqc.modules.quantms.quantms",
        ],
    )
    def test_imports_without_pymzqc(self, module_name):
        """Nothing that exports mzQC may import mzqc at module level."""
        if module_name.startswith("pmultiqc.modules."):
            pytest.importorskip("multiqc")
        result = subprocess.run(
            [sys.executable, "-c", _IMPORT_WITHOUT_PYMZQC, module_name],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr[-2000:]


class TestWithPymzqc:
    def test_exporter_writes_an_mzqc_document(self, tmp_path):
        pytest.importorskip("mzqc")
        from pmultiqc.export.mzqc_exporter import MzQcExporter

        exporter = MzQcExporter(
            pipeline_name="MaxQuant", raw_data={}, output_dir=str(tmp_path)
        )
        metrics = exporter._execute_extraction_loop(
            {"missed_cleavages": {"0": 0.8, "1": 0.2}, "absent_metric": None}
        )
        # the None-valued metric is skipped rather than exported empty
        assert len(metrics) == 1

        path = exporter.export_to_file(metrics, filename="test_qc.mzQC")

        with open(path, encoding="utf-8") as handle:
            document = json.load(handle)
        exported = document["mzQC"]["runQualities"][0]["qualityMetrics"]
        assert len(exported) == len(metrics)
        assert exported[0]["value"] == {"0": 0.8, "1": 0.2}


def _strict_json_load(path):
    """Parse as strict JSON: NaN and Infinity are not valid JSON values."""

    def reject(constant):
        raise ValueError(f"non-standard JSON constant {constant}")

    with open(path, encoding="utf-8") as handle:
        return json.load(handle, parse_constant=reject)


class TestMzqcFileIsValidJson:
    """The exporter must write strict JSON or nothing at all.

    In CI the quantms exports hit a numpy int64 half-way through json.dump,
    which left a truncated file on disk, and the MaxQuant DIA export wrote
    bare NaN values, which strict JSON parsers reject.
    """

    @pytest.fixture
    def exporter(self, tmp_path):
        pytest.importorskip("mzqc")
        from pmultiqc.export.mzqc_exporter import MzQcExporter

        return MzQcExporter(pipeline_name="quantms", raw_data={}, output_dir=str(tmp_path))

    @staticmethod
    def _export(exporter, registry):
        metrics = exporter._execute_extraction_loop(registry)
        return exporter.export_to_file(metrics, filename="qc.mzQC")

    def test_numpy_scalars_nested_in_values(self, exporter):
        np = pytest.importorskip("numpy")
        path = self._export(
            exporter,
            {"identified_spectra": {"0 ~ 100": np.int64(42), "100 ~ 200": np.float32(1.5)}},
        )
        value = _strict_json_load(path)["mzQC"]["runQualities"][0]["qualityMetrics"][0]["value"]
        assert value == {"0 ~ 100": 42, "100 ~ 200": 1.5}

    def test_numpy_integer_keys(self, exporter):
        np = pytest.importorskip("numpy")
        path = self._export(exporter, {"charge_states": {np.int64(2): 10, np.int64(3): 5}})
        value = _strict_json_load(path)["mzQC"]["runQualities"][0]["qualityMetrics"][0]["value"]
        assert value == {"2": 10, "3": 5}

    def test_non_finite_numbers_become_null(self, exporter):
        np = pytest.importorskip("numpy")
        path = self._export(
            exporter,
            {"mass_error": [1.0, float("nan"), np.float64("inf"), -np.inf]},
        )
        value = _strict_json_load(path)["mzQC"]["runQualities"][0]["qualityMetrics"][0]["value"]
        assert value == [1.0, None, None, None]

    def test_pandas_objects_nested_in_values(self, exporter):
        pd = pytest.importorskip("pandas")
        frame = pd.DataFrame({"run": ["a", "b"], "count": [3, None]})
        path = self._export(exporter, {"peptide_intensity": {"table": frame}})
        value = _strict_json_load(path)["mzQC"]["runQualities"][0]["qualityMetrics"][0]["value"]
        assert value == {"table": [{"run": "a", "count": 3.0}, {"run": "b", "count": None}]}

    def test_failed_write_leaves_no_partial_file(self, exporter, tmp_path, monkeypatch):
        target = tmp_path / "qc.mzQC"
        metrics = exporter._execute_extraction_loop({"missed_cleavages": {"0": 1}})

        def broken_dumps(*args, **kwargs):
            raise TypeError("simulated serialisation failure")

        exporter_module = sys.modules[type(exporter).__module__]
        monkeypatch.setattr(exporter_module.json, "dumps", broken_dumps)
        with pytest.raises(OSError):
            exporter.export_to_file(metrics, filename="qc.mzQC")
        assert not target.exists()
        assert list(tmp_path.iterdir()) == []

    def test_no_metrics_returns_empty_path(self, exporter, tmp_path):
        assert exporter.export_to_file([], filename="qc.mzQC") == ""
        assert list(tmp_path.iterdir()) == []
