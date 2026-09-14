"""mzQC export is optional: it runs only when pymzqc is installed.

These tests pass with and without the ``mzqc`` extra installed. The checks that
need pymzqc itself are skipped when it is absent, and the checks for the absent
case hide the package so they also run on installations that have it.
"""

import importlib
import json
import logging
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
        monkeypatch.setattr(
            importlib.util, "find_spec", lambda name: looked_up.append(name)
        )
        mzqc_export.is_mzqc_available()
        assert looked_up == ["mzqc"]

    def test_result_is_computed_once(self, monkeypatch):
        calls = []
        monkeypatch.setattr(
            importlib.util, "find_spec", lambda name: calls.append(name)
        )
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
        assert 'pip install "pmultiqc[mzqc]"' in messages[0]
        # skipping is the expected state without the extra, not a problem
        assert caplog.records[0].levelno == logging.INFO

    def test_nothing_is_logged_when_pymzqc_is_installed(self, monkeypatch, caplog):
        monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
        with caplog.at_level(logging.DEBUG, logger="pmultiqc.export"):
            mzqc_export.is_mzqc_available()
        assert caplog.records == []


class TestWithoutPymzqc:
    @pytest.fixture
    def hide_pymzqc(self, monkeypatch):
        """Make ``import mzqc`` fail even where the package is installed."""
        for name in [m for m in sys.modules if m == "mzqc" or m.startswith("mzqc.")]:
            monkeypatch.delitem(sys.modules, name)
        # a None entry in sys.modules makes ``import mzqc`` raise ImportError and
        # find_spec("mzqc") return None, without touching any other import
        monkeypatch.setitem(sys.modules, "mzqc", None)

    def test_export_package_imports_without_pymzqc(self, hide_pymzqc, monkeypatch):
        monkeypatch.delitem(sys.modules, "pmultiqc.export")
        module = importlib.import_module("pmultiqc.export")
        assert module.is_mzqc_available() is False

    @pytest.mark.parametrize(
        "module_name",
        [
            "pmultiqc.modules.maxquant.maxquant",
            "pmultiqc.modules.diann.diann",
            "pmultiqc.modules.quantms.quantms",
        ],
    )
    def test_plugin_modules_import_without_pymzqc(self, hide_pymzqc, monkeypatch, module_name):
        """The modules that export mzQC must not import mzqc at module level."""
        pytest.importorskip("multiqc")
        # re-import from scratch; monkeypatch puts the original module back afterwards
        monkeypatch.delitem(sys.modules, module_name, raising=False)
        importlib.import_module(module_name)


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
        exported = document["mzQC"]["runQuality"][0]["qualityMetrics"]
        assert len(exported) == len(metrics)
        assert exported[0]["value"] == {"0": 0.8, "1": 0.2}
