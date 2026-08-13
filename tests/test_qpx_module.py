"""End-to-end tests for QpxModule.

The helper modules are unit-tested elsewhere; this drives the module itself -- data
discovery, the section orchestration in draw_plots, and the degradation paths -- against
real (subsampled) quantms.io parquet files, so the wiring is covered rather than just
the pieces.
"""

from pathlib import Path

import pytest

from pmultiqc.modules.qpx.qpx import QpxModule


QPX_DATA_DIR = Path("tests/resources/qpx")

SUB_SECTION_KEYS = [
    "experiment", "summary", "identification", "search_engine", "contaminants",
    "quantification", "ms1", "ms2", "mass_error", "rt_qc",
]

# Search-pattern key -> the file suffix the module expects for it.
PATTERN_SUFFIX = {
    "pmultiqc/qpx_psm": ".psm.parquet",
    "pmultiqc/qpx_pg": ".pg.parquet",
    "pmultiqc/qpx_feature": ".feature.parquet",
    "pmultiqc/qpx_run": ".run.parquet",
    "pmultiqc/qpx_sample": ".sample.parquet",
}


def _heatmap_colors():
    return [[round(i * 0.1, 1), "#ff0000"] for i in range(11)]


def _find_log_files(directory, exclude=()):
    """Stand in for MultiQC's find_log_files over a directory of parquet files."""
    directory = Path(directory)

    def finder(search_pattern, filecontents=False):
        suffix = PATTERN_SUFFIX.get(search_pattern)
        if suffix is None:
            return
        for path in sorted(directory.glob(f"*{suffix}")):
            if path.name in exclude:
                continue
            yield {"root": str(path.parent), "fn": path.name}

    return finder


def _module(directory=QPX_DATA_DIR, exclude=()):
    return QpxModule(
        _find_log_files(directory, exclude),
        {key: [] for key in SUB_SECTION_KEYS},
        _heatmap_colors(),
    )


@pytest.fixture(autouse=True)
def _multiqc_config(monkeypatch):
    """MultiQC's global kwargs are read directly by the module."""
    from multiqc import config

    defaults = {
        "remove_decoy": False,
        "contaminant_affix": "CONT",
        "keep_raw": False,
        "condition": None,
        "disable_hoverinfo": False,
    }
    monkeypatch.setattr(config, "kwargs", {**getattr(config, "kwargs", {}), **defaults})


@pytest.mark.skipif(not QPX_DATA_DIR.exists(), reason="QPX fixtures not available")
class TestQpxModuleDataLoading:
    def test_loads_all_tables(self):
        module = _module()

        assert module.get_data() is True
        assert module.psm_df_valid and module.pg_df_valid and module.feature_df_valid
        assert module.run_df is not None and module.sample_parquet_df is not None

    def test_identification_frame_prefers_psm(self):
        module = _module()
        module.get_data()

        assert module.id_source == "psm"
        assert module.id_df_valid

    def test_falls_back_to_feature_without_psm(self):
        """DIA-NN projects have no psm.parquet at all."""
        module = _module(exclude={"test.psm.parquet"})
        module.get_data()

        assert module.id_source == "feature"
        assert module.id_df_valid

    def test_reports_no_data_for_an_empty_directory(self, tmp_path):
        module = _module(tmp_path)

        assert module.get_data() is False


@pytest.mark.skipif(not QPX_DATA_DIR.exists(), reason="QPX fixtures not available")
class TestQpxModuleSections:
    @staticmethod
    def _run(exclude=()):
        module = _module(exclude=exclude)
        assert module.get_data() is True
        module.draw_plots()
        return module

    @staticmethod
    def _plot_ids(module):
        ids = set()
        for section in module.sub_sections.values():
            for entry in section:
                plot = entry.get("plot") if isinstance(entry, dict) else None
                config = getattr(plot, "pconfig", None)
                anchor = getattr(config, "id", None) or getattr(plot, "anchor", None)
                if anchor:
                    ids.add(str(anchor))
        return ids

    def test_draw_plots_completes(self):
        module = self._run()

        assert any(module.sub_sections[key] for key in SUB_SECTION_KEYS)

    def test_core_sections_are_populated(self):
        module = self._run()

        for key in ["summary", "identification", "quantification"]:
            assert module.sub_sections[key], f"{key} produced no plots"

    def test_design_is_derived_without_an_external_file(self):
        module = self._run()

        assert module.enable_exp, "derived design must be marked enabled"
        assert not module.sample_df.empty and not module.file_df.empty
        assert module.sub_sections["experiment"], "no experimental design section"

    def test_missed_cleavages_feed_the_heatmap(self):
        module = self._run()

        assert module.missed_cleavages_by_run, "heatmap needs the per-run counts"

    def test_unambiguous_peptides_are_counted(self):
        module = self._run()

        runs = module.cal_num_table_data.get("ms_runs", {})
        assert runs, "no per-run statistics"
        assert any(
            stats.get("unique_peptide_num") not in (None, "")
            for stats in runs.values()
        ), "#Unambiguous Peptide IDs must not be blank"

    def test_survives_without_pg(self):
        module = self._run(exclude={"test.pg.parquet"})

        assert module.sub_sections["identification"], "identification needs no pg table"

    def test_survives_without_feature(self):
        module = self._run(exclude={"test.feature.parquet"})

        assert module.sub_sections["summary"]

    def test_survives_with_only_psm(self):
        module = self._run(
            exclude={"test.pg.parquet", "test.feature.parquet",
                     "test.run.parquet", "test.sample.parquet"}
        )

        assert module.sub_sections["identification"]


class TestHostDelegation:
    """QpxModule can be hosted by another module (quantms) rather than run standalone."""

    def test_section_group_registration_can_be_suppressed(self):
        module = _module()
        module.register_section_groups = False
        assert module.get_data() is True

        module.draw_plots()

        # Sections are still populated; only the group registration is the host's job.
        assert module.sub_sections["identification"]

    def test_design_is_derived_even_when_not_rendered(self):
        """The host renders the design, but the sample-level plots still need it."""
        module = _module()
        module.draw_experimental_design = False
        module.get_data()
        module.draw_plots()

        assert not module.file_df.empty, "run -> sample mapping must still be derived"
        assert "Sample" in module.file_df.columns
        assert not module.sub_sections["experiment"], "host owns this section"
