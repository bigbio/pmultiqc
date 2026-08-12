"""Tests for the adaptive QPX sections: mass error, oversampling, contaminants, scores."""

import numpy as np
import pandas as pd

from pmultiqc.modules.qpx.qpx_io import has_data
from pmultiqc.modules.qpx.qpx_mass_error import calculate_mass_error
from pmultiqc.modules.qpx.qpx_quant import (
    calculate_intensity_std,
    create_qpx_peptide_table,
    protein_group_key,
    protein_intensity_pca,
)
from pmultiqc.modules.qpx.qpx_sections import (
    calculate_contaminants,
    calculate_oversampling,
    calculate_search_engine_scores,
)


def _scores(*pairs):
    return [{"score_name": n, "score_value": v, "higher_better": False} for n, v in pairs]


class TestHasData:
    """Presence is not usability: the spec makes many fields nullable."""

    def test_true_only_when_populated(self):
        df = pd.DataFrame({"a": [1, None], "empty": [None, None]})

        assert has_data(df, "a")
        assert not has_data(df, "empty"), "present but all-null is not usable"
        assert not has_data(df, "absent")
        assert not has_data(df, "a", "empty"), "all named columns must be usable"

    def test_false_for_empty_or_none_frame(self):
        assert not has_data(None, "a")
        assert not has_data(pd.DataFrame(), "a")


class TestMassError:
    @staticmethod
    def _psm(with_column=False, n=200):
        rng = np.random.default_rng(0)
        calculated = pd.Series(np.linspace(400.0, 900.0, n))
        # ~2 ppm systematic offset plus noise
        observed = calculated * (1 + (2.0 + rng.normal(0, 0.5, n)) / 1e6)
        df = pd.DataFrame(
            {"calculated_mz": calculated, "observed_mz": observed, "charge": [2] * n}
        )
        df["mass_error_ppm"] = (
            pd.Series(np.full(n, 5.0)) if with_column else pd.Series([None] * n)
        )
        return df

    def test_derives_when_column_is_empty(self):
        ppm, da, source = calculate_mass_error(self._psm())

        assert source == "derived"
        assert ppm is not None and da is not None
        centre = max(ppm["count"], key=ppm["count"].get)
        assert 1.0 < centre < 3.0, "derived ppm should centre near the 2 ppm offset"

    def test_prefers_the_column_when_populated(self):
        ppm, _, source = calculate_mass_error(self._psm(with_column=True))

        assert source == "column"
        assert ppm is None, "a constant column has no spread to bin"

    def test_no_mz_and_no_column(self):
        df = pd.DataFrame({"charge": [2, 3]})
        assert calculate_mass_error(df) == (None, None, None)

    def test_da_needs_charge(self):
        df = self._psm().drop(columns=["charge"])
        ppm, da, _ = calculate_mass_error(df)

        assert ppm is not None
        assert da is None, "Da conversion is impossible without charge"


class TestOversampling:
    def test_bins_precursor_repeat_counts(self):
        df = pd.DataFrame(
            {
                "run": ["r1"] * 6,
                # A seen 1x, B 2x, C 3x
                "peptidoform": ["A", "B", "B", "C", "C", "C"],
                "charge": [2, 2, 2, 2, 2, 2],
            }
        )
        data, cats = calculate_oversampling(df)

        assert set(cats) == {"1", "2", ">=3"}
        assert data["r1"] == {"1": 1 / 3 * 100, "2": 1 / 3 * 100, ">=3": 1 / 3 * 100}

    def test_skipped_without_required_columns(self):
        assert calculate_oversampling(pd.DataFrame({"run": ["r1"]})) == (None, None)
        assert calculate_oversampling(None) == (None, None)


class TestContaminants:
    @staticmethod
    def _pg(contaminant):
        return pd.DataFrame(
            {
                "run": ["r1", "r1"],
                "anchor_protein": ["P1", "CONT_P2"],
                "pg_accessions": [["P1"], ["CONT_P2"]],
                "intensity": [900.0, 100.0],
                "contaminant": contaminant,
            }
        )

    def test_uses_the_flag(self):
        per_run, top_n = calculate_contaminants(self._pg([False, True]), "CONT")

        assert per_run["r1"]["Potential Contaminants"] == 10.0
        assert top_n["plot_data"]["r1"]["CONT_P2"] == 10.0
        assert "CONT_P2" in top_n["cats"]

    def test_all_null_flag_falls_back_to_the_affix(self):
        """Null means unknown, so fall through rather than claiming 0%."""
        per_run, _ = calculate_contaminants(self._pg([None, None]), "CONT")

        assert per_run["r1"]["Potential Contaminants"] == 10.0

    def test_all_null_and_no_affix_is_skipped(self):
        assert calculate_contaminants(self._pg([None, None]), "") == (None, None)

    def test_known_clean_is_skipped(self):
        """Known-clean is a real answer, but an all-zero chart conveys nothing."""
        clean = pd.DataFrame(
            {
                "run": ["r1", "r1"],
                "anchor_protein": ["P1", "P2"],
                "pg_accessions": [["P1"], ["P2"]],
                "intensity": [900.0, 100.0],
                "contaminant": [False, False],
            }
        )
        assert calculate_contaminants(clean, "CONT") == (None, None)

    def test_affix_match_overrides_a_false_flag(self):
        """Converters disagree: OpenMS-consensus reports False for Cont_-prefixed
        accessions, so an explicit affix must still find them."""
        mislabelled = pd.DataFrame(
            {
                "run": ["r1", "r1"],
                "anchor_protein": ["sp|P1|A", "sp|Cont_P00330|ADH1_YEAST"],
                "pg_accessions": [["sp|P1|A"], ["sp|Cont_P00330|ADH1_YEAST"]],
                "intensity": [900.0, 100.0],
                "contaminant": [False, False],
            }
        )
        per_run, top_n = calculate_contaminants(mislabelled, "CONT")

        assert per_run["r1"]["Potential Contaminants"] == 10.0, "case-insensitive match"
        assert top_n is not None

    def test_partial_contamination_still_reports_top_n(self):
        per_run, top_n = calculate_contaminants(self._pg([False, True]), "CONT")

        assert per_run["r1"]["Potential Contaminants"] == 10.0
        assert set(top_n["plot_data"]["r1"]) == {"CONT_P2"}


class TestSearchEngineScores:
    def test_prefers_a_named_engine_score(self):
        df = pd.DataFrame(
            {
                "run": ["r1"] * 4,
                "additional_scores": [
                    _scores(("andromeda_score", v), ("q-value", 0.01))
                    for v in (10.0, 50.0, 120.0, 200.0)
                ],
                "posterior_error_probability": [0.1, 0.2, 0.3, 0.4],
            }
        )
        _, name = calculate_search_engine_scores(df)
        assert name == "andromeda_score"

    def test_falls_back_to_pep_then_qvalue(self):
        pep_only = pd.DataFrame(
            {"run": ["r1"] * 4, "posterior_error_probability": [0.1, 0.2, 0.3, 0.4]}
        )
        _, name = calculate_search_engine_scores(pep_only)
        assert name == "posterior_error_probability"

        qvalue_only = pd.DataFrame(
            {
                "run": ["r1"] * 4,
                "additional_scores": [_scores(("q-value", v)) for v in (0.001, 0.01, 0.02, 0.05)],
            }
        )
        _, name = calculate_search_engine_scores(qvalue_only)
        assert name == "q-value"

    def test_bins_adapt_to_the_observed_range(self):
        """A 0-1 PEP and a 0-300 engine score must not share fixed bins."""
        small = pd.DataFrame(
            {"run": ["r1"] * 50, "posterior_error_probability": np.linspace(0, 0.5, 50)}
        )
        large = pd.DataFrame(
            {
                "run": ["r1"] * 50,
                "additional_scores": [_scores(("hyperscore", v)) for v in np.linspace(0, 300, 50)],
            }
        )
        small_data, _ = calculate_search_engine_scores(small)
        large_data, _ = calculate_search_engine_scores(large)

        assert sum(small_data["r1"].values()) > 0
        assert sum(large_data["r1"].values()) > 0
        # Bin labels differ by orders of magnitude between the two.
        assert list(small_data["r1"])[-1] != list(large_data["r1"])[-1]

    def test_no_score_available(self):
        assert calculate_search_engine_scores(pd.DataFrame({"run": ["r1"]})) == (None, None)

    def test_constant_score_has_no_spread(self):
        df = pd.DataFrame({"run": ["r1"] * 5, "posterior_error_probability": [0.1] * 5})
        assert calculate_search_engine_scores(df) == (None, None)


class TestProteinGroupKey:
    def test_keys_on_the_accession_set_not_the_anchor(self):
        """Two distinct groups may share an anchor; they must stay distinct."""
        df = pd.DataFrame(
            {
                "anchor_protein": ["P1", "P1"],
                "pg_accessions": [["P1", "P2"], ["P1", "P3"]],
            }
        )
        keys = protein_group_key(df)

        assert keys.nunique() == 2, "shared anchor must not collapse distinct groups"

    def test_order_independent(self):
        df = pd.DataFrame({"pg_accessions": [["P2", "P1"], ["P1", "P2"]]})
        assert protein_group_key(df).nunique() == 1

    def test_falls_back_to_anchor(self):
        df = pd.DataFrame({"anchor_protein": ["P1", "P2"]})
        assert protein_group_key(df).nunique() == 2


class TestPca:
    def test_returns_a_point_per_run(self):
        rng = np.random.default_rng(1)
        rows = []
        for protein in range(30):
            for run in ["r1", "r2", "r3"]:
                rows.append(
                    {
                        "run": run,
                        "pg_accessions": [f"P{protein}"],
                        "intensity": float(rng.uniform(1e5, 1e7)),
                    }
                )
        result = protein_intensity_pca(pd.DataFrame(rows))

        assert set(result) == {"r1", "r2", "r3"}
        assert all({"x", "y"} == set(v) for v in result.values())

    def test_needs_at_least_two_runs(self):
        df = pd.DataFrame(
            {"run": ["r1"] * 3, "pg_accessions": [["P1"], ["P2"], ["P3"]],
             "intensity": [1.0, 2.0, 3.0]}
        )
        assert protein_intensity_pca(df) is None


class TestPeptideTableAndStd:
    @staticmethod
    def _feature():
        return pd.DataFrame(
            {
                "run": ["r1", "r2", "r1", "r2"],
                "peptidoform": ["A", "A", "B", "B"],
                "charge": [2, 2, 3, 3],
                "anchor_protein": ["P1", "P1", "P2", "P2"],
                "intensities": [
                    [{"label": "LFQ", "intensity": v}] for v in (100.0, 200.0, 300.0, 400.0)
                ],
            }
        )

    def test_peptide_table(self):
        table, headers = create_qpx_peptide_table(self._feature(), pd.DataFrame(), pd.DataFrame())

        rows = {e["PeptideID"]: e for e in table.values()}
        assert set(rows) == {"A", "B"}
        assert rows["A"]["Average Intensity"] == float(np.log10(150.0))
        assert "Charge" in headers and "ProteinName" in headers

    def test_peptide_table_without_intensities(self):
        df = self._feature().drop(columns=["intensities"])
        assert create_qpx_peptide_table(df, None, None) == (None, None)

    def test_intensity_std_needs_a_design(self):
        assert calculate_intensity_std(self._feature(), pd.DataFrame(), pd.DataFrame()) == {}

    def test_intensity_std_per_condition(self):
        sample_df = pd.DataFrame({"Sample": [1, 2], "Condition": ["A", "A"]})
        file_df = pd.DataFrame({"Sample": [1, 2], "Run": ["r1", "r2"]})

        result = calculate_intensity_std(self._feature(), sample_df, file_df)

        assert set(result) == {"A"}
        assert len(result["A"]) == 2, "one std per peptidoform"


class TestContaminantsZeroSignal:
    def test_all_zero_is_skipped_not_drawn(self):
        """0% everywhere is a real answer, but an all-zero bar chart shows nothing."""
        pg = pd.DataFrame(
            {
                "run": ["r1", "r2"],
                "anchor_protein": ["P1", "P2"],
                "pg_accessions": [["P1"], ["P2"]],
                "intensity": [900.0, 100.0],
                "contaminant": [False, False],
            }
        )
        assert calculate_contaminants(pg, "CONT") == (None, None)


class TestSummariseBoxData:
    """MultiQC flips a plot to a static image past 100k points; that image does not
    fill the panel. Summarising to box statistics keeps it interactive."""

    @staticmethod
    def _raw():
        rng = np.random.default_rng(0)
        return [{"r1": list(rng.normal(26, 1.2, 50000)) + [36.1, 20.3]}]

    def test_collapses_to_six_numbers(self):
        from pmultiqc.modules.common.plots.general import summarise_box_data

        stats = summarise_box_data(self._raw())[0]["r1"]

        assert set(stats) == {"min", "q1", "median", "q3", "max", "mean"}

    def test_quartiles_match_the_raw_data(self):
        from pmultiqc.modules.common.plots.general import summarise_box_data

        raw = self._raw()
        expected = np.percentile(raw[0]["r1"], [25, 50, 75])
        stats = summarise_box_data(raw)[0]["r1"]

        assert np.allclose(
            [stats["q1"], stats["median"], stats["q3"]], expected, atol=1e-6
        )

    def test_fences_exclude_outliers(self):
        from pmultiqc.modules.common.plots.general import summarise_box_data

        stats = summarise_box_data(self._raw())[0]["r1"]

        assert stats["max"] < 36.1, "upper fence must exclude the extreme point"
        assert stats["min"] > 20.3, "lower fence must exclude the extreme point"

    def test_point_count_drops_below_the_flat_threshold(self):
        from pmultiqc.modules.common.plots.general import FLAT_THRESHOLD, summarise_box_data

        raw = [{f"r{i}": list(np.random.default_rng(i).normal(26, 1, 40000)) for i in range(6)}]
        assert sum(len(v) for v in raw[0].values()) > FLAT_THRESHOLD

        stats = summarise_box_data(raw)
        assert sum(len(v) for v in stats[0].values()) < FLAT_THRESHOLD

    def test_already_summarised_is_passed_through(self):
        from pmultiqc.modules.common.plots.general import summarise_box_data

        stats = {"r1": {"min": 1, "q1": 2, "median": 3, "q3": 4, "max": 5}}
        assert summarise_box_data([stats])[0] == stats

    def test_preserves_shape(self):
        from pmultiqc.modules.common.plots.general import summarise_box_data

        as_dict = {"r1": list(np.random.default_rng(2).normal(0, 1, 500))}
        assert isinstance(summarise_box_data(as_dict), dict)
        assert isinstance(summarise_box_data([as_dict]), list)


class TestReviewRegressions:
    """Guards for the defects the code review surfaced."""

    def test_feature_columns_cover_the_psm_fallback(self):
        """DIA-NN has no psm.parquet, so feature must carry what those plots read."""
        from pmultiqc.modules.qpx.qpx_io import QPX_COLUMNS

        for column in ["scan", "posterior_error_probability", "additional_scores"]:
            assert column in QPX_COLUMNS["feature"], column

    def test_psm_reads_additional_scores(self):
        """Search-engine scores live here; without it the whole path is dead code."""
        from pmultiqc.modules.qpx.qpx_io import QPX_COLUMNS

        assert "additional_scores" in QPX_COLUMNS["psm"]

    def test_display_name_falls_back_to_the_key_not_an_index(self):
        from pmultiqc.modules.qpx.qpx_quant import create_qpx_protein_table

        pg = pd.DataFrame(
            {
                "run": ["r1", "r1"],
                "pg_accessions": [["P1"], ["P2"]],
                "intensity": [10.0, 20.0],
            }
        )
        table, _ = create_qpx_protein_table(pg, None, None, None)
        names = {e["ProteinName"] for e in table.values()}

        assert names == {"P1", "P2"}, "must not label proteins '0'/'1'"

    def test_peptides_number_header_only_when_populated(self):
        from pmultiqc.modules.qpx.qpx_quant import create_qpx_protein_table

        pg = pd.DataFrame(
            {"run": ["r1"], "pg_accessions": [["P1"]], "intensity": [10.0]}
        )
        # feature keys on a different protein group -> no counts join
        feature = pd.DataFrame(
            {"sequence": ["PEP"], "pg_accessions": [[{"accession": "OTHER"}]]}
        )
        _, headers = create_qpx_protein_table(pg, feature, None, None)

        assert "Peptides_Number" not in headers

    def test_explicit_contaminant_flag_is_not_overridden(self):
        """The affix may only supplement a flag that marks nothing, never flip a True/False."""
        from pmultiqc.modules.qpx.qpx_sections import resolve_contaminant_flags

        df = pd.DataFrame(
            {
                "anchor_protein": ["sp|Cont_A|X", "sp|B|Y"],
                "pg_accessions": [["sp|Cont_A|X"], ["sp|B|Y"]],
                "contaminant": [False, True],
            }
        )
        flags = resolve_contaminant_flags(df, "CONT")

        assert list(flags) == [False, True], "producer flag wins when it marks anything"

    def test_affix_supplements_an_empty_flag(self):
        from pmultiqc.modules.qpx.qpx_sections import resolve_contaminant_flags

        df = pd.DataFrame(
            {
                "anchor_protein": ["sp|Cont_A|X", "sp|B|Y"],
                "pg_accessions": [["sp|Cont_A|X"], ["sp|B|Y"]],
                "contaminant": [False, False],
            }
        )
        assert list(resolve_contaminant_flags(df, "CONT")) == [True, False]

    def test_sample_merge_survives_mixed_dtypes(self):
        """Sample is int from the derived design and str from an OpenMS design file."""
        from pmultiqc.modules.qpx.qpx_quant import create_qpx_protein_table

        pg = pd.DataFrame(
            {"run": ["r1"], "pg_accessions": [["P1"]], "intensity": [10.0]}
        )
        sample_df = pd.DataFrame({"Sample": ["1"], "Condition": ["A"]})   # str
        file_df = pd.DataFrame({"Sample": [1], "Run": ["r1"]})            # int

        table, headers = create_qpx_protein_table(pg, None, sample_df, file_df)

        assert "A" in headers, "condition column must survive the dtype mismatch"
