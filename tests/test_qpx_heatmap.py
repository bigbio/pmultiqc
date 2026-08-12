"""Tests for the QPX QC heatmap and protein-level quantification."""

import numpy as np
import pandas as pd

from pmultiqc.modules.qpx.qpx_heatmap import (
    INTENSITY_REFERENCE,
    METRIC_ORDER,
    _contaminant_scores,
    _oversampling_scores,
    _peptide_intensity_scores,
    _pep_missing_value_scores,
    _variance_scores,
    calculate_qpx_heatmap,
)
from pmultiqc.modules.qpx.qpx_quant import (
    create_qpx_protein_table,
    get_protein_intensity,
)


def _psm_df():
    return pd.DataFrame(
        {
            "run": ["r1"] * 4 + ["r2"] * 4,
            "peptidoform": ["A", "A", "B", "C", "A", "B", "D", "D"],
            "charge": [2, 2, 3, 2, 2, 2, 3, 3],
            "rt": [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0],
            "sequence": ["PEPA", "PEPA", "PEPB", "PEPC", "PEPA", "PEPB", "PEPD", "PEPD"],
        }
    )


def _pg_df(contaminant=None):
    df = pd.DataFrame(
        {
            "run": ["r1", "r1", "r2", "r2"],
            "anchor_protein": ["P1", "CONT_P2", "P1", "CONT_P2"],
            "intensity": [900.0, 100.0, 800.0, 200.0],
            "global_qvalue": [0.001, 0.02, 0.002, 0.03],
        }
    )
    if contaminant is not None:
        df["contaminant"] = contaminant
    return df


def _feature_df(intensity=1000.0):
    return pd.DataFrame(
        {
            "run": ["r1", "r1", "r2", "r2"],
            "sequence": ["PEPA", "PEPB", "PEPA", "PEPD"],
            "anchor_protein": ["P1", "P1", "P1", "P2"],
            "intensities": [
                [{"label": "LFQ", "intensity": intensity}] for _ in range(4)
            ],
        }
    )


class TestHeatmapAssembly:
    def test_full_metric_set(self):
        scores, xnames, ynames = calculate_qpx_heatmap(
            _psm_df(), _pg_df(), _feature_df(), {"r1": {0: 8, 1: 2}, "r2": {0: 5, 1: 5}}, "CONT"
        )

        assert xnames == METRIC_ORDER
        assert ynames == ["r1", "r2"]
        assert len(scores) == 2 and all(len(row) == len(METRIC_ORDER) for row in scores)
        assert all(0.0 <= v <= 1.0 for row in scores for v in row)

    def test_metrics_without_data_are_dropped(self):
        # No pg -> no Contaminants; no feature -> no Peptide Intensity;
        # no missed-cleavage counts -> no Missed Cleavages or its variance.
        _, xnames, _ = calculate_qpx_heatmap(_psm_df(), None, None, {}, "CONT")

        assert xnames == ["Charge", "ID rate over RT", "MS2 OverSampling", "Pep Missing Values"]

    def test_metric_covering_only_some_runs_is_dropped(self):
        """A partially-covered metric must not be padded with invented cells."""
        _, xnames, ynames = calculate_qpx_heatmap(
            _psm_df(), None, None, {"r1": {0: 8, 1: 2}}, "CONT"
        )

        assert ynames == ["r1", "r2"]
        assert "Missed Cleavages" not in xnames

    def test_no_psm_data_returns_empty(self):
        assert calculate_qpx_heatmap(None, _pg_df(), _feature_df(), {}, "CONT") == ([], [], [])
        assert calculate_qpx_heatmap(
            pd.DataFrame(), _pg_df(), _feature_df(), {}, "CONT"
        ) == ([], [], [])

    def test_score_ordering_matches_xnames(self):
        scores, xnames, ynames = calculate_qpx_heatmap(
            _psm_df(), _pg_df(), _feature_df(), {"r1": {0: 10}, "r2": {0: 5, 1: 5}}, "CONT"
        )
        mc = scores[ynames.index("r1")][xnames.index("Missed Cleavages")]
        assert mc == 1.0, "r1 has only zero-missed-cleavage PSMs"
        assert scores[ynames.index("r2")][xnames.index("Missed Cleavages")] == 0.5


class TestContaminants:
    def test_uses_contaminant_flag_when_present(self):
        scores = _contaminant_scores(_pg_df(contaminant=[False, True, False, True]), "CONT")
        # r1: 100 of 1000 contaminant -> 0.9 ; r2: 200 of 1000 -> 0.8
        assert scores == {"r1": 0.9, "r2": 0.8}

    def test_falls_back_to_accession_affix(self):
        scores = _contaminant_scores(_pg_df(), "CONT")
        assert scores == {"r1": 0.9, "r2": 0.8}

    def test_affix_match_is_unioned_with_the_flag(self):
        """A False flag must not hide accessions the affix matches.

        The OpenMS-consensus converter reports contaminant=False for the "Cont_"-prefixed
        accessions these projects actually use, so trusting the flag alone reports a
        confident 0% contamination for data that is 5-10% contaminant.
        """
        scores = _contaminant_scores(_pg_df(contaminant=[False] * 4), "CONT")
        assert scores == {"r1": 0.9, "r2": 0.8}

    def test_case_insensitive_affix(self):
        df = _pg_df()
        df["anchor_protein"] = ["P1", "sp|Cont_P00330|ADH1", "P1", "sp|Cont_P00330|ADH1"]
        df["pg_accessions"] = [["P1"], ["sp|Cont_P00330|ADH1"], ["P1"], ["sp|Cont_P00330|ADH1"]]
        scores = _contaminant_scores(df, "CONT")
        assert scores == {"r1": 0.9, "r2": 0.8}, "'Cont_' must match affix 'CONT'"

    def test_no_intensity_column(self):
        assert _contaminant_scores(pd.DataFrame({"run": ["r1"]}), "CONT") == {}


class TestPeptideIntensity:
    def test_scaled_against_reference(self):
        scores = _peptide_intensity_scores(_feature_df(intensity=INTENSITY_REFERENCE / 4))
        assert scores == {"r1": 0.25, "r2": 0.25}

    def test_clamped_at_one(self):
        scores = _peptide_intensity_scores(_feature_df(intensity=INTENSITY_REFERENCE * 10))
        assert scores == {"r1": 1.0, "r2": 1.0}

    def test_zero_intensities_ignored(self):
        assert _peptide_intensity_scores(_feature_df(intensity=0.0)) == {}


class TestOtherMetrics:
    def test_oversampling_counts_single_hit_precursors(self):
        scores = _oversampling_scores(_psm_df())
        # r1 precursors: A/2 seen twice, B/3 once, C/2 once -> 2 of 3 single
        assert scores["r1"] == 2 / 3
        # r2: A/2 once, B/2 once, D/3 twice -> 2 of 3
        assert scores["r2"] == 2 / 3

    def test_pep_missing_values_is_share_of_global(self):
        scores = _pep_missing_value_scores(_psm_df())
        # global peptidoforms = {A,B,C,D}; r1 has {A,B,C}, r2 has {A,B,D}
        assert scores == {"r1": 0.75, "r2": 0.75}

    def test_variance_scores(self):
        assert _variance_scores({"r1": 1.0, "r2": 0.5}) == {"r1": 0.75, "r2": 0.75}
        assert _variance_scores({}) == {}


class TestProteinQuantTable:
    def test_table_contents(self):
        table, headers = create_qpx_protein_table(
            _pg_df(), _feature_df(), pd.DataFrame(), pd.DataFrame()
        )

        assert set(headers) >= {"ProteinName", "Average Intensity", "Peptides_Number"}
        rows = {entry["ProteinName"]: entry for entry in table.values()}
        assert set(rows) == {"P1", "CONT_P2"}
        # P1: mean(900, 800) = 850
        assert rows["P1"]["Average Intensity"] == float(np.log10(850.0))
        assert rows["P1"]["Peptides_Number"] == 2, "PEPA and PEPB map to P1"
        assert rows["P1"]["Global Q-value"] == 0.001, "best (lowest) q-value"

    def test_sorted_by_intensity_descending(self):
        table, _ = create_qpx_protein_table(
            _pg_df(), _feature_df(), pd.DataFrame(), pd.DataFrame()
        )
        intensities = [entry["Average Intensity"] for entry in table.values()]
        assert intensities == sorted(intensities, reverse=True)

    def test_condition_columns_added_when_design_available(self):
        sample_df = pd.DataFrame({"Sample": [1, 2], "Condition": ["A", "B"]})
        file_df = pd.DataFrame({"Sample": [1, 2], "Run": ["r1", "r2"]})

        table, headers = create_qpx_protein_table(_pg_df(), _feature_df(), sample_df, file_df)

        assert "A" in headers and "B" in headers
        rows = {entry["ProteinName"]: entry for entry in table.values()}
        assert rows["P1"]["A"] == float(np.log10(900.0))
        assert rows["P1"]["B"] == float(np.log10(800.0))

    def test_returns_none_without_usable_data(self):
        assert create_qpx_protein_table(None, None, None, None) == (None, None)
        assert create_qpx_protein_table(pd.DataFrame(), None, None, None) == (None, None)

        no_intensity = pd.DataFrame({"run": ["r1"], "anchor_protein": ["P1"], "intensity": [0.0]})
        assert create_qpx_protein_table(no_intensity, None, None, None) == (None, None)


class TestProteinIntensity:
    def test_by_run_and_by_sample(self):
        file_df = pd.DataFrame({"Sample": [1, 2], "Run": ["r1", "r2"]})
        by_run, by_sample = get_protein_intensity(_pg_df(), file_df)

        assert set(by_run) == {"r1", "r2"}
        assert set(by_sample) == {"Sample 1", "Sample 2"}

    def test_by_run_only_without_design(self):
        by_run, by_sample = get_protein_intensity(_pg_df(), pd.DataFrame())
        assert set(by_run) == {"r1", "r2"}
        assert by_sample == {}

    def test_empty_input(self):
        assert get_protein_intensity(None, None) == [{}, {}]


class TestSelectColumns:
    """The reader omits absent columns, so downstream selection must be guarded."""

    def test_returns_frame_when_all_present(self):
        from pmultiqc.modules.qpx.qpx_io import select_columns

        df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
        out = select_columns(df, ["a", "b"], "ctx")

        assert list(out.columns) == ["a", "b"]

    def test_returns_none_when_a_column_is_missing(self):
        from pmultiqc.modules.qpx.qpx_io import select_columns

        df = pd.DataFrame({"a": [1]})
        assert select_columns(df, ["a", "missing"], "ctx") is None

    def test_returns_none_for_empty_or_absent_frame(self):
        from pmultiqc.modules.qpx.qpx_io import select_columns

        assert select_columns(None, ["a"], "ctx") is None
        assert select_columns(pd.DataFrame(), ["a"], "ctx") is None

    def test_result_is_a_copy(self):
        from pmultiqc.modules.qpx.qpx_io import select_columns

        df = pd.DataFrame({"a": [1]})
        out = select_columns(df, ["a"], "ctx")
        out.loc[0, "a"] = 99

        assert df.loc[0, "a"] == 1, "caller mutations must not touch the source frame"


class TestRunStatWithoutModifications:
    def test_missing_modifications_column_is_tolerated(self):
        from pmultiqc.modules.qpx.qpx_utils import calculate_run_stat

        df = pd.DataFrame({"peptidoform": ["A", "B"]})
        stat, data = calculate_run_stat(df, {"P1"})

        assert stat["peptide_num"] == 2
        assert stat["modified_peptide_num"] == 0
        assert data["modified_peps"] == set()
