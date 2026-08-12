"""Tests for deriving an experimental design from quantms.io run/sample parquet."""

import numpy as np
import pandas as pd

from pmultiqc.modules.qpx.qpx_design import (
    _conditions,
    _label_number,
    _number_samples,
    build_design_from_parquet,
)


def _run_df(rows=None):
    """A run.parquet-shaped frame: 2 samples x 3 runs, LFQ, single fraction."""
    if rows is None:
        rows = [
            ("run_A_01", "Sample 1", 1, 1),
            ("run_A_02", "Sample 1", 1, 1),
            ("run_B_01", "Sample 2", 1, 1),
        ]
    return pd.DataFrame(
        {
            "run_accession": [r[0] for r in rows],
            "run_file_name": [r[0] for r in rows],
            "file_name": [f"{r[0]}.raw" for r in rows],
            "samples": [
                np.array(
                    [
                        {
                            "sample_accession": r[1],
                            "label": "LFQ",
                            "biological_replicate": r[2],
                            "technical_replicate": 1,
                        }
                    ],
                    dtype=object,
                )
                for r in rows
            ],
            "fraction": [r[3] for r in rows],
            "enzymes": [np.array(["Trypsin"], dtype=object) for _ in rows],
        }
    )


def _sample_df(spiked=("12500 amol", "125 amol")):
    return pd.DataFrame(
        {
            "sample_accession": ["Sample 1", "Sample 2"],
            "organism": ["Homo sapiens", "Homo sapiens"],
            "cell_line": ["K562", "K562"],
            "spiked_compound": list(spiked),
        }
    )


class TestBuildDesign:
    def test_produces_read_openms_design_schema(self):
        sample_df, file_df = build_design_from_parquet(_run_df(), _sample_df())

        assert list(file_df.columns) == [
            "FractionGroup", "Fraction", "Label", "Sample", "Filename", "Run",
        ]
        assert list(sample_df.columns) == ["Sample", "Condition", "BioReplicate"]

    def test_run_to_sample_mapping(self):
        _, file_df = build_design_from_parquet(_run_df(), _sample_df())

        assert len(file_df) == 3
        assert dict(zip(file_df["Run"], file_df["Sample"])) == {
            "run_A_01": 1,
            "run_A_02": 1,
            "run_B_01": 2,
        }

    def test_sample_and_label_are_integers(self):
        """Downstream code does file_df["Sample"].astype(int)."""
        sample_df, file_df = build_design_from_parquet(_run_df(), _sample_df())

        file_df["Sample"].astype(int)  # must not raise
        assert all(isinstance(v, (int, np.integer)) for v in file_df["Sample"])
        assert all(isinstance(v, (int, np.integer)) for v in file_df["Label"])
        assert all(int(v) == 1 for v in file_df["Label"]), "LFQ is a single channel"
        assert all(isinstance(int(v), int) for v in sample_df["BioReplicate"])

    def test_condition_uses_the_varying_characteristic(self):
        sample_df, _ = build_design_from_parquet(_run_df(), _sample_df())

        # cell_line is constant, spiked_compound varies -> spiked_compound wins.
        assert dict(zip(sample_df["Sample"], sample_df["Condition"])) == {
            1: "12500 amol",
            2: "125 amol",
        }

    def test_multiple_varying_characteristics_use_key_value_form(self):
        samples = _sample_df()
        samples["cell_line"] = ["K562", "HeLa"]
        sample_df, _ = build_design_from_parquet(_run_df(), samples)

        for condition in sample_df["Condition"]:
            assert "=" in condition and ";" in condition
            assert "cell_line=" in condition and "spiked_compound=" in condition

    def test_no_varying_characteristic_falls_back_to_a_constant(self):
        sample_df, _ = build_design_from_parquet(
            _run_df(), _sample_df(spiked=("125 amol", "125 amol"))
        )
        assert set(sample_df["Condition"]) == {"Homo sapiens"}

    def test_works_without_sample_parquet(self):
        sample_df, file_df = build_design_from_parquet(_run_df(), None)

        assert len(file_df) == 3
        assert set(sample_df["Sample"]) == {1, 2}
        assert set(sample_df["Condition"]) == {"not available"}

    def test_tmt_labels_map_to_channels(self):
        run_df = _run_df()
        # One multiplexed run carrying both samples in two TMT channels.
        multiplexed = np.array(
            [
                {"sample_accession": "Sample 1", "label": "TMT126",
                 "biological_replicate": 1, "technical_replicate": 1},
                {"sample_accession": "Sample 2", "label": "TMT127N",
                 "biological_replicate": 1, "technical_replicate": 1},
            ],
            dtype=object,
        )
        samples = list(run_df["samples"])
        samples[0] = multiplexed
        run_df["samples"] = samples
        _, file_df = build_design_from_parquet(run_df, _sample_df())

        first_run = file_df[file_df["Run"] == "run_A_01"]
        assert sorted(first_run["Label"]) == [1, 2]
        assert sorted(first_run["Sample"]) == [1, 2]

    def test_multi_fraction_runs_are_kept_separate(self):
        run_df = _run_df(
            [("run_A_f1", "Sample 1", 1, 1), ("run_A_f2", "Sample 1", 1, 2)]
        )
        _, file_df = build_design_from_parquet(run_df, _sample_df())

        assert sorted(file_df["Fraction"]) == [1, 2]
        assert set(file_df["Sample"]) == {1}


class TestGracefulDegradation:
    """A project that cannot supply a design must fall back, never emit a broken one."""

    def test_none_run_df(self):
        assert build_design_from_parquet(None, _sample_df()) == (None, None)

    def test_empty_run_df(self):
        assert build_design_from_parquet(pd.DataFrame(), _sample_df()) == (None, None)

    def test_run_df_without_samples_column(self):
        df = pd.DataFrame({"run_file_name": ["a"], "fraction": [1]})
        assert build_design_from_parquet(df, _sample_df()) == (None, None)

    def test_run_df_with_empty_sample_lists(self):
        df = _run_df()
        df["samples"] = [np.array([], dtype=object) for _ in range(len(df))]
        assert build_design_from_parquet(df, _sample_df()) == (None, None)


class TestHelpers:
    def test_number_samples_prefers_trailing_digits(self):
        assert _number_samples(["Sample 1", "Sample 2"]) == {"Sample 1": 1, "Sample 2": 2}

    def test_number_samples_falls_back_to_order_when_ambiguous(self):
        # Duplicate trailing numbers cannot be used as ids.
        assert _number_samples(["A_1", "B_1"]) == {"A_1": 1, "B_1": 2}

    def test_number_samples_falls_back_when_no_digits(self):
        assert _number_samples(["alpha", "beta"]) == {"alpha": 1, "beta": 2}

    def test_label_number(self):
        assert _label_number("LFQ") == 1
        assert _label_number(None) == 1
        assert _label_number("TMT126") == 1
        assert _label_number("TMT131C") == 11
        assert _label_number("3") == 3
        assert _label_number("something-unknown") == 1

    def test_conditions_ignores_accession_columns(self):
        # sample_accession varies by definition; it must not become the Condition.
        df = pd.DataFrame(
            {"sample_accession": ["S1", "S2"], "organism": ["Homo sapiens"] * 2}
        )
        assert set(_conditions(df, ["S1", "S2"]).values()) == {"Homo sapiens"}
