"""Tests for LFQ groups written by quantms --lfq_group_by (lfq_groups.tsv)."""

# Each group is quantified by its own ProteomicsLFQ run, so summaries that compare runs
# must stay within a group.

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pmultiqc.modules.qpx.qpx_heatmap import _variance_scores, calculate_qpx_heatmap
from pmultiqc.modules.qpx.qpx_io import read_lfq_groups, split_by_lfq_group
from pmultiqc.modules.qpx.qpx_quant import calculate_intensity_std

HEADER = ["group_id", "run_id", "sdrf_file", "design_file"]


def _write_manifest(path, header, rows):
    path.write_text("\n".join("\t".join(row) for row in [header, *rows]) + "\n")
    return path


class TestReadLfqGroups(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp_path = Path(tmp.name)

    def test_labels_runs_by_their_sdrf_values(self):
        manifest = _write_manifest(
            self.tmp_path / "lfq_groups.tsv",
            [*HEADER, "source name", "instrument"],
            [
                ["lfq_group_a", "r1", "a.sdrf.tsv", "a_design.tsv", "S1", "QE"],
                ["lfq_group_a", "r2", "a.sdrf.tsv", "a_design.tsv", "S1", "QE"],
                ["lfq_group_b", "r3", "b.sdrf.tsv", "b_design.tsv", "S2", "QE"],
            ],
        )

        self.assertEqual(
            read_lfq_groups([manifest]), {"r1": "S1 | QE", "r2": "S1 | QE", "r3": "S2 | QE"}
        )

    def test_falls_back_to_the_group_id(self):
        manifest = _write_manifest(
            self.tmp_path / "lfq_groups.tsv", HEADER, [["lfq_group_a", "r1", "a.sdrf.tsv", "a.tsv"]]
        )

        self.assertEqual(read_lfq_groups([manifest]), {"r1": "lfq_group_a"})

    def test_ignores_a_table_without_run_ids(self):
        other = _write_manifest(self.tmp_path / "lfq_groups.tsv", ["sample", "group"], [["s1", "g1"]])

        self.assertEqual(read_lfq_groups([other]), {})


class TestSplitByLfqGroup(unittest.TestCase):
    def test_without_groups_the_table_stays_whole(self):
        df = pd.DataFrame({"run": ["r1", "r2"]})

        [(group, rows)] = list(split_by_lfq_group(df, {}))

        self.assertIsNone(group)
        self.assertIs(rows, df)

    def test_runs_without_a_group_are_left_out(self):
        df = pd.DataFrame({"run": ["r1", "r2", "r3"], "x": [1, 2, 3]})

        parts = {
            group: rows["x"].tolist()
            for group, rows in split_by_lfq_group(df, {"r1": "g1", "r3": "g2"})
        }

        self.assertEqual(parts, {"g1": [1], "g2": [3]})


class TestHeatmapWithinGroups(unittest.TestCase):
    @staticmethod
    def _metric(result, metric):
        scores, xnames, runs = result
        column = xnames.index(metric)
        return {run: row[column] for run, row in zip(runs, scores)}

    def test_runs_are_compared_only_with_their_group(self):
        # r1 and r2 share no peptidoform and no charge state.
        psm = pd.DataFrame(
            {
                "run": ["r1", "r1", "r2", "r2"],
                "peptidoform": ["A", "B", "C", "D"],
                "charge": [2, 2, 3, 3],
            }
        )

        project = calculate_qpx_heatmap(psm, None, None, {}, "CONT")
        grouped = calculate_qpx_heatmap(
            psm, None, None, {}, "CONT", run_groups={"r1": "g1", "r2": "g2"}
        )

        for metric in ["Charge", "Pep Missing Values"]:
            with self.subTest(metric=metric):
                self.assertEqual(self._metric(project, metric), {"r1": 0.5, "r2": 0.5})
                self.assertEqual(self._metric(grouped, metric), {"r1": 1.0, "r2": 1.0})

    def test_variance_uses_the_group_median(self):
        groups = {"r1": "g1", "r2": "g2"}

        self.assertEqual(
            _variance_scores({"r1": 1.0, "r2": 0.5}, groups), {"r1": 1.0, "r2": 1.0}
        )


class TestIntensityStdWithinGroups(unittest.TestCase):
    def test_only_conditions_in_several_groups_are_split(self):
        runs = ["r1", "r2", "r3", "r4", "r5", "r6"]
        feature = pd.DataFrame(
            {
                "run": runs,
                "peptidoform": ["A"] * 6,
                "intensities": [
                    [{"label": "LFQ", "intensity": v}]
                    for v in (100.0, 200.0, 400.0, 800.0, 100.0, 300.0)
                ],
            }
        )
        sample_df = pd.DataFrame({"Sample": range(1, 7), "Condition": ["A"] * 4 + ["B"] * 2})
        file_df = pd.DataFrame({"Sample": range(1, 7), "Run": runs})
        groups = {"r1": "g1", "r2": "g1", "r3": "g2", "r4": "g2", "r5": "g2", "r6": "g2"}

        self.assertEqual(set(calculate_intensity_std(feature, sample_df, file_df)), {"A", "B"})
        self.assertEqual(
            set(calculate_intensity_std(feature, sample_df, file_df, groups)),
            {"g1: A", "g2: A", "B"},
        )
