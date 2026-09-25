import unittest
from pathlib import Path

from pmultiqc.modules.mzqc.mzqc import (
    MzQCModule,
    MzQCRun,
    _without_run_level_tolerance_general_stats,
)


class TestMzQCGeneralStatsRegression(unittest.TestCase):
    def test_general_stats_removes_run_level_tolerance_columns(self) -> None:
        data = {
            "run-a": {
                "raw_precursor": 16.42,
                "raw_fragment": 22.58,
                "ms1": 4953,
            },
            "run-b": {
                "raw_precursor": 23.13,
                "raw_fragment": 26.67,
                "ms1": 4347,
            },
        }
        headers = {
            "raw_precursor": {
                "title": "Precursor tol.",
                "description": "Suggested precursor search tolerance",
            },
            "raw_fragment": {
                "title": "Fragment tol.",
                "description": "Suggested fragment search tolerance",
            },
            "ms1": {"title": "number of MS1 spectra"},
        }

        filtered_data, filtered_headers = _without_run_level_tolerance_general_stats(
            data,
            headers,
        )

        self.assertEqual(set(filtered_headers), {"ms1"})
        self.assertEqual(
            filtered_data,
            {
                "run-a": {"ms1": 4953},
                "run-b": {"ms1": 4347},
            },
        )

    def test_group_general_stats_payload_is_complete_and_numeric(self) -> None:
        runs = [
            MzQCRun("run-a", Path("a.mzQC"), {}, ()),
            MzQCRun("run-b", Path("b.mzQC"), {}, ()),
            MzQCRun("run-c", Path("c.mzQC"), {}, ()),
        ]
        module = object.__new__(MzQCModule)
        module.runs = runs
        module._run_groups = {
            "run-a": "Experiment group 1",
            "run-b": "Experiment group 1",
            "run-c": "Experiment group 2",
        }
        module._run_group_summaries = {
            "Experiment group 1": {
                "files": 2,
                "precursor_median_ppm": 7.74,
                "fragment_median": 22.307,
                "fragment_unit": "ppm",
            },
            "Experiment group 2": {
                "files": 1,
                "precursor_median_ppm": 13.99,
                "fragment_median": 7.173,
                "fragment_unit": "ppm",
            },
        }

        data, headers = module._mass_accuracy_general_stats_payload()

        self.assertEqual(
            data["run-a"],
            {
                "prideqc_group_precursor_tolerance_ppm": 8,
                "prideqc_group_fragment_tolerance_ppm": 22,
                "prideqc_experiment_group": 1,
            },
        )
        self.assertEqual(data["run-b"], data["run-a"])
        self.assertEqual(
            data["run-c"],
            {
                "prideqc_group_precursor_tolerance_ppm": 14,
                "prideqc_group_fragment_tolerance_ppm": 7,
                "prideqc_experiment_group": 2,
            },
        )
        self.assertEqual(
            headers["prideqc_group_precursor_tolerance_ppm"]["title"],
            "Precursor median tol.",
        )
        self.assertEqual(
            headers["prideqc_group_fragment_tolerance_ppm"]["title"],
            "Fragment median tol.",
        )
        self.assertEqual(headers["prideqc_experiment_group"]["format"], "{:,.0f}")
        self.assertFalse(headers["prideqc_experiment_group"]["scale"])

    def test_group_general_stats_keeps_group_id_when_tolerance_is_unavailable(self) -> None:
        run = MzQCRun("run-a", Path("a.mzQC"), {}, ())
        module = object.__new__(MzQCModule)
        module.runs = [run]
        module._run_groups = {"run-a": "Experiment group 2"}
        module._run_group_summaries = {
            "Experiment group 1": {
                "files": 3,
                "precursor_median_ppm": 9.01,
                "fragment_median": 15.434,
                "fragment_unit": "ppm",
            },
            "Experiment group 2": {"files": 1},
        }

        data, headers = module._mass_accuracy_general_stats_payload()

        self.assertEqual(data["run-a"], {"prideqc_experiment_group": 2})
        self.assertNotIn("prideqc_group_precursor_tolerance_ppm", headers)
        self.assertNotIn("prideqc_group_fragment_tolerance_ppm", headers)
        self.assertIn("prideqc_experiment_group", headers)


if __name__ == "__main__":
    unittest.main()
