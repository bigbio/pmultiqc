"""The DIA-NN parquet reader must not materialise the whole report (bigbio/pmultiqc#717)."""

import ast
import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pmultiqc.modules.common.ms import diann as diann_reader

PKG = Path(diann_reader.__file__).resolve().parents[3]  # .../pmultiqc
# Every column of a DIA-NN 2.5.1 report.parquet (PXD030304). A string literal in the
# code that equals one of these is a report column being read.
DIANN_REPORT_VOCABULARY = frozenset(
    (
        "Run.Index",
        "Run",
        "Channel",
        "Precursor.Id",
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Precursor.Lib.Index",
        "Decoy",
        "Proteotypic",
        "Precursor.Mz",
        "Protein.Ids",
        "Protein.Group",
        "Protein.Names",
        "Genes",
        "RT",
        "iRT",
        "Predicted.RT",
        "Predicted.iRT",
        "IM",
        "iIM",
        "Predicted.IM",
        "Predicted.iIM",
        "Precursor.Quantity",
        "Precursor.Normalised",
        "Ms1.Area",
        "Ms1.Normalised",
        "Ms1.Apex.Area",
        "Ms1.Apex.Mz.Delta",
        "Normalisation.Factor",
        "Quantity.Quality",
        "Empirical.Quality",
        "Normalisation.Noise",
        "Ms1.Profile.Corr",
        "Evidence",
        "Mass.Evidence",
        "Channel.Evidence",
        "Ms1.Total.Signal.Before",
        "Ms1.Total.Signal.After",
        "RT.Start",
        "RT.Stop",
        "FWHM",
        "PG.TopN",
        "PG.MaxLFQ",
        "Genes.TopN",
        "Genes.MaxLFQ",
        "Genes.MaxLFQ.Unique",
        "PG.MaxLFQ.Quality",
        "Genes.MaxLFQ.Quality",
        "Genes.MaxLFQ.Unique.Quality",
        "Q.Value",
        "PEP",
        "Global.Q.Value",
        "Lib.Q.Value",
        "Peptidoform.Q.Value",
        "Global.Peptidoform.Q.Value",
        "Lib.Peptidoform.Q.Value",
        "PTM.Site.Confidence",
        "Site.Occupancy.Probabilities",
        "Protein.Sites",
        "Lib.PTM.Site.Confidence",
        "Translated.Q.Value",
        "Channel.Q.Value",
        "PG.Q.Value",
        "PG.PEP",
        "GG.Q.Value",
        "Protein.Q.Value",
        "Global.PG.Q.Value",
        "Lib.PG.Q.Value",
        "Best.Fr.Mz",
        "Best.Fr.Mz.Delta",
    )
)


def _report(tmp_path, extra_cols=True, drop=()):
    data = {
        "Run": ["r1", "r1", "r2"],
        "Modified.Sequence": ["PEPT(SILAC)IDEK", "PEPTIDER", "PEPT(SILAC)IDEK"],
        "Stripped.Sequence": ["PEPTIDEK", "PEPTIDER", "PEPTIDEK"],
        "Protein.Group": ["P1", "P2", "P1"],
        "Protein.Names": ["N1", "N2", "N1"],
        "Precursor.Quantity": [1.0, 2.0, 3.0],
        "Precursor.Charge": [2, 2, 3],
        "RT": [10.0, 20.0, 11.0],
        "Q.Value": [0.001, 0.002, 0.003],
    }
    if extra_cols:
        data["Ms1.Total.Signal.Before"] = [0.1, 0.2, 0.3]
        data["Evidence"] = [5.0, 6.0, 7.0]
    for c in drop:
        data.pop(c)
    path = tmp_path / "diann_report.parquet"
    pq.write_table(pa.table(data), path)
    return path


def test_reader_projects_to_used_columns(tmp_path):
    reader = diann_reader.DiannReader(_report(tmp_path))
    reader.parse()
    df = reader.report_data
    assert "Evidence" not in df.columns
    assert "Ms1.Total.Signal.Before" not in df.columns
    assert set(df.columns) <= set(diann_reader.DIANN_REPORT_COLUMNS)
    assert len(df) == 3


def test_reader_strings_are_categorical_and_sanitised(tmp_path):
    reader = diann_reader.DiannReader(_report(tmp_path))
    reader.parse()
    df = reader.report_data
    for c in diann_reader.DIANN_REPORT_STRING_COLUMNS:
        assert isinstance(df[c].dtype, pd.CategoricalDtype), c
    assert list(df["Modified.Sequence"]) == ["PEPTIDEK", "PEPTIDER", "PEPTIDEK"]
    assert "(SILAC)" not in "".join(df["Modified.Sequence"].cat.categories)


def test_reader_tolerates_missing_optional_columns(tmp_path):
    reader = diann_reader.DiannReader(_report(tmp_path, drop=("Q.Value", "Protein.Names")))
    reader.parse()
    assert "Q.Value" not in reader.report_data.columns
    assert len(reader.report_data) == 3


def test_diann_report_columns_are_declared():
    """Every DIA-NN report column referenced in pmultiqc must be in the allowlist.

    The reader projects to DIANN_REPORT_COLUMNS, so a column read anywhere else
    but not declared here is silently absent at runtime (the first version of
    this test matched dotted names only and missed "Decoy": decoy filtering
    stopped without any error). Add the column to the list rather than widening
    the read.
    """
    declared = set(diann_reader.DIANN_REPORT_COLUMNS)
    used = set()
    for py in (PKG / "modules").rglob("*.py"):
        for node in ast.walk(ast.parse(py.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value in DIANN_REPORT_VOCABULARY:
                used.add(node.value)
    undeclared = sorted(used - declared)
    assert not undeclared, f"DIA-NN columns used but not declared in DIANN_REPORT_COLUMNS: {undeclared}"


def test_statistics_work_on_categorical_columns():
    """The reader yields Categoricals; the statistics must not hash lists (#717 follow-up)."""
    import numpy as np
    from pmultiqc.modules.common import dia_utils

    rng = np.random.default_rng(1)
    n = 500
    df = pd.DataFrame({
        "Run": pd.Categorical(rng.choice(["r1", "r2", "r3"], n)),
        "Modified.Sequence": pd.Categorical(rng.choice([f"PEP{i}K" for i in range(40)], n)),
        "Protein.Group": pd.Categorical(rng.choice([f"P{i}" for i in range(8)], n)),
        "Precursor.Quantity": rng.random(n),
        "Precursor.Normalised": rng.random(n),
        "Q.Value": rng.random(n) * 0.01,
    })
    total_protein, total_peptide, pep_plot = dia_utils._process_diann_statistics(df)
    assert total_peptide == df["Modified.Sequence"].nunique()
    assert pep_plot is not None


def test_reader_drops_decoys_while_reading(tmp_path):
    data = {
        "Run": ["r1", "r1", "r2", "r2"],
        "Modified.Sequence": ["PEPTIDEK", "DECOYK", "PEPTIDER", "DECOYR"],
        "Protein.Group": ["P1", "rev_P1", "P2", "rev_P2"],
        "Precursor.Quantity": [1.0, 2.0, 3.0, 4.0],
        "Decoy": [0, 1, 0, 1],
    }
    path = tmp_path / "diann_report.parquet"
    pq.write_table(pa.table(data), path)
    reader = diann_reader.DiannReader(path)
    reader.parse()
    df = reader.report_data
    assert len(df) == 2 and set(df["Modified.Sequence"]) == {"PEPTIDEK", "PEPTIDER"}
    assert "Decoy" not in df.columns, "Decoy must not be carried into pandas"


def test_modifications_computed_per_category():
    from pmultiqc.modules.common import dia_utils

    seqs = pd.Categorical(["PEP(UniMod:4)TIDEK", "PEPTIDEK", "PEP(UniMod:4)TIDEK", "PEP(UniMod:35)K"])
    df = pd.DataFrame({"Modified.Sequence": seqs})
    assert dia_utils._process_modifications(df) is True
    mods = df["Modifications"]
    assert isinstance(mods.dtype, pd.CategoricalDtype)
    assert list(mods.astype(str))[1] == "Unmodified"
    assert list(mods.astype(str))[0] == list(mods.astype(str))[2] != "Unmodified"


def test_identification_counts_match_the_set_based_logic():
    """Per-run/per-sample counts must equal what the old per-run sets produced (#717)."""
    import numpy as np
    from pmultiqc.modules.common import dia_utils
    from pmultiqc.modules.common.common_utils import cal_num_table_at_sample

    rng = np.random.default_rng(3)
    n = 3000
    runs = [f"run{i}" for i in range(6)]
    df = pd.DataFrame({
        "Run": pd.Categorical(rng.choice(runs, n)),
        "Modified.Sequence": pd.Categorical(rng.choice([f"PEP{i}K" for i in range(150)] + [f"PEP{i}(UniMod:4)K" for i in range(50)], n)),
        "Protein.Group": pd.Categorical(rng.choice([f"P{i}" for i in range(30)], n)),
        "Proteotypic": rng.integers(0, 2, n),
    })
    file_df = pd.DataFrame({"Run": runs, "Sample": [1, 1, 2, 2, 3, 3]})
    old_run, data_per_run = {}, {}
    for run, group in df.groupby("Run", observed=True):
        old_run[str(run)], data_per_run[str(run)] = dia_utils._calculate_run_statistics(group)
    old_sample = cal_num_table_at_sample(file_df, data_per_run)
    assert dia_utils._run_identification_counts(df) == old_run
    assert dia_utils._sample_identification_counts(df, file_df) == old_sample


def test_sample_level_modifications_without_merge_matches_merge_reference():
    import numpy as np
    from pmultiqc.modules.common import dia_utils
    from pmultiqc.modules.common.common_utils import summarize_modifications

    rng = np.random.default_rng(5)
    n = 4000
    runs = [f"run{i}" for i in range(6)]
    df = pd.DataFrame({
        "Run": pd.Categorical(rng.choice(runs, n)),
        "Modified.Sequence": pd.Categorical(rng.choice([f"PEP{i}K" for i in range(80)], n)),
        "Modifications": pd.Categorical(rng.choice(["Unmodified", "Oxidation", "Carbamidomethyl,Oxidation"], n)),
        "Protein.Group": pd.Categorical(rng.choice([f"P{i}" for i in range(10)], n)),
    })
    file_df = pd.DataFrame({"Run": runs[:5], "Sample": [1, 1, 2, 2, 3]})  # run5 unmapped
    ref = df.merge(file_df[["Sample", "Run"]].drop_duplicates(), on="Run")
    ref["Sample"] = ref["Sample"].astype(int)
    expected = {f"Sample {s}": summarize_modifications(g.drop_duplicates())[0] for s, g in ref.groupby("Sample", sort=True)}
    got = dia_utils.dia_sample_level_modifications(df, file_df)
    assert got.keys() == expected.keys()
    for k in expected:
        assert got[k].keys() == expected[k].keys()
        for m in expected[k]:
            assert np.isclose(got[k][m], expected[k][m]), (k, m)
    codes = dia_utils.run_to_sample_codes(df["Run"], file_df)
    assert codes.isna().sum() == (df["Run"].astype(str) == "run5").sum()


def test_dia_utils_imports_cleanly_first():
    """multiqc imports dia_utils before plots.dia; a cycle there skips the whole module (#717)."""
    import subprocess, sys
    code = "import pmultiqc.modules.common.dia_utils as d; import pmultiqc.modules.common.plots.dia; print(callable(d.run_to_sample_codes))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0 and out.stdout.strip() == "True", out.stderr[-800:]
