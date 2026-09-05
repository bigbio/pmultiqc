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
