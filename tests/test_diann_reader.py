"""The DIA-NN parquet reader must not materialise the whole report (bigbio/pmultiqc#717)."""

import ast
import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from pmultiqc.modules.common.ms import diann as diann_reader

PKG = Path(diann_reader.__file__).resolve().parents[3]  # .../pmultiqc
DIANN_LITERAL = re.compile(r"^(Run|RT|FWHM|iRT|[A-Z][A-Za-z0-9]*(\.[A-Za-z0-9]+)+)$")


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
    """Every DIA-NN column literal used anywhere in pmultiqc must be in the allowlist.

    If this fails, a new column is being read from the report: add it to
    DIANN_REPORT_COLUMNS (and DIANN_REPORT_STRING_COLUMNS if it is a string) rather
    than widening the read.
    """
    declared = set(diann_reader.DIANN_REPORT_COLUMNS)
    used = set()
    for py in (PKG / "modules").rglob("*.py"):
        for node in ast.walk(ast.parse(py.read_text())):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and DIANN_LITERAL.match(node.value):
                used.add(node.value)
    # Literals that look like DIA-NN columns but are not report columns pmultiqc reads.
    not_report_columns = {c for c in used if c not in declared}
    unexpected = sorted(c for c in not_report_columns if c.count(".") <= 3 and c not in {"Q.Value"} | declared)
    assert not unexpected, f"DIA-NN columns used but not declared in DIANN_REPORT_COLUMNS: {unexpected}"
