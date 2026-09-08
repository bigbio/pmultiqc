"""The MSstats input parse only feeds the quantification tables; skip it when tables are disabled (#717)."""

import ast
from pathlib import Path

import pmultiqc.modules.quantms.quantms as qm


def test_msstats_parse_is_gated_on_disable_table():
    src = Path(qm.__file__).read_text(encoding="utf-8")
    tree = ast.parse(src)
    gate_found = False
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            cond = ast.unparse(node.test)
            body = ast.unparse(node)
            if "msstats_input_valid" in cond and "parse_msstats_input()" in body:
                assert "disable_table" in cond, f"MSstats parse not gated on disable_table: {cond}"
                gate_found = True
    assert gate_found, "could not find the msstats parse call site"
