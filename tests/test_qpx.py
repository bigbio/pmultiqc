"""Regression tests for the QPX module and for plugin option wiring.

These cover the defects found while reviewing the QPX pull request:
  * ``--qpx-plugin`` was bound to the ``mhcquant_plugin`` variable, silently
    dropping ``--mhcquant-plugin``.
  * ``quantms.py`` still passed ``quantms_modified=`` / ``quantms_missed_cleavages=``
    after those parameters were renamed.
  * ``#Identified MS2 Spectra`` grouped on the first *character* of the scan id.
  * ``_get_enzyme_name`` indexed ``[0]`` of a cell that may be a plain string.
  * Only the first parquet file per QPX table was read.
  * The SDRF subset was merged without de-duplication.
"""

import ast
import re
from pathlib import Path

import pandas as pd
import pytest

from pmultiqc.modules.qpx.qpx_utils import _get_enzyme_name, get_pep_intensity


REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_DIR = REPO_ROOT / "pmultiqc"

# Wrappers that take the real callee as their first positional argument and forward
# **kwargs to it, mapped to the keyword names they consume themselves. These swallow
# exceptions, so a stale keyword here fails silently at runtime.
FORWARDING_WRAPPERS = {"_safe_draw": {"name"}}


class TestEnzymeName:
    """``_get_enzyme_name`` must not slice a scalar string into a single letter."""

    def test_scalar_string_cell(self):
        df = pd.DataFrame({"enzymes": ["Lys-C"]})
        assert _get_enzyme_name(df) == "Lys-C"

    def test_list_cell(self):
        df = pd.DataFrame({"enzymes": [["Lys-C", "Trypsin"]]})
        assert _get_enzyme_name(df) == "Lys-C"

    def test_struct_cell(self):
        df = pd.DataFrame({"enzymes": [[{"name": "Chymotrypsin"}]]})
        assert _get_enzyme_name(df) == "Chymotrypsin"

    def test_defaults_to_trypsin(self):
        assert _get_enzyme_name(None) == "Trypsin"
        assert _get_enzyme_name(pd.DataFrame()) == "Trypsin"
        assert _get_enzyme_name(pd.DataFrame({"other": [1]})) == "Trypsin"
        assert _get_enzyme_name(pd.DataFrame({"enzymes": [None]})) == "Trypsin"

    def test_empty_collection_cell(self):
        from pmultiqc.modules.qpx.qpx_utils import _first_enzyme_value

        assert _first_enzyme_value([]) is None
        assert _first_enzyme_value("") is None


class TestPepIntensitySdrfMerge:
    """Duplicated SDRF rows must not inflate the run-level intensity distribution."""

    @staticmethod
    def _pep_table():
        return pd.DataFrame(
            {
                "run": ["run1", "run1", "run1"],
                "intensities": [
                    [{"label": "TMT126", "intensity": 10.0}],
                    [{"label": "TMT126", "intensity": 20.0}],
                    [{"label": "TMT126", "intensity": 30.0}],
                ],
            }
        )

    def test_duplicate_sdrf_rows_do_not_duplicate_intensities(self):
        unique_sdrf = pd.DataFrame(
            {"Label": ["1"], "Run": ["run1"], "Sample": [1], "Extra": ["a"]}
        )
        # Same Label/Run/Sample repeated, as happens when the SDRF has several
        # rows per run (e.g. one per fraction or per technical replicate).
        duplicated_sdrf = pd.DataFrame(
            {
                "Label": ["1", "1", "1"],
                "Run": ["run1", "run1", "run1"],
                "Sample": [1, 1, 1],
                "Extra": ["a", "b", "c"],
            }
        )

        by_run_unique, _ = get_pep_intensity(self._pep_table(), unique_sdrf)
        by_run_dup, _ = get_pep_intensity(self._pep_table(), duplicated_sdrf)

        assert by_run_unique == by_run_dup


class TestIdentifiedSpectraCount:
    """The identified-MS2 count is over distinct (run, scan) pairs.

    ``scan`` is a ``list<int32>`` column in psm.parquet, so the count unwraps it with
    ``.str[0]``. The raw column holds numpy arrays and is unhashable, so grouping on
    the column directly raises TypeError -- this pins that down.
    """

    @staticmethod
    def _psm_df():
        import numpy as np

        return pd.DataFrame(
            {
                "run": ["run1", "run1", "run1", "run2"],
                "scan": [
                    np.array([1], dtype="int32"),
                    np.array([2], dtype="int32"),
                    np.array([2], dtype="int32"),
                    np.array([1], dtype="int32"),
                ],
                "sequence": ["PEPTIDEA", "PEPTIDEB", "PEPTIDEB", "PEPTIDEA"],
            }
        )

    def test_counts_distinct_run_scan_pairs(self):
        psm_df = self._psm_df()
        # run1 has scans {1, 2} (2 rows share scan 2), run2 has scan {1} -> 3 pairs.
        assert psm_df.groupby(["run", psm_df["scan"].str[0]]).ngroups == 3

    def test_scan_column_is_not_directly_groupable(self):
        psm_df = self._psm_df()
        with pytest.raises(TypeError):
            psm_df.groupby(["run", "scan"]).ngroups


class TestQpxParquetDiscovery:
    """Every matching parquet file is read, not just the first one."""

    class _StubModule:
        """Minimal stand-in exposing only what ``_read_qpx_parquets`` uses."""

        def __init__(self, paths):
            self._paths = paths

        def find_log_files(self, search_pattern, filecontents=False):
            for path in self._paths:
                yield {"root": str(path.parent), "fn": path.name}

    def _read(self, paths):
        from pmultiqc.modules.qpx.qpx import QpxModule

        stub = self._StubModule(paths)
        return QpxModule._read_qpx_parquets(stub, "pmultiqc/qpx_run", None)

    def test_concatenates_multiple_files(self, tmp_path):
        paths = []
        for i in range(3):
            path = tmp_path / f"part{i}.parquet"
            pd.DataFrame({"run": [f"run{i}"], "value": [i]}).to_parquet(path)
            paths.append(path)

        result = self._read(paths)

        assert len(result) == 3
        assert sorted(result["run"]) == ["run0", "run1", "run2"]

    def test_single_file(self, tmp_path):
        path = tmp_path / "only.parquet"
        pd.DataFrame({"run": ["run0"], "value": [0]}).to_parquet(path)

        assert len(self._read([path])) == 1

    def test_no_files_returns_none(self):
        assert self._read([]) is None


def _click_option_params():
    """Map each ``pmultiqc.cli`` module attribute to the click parameter it creates."""
    import click

    from pmultiqc import cli

    params = {}
    for attr_name, value in vars(cli).items():
        if attr_name.startswith("_") or not callable(value):
            continue

        def _target():
            pass

        try:
            value(_target)
        except (TypeError, AttributeError):
            continue

        click_params = getattr(_target, "__click_params__", [])
        if len(click_params) == 1 and isinstance(click_params[0], click.Option):
            params[attr_name] = click_params[0]

    return params


class TestCliPluginOptions:
    """Each CLI option must be reachable: right variable name, and an entry point."""

    def test_option_variable_matches_parameter_name(self):
        mismatched = {
            attr: param.name
            for attr, param in _click_option_params().items()
            if attr != param.name
        }
        assert not mismatched, (
            f"CLI option assigned to a variable that does not match its parameter "
            f"name (the later assignment silently clobbers the earlier one): {mismatched}"
        )

    def test_expected_plugin_flags_exist(self):
        params = _click_option_params()
        for plugin in [
            "mzid",
            "proteobench",
            "maxquant",
            "diann",
            "quantms",
            "fragpipe",
            "mhcquant",
            "qpx",
        ]:
            attr = f"{plugin}_plugin"
            assert attr in params, f"missing CLI option variable {attr}"
            assert f"--{plugin.replace('_', '-')}-plugin" in params[attr].opts

    def test_plugin_options_are_registered_as_entry_points(self):
        pyproject = REPO_ROOT / "pyproject.toml"
        if not pyproject.exists():
            pytest.skip("pyproject.toml not available (installed package)")

        section = re.search(
            r'\[project\.entry-points\."multiqc\.cli_options\.v1"\](.*?)(?=\n\[|\Z)',
            pyproject.read_text(encoding="utf-8"),
            re.DOTALL,
        )
        assert section, "cli_options entry-point section not found"
        registered = set(re.findall(r"^(\w+)\s*=", section.group(1), re.MULTILINE))

        plugin_options = {a for a in _click_option_params() if a.endswith("_plugin")}
        missing = plugin_options - registered
        assert not missing, (
            f"plugin CLI options with no multiqc.cli_options.v1 entry point "
            f"(the flag would not be exposed by multiqc): {sorted(missing)}"
        )


class TestInternalCallSiteKeywords:
    """Guard against call sites left behind by a parameter rename."""

    def test_no_unknown_keyword_arguments(self):
        trees = {}
        for path in sorted(PACKAGE_DIR.rglob("*.py")):
            # Explicit encoding: the runner's default may be ASCII, and several
            # sources contain non-ASCII characters.
            trees[path] = ast.parse(path.read_text(encoding="utf-8"))

        # Collect unambiguous, **kwargs-free function definitions.
        definitions = {}
        for tree in trees.values():
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    args = node.args
                    params = {
                        a.arg for a in args.posonlyargs + args.args + args.kwonlyargs
                    }
                    definitions.setdefault(node.name, []).append(
                        (params, args.kwarg is not None)
                    )

        def resolve(node):
            """Return (callee_name, keywords_to_check) for a call we can attribute.

            Handles both plain ``helper(...)`` calls and ``self._safe_draw(helper, ...)``
            forwarding wrappers. Only ``ast.Name`` callees are considered, so pandas or
            stdlib methods sharing a name with a local helper (e.g. ``.to_dict``) are
            never misattributed.
            """
            if isinstance(node.func, ast.Name):
                return node.func.id, node.keywords
            if (
                isinstance(node.func, ast.Attribute)
                and node.func.attr in FORWARDING_WRAPPERS
                and node.args
                and isinstance(node.args[0], ast.Name)
            ):
                skip = FORWARDING_WRAPPERS[node.func.attr]
                return node.args[0].id, [k for k in node.keywords if k.arg not in skip]
            return None, []

        problems = []
        for path, tree in trees.items():
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call) or not node.keywords:
                    continue
                callee, keywords = resolve(node)
                defs = definitions.get(callee)
                if not defs or any(has_kwargs for _, has_kwargs in defs):
                    continue
                # Several modules define same-named drawers (e.g. draw_peptides_table in
                # both dia_utils and mzidentml_utils). That is still checkable as long as
                # every definition accepts the same parameters -- only genuinely
                # divergent signatures are ambiguous enough to skip.
                param_sets = {frozenset(params) for params, _ in defs}
                if len(param_sets) > 1:
                    continue
                params = defs[0][0]
                for keyword in keywords:
                    if keyword.arg is not None and keyword.arg not in params:
                        problems.append(
                            f"{path.relative_to(REPO_ROOT)}:{node.lineno}: "
                            f"{callee}() has no parameter '{keyword.arg}' "
                            f"(accepts {sorted(params)})"
                        )

        assert not problems, "Unknown keyword arguments:\n" + "\n".join(problems)
