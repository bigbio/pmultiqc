from __future__ import annotations

from pathlib import Path
import os
import pandas as pd
from datetime import datetime

from pmultiqc.modules.common.ms.base import BaseParser
from pmultiqc.modules.common.logging import get_logger


class DiannReader(BaseParser):
    def __init__(
        self,
        file_path: Path | str,
    ) -> None:

        super().__init__([file_path])

        self.file_path = file_path

        # Outputs populated by parse()
        self.report_data: pd.DataFrame = pd.DataFrame()

        self.log = get_logger("pmultiqc.modules.common.ms.diann")

    def parse(self, **_kwargs) -> None:

        # parse DIA-NN report data
        if os.path.splitext(self.file_path)[1] == ".tsv":
            report_data = pd.read_csv(
                self.file_path, header=0, sep="\t", on_bad_lines="warn"
            )

        else:
            report_data = _read_report_parquet(self.file_path, self.log)
            report_data["Modified.Sequence"] = _sanitize_column(report_data["Modified.Sequence"])

        self.log.info(
            "{}: Done parsing DIANN file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.file_path
            )
        )

        self.report_data = report_data

        return None

def _sanitize_sequence(seq):
    seq = seq.replace("(SILAC)", "")
    return seq


def _sanitize_column(col: pd.Series) -> pd.Series:
    """Apply ``_sanitize_sequence`` without densifying a categorical column.

    On a 231 M-row report ``Series.apply`` on an object column is one Python call
    per row and, on a categorical, would first expand it back to objects.
    Rewriting the categories touches each distinct sequence once instead.
    """
    if isinstance(col.dtype, pd.CategoricalDtype):
        return col.cat.rename_categories([_sanitize_sequence(c) for c in col.cat.categories])
    return col.apply(_sanitize_sequence)


# Every DIA-NN report column pmultiqc reads. Keep in step with the code: the
# test ``test_diann_report_columns_are_declared`` scans pmultiqc/modules for
# DIA-NN column literals and fails if one is used that is not listed here.
DIANN_REPORT_COLUMNS = (
    "Run",
    "RT",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "Precursor.Charge",
    "Protein.Group",
    "Protein.Names",
    "Q.Value",
    "iRT",
    "Predicted.RT",
    "RT.Start",
    "RT.Stop",
    "FWHM",
    "Normalisation.Factor",
    "Ms1.Area",
    "Ms1.Apex.Mz.Delta",
    # Undotted names are easy to miss; the test checks against the full
    # DIA-NN column vocabulary, not a dotted-name regex.
    "Decoy",
    "Proteotypic",
)

# String columns; read dictionary-encoded so pandas gets a Categorical instead of
# one Python object per row.
DIANN_REPORT_STRING_COLUMNS = (
    "Run",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Protein.Group",
    "Protein.Names",
)


def _read_report_parquet(path, log) -> pd.DataFrame:
    """Read only the columns pmultiqc uses, with strings as categoricals.

    A DIA-NN report has ~70 columns; pmultiqc reads 18. Reading the whole table
    into pandas materialises every column, with the five string columns as
    Python objects, which is why the summary step needed 234 GB on a 231 M-row
    report (bigbio/pmultiqc#717). Projecting to the used columns and reading the
    string columns dictionary-encoded brought the measured load to 52 GB.

    Columns missing from the file are simply not requested, so older or
    differently configured DIA-NN reports still load; downstream code already
    guards the optional ones.
    """
    import pyarrow.parquet as pq

    # pyarrow raises KeyError for a read_dictionary column the file lacks, so
    # look at the schema first (footer only, cheap) and request what exists.
    present = set(pq.ParquetFile(path).schema_arrow.names)
    dictionary = [c for c in DIANN_REPORT_STRING_COLUMNS if c in present]
    columns = [c for c in DIANN_REPORT_COLUMNS if c in present]
    missing = [c for c in DIANN_REPORT_COLUMNS if c not in present]
    if missing:
        log.debug("DIA-NN report lacks %s; reading without them.", ", ".join(missing))
    # Drop decoys while reading, in Arrow. Doing it in pandas afterwards costs two
    # more copies of the whole frame (the boolean index and the .copy()) at the
    # point where memory is already highest; on PXD030304 that alone crossed the
    # 72 GB budget. The Decoy column is then not carried into the frame, so the
    # pandas-side guard in _load_and_preprocess_diann_data does nothing.
    filters = [("Decoy", "==", 0)] if "Decoy" in present else None
    columns = [c for c in columns if c != "Decoy"]
    log.info("Reading %d of %d columns from %s%s", len(columns), len(present), path,
             " (decoys filtered while reading)" if filters else "")
    table = pq.read_table(path, columns=columns, filters=filters, read_dictionary=dictionary)
    # self_destruct releases each Arrow buffer as its column is converted, and
    # split_blocks avoids consolidating into one big block first. Without them
    # ~33 GB of Arrow memory stayed resident next to the 16.8 GB pandas frame
    # on PXD030304 even after the table was deleted (measured: 49.9 -> 25.2 GiB
    # after conversion), and that baseline is what pushed the first plot stage
    # over the pipeline's 72 GB (bigbio/pmultiqc#717).
    return table.to_pandas(self_destruct=True, split_blocks=True)
