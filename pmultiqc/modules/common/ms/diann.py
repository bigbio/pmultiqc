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
            report_data = pd.read_parquet(self.file_path)

        self.log.info(
            "{}: Done parsing DIANN file {}...".format(
                datetime.now().strftime("%H:%M:%S"), self.file_path
            )
        )

        self.report_data = report_data

        return None