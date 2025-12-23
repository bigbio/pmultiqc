from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.ms.base import BaseParser
from pmultiqc.modules.common.ms_io import (
    get_ms_qc_info,
    add_ms_values_df,
    spectra_ref_check,
)


class MsInfoReader(BaseParser):
    def __init__(
        self,
        file_paths: list[str | Path],
        ms_with_psm,
        identified_spectrum,
        mzml_charge_plot,
        mzml_peak_distribution_plot,
        mzml_peaks_ms2_plot,
        mzml_charge_plot_1,
        mzml_peak_distribution_plot_1,
        mzml_peaks_ms2_plot_1,
        ms_without_psm,
        enable_dia: bool = False,
    ) -> None:
        super().__init__(file_paths)
        self.ms_with_psm = ms_with_psm
        self.identified_spectrum = identified_spectrum
        self.mzml_charge_plot = mzml_charge_plot
        self.mzml_peak_distribution_plot = mzml_peak_distribution_plot
        self.mzml_peaks_ms2_plot = mzml_peaks_ms2_plot
        self.mzml_charge_plot_1 = mzml_charge_plot_1
        self.mzml_peak_distribution_plot_1 = mzml_peak_distribution_plot_1
        self.mzml_peaks_ms2_plot_1 = mzml_peaks_ms2_plot_1
        self.ms_without_psm = ms_without_psm
        self.enable_dia = enable_dia

        # Outputs populated by parse()
        self.mzml_table: dict = {}
        self.heatmap_charge: dict = {}
        self.total_ms2_spectra: int = 0
        self.ms1_tic: dict = {}
        self.ms1_bpc: dict = {}
        self.ms1_peaks: dict = {}
        self.ms1_general_stats: dict = {}

        self.log = get_logger("pmultiqc.modules.common.ms.msinfo")

    def parse(self, **_kwargs) -> None:
        mzml_table = {}
        heatmap_charge = {}
        total_ms2_spectra = 0
        ms1_tic = {}
        ms1_bpc = {}
        ms1_peaks = {}
        ms1_general_stats = {}

        for file in self.file_paths:
            self.log.info(
                "{}: Parsing ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )
            mzml_df = pd.read_parquet(file)
            m_name = file_prefix(file).replace("_ms_info", "")
            if m_name not in mzml_table:
                mzml_table[m_name] = dict.fromkeys(["MS1_Num", "MS2_Num", "Charge_2"], 0)

            charge_group = mzml_df.groupby("precursor_charge").size()
            ms_level_group = mzml_df.groupby("ms_level").size()
            charge_2 = charge_group[2] if 2 in charge_group else 0
            ms1_number = int(ms_level_group[1]) if 1 in ms_level_group else 0
            ms2_number = int(ms_level_group[2]) if 2 in ms_level_group else 0
            total_ms2_spectra += ms2_number
            mzml_table[m_name].update({"MS1_Num": mzml_table[m_name]["MS1_Num"] + ms1_number})
            mzml_table[m_name].update({"MS2_Num": mzml_table[m_name]["MS2_Num"] + ms2_number})
            mzml_table[m_name].update({"Charge_2": mzml_table[m_name]["Charge_2"] + charge_2})

            (
                ms1_tic[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_bpc[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_peaks[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_general_stats[os.path.basename(file).replace("_ms_info.parquet", "")],
            ) = get_ms_qc_info(mzml_df)

            group = mzml_df[mzml_df["ms_level"] == 2]
            del mzml_df

            if self.enable_dia:
                identified_spectrum_scan_id = []
            else:
                if not self.identified_spectrum:
                    raise ValueError(
                        "ms_io: The identified_spectrum is missing. Please check your mzTab file!"
                    )

                if m_name not in self.identified_spectrum:
                    raise ValueError(
                        f"identified_spectrum missing entries for '{m_name}'. Check your mzTab file."
                    )
                identified_spectrum_scan_id = [
                    spectra_ref_check(spectrum_id)
                    for spectrum_id in self.identified_spectrum[m_name]
                ]

            add_ms_values_df(
                group,
                m_name,
                self.ms_with_psm,
                identified_spectrum_scan_id,
                self.mzml_charge_plot,
                self.mzml_peak_distribution_plot,
                self.mzml_peaks_ms2_plot,
                self.mzml_charge_plot_1,
                self.mzml_peak_distribution_plot_1,
                self.mzml_peaks_ms2_plot_1,
                self.ms_without_psm,
                self.enable_dia,
            )

            for m in mzml_table.keys():
                if mzml_table[m]["MS2_Num"] > 0:
                    heatmap_charge[m] = mzml_table[m]["Charge_2"] / mzml_table[m]["MS2_Num"]
                else:
                    heatmap_charge[m] = 0

            self.log.info(
                "{}: Done aggregating ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )

        self.mzml_table = mzml_table
        self.heatmap_charge = heatmap_charge
        self.total_ms2_spectra = total_ms2_spectra

        self.ms1_tic = {
            k: v for k, v in ms1_tic.items() if v is not None
        }
        self.ms1_bpc = {
            k: v for k, v in ms1_bpc.items() if v is not None
        }
        self.ms1_peaks = {
            k: v for k, v in ms1_peaks.items() if v is not None
        }
        self.ms1_general_stats = {
            k: v for k, v in ms1_general_stats.items() if v is not None
        }

        return None
