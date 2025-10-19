from __future__ import annotations

from pathlib import Path
from datetime import datetime
from pyopenms import MzMLFile, MSExperiment
import pandas as pd
import numpy as np

from pmultiqc.modules.common.ms.base import BaseParser
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.ms_io import get_ms_qc_info


class MzMLReader(BaseParser):
    def __init__(
        self,
        file_paths: list[str, Path],
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
        enable_mzid: bool = False
    ) -> None:

        """
        Read mzML files and extract information

        Args:
            ms_paths: List of paths to mzML files
            ms_with_psm: List of MS files with PSMs
            identified_spectrum: Dictionary of identified spectra by MS file
            mzml_charge_plot: Histogram for charge distribution of identified spectra
            mzml_peak_distribution_plot: Histogram for peak distribution of identified spectra
            mzml_peaks_ms2_plot: Histogram for peaks per MS2 of identified spectra
            mzml_charge_plot_1: Histogram for charge distribution of unidentified spectra
            mzml_peak_distribution_plot_1: Histogram for peak distribution of unidentified spectra
            mzml_peaks_ms2_plot_1: Histogram for peaks per MS2 of unidentified spectra
            ms_without_psm: List of MS files without PSMs
            enable_dia: Whether DIA mode is enabled
            enable_mzid: Whether mzid_plugin mode is enabled
        """

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
        self.enable_mzid = enable_mzid

        # Outputs populated by parse()
        self.mzml_table: dict = {}
        self.heatmap_charge: dict = {}
        self.total_ms2_spectra: int = 0
        self.ms1_tic: dict = {}
        self.ms1_bpc: dict = {}
        self.ms1_peaks: dict = {}
        self.ms1_general_stats: dict = {}
        self.mzml_ms_df: pd.DataFrame= pd.DataFrame()

        self.log = get_logger("pmultiqc.modules.common.ms.mzml")

    def parse(self, **kwargs) -> None:
        mzml_table = {}
        heatmap_charge = {}
        total_ms2_spectra = 0

        mzml_ms_dicts = list()
        ms1_tic = dict()
        ms1_bpc = dict()
        ms1_peaks = dict()
        ms1_general_stats = dict()

        for file_name in self.file_paths:
            ms1_number = 0
            ms2_number = 0

            self.log.info(
                "{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), file_name)
            )

            mzml_exp = MSExperiment()
            MzMLFile().load(file_name, mzml_exp)
            self.log.info(
                "{}: Done parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), file_name)
            )

            m_name = file_prefix(file_name)
            self.log.info(
                "{}: Aggregating mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m_name)
            )

            spectrums_data = list()

            if mzml_exp.getMetaValue("acquisition_date_time"):
                acquisition_datetime = str(mzml_exp.getMetaValue("acquisition_date_time"))
            elif mzml_exp.getDateTime().get():
                acquisition_datetime = str(mzml_exp.getDateTime().get())
            else:
                acquisition_datetime = ""

            charge_2 = 0
            for spectrum in mzml_exp:

                mz_array, intensity_array = spectrum.get_peaks()
                num_peaks = len(mz_array)

                if num_peaks == 0:
                    continue

                summed_peak_intensities = float(np.sum(intensity_array))

                ms_level = spectrum.getMSLevel()
                rt = spectrum.getRT()
                scan_id = spectrum.getNativeID()

                spectrums_data.append(
                    {
                        "scan_id": scan_id,
                        "ms_level": ms_level,
                        "summed_peak_intensities": summed_peak_intensities,
                        "rt": rt,
                        "num_peaks": num_peaks,
                        "acquisition_datetime": acquisition_datetime,
                    }
                )

                if ms_level == 1:
                    ms1_number += 1

                elif ms_level == 2:
                    ms2_number += 1

                    charge_state = spectrum.getPrecursors()[0].getCharge()

                    if self.enable_mzid:
                        # retention_time: minute
                        mzml_ms_dicts.append(
                            {
                                "spectrumID": scan_id,
                                "intensity": float(intensity_array.max()),
                                "retention_time": rt / 60,
                                "filename": m_name,
                            }
                        )

                    peak_per_ms2 = len(mz_array)
                    if spectrum.getMetaValue("base peak intensity"):
                        base_peak_intensity = spectrum.getMetaValue("base peak intensity")
                    else:
                        base_peak_intensity = max(intensity_array) if len(intensity_array) > 0 else None

                    if charge_state == 2:
                        charge_2 += 1

                    if self.enable_dia:
                        self.mzml_charge_plot.add_value(charge_state)
                        self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                        self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                        continue

                    if m_name in self.ms_with_psm:
                        if scan_id in self.identified_spectrum[m_name]:
                            self.mzml_charge_plot.add_value(charge_state)
                            self.mzml_peak_distribution_plot.add_value(base_peak_intensity)
                            self.mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                        else:
                            self.mzml_charge_plot_1.add_value(charge_state)
                            self.mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                            self.mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
                    else:
                        if m_name not in self.ms_without_psm:
                            self.ms_without_psm.append(m_name)

            heatmap_charge[m_name] = charge_2 / ms2_number if ms2_number > 0 else 0
            total_ms2_spectra += ms2_number
            mzml_table[m_name] = {"MS1_Num": ms1_number}
            mzml_table[m_name]["MS2_Num"] = ms2_number
            self.log.info(
                "{}: Done aggregating mzML file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m_name
                )
            )

            spectrums_df = pd.DataFrame(spectrums_data)

            (
                ms1_tic[m_name],
                ms1_bpc[m_name],
                ms1_peaks[m_name],
                ms1_general_stats[m_name],
            ) = get_ms_qc_info(spectrums_df)

        self.mzml_table = mzml_table
        self.heatmap_charge = heatmap_charge
        self.total_ms2_spectra = total_ms2_spectra
        self.mzml_ms_df = pd.DataFrame(mzml_ms_dicts)

        ms1_tic = {k: v for k, v in ms1_tic.items() if v is not None}
        ms1_bpc = {k: v for k, v in ms1_bpc.items() if v is not None}
        ms1_peaks = {k: v for k, v in ms1_peaks.items() if v is not None}
        ms1_general_stats = {k: v for k, v in ms1_general_stats.items() if v is not None}

        self.ms1_tic = ms1_tic
        self.ms1_bpc = ms1_bpc
        self.ms1_peaks = ms1_peaks
        self.ms1_general_stats = ms1_general_stats

        return None
