

"""
MzML file reading functionality
"""

from __future__ import absolute_import
import logging
import pandas as pd
from datetime import datetime
from pyopenms import MzMLFile, MSExperiment
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.ms.msreader import MSReader

# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

class MzMLReader(MSReader):
    """Class for reading and processing mzML files"""
    def read(
        self,
        file_paths=None,
        ms_with_psm=None,
        identified_spectrum=None,
        mzml_charge_plot=None,
        mzml_peak_distribution_plot=None,
        mzml_peaks_ms2_plot=None,
        mzml_charge_plot_1=None,
        mzml_peak_distribution_plot_1=None,
        mzml_peaks_ms2_plot_1=None,
        ms_without_psm=None,
        enable_dia=False,
        enable_mzid=False,
    ):
        """
        Read mzML files and extract information

        Args:
            file_paths: List of paths to mzML files (optional, uses self.file_paths if not provided)
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

        Returns:
            tuple: (mzml_table, heatmap_charge, total_ms2_spectra)
        """
        # Use provided file_paths or fall back to instance attribute
        ms_paths = file_paths or self.file_paths
        mzml_table = {}
        heatmap_charge = {}
        # Reset spectrum counts for this reading session
        self._reset_spectrum_counts()

        mzml_ms_dicts = list()
        for m in ms_paths:
            file_ms1_number = 0
            file_ms2_number = 0
            log.info("{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
            exp = MSExperiment()
            MzMLFile().load(m, exp)
            log.info("{}: Done parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
            m_name = file_prefix(m)
            log.info(
                "{}: Aggregating mzML file {}...".format(datetime.now().strftime("%H:%M:%S"), m_name)
            )

            charge_2 = 0
            for i in exp:
                if i.getMSLevel() == 1:
                    file_ms1_number += 1
                    self._increment_ms1_count()
                elif i.getMSLevel() == 2:
                    file_ms2_number += 1
                    self._increment_ms2_count()
                    charge_state = i.getPrecursors()[0].getCharge()
                    peaks_tuple = i.get_peaks()

                    if enable_mzid:
                        # retention_time: minute
                        mzml_ms_dicts.append(
                            {
                                "spectrumID": i.getNativeID(),
                                "intensity": float(peaks_tuple[1].max()),
                                "retention_time": i.getRT() / 60,
                                "filename": m_name,
                            }
                        )

                    peak_per_ms2 = len(peaks_tuple[0])
                    if i.getMetaValue("base peak intensity"):
                        base_peak_intensity = i.getMetaValue("base peak intensity")
                    else:
                        base_peak_intensity = max(peaks_tuple[1]) if len(peaks_tuple[1]) > 0 else None

                    if charge_state == 2:
                        charge_2 += 1

                    if enable_dia:
                        mzml_charge_plot.add_value(charge_state)
                        mzml_peak_distribution_plot.add_value(base_peak_intensity)
                        mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                        continue

                    if m_name in ms_with_psm:
                        if i.getNativeID() in identified_spectrum[m_name]:
                            mzml_charge_plot.add_value(charge_state)
                            mzml_peak_distribution_plot.add_value(base_peak_intensity)
                            mzml_peaks_ms2_plot.add_value(peak_per_ms2)
                        else:
                            mzml_charge_plot_1.add_value(charge_state)
                            mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                            mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
                    else:
                        if m_name not in ms_without_psm:
                            ms_without_psm.append(m_name)

            heatmap_charge[m_name] = charge_2 / file_ms2_number if file_ms2_number > 0 else 0
            mzml_table[m_name] = {"MS1_Num": file_ms1_number}
            mzml_table[m_name]["MS2_Num"] = file_ms2_number
            log.info(
                "{}: Done aggregating mzML file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m_name
                )
            )

        if enable_mzid:
            return mzml_table, heatmap_charge, self.total_ms2_number, pd.DataFrame(mzml_ms_dicts)
        else:
            return mzml_table, heatmap_charge, self.total_ms2_number