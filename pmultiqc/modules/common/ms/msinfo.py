"""
MS info file reading functionality
"""

from __future__ import absolute_import
import os
import logging
import pandas as pd
import re
from datetime import datetime

from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.utils import get_ms_qc_info
from pmultiqc.modules.common.ms.msreader import MSReader

# Initialise the main MultiQC logger
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class MSInfoReader(MSReader):
    """Class for reading and processing MS info files"""

    def _spectra_ref_check(self, spectra_ref):
        """
        Extract scan number from spectra reference string
        
        Args:
            spectra_ref: String containing spectra reference (e.g., "scan=123", "spectrum=456", or "123")
            
        Returns:
            str: Extracted scan number
            
        Raises:
            ValueError: If spectra_ref format is not recognized
        """
        match_scan = re.search(r"scan=(\d+)", spectra_ref)
        if match_scan:
            return match_scan.group(1)

        match_spectrum = re.search(r"spectrum=(\d+)", spectra_ref)
        if match_spectrum:
            return match_spectrum.group(1)

        try:
            if int(spectra_ref):
                return spectra_ref
        except ValueError:
            raise ValueError("Please check the 'spectra_ref' field in your mzTab file.")

    def _add_ms_values(
        self,
        info_df,
        ms_name,
        ms_with_psm,
        identified_spectrum_scan_id,
        mzml_charge_plot,
        mzml_peak_distribution_plot,
        mzml_peaks_ms2_plot,
        mzml_charge_plot_1,
        mzml_peak_distribution_plot_1,
        mzml_peaks_ms2_plot_1,
        ms_without_psm,
        enable_dia=False,
    ):
        """
        Process MS values from a dataframe row and add them to the appropriate histograms

        Args:
            info_df: Series row containing MS information
            ms_name: Name of the MS file
            ms_with_psm: List of MS files with PSMs
            identified_spectrum_scan_id: List of identified spectra by MS file
            mzml_charge_plot: Histogram for charge distribution of identified spectra
            mzml_peak_distribution_plot: Histogram for peak distribution of identified spectra
            mzml_peaks_ms2_plot: Histogram for peaks per MS2 of identified spectra
            mzml_charge_plot_1: Histogram for charge distribution of unidentified spectra
            mzml_peak_distribution_plot_1: Histogram for peak distribution of unidentified spectra
            mzml_peaks_ms2_plot_1: Histogram for peaks per MS2 of unidentified spectra
            ms_without_psm: List of MS files without PSMs
            enable_dia: Whether DIA mode is enabled
        """
        # info_df is a Pandas.Series not a DataFrame
        # "precursor_charge" --> "Charge",
        # "base_peak_intensity" --> "Base_Peak_Intensity",
        # "num_peaks": --> "MS_peaks"

        charge_state = (
            int(info_df["precursor_charge"]) if pd.notna(info_df["precursor_charge"]) else None
        )

        base_peak_intensity = (
            float(info_df["base_peak_intensity"]) if pd.notna(info_df["base_peak_intensity"]) else None
        )

        peak_per_ms2 = int(info_df["num_peaks"]) if pd.notna(info_df["num_peaks"]) else None

        if enable_dia:
            mzml_charge_plot.add_value(charge_state)
            mzml_peak_distribution_plot.add_value(base_peak_intensity)
            mzml_peaks_ms2_plot.add_value(peak_per_ms2)
            return

        if ms_name in ms_with_psm:
            # only "scan" in info_df not "SpectrumID"
            if info_df["scan"] in identified_spectrum_scan_id:
                mzml_charge_plot.add_value(charge_state)
                mzml_peak_distribution_plot.add_value(base_peak_intensity)
                mzml_peaks_ms2_plot.add_value(peak_per_ms2)
            else:
                mzml_charge_plot_1.add_value(charge_state)
                mzml_peak_distribution_plot_1.add_value(base_peak_intensity)
                mzml_peaks_ms2_plot_1.add_value(peak_per_ms2)
        else:
            ms_without_psm.append(ms_name)

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
    ):
        """
        Read MS info files and extract information

        Args:
            file_paths: List of paths to MS info files (optional, uses self.file_paths if not provided)
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

        Returns:
            tuple: (mzml_table, heatmap_charge, total_ms2_spectra, ms1_tic, ms1_bpc, ms1_peaks, ms1_general_stats)
        """
        # Use provided file_paths or fall back to instance attribute
        ms_info_path = file_paths or self.file_paths
        mzml_table = {}
        heatmap_charge = {}
        total_ms2_spectra = 0
        ms1_tic = {}
        ms1_bpc = {}
        ms1_peaks = {}
        ms1_general_stats = {}

        for file in ms_info_path:
            log.info(
                "{}: Parsing ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )
            mzml_df = pd.read_parquet(file)
            m = file_prefix(file).replace("_ms_info", "")
            if m not in mzml_table:
                mzml_table[m] = dict.fromkeys(["MS1_Num", "MS2_Num", "Charge_2"], 0)
            charge_group = mzml_df.groupby("precursor_charge").size()
            ms_level_group = mzml_df.groupby("ms_level").size()
            charge_2 = charge_group[2] if 2 in charge_group else 0
            ms1_number = int(ms_level_group[1]) if 1 in ms_level_group else 0
            ms2_number = int(ms_level_group[2]) if 2 in ms_level_group else 0
            total_ms2_spectra += ms2_number
            mzml_table[m].update({"MS1_Num": mzml_table[m]["MS1_Num"] + ms1_number})
            mzml_table[m].update({"MS2_Num": mzml_table[m]["MS2_Num"] + ms2_number})
            mzml_table[m].update({"Charge_2": mzml_table[m]["Charge_2"] + charge_2})

            (
                ms1_tic[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_bpc[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_peaks[os.path.basename(file).replace("_ms_info.parquet", "")],
                ms1_general_stats[os.path.basename(file).replace("_ms_info.parquet", "")],
            ) = get_ms_qc_info(mzml_df)

            group = mzml_df[mzml_df["ms_level"] == 2]
            del mzml_df

            if enable_dia:
                identified_spectrum_scan_id = []
            else:
                if not identified_spectrum:
                    raise ValueError(
                        "ms_io: The identified_spectrum is missing. Please check your mzTab file!"
                    )
                identified_spectrum_scan_id = [
                    self._spectra_ref_check(spectrum_id) for spectrum_id in identified_spectrum[m]
                ]

            # Apply _add_ms_values to each row
            for _, row in group.iterrows():
                self._add_ms_values(
                    row,
                    m,
                    ms_with_psm,
                    identified_spectrum_scan_id,
                    mzml_charge_plot,
                    mzml_peak_distribution_plot,
                    mzml_peaks_ms2_plot,
                    mzml_charge_plot_1,
                    mzml_peak_distribution_plot_1,
                    mzml_peaks_ms2_plot_1,
                    ms_without_psm,
                    enable_dia,
                )

            for m in mzml_table.keys():
                if mzml_table[m]["MS2_Num"] > 0:
                    heatmap_charge[m] = mzml_table[m]["Charge_2"] / mzml_table[m]["MS2_Num"]
                else:
                    heatmap_charge[m] = 0

            log.info(
                "{}: Done aggregating ms_statistics dataframe {}...".format(
                    datetime.now().strftime("%H:%M:%S"), file
                )
            )

        return (
            mzml_table,
            heatmap_charge,
            total_ms2_spectra,
            ms1_tic,
            ms1_bpc,
            ms1_peaks,
            ms1_general_stats,
        )