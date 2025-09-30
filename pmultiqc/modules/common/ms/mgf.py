


import logging
import copy
import pandas as pd
import numpy as np
from collections import OrderedDict
from datetime import datetime
from pyteomics import mgf
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.ms.msreader import MSReader

# Initialise the logger
log = logging.getLogger(__name__)


class MGFReader(MSReader):
    """Class for reading and processing MGF files"""
    def read(self, file_paths=None, ms_with_psm=None, identified_spectrum=None, ms_without_psm=None):
        """
        Parse MGF files and extract information
        
        Args:
            file_paths: List of MGF file paths (optional, uses self.file_paths if not provided)
            ms_with_psm: List of files with PSM data
            identified_spectrum: Dictionary of identified spectra per file
            ms_without_psm: List of files without PSM data
            
        Returns:
            Dictionary containing parsed MGF data
        """
        def get_spectrum_id(spectrum_title, index):
            if "scan=" in spectrum_title:
                spectrum_id = spectrum_title
            else:
                spectrum_id = "index=" + str(index)
            return spectrum_id

        mgf_peak_distribution_plot = Histogram(
            "Peak Intensity",
            plot_category="range",
            breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000],
        )
        mgf_charge_plot = Histogram("Precursor Charge", plot_category="frequency")
        mgf_peaks_ms2_plot = Histogram(
            "#Peaks per MS/MS spectrum",
            plot_category="range",
            breaks=[i for i in range(0, 1001, 100)],
        )

        mgf_peak_distribution_plot_1 = copy.deepcopy(mgf_peak_distribution_plot)
        mgf_charge_plot_1 = copy.deepcopy(mgf_charge_plot)
        mgf_peaks_ms2_plot_1 = copy.deepcopy(mgf_peaks_ms2_plot)

        # Use provided file_paths or fall back to instance attribute
        mgf_paths = file_paths or self.file_paths
        heatmap_charge = dict()
        mgf_rtinseconds = {"spectrumID": [], "title": [], "filename": [], "retention_time": []}
        # Reset spectrum counts for this reading session
        self._reset_spectrum_counts()

        for m in mgf_paths:
            log.info("{}: Parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m))
            mgf_data = mgf.MGF(m)
            log.info(
                "{}: Done parsing MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )
            m = file_prefix(m)
            log.info(
                "{}: Aggregating MGF file {}...".format(datetime.now().strftime("%H:%M:%S"), m)
            )

            charge_2 = 0

            for i, spectrum in enumerate(mgf_data):
                charge_state = int(spectrum.get("params", {}).get("charge", [])[0])
                if charge_state == 2:
                    charge_2 += 1

                peak_per_ms2 = len(spectrum["m/z array"])
                base_peak_intensity = (
                    max(spectrum["intensity array"])
                    if len(spectrum["intensity array"]) > 0
                    else None
                )

                raw_title = spectrum.get("params", {}).get("title", [])
                mgf_rtinseconds["title"].append(raw_title)
                title = get_spectrum_id(raw_title, i)
                mgf_rtinseconds["spectrumID"].append(title)
                mgf_rtinseconds["filename"].append(m)

                rtinseconds = float(spectrum.get("params", {}).get("rtinseconds", None))
                mgf_rtinseconds["retention_time"].append(rtinseconds)

                if m in ms_with_psm:
                    if title in identified_spectrum[m]:
                        mgf_charge_plot.add_value(charge_state)
                        mgf_peak_distribution_plot.add_value(base_peak_intensity)
                        mgf_peaks_ms2_plot.add_value(peak_per_ms2)
                    else:
                        mgf_charge_plot_1.add_value(charge_state)
                        mgf_peak_distribution_plot_1.add_value(base_peak_intensity)
                        mgf_peaks_ms2_plot_1.add_value(peak_per_ms2)
                else:
                    if m not in ms_without_psm:
                        ms_without_psm.append(m)
            ms2_number = i + 1
            # Increment the total MS2 count for this file
            self._increment_ms2_count(ms2_number)

            heatmap_charge[m] = charge_2 / ms2_number
            log.info(
                "{}: Done aggregating MGF file {}...".format(
                    datetime.now().strftime("%H:%M:%S"), m
                )
            )

        mgf_rtinseconds_df = pd.DataFrame(mgf_rtinseconds)

        for i in ms_without_psm:
            log.warning("No PSM found in '{}'!".format(i))

        mgf_peaks_ms2_plot.to_dict()
        mgf_peak_distribution_plot.to_dict()
        mgf_peaks_ms2_plot_1.to_dict()
        mgf_peak_distribution_plot_1.to_dict()
        mgf_charge_plot.to_dict()
        mgf_charge_plot_1.to_dict()

        mgf_charge_plot.dict["cats"].update(mgf_charge_plot_1.dict["cats"])
        charge_cats_keys = [int(i) for i in mgf_charge_plot.dict["cats"]]
        charge_cats_keys.sort()
        mgf_charge_plot.dict["cats"] = OrderedDict(
            {str(i): mgf_charge_plot.dict["cats"][str(i)] for i in charge_cats_keys}
        )

        ms_info = {
            "charge_distribution": {
                "identified_spectra": mgf_charge_plot.dict["data"],
                "unidentified_spectra": mgf_charge_plot_1.dict["data"],
            },
            "peaks_per_ms2": {
                "identified_spectra": mgf_peaks_ms2_plot.dict["data"],
                "unidentified_spectra": mgf_peaks_ms2_plot_1.dict["data"],
            },
            "peak_distribution": {
                "identified_spectra": mgf_peak_distribution_plot.dict["data"],
                "unidentified_spectra": mgf_peak_distribution_plot_1.dict["data"],
            }
        }

        median = np.median(list(heatmap_charge.values()))
        heatmap_charge_score = dict(
            zip(
                heatmap_charge.keys(),
                list(map(lambda v: 1 - np.abs(v - median), heatmap_charge.values())),
            )
        )

        return {
            "mgf_rtinseconds": mgf_rtinseconds_df,
            "mgf_charge_plot": mgf_charge_plot,
            "mgf_peak_distribution_plot": mgf_peak_distribution_plot,
            "mgf_peaks_ms2_plot": mgf_peaks_ms2_plot,
            "mgf_charge_plot_1": mgf_charge_plot_1,
            "mgf_peak_distribution_plot_1": mgf_peak_distribution_plot_1,
            "mgf_peaks_ms2_plot_1": mgf_peaks_ms2_plot_1,
            "ms_info": ms_info,
            "heatmap_charge_score": heatmap_charge_score,
            "total_ms2_spectra": self.total_ms2_number,
            "oversampling_plot": None  # This will be set by the calling code
        }