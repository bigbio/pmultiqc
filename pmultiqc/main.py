#!/usr/bin/env python
"""
MultiQC pmultiqc plugin functions
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging
from multiqc import config

from pathlib import Path
import os
from pmultiqc.modules.common.file_utils import is_archive_file, extract_archive_file, get_clean_stem

# Initialise the main MultiQC logger
log = logging.getLogger("pmultiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.pmultiqc_version = get_distribution("pmultiqc").version


# Add default config options for the things that are used in MultiQC_NGI
def pmultiqc_plugin_execution_start():
    """Code to execute after the config files and
    command line flags have been parsed self.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get("disable_plugin", True):
        return None

    log.warning("Running pmultiqc Plugin v{}".format(config.pmultiqc_version))

    # If a compressed file is submitted, extract it first.
    analysis_dir_new = list()
    for anal_dir in config.analysis_dir:
        if is_archive_file(anal_dir):

            root = Path(anal_dir).parent
            file = Path(anal_dir).name
            extract_archive_file(root, file)

            output_dir = os.path.join(root, get_clean_stem(anal_dir))

            if os.path.exists(output_dir):
                analysis_dir_new.append(output_dir)
            else:
                raise SystemExit(f"Illegal file path: {output_dir}")
        else:
            analysis_dir_new.append(anal_dir)
    config.analysis_dir = analysis_dir_new
    
    # Module filename search patterns
    if "quantms/exp_design" not in config.sp:
        config.update_dict(config.sp, {"quantms/exp_design": {"fn": "experimental_design.tsv", "num_lines": 0}})

    if "quantms/sdrf" not in config.sp:
        config.update_dict(config.sp, {"quantms/sdrf": {"fn": "*.sdrf.tsv", "num_lines": 0}})

    if "quantms/mztab" not in config.sp:
        config.update_dict(config.sp, {"quantms/mztab": {"fn": "*.mzTab", "num_lines": 0}})

    if "quantms/mzML" not in config.sp:
        config.update_dict(config.sp, {"quantms/mzML": {"fn": "*.mzML", "num_lines": 0}})

    if "quantms/mgf" not in config.sp:
        config.update_dict(config.sp, {"quantms/mgf": {"fn": "*.mgf", "num_lines": 0}})

    if "quantms/mzid" not in config.sp:
        config.update_dict(config.sp, {"quantms/mzid": {"fn": "*.mzid", "num_lines": 0}})

    if "quantms/ms_info" not in config.sp:
        config.update_dict(config.sp, {"quantms/ms_info": {"fn": "*_ms_info.parquet", "num_lines": 0}})

    if "quantms/idXML" not in config.sp:
        config.update_dict(config.sp, {"quantms/idXML": {"fn": "*.idXML", "num_lines": 0}})

    if "quantms/msstats" not in config.sp:
        config.update_dict(config.sp, {"quantms/msstats": {"fn": "*msstats_in.csv", "num_lines": 0}})

    if "quantms/diann_report" not in config.sp:
        config.update_dict(config.sp, {"quantms/diann_report": {"fn": "*report.tsv", "num_lines": 0}})

    if "quantms/maxquant_result" not in config.sp:
        config.update_dict(config.sp, {"quantms/maxquant_result": {"fn": "*.txt", "num_lines": 0}})

    if "quantms/proteobench_result" not in config.sp:
        config.update_dict(config.sp, {"quantms/proteobench_result": {"fn": "result_performance.*", "num_lines": 0}})

    config.update({"log_filesize_limit": 200 * pow(1024, 3), "thousandsSep_format": ""})
