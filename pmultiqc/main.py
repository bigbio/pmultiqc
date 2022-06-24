#!/usr/bin/env python
""" MultiQC pmultiqc plugin functions

"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging
from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger('pmultiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.pmultiqc_version = get_distribution("pmultiqc").version


# Add default config options for the things that are used in MultiQC_NGI
def pmultiqc_plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsed self.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
# 
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None

    log.warning("Running pmultiqc Plugin v{}".format(config.pmultiqc_version))

    if 'quantms/exp_design' not in config.sp:
        config.update_dict(config.sp, {'quantms/exp_design': {'fn': 'experimental_design.tsv', 'num_lines': 0}, 'shared': False})

    if 'quantms/sdrf' not in config.sp:
        config.update_dict(config.sp, {'quantms/sdrf': {'fn': '*.sdrf.tsv', 'num_lines': 0}, 'shared': False})

    # Add to the search patterns used by modules
    if 'quantms/mztab' not in config.sp:
        config.update_dict(config.sp, {'quantms/mztab': {'fn': '*.mzTab', 'num_lines': 0}})

    if 'quantms/mzML' not in config.sp:
        config.update_dict(config.sp, {'quantms/mzML': {'fn': '*.mzML', 'num_lines': 0}})

    if 'quantms/idXML' not in config.sp:
        config.update_dict(config.sp, {'quantms/idXML': {'fn': '*.idXML', 'num_lines': 0}})

    if 'quantms/msstats' not in config.sp:
        config.update_dict(config.sp, {'quantms/msstats': {'fn': '*msstats_in.csv', 'num_lines': 0}})

    if 'quantms/diann_report' not in config.sp:
        # TODO Why is the 'shared' suddenly inside the dict???
        config.update_dict(config.sp, {'quantms/diann_report': {'fn': '*report.tsv', 'num_lines': 0, 'shared': False}})

    config.update({'log_filesize_limit': 20000000000, 'thousandsSep_format': ''})