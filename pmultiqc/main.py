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

    log.info("Running pmultiqc Plugin v{}".format(config.pmultiqc_version))

    # Add to the search patterns used by modules
    if 'quantms/out_mzTab' not in config.sp:
        config.update_dict(config.sp, {'quantms/out_mzTab': {'fn': 'out.mzTab'}})
