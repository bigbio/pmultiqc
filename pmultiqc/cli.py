#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.

See the Click documentation for more command line flag types:
http://click.pocoo.org/5/
"""

import click


exp_design = click.option('--exp_design', 'exp_design', required=True)
mzMLs = click.option('--mzMLs', 'mzMLs', required=True)
raw_ids = click.option('--raw_ids', 'raw_ids', required=True)


disable_plugin = click.option('--disable-example-plugin', 'disable_plugin',
                              is_flag=True,
                              help="Disable the Example MultiQC plugin on this run"
                              )
