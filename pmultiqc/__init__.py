"""
pmultiqc - A MultiQC plugin for proteomics data analysis
"""

from importlib import metadata

try:
    __version__ = metadata.version("pmultiqc")
except metadata.PackageNotFoundError:
    __version__ = "unknown"
