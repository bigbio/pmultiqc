"""MultiQC support for mzQC quality-control documents."""

from .mzqc import MzQCModule, MzQCRun, parse_mzqc_document

__all__ = ["MzQCModule", "MzQCRun", "parse_mzqc_document"]
