"""mzQC export.

mzQC output is optional. It needs the ``pymzqc`` package, which is not a core
dependency of pmultiqc; install it with ``pip install "pmultiqc[mzqc]"``. When
the package is present the plugin modules compute and write an mzQC file next
to the MultiQC report, and when it is absent they skip that step entirely.

Nothing in this package module imports ``mzqc``, so it is safe to import on
any installation.
"""

import importlib.util
import logging

log = logging.getLogger("pmultiqc.export")

_availability = None


def is_mzqc_available() -> bool:
    """Return True when ``pymzqc`` is installed and mzQC export can run.

    The result is computed once per process. The first time the package is
    found missing an informational message explains how to enable the export;
    later calls stay quiet.
    """
    global _availability

    if _availability is None:
        # find_spec checks installability without importing the package
        _availability = importlib.util.find_spec("mzqc") is not None
        if not _availability:
            log.info(
                'mzQC export skipped: pymzqc is not installed '
                '(install it with: pip install "pmultiqc[mzqc]")'
            )

    return _availability
