"""Software identities for mzQC exports; input provenance is handled separately."""

import logging

log = logging.getLogger("pmultiqc.mzqc_generator")

# PSI-MS terms: https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo
SOFTWARE = {
    "diann": ("MS:1003253", "DIA-NN", "https://github.com/vdemichev/DiaNN"),
    "maxquant": ("MS:1001583", "MaxQuant", "https://www.maxquant.org"),
    "quantms": ("MS:1003425", "quantms", "https://quantms.org"),
}


def build_analysis_software(
    pipeline_name: str, software_version: str | None = None
) -> dict[str, str]:
    """
    Describe a supported workflow without inventing a software version.

    mzQC requires a version string. Missing versions are explicitly recorded as
    unknown. Unknown workflow names raise rather than claiming a different tool.
    """
    pipeline = pipeline_name.lower().replace("-", "")
    if pipeline not in SOFTWARE:
        raise ValueError(f"Unsupported mzQC pipeline: {pipeline_name}")
    accession, name, uri = SOFTWARE[pipeline]
    version = str(software_version).strip() if software_version is not None else ""
    if not version:
        log.warning("mzQC: %s software version is unavailable; recording 'unknown'.", name)
        version = "unknown"
    return {"accession": accession, "name": name, "version": version, "uri": uri}
