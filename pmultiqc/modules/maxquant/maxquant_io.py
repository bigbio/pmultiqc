import os
from ..common.file_utils import file_prefix


# MaxQuant File Paths
def maxquant_file_path(find_log_files):

    required_files = [
        "parameters.txt",
        "summary.txt",
        "proteinGroups.txt",
        "evidence.txt",
        "msms.txt",
        "msScans.txt",
        "msmsScans.txt",
    ]

    maxquant_paths = {}
    for maxquant_file in find_log_files("pmultiqc/maxquant_result", filecontents=False):
        if maxquant_file["fn"] in required_files:
            f_path = os.path.join(maxquant_file["root"], maxquant_file["fn"])
            maxquant_paths[file_prefix(f_path)] = f_path

    return maxquant_paths
