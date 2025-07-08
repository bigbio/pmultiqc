from typing import Union, Optional
from pathlib import Path
import io
import zipfile
import gzip
import tarfile
import os
import shutil

import logging

log = logging.getLogger(__name__)


def get_filename(file: Union[Path, io.BufferedReader, io.StringIO, str]) -> Optional[str]:
    """
    Extract the filename from different file input types.

    Parameters:
    -----------
    file : Path | ReadCsvBuffer[bytes] | ReadCsvBuffer[str] | str
        The file input, which could be:
        - A pathlib.Path object
        - A pandas-compatible buffer object (StringIO or BytesIO)
        - A string filepath

    Returns:
    --------
    Optional[str]
        The extracted filename, or None if it couldn't be determined
    """
    # Check if it's a Path object
    if isinstance(file, Path):
        return file.name

    # Check if it's a string (filepath)
    if isinstance(file, str):
        return Path(file).name if "/" in file or "\\" in file else file

    # Check if it's a buffer with a name attribute (some file objects have this)
    if hasattr(file, "name") and isinstance(file.name, str):
        return Path(file.name).name

    # For StringIO/BytesIO objects without names or other cases
    return None


def extract_zip(file_path, extract_to):
    with zipfile.ZipFile(file_path, "r") as zip_ref:
        zip_ref.extractall(extract_to)
        log.info(f"Extracted {file_path} to {extract_to}")


def extract_gz(file_path, extract_to):
    with gzip.open(file_path, "rb") as f_in:
        out_path = os.path.join(extract_to, os.path.basename(file_path).replace(".gz", ""))
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
            log.info(f"Extracted {file_path} to {out_path}")


def extract_tar(file_path, extract_to):
    with tarfile.open(file_path, "r:*") as tar_ref:
        tar_ref.extractall(extract_to)
        log.info(f"Extracted {file_path} to {extract_to}")


def extract_archive_file(root_dir, file_name):
    file_path = os.path.join(root_dir, file_name)
    
    try:
        # *.zip
        if file_name.endswith(".zip"):
            extract_zip(file_path, root_dir)
        # *.gz
        elif file_name.endswith(".gz") and not file_name.endswith(".tar.gz"):
            extract_gz(file_path, root_dir)
        # .tar, .tar.gz, .tar.bz2
        elif (
            file_name.endswith(".tar")
            or file_name.endswith(".tar.gz")
            or file_name.endswith(".tgz")
            or file_name.endswith(".tar.bz2")
        ):
            extract_tar(file_path, root_dir)
    except Exception:
            raise SystemExit(f"Failed to extract: {file_path}")


def extract_files(folder_path):
    for root, _, files in os.walk(folder_path):
        for file in files:
            extract_archive_file(root, file)


def is_archive_file(file_path):
    archive_types = [".zip", ".gz", ".tar", ".tar.gz", ".tgz", ".tar.bz2"]
    return any(file_path.lower().endswith(arch_type) for arch_type in archive_types)


def get_clean_stem(path):
    name = Path(path).name
    for ext in [".tar.gz", ".tar.bz2"]:
        if name.endswith(ext):
            return name[: -len(ext)]
    return Path(path).stem


def file_prefix(path):
    try:
        path = os.path.normpath(path)
        if "\\" in path:
            path = path.replace("\\", "/")
        return Path(path).stem
    except:
        raise SystemExit(f"Illegal file path: {path}")
