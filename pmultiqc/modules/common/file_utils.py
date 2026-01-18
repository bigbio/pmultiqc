import gzip
import io
import logging
import os
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import Union, Optional

import pandas as pd

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


def extract_zip(file_path: str, extract_to: str) -> None:
    """Extract a zip file with security checks against path traversal and zip bombs."""
    # Security: Check file size first
    file_size = os.path.getsize(file_path)
    max_size = 10 * 1024 * 1024 * 1024  # 10GB limit

    with zipfile.ZipFile(file_path, "r") as zip_ref:
        # Get absolute path of extract_to for proper comparison
        extract_to_abs = os.path.abspath(extract_to)

        # Security: Check for path traversal in zip entries
        for member in zip_ref.namelist():
            # Check for explicit path traversal attempts (../ or ..\ patterns)
            # Allow filenames containing ".." as a substring (e.g., "file..name.txt")
            if "../" in member or "..\\" in member or member.startswith("../") or member.startswith("..\\"):
                raise ValueError(f"Invalid path in zip file: {member}")
            # Check for absolute paths
            if member.startswith("/") or (os.name == "nt" and len(member) > 1 and member[1] == ":"):
                raise ValueError(f"Invalid path in zip file: {member}")

            # Normalize the member path (remove leading slashes and normalize)
            member_normalized = os.path.normpath(member.lstrip("/"))

            # Join with extract_to and normalize to get the full path
            member_path = os.path.normpath(os.path.join(extract_to_abs, member_normalized))

            # Check if the resolved path is within extract_to directory
            if not member_path.startswith(extract_to_abs + os.sep) and member_path != extract_to_abs:
                raise ValueError(f"Attempted path traversal in zip file: {member}")

        # Security: Check for zip bombs
        total_uncompressed_size = sum(info.file_size for info in zip_ref.infolist())
        compression_ratio = total_uncompressed_size / file_size if file_size > 0 else 0
        if compression_ratio > 100 or total_uncompressed_size > max_size:
            raise ValueError(f"Suspicious zip file detected (compression ratio: {compression_ratio:.1f}:1)")

        zip_ref.extractall(extract_to)
        log.info(f"Extracted {file_path} to {extract_to}")


def extract_gz(file_path: str, extract_to: str) -> None:
    with gzip.open(file_path, "rb") as f_in:
        out_path = os.path.join(extract_to, os.path.basename(file_path).replace(".gz", ""))
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
            log.info(f"Extracted {file_path} to {out_path}")


def extract_tar(file_path: str, extract_to: str) -> None:
    """Extract a tar file with security checks against path traversal."""
    with tarfile.open(file_path, "r:*") as tar_ref:
        # Get absolute path of extract_to for proper comparison
        extract_to_abs = os.path.abspath(extract_to)

        # Security: Check for path traversal in tar entries
        for member in tar_ref.getmembers():
            # Check for explicit path traversal attempts (../ or ..\ patterns)
            # Allow filenames containing ".." as a substring (e.g., "file..name.txt")
            if "../" in member.name or "..\\" in member.name or member.name.startswith("../") or member.name.startswith("..\\"):
                raise ValueError(f"Invalid path in tar file: {member.name}")
            # Check for absolute paths
            if member.name.startswith("/") or (os.name == "nt" and len(member.name) > 1 and member.name[1] == ":"):
                raise ValueError(f"Invalid path in tar file: {member.name}")

            # Normalize the member path (remove leading slashes and normalize)
            member_normalized = os.path.normpath(member.name.lstrip("/"))

            # Join with extract_to and normalize to get the full path
            member_path = os.path.normpath(os.path.join(extract_to_abs, member_normalized))

            # Check if the resolved path is within extract_to directory
            if not member_path.startswith(extract_to_abs + os.sep) and member_path != extract_to_abs:
                raise ValueError(f"Attempted path traversal in tar file: {member.name}")

        tar_ref.extractall(extract_to)
        log.info(f"Extracted {file_path} to {extract_to}")


def extract_archive_file(root_dir: str, file_name: str) -> None:
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
    except (zipfile.BadZipFile, tarfile.TarError, gzip.BadGzipFile, OSError, ValueError) as e:
        raise ValueError(f"Failed to extract {file_path}: {e}") from e


def extract_files(folder_path: str) -> None:
    for root, _, files in os.walk(folder_path):
        for file in files:
            extract_archive_file(root, file)


def is_archive_file(file_path: str) -> bool:
    archive_types = [".zip", ".gz", ".tar", ".tar.gz", ".tgz", ".tar.bz2"]
    return any(file_path.lower().endswith(arch_type) for arch_type in archive_types)


def get_clean_stem(path: str) -> str:
    name = Path(path).name
    for ext in [".tar.gz", ".tar.bz2"]:
        if name.endswith(ext):
            return name[: -len(ext)]
    return Path(path).stem


def file_prefix(path: str) -> str:
    try:
        path = os.path.normpath(path)
        if "\\" in path:
            path = path.replace("\\", "/")
        return Path(path).stem
    except (OSError, ValueError) as e:
        raise ValueError(f"Illegal file path: {path}") from e


def drop_empty_row(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """
    Remove rows from a DataFrame where any of the specified columns are empty or NaN.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to filter.
    cols : list or iterable
        List of column names to check for empty or NaN values.

    Returns
    -------
    pandas.DataFrame
        A copy of the DataFrame with rows removed where any of the specified columns are empty or NaN.
    """
    mask = pd.Series(True, index=df.index)
    for col in cols:
        mask &= df[col].notna() & (df[col] != "")
    return df[mask].copy()