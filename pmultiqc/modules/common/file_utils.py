from typing import Union, Optional
from pathlib import Path
import io


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
