"""
Abstract base class for reading and processing identification files
"""

from __future__ import absolute_import
from abc import ABC, abstractmethod
import logging
from typing import Any

# Initialise the logger
log = logging.getLogger(__name__)


class IDReader(ABC):
    """Abstract base class for reading and processing identification files"""

    def __init__(self, file_paths=None):
        """
        Initialize the ID reader with common attributes

        Args:
            file_paths: List of file paths to read data from
        """
        self.logger = log
        self.file_paths = file_paths or []

    @abstractmethod
    def read(self, *args, **kwargs) -> Any:
        """
        Abstract method for reading identification files.

        This method must be implemented by subclasses to read and process
        identification data from files.

        Args:
            *args: Variable length argument list
            **kwargs: Arbitrary keyword arguments

        Returns:
            Any: Processed identification data (format depends on subclass implementation)

        Raises:
            NotImplementedError: If the subclass hasn't implemented this method
        """
        raise NotImplementedError("Subclasses must implement the read method")

    def _log_info(self, message: str) -> None:
        """Log an info message"""
        self.logger.info(message)

    def _log_warning(self, message: str) -> None:
        """Log a warning message"""
        self.logger.warning(message)

    def _log_error(self, message: str) -> None:
        """Log an error message"""
        self.logger.error(message)

    def _validate_paths(self, paths) -> None:
        """
        Validate if all provided paths exist.

        Args:
            paths: List of file paths to validate

        Raises:
            FileNotFoundError: If any path doesn't exist
        """
        import os
        for path in paths:
            if not os.path.exists(path):
                self._log_error(f"File not found: {path}")
                raise FileNotFoundError(f"File not found: {path}")
