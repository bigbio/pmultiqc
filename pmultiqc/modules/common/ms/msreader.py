"""
Abstract base class for MS file readers
"""

from __future__ import absolute_import
from abc import ABC, abstractmethod
import logging
from typing import List, Dict, Any

# Initialise the logger
log = logging.getLogger(__name__)


class MSReader(ABC):
    """Abstract base class for reading and processing MS files"""

    def __init__(self, file_paths: List[str] = None):
        """
        Initialize the MS reader with common attributes
        
        Args:
            file_paths: List of file paths to read data from
        """
        self.logger = log
        self.total_ms1_number = 0
        self.total_ms2_number = 0
        self.file_paths = file_paths or []

    @abstractmethod
    def read(self, *args, **kwargs) -> Dict[str, Any]:
        """
        Abstract method for reading MS files.

        Each subclass must implement this method with the appropriate
        parameters for their specific file type.

        Returns:
            Dict[str, Any]: Dictionary containing parsed data
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

    def _validate_paths(self, paths: List[str]) -> List[str]:
        """
        Validate that file paths exist and are accessible.

        Args:
            paths: List of file paths to validate

        Returns:
            List[str]: Valid file paths

        Raises:
            FileNotFoundError: If any path doesn't exist
        """
        import os
        valid_paths = []
        for path in paths:
            if os.path.exists(path):
                valid_paths.append(path)
            else:
                self._log_warning(f"File not found: {path}")
        return valid_paths

    def _get_common_parameters(self) -> Dict[str, Any]:
        """
        Get common parameters that are shared across all MS readers.

        Returns:
            Dict[str, Any]: Common parameters dictionary
        """
        return {
            "ms_with_psm": [],
            "identified_spectrum": {},
            "ms_without_psm": [],
            "enable_dia": False
        }

    def _increment_ms1_count(self, count: int = 1) -> None:
        """
        Increment the total MS1 spectrum count.

        Args:
            count: Number to add to total MS1 count (default: 1)
        """
        self.total_ms1_number += count

    def _increment_ms2_count(self, count: int = 1) -> None:
        """
        Increment the total MS2 spectrum count.

        Args:
            count: Number to add to total MS2 count (default: 1)
        """
        self.total_ms2_number += count

    def _reset_spectrum_counts(self) -> None:
        """Reset both total MS1 and MS2 spectrum counts to zero."""
        self.total_ms1_number = 0
        self.total_ms2_number = 0

    def _get_spectrum_counts(self) -> Dict[str, int]:
        """
        Get current total spectrum counts.

        Returns:
            Dict[str, int]: Dictionary with total MS1 and MS2 counts
        """
        return {
            "total_ms1_number": self.total_ms1_number,
            "total_ms2_number": self.total_ms2_number
        }
