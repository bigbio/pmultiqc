"""
Test cases for the logging module.
"""

import logging
import sys
import unittest.mock as mock

import pytest

from pmultiqc.modules.common.logging import get_logger, log_system_info


class TestLogging:
    """Test logging functionality."""

    def test_get_logger(self):
        """Test that get_logger creates a logger."""
        logger = get_logger("test_logger")
        assert logger is not None
        assert isinstance(logger, logging.Logger)
        assert logger.name == "test_logger"

    def test_log_system_info_with_psutil(self):
        """Test log_system_info when psutil is available."""
        logger = get_logger("test_logger")
        # This should not raise an error whether psutil is installed or not
        log_system_info(logger)

    def test_log_system_info_without_psutil(self):
        """Test log_system_info when psutil is not available."""
        logger = get_logger("test_logger")
        
        # Mock the import to simulate psutil not being available
        with mock.patch.dict(sys.modules, {'psutil': None}):
            # Remove psutil from sys.modules if it exists
            original_psutil = sys.modules.pop('psutil', None)
            try:
                # This should not raise an error
                log_system_info(logger)
            finally:
                # Restore psutil if it was available
                if original_psutil is not None:
                    sys.modules['psutil'] = original_psutil
