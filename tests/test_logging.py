"""
Test cases for the logging module.
"""

import logging

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

    def test_log_system_info(self):
        """Test log_system_info logs system information."""
        logger = get_logger("test_logger")
        # This should not raise an error
        log_system_info(logger)
