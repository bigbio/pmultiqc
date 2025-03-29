#!/usr/bin/env python
"""
Logging system for pmultiqc package.

This module provides a centralized logging system with multiple verbosity levels
to track the execution flow and provide detailed information about the package's operations.
"""

import logging
import os
import sys
import time
from datetime import datetime

# Define log levels
LOG_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}

# Default format for log messages
DEFAULT_LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
DETAILED_LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s (%(filename)s:%(lineno)d): %(message)s"

# Global logger dictionary to keep track of created loggers
_loggers = {}


def get_logger(name="pmultiqc", level=None, log_file=None, detailed=False):
    """
    Get or create a logger with the specified name and configuration.

    Args:
        name (str): Name of the logger (default: "pmultiqc")
        level (str or int): Log level (default: from environment or INFO)
        log_file (str): Path to log file (default: None, logs to console only)
        detailed (bool): Whether to use detailed log format with file and line info

    Returns:
        logging.Logger: Configured logger instance
    """
    # Check if logger already exists
    if name in _loggers:
        return _loggers[name]

    # Create new logger
    logger = logging.getLogger(name)

    # Set level from parameter, environment, or default to INFO
    if level is None:
        level = os.environ.get("PMULTIQC_LOG_LEVEL", "INFO")

    if isinstance(level, str):
        level = LOG_LEVELS.get(level.upper(), logging.INFO)

    logger.setLevel(level)

    # Remove existing handlers if any
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)

    # Set formatter
    log_format = DETAILED_LOG_FORMAT if detailed else DEFAULT_LOG_FORMAT
    formatter = logging.Formatter(log_format)
    console_handler.setFormatter(formatter)

    # Add console handler to logger
    logger.addHandler(console_handler)

    # Add file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Store logger in global dictionary
    _loggers[name] = logger

    return logger


class Timer:
    """
    Timer class for measuring and logging execution time of operations.

    Usage:
        with Timer(logger, "Operation name"):
            # code to time
    """

    def __init__(self, logger, operation_name):
        """
        Initialize timer with logger and operation name.

        Args:
            logger (logging.Logger): Logger to use for timing messages
            operation_name (str): Name of the operation being timed
        """
        self.logger = logger
        self.operation_name = operation_name
        self.start_time = None

    def __enter__(self):
        """Start timing when entering context."""
        self.start_time = time.time()
        self.logger.info(f"Starting {self.operation_name}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Log elapsed time when exiting context."""
        elapsed_time = time.time() - self.start_time
        if exc_type:
            self.logger.error(f"{self.operation_name} failed after {elapsed_time:.2f} seconds")
        else:
            self.logger.info(f"Completed {self.operation_name} in {elapsed_time:.2f} seconds")


def log_system_info(logger):
    """
    Log system information at the start of execution.

    Args:
        logger (logging.Logger): Logger to use
    """
    import platform
    import multiprocessing

    logger.info("=" * 50)
    logger.info(f"pmultiqc execution started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"System: {platform.system()} {platform.release()} ({platform.machine()})")
    logger.info(f"Python: {platform.python_version()}")
    logger.info(f"CPU cores: {multiprocessing.cpu_count()}")
    try:
        import psutil

        memory = psutil.virtual_memory()
        logger.info(
            f"Memory: {memory.total / (1024**3):.1f} GB total, {memory.available / (1024**3):.1f} GB available"
        )
    except ImportError:
        logger.debug("psutil not available, skipping memory info")
    logger.info("=" * 50)


def configure_package_logging(level=None, log_file=None, detailed=False):
    """
    Configure logging for the entire pmultiqc package.

    Args:
        level (str or int): Log level (default: from environment or INFO)
        log_file (str): Path to log file (default: None, logs to console only)
        detailed (bool): Whether to use detailed log format with file and line info

    Returns:
        logging.Logger: Root logger for the package
    """
    # Configure root logger
    root_logger = get_logger("pmultiqc", level, log_file, detailed)

    # Configure module-specific loggers
    modules = ["main", "cli", "modules.quantms"]
    for module in modules:
        get_logger(f"pmultiqc.{module}", level, log_file, detailed)

    # Log system information
    log_system_info(root_logger)

    return root_logger
