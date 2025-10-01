"""
Base module class for all PMultiQC modules.

This abstract base class provides common functionality shared across
all analysis modules (quantms, maxquant, diann, mzidentml, proteobench).
"""

import logging
import os
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Callable

from multiqc import config
from pmultiqc.modules.common.plots import HEATMAP_COLOR_LIST
from pmultiqc.modules.core.section_groups import add_group_modules


class BaseModule(ABC):
    """
    Abstract base class for all PMultiQC analysis modules.

    Provides common initialization, file handling, and plotting infrastructure
    that all analysis modules can extend.
    """

    def __init__(
        self,
        find_log_files_func: Callable,
        sub_sections: Dict[str, Any],
        module_name: str
    ):
        """
        Initialize base module with common functionality.

        Args:
            find_log_files_func: Function to discover input files
            sub_sections: Dictionary of report section containers
            module_name: Name of the specific module (for logging)
        """
        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections
        self.module_name = module_name

        # Initialize logger
        self.log = logging.getLogger(f"{__name__}.{module_name}")

        # Common attributes
        self.heatmap_color_list = HEATMAP_COLOR_LIST
        self.results = {}
        self.file_paths = {}

        # Initialize module-specific data
        self._initialize_module()

        self.log.info(f"Initialized {module_name} module")

    def _initialize_module(self):
        """Initialize module-specific attributes. Override in subclasses."""

    def _setup_output_directory(self, output_dir: Optional[str] = None) -> str:
        """
        Setup output directory for module results.

        Args:
            output_dir: Optional custom output directory

        Returns:
            Path to the output directory
        """
        if output_dir:
            base_dir = output_dir
        elif config.output_dir:
            base_dir = config.output_dir
        else:
            base_dir = os.getcwd()

        if not os.path.exists(base_dir):
            os.makedirs(base_dir)

        return base_dir

    def _discover_files(self, file_patterns: Dict[str, str]) -> Dict[str, str]:
        """
        Discover input files using specified patterns.

        Args:
            file_patterns: Dictionary mapping file types to search patterns

        Returns:
            Dictionary mapping file types to discovered file paths
        """
        discovered_files = {}

        for file_type, pattern in file_patterns.items():
            for file_info in self.find_log_files(pattern, filecontents=False):
                file_path = os.path.join(file_info["root"], file_info["fn"])
                discovered_files[file_type] = file_path
                self.log.info(f"Found {file_type} file: {file_path}")

        return discovered_files

    def _validate_required_files(self, required_files: Dict[str, str]) -> None:
        """
        Validate that required files are present.

        Args:
            required_files: Dictionary of required file types and their paths

        Raises:
            ValueError: If required files are missing
        """
        missing_files = []

        for file_type, file_path in required_files.items():
            if file_path is None:
                missing_files.append(file_type)
            elif not os.path.exists(file_path):
                missing_files.append(f"{file_type} (not found: {file_path})")

        if missing_files:
            raise ValueError(f"Missing required files for {self.module_name}: {', '.join(missing_files)}")

    def _create_section_groups(self, section_mapping: Dict[str, str]) -> None:
        """
        Create section groups for the report.

        Args:
            section_mapping: Dictionary mapping internal section names to report sections
        """
        section_group_dict = {}

        for internal_name, section_key in section_mapping.items():
            if section_key in self.sub_sections:
                section_group_dict[f"{internal_name}_sub_section"] = self.sub_sections[section_key]

        add_group_modules(section_group_dict, self.module_name)

    @abstractmethod
    def get_data(self) -> Dict[str, Any]:
        """
        Extract and process data from input files.

        Returns:
            Dictionary containing processed data
        """

    @abstractmethod
    def draw_plots(self) -> None:
        """
        Generate plots and add them to report sections.
        """

    def get_results(self) -> Dict[str, Any]:
        """
        Get processed results from the module.

        Returns:
            Dictionary containing module results
        """
        return self.results

    def cleanup(self) -> None:
        """Clean up temporary files and resources."""

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup()