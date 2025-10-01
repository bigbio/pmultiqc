import logging

from .proteobench_utils import get_pb_data
from ..core.section_groups import add_group_modules, add_sub_section
from ..common.base_module import BaseModule


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class ProteoBenchModule(BaseModule):

    def __init__(self, find_log_files_func, sub_sections):
        
        # Define file patterns for this module
        self.file_patterns = {
            "proteobench_result": "pmultiqc/proteobench_result"
        }

        # Call parent constructor
        super().__init__(
            find_log_files_func=find_log_files_func,
            sub_sections=sub_sections,
            module_name="proteobench"
        )

        # Initialize module-specific attributes
        self.pb_results = None

    def _initialize_module(self):
        
        # Discover required files
        self.file_paths = self._discover_files(self.file_patterns)

        # Validate required files
        self._validate_required_files(self.file_paths)

    def get_data(self):
        
        return self.get_proteobench()

    def get_proteobench(self):
        
        pb_file_path = self.file_paths.get("proteobench_result")

        if not pb_file_path:
            raise ValueError("ProteoBench result file not found")

        self.pb_results = get_pb_data(pb_file_path)
        self.results = self.pb_results

        return self.pb_results

    def draw_plots(self):
        
        self.log.info("Start plotting the ProteoBench results...")

        # 1. precursor ion
        precursor_sub_section = list()

        add_sub_section(
            sub_section=precursor_sub_section,
            plot=self.pb_results["charge_html"],
            order=1,
            description="""
                This bar chart shows the distribution of precursor ion charges.
                """,
            helptext="""
                [result_performance.csv] The precursor ion charge extracted from the 'precursor ion' column.
                """,
        )

        # 2. Intensity
        log_mean_sub_section = list()

        add_sub_section(
            sub_section=log_mean_sub_section,
            plot=self.pb_results["log_mean_html"]["linegraph_html"],
            order=1,
            description="""
                Distribution of 'log_Intensity_mean' under Condition A and B.
                """,
            helptext="""
                [result_performance.csv] This plot visualizes the distribution of mean intensity 
                values after log2 transformation for both Condition A and Condition B.
                """,
        )
        add_sub_section(
            sub_section=log_mean_sub_section,
            plot=self.pb_results["log_mean_html"]["bar_html"],
            order=2,
            description="""
                Number of missing (NA) values in 'log_Intensity_mean' for Condition A and B.
                """,
            helptext="""
                [result_performance.csv] This plot shows the number of missing (NA) values in the 
                mean log2-transformed intensities for Condition A and Condition B.
                """,
        )
        add_sub_section(
            sub_section=log_mean_sub_section,
            plot=self.pb_results["log_intensity_html"]["linegraph_html"],
            order=3,
            description="""
                Distribution of intensity for each run.
                """,
            helptext="""
                [result_performance.csv] Distribution of intensity for each run.
                """,
        )
        add_sub_section(
            sub_section=log_mean_sub_section,
            plot=self.pb_results["log_intensity_html"]["bar_html"],
            order=4,
            description="""
                Number of missing (NA) values for each run.
                """,
            helptext="""
                [result_performance.csv] Number of missing (NA) values for each run.
                """,
        )
        add_sub_section(
            sub_section=log_mean_sub_section,
            plot=self.pb_results["num_inten_per_file_html"],
            order=5,
            description="""
                Number of detected (Non-NA) and missing (NA) values per sample file.
                """,
            helptext="""
                [result_performance.csv] This plot shows the number of missing (NA) values in the 
                mean log2-transformed intensities for Condition A and Condition B.
                """,
        )

        # 3. Standard Deviation of Intensity
        log_std_sub_section = list()

        add_sub_section(
            sub_section=log_std_sub_section,
            plot=self.pb_results["log_std_html"]["box_html"],
            order=1,
            description="""
                This plot shows the distribution of 'log_Intensity_std' values for Condition A and Condition B.
                """,
            helptext="""
                [result_performance.csv] This plot shows the distribution of standard deviations 
                calculated from log2-transformed intensity values for Condition A and Condition B.
                """,
        )

        # 4. CV
        cv_sub_section = list()

        add_sub_section(
            sub_section=cv_sub_section,
            plot=self.pb_results["cv_html"]["linegraph_html"],
            order=1,
            description="""
                This plot shows the distribution of coefficient of variation (CV) values for Condition A and Condition B.
                """,
            helptext="""
                [result_performance.csv] This plot shows the distribution of coefficient of variation (CV) values 
                for Condition A and Condition B.
                """,
        )
        add_sub_section(
            sub_section=cv_sub_section,
            plot=self.pb_results["cv_html"]["bar_html"],
            order=2,
            description="""
                This plot shows the number of missing values (NAs) in the coefficient of variation (CV) for condition A and B.
                """,
            helptext="""
                [result_performance.csv] This plot shows the number of missing values (NAs) in the coefficient of 
                variation (CV) for condition A and B.
                """,
        )

        # 5. log2_A_vs_B
        log_vs_sub_section = list()

        add_sub_section(
            sub_section=log_vs_sub_section,
            plot=self.pb_results["log_vs_html"]["linegraph_html"],
            order=1,
            description="""
                This plot shows the distribution of 'log2_A_vs_B' values for Condition A and Condition B.
                """,
            helptext="""
                [result_performance.csv] Distribution of log₂ fold changes (log2_A_vs_B) between Condition A and B 
                based on mean log2-transformed intensities.
                """,
        )
        add_sub_section(
            sub_section=log_vs_sub_section,
            plot=self.pb_results["logfc_logmean_html"],
            order=2,
            description="""
                Distribution of mean intensity across all runs and log2 fold change. Legend: ECOLI (blue), HUMAN (green), YEAST (red).
                """,
            helptext="""
                [result_performance.csv] Distribution of mean intensity across all runs and log2 fold change (log2FC). 
                Legend: ECOLI (blue), HUMAN (green), YEAST (red).
                """,
        )
        add_sub_section(
            sub_section=log_vs_sub_section,
            plot=self.pb_results["epsilon_html"]["linegraph_html"],
            order=3,
            description="""
                Distribution of 'epsilon' values (difference between observed and expected log2 fold changes).
                """,
            helptext="""
                [result_performance.csv] 'Epsilon' measures the deviation between observed and expected log2 fold changes, 
                indicating agreement between data and expectations.
                """,
        )

        # Create section groups manually (legacy approach for now)
        self.section_group_dict = {
            "precursor_sub_section": precursor_sub_section,
            "log_mean_sub_section": log_mean_sub_section,
            "log_std_sub_section": log_std_sub_section,
            "cv_sub_section": cv_sub_section,
            "log_vs_sub_section": log_vs_sub_section,
        }

        add_group_modules(self.section_group_dict, "proteobench")
