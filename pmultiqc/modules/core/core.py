from importlib import import_module

from matplotlib.colors import LinearSegmentedColormap, to_hex
from multiqc import BaseMultiqcModule, config

from pmultiqc.modules.common.logging import get_logger


# Initialize logger for this module
logger = get_logger("pmultiqc.modules.core")


PLUGIN_MAP = {
    "quantms_plugin": ("quantms", "QuantMSModule"),
    "diann_plugin": ("diann", "DiannModule"),
    "maxquant_plugin": ("maxquant", "MaxQuantModule"),
    "mzid_plugin": ("mzidentml", "MzIdentMLModule"),
    "proteobench_plugin": ("proteobench", "ProteoBenchModule"),
    "fragpipe_plugin": ("fragpipe", "FragPipeModule"),
}

class PMultiQC(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent module Class object
        super().__init__(
            name="pmultiqc",
            target="pmultiqc",
            anchor="pmultiqc",
            href="https://github.com/bigbio/pmultiqc",
            info="""
                is a MultiQC module to show the pipeline performance of mass spectrometry based quantification
                pipelines such as <a href='https://nf-co.re/quantms'>nf-core/quantms</a>,
                <a href='https://www.maxquant.org'>MaxQuant</a>,
                <a href='https://aptila.bio'>DIA-NN</a>,
                and <a href='https://fragpipe.nesvilab.org/'>FragPipe</a>
                """,
        )

        self.add_section(
            name="pmultiqc",
            anchor="pmultiqc",
            description="",
            content="",
        )

        # HeatMap color list
        color_map = LinearSegmentedColormap.from_list("red_green", ["#ff0000", "#00ff00"])
        self.heatmap_color_list = [
            [s, to_hex(color_map(s))] for s in [round(i * 0.1, 1) for i in range(11)]
        ]

        # section groups
        self.section_group_dict = {}
        self.sub_sections = {
            "experiment": [],
            "summary": [],
            "identification": [],
            "search_engine": [],
            "contaminants": [],
            "quantification": [],
            "ms1": [],
            "ms2": [],
            "mass_error": [],
            "rt_qc": [],
            "long_trends": [],
        }

        if config.kwargs.get("disable_hoverinfo", False):
            logger.info("disable_hoverinfo has been enabled by the user; hoverinfo will no longer be displayed.")

        # Run
        self.load_and_run_plugin()

    def load_and_run_plugin(self):

        plugin_loaded = False

        for flag, (module_name, class_name) in PLUGIN_MAP.items():

            if config.kwargs.get(flag, False):

                ModuleClass = get_module(module_name, class_name)

                if "proteobench_plugin" == flag:
                    plugin = ModuleClass(self.find_log_files, None, None)
                else:
                    plugin = ModuleClass(
                        self.find_log_files,
                        self.sub_sections,
                        self.heatmap_color_list
                    )

                if plugin.get_data():
                    plugin.draw_plots()

                plugin_loaded = True
                return

        if not plugin_loaded:
            raise ValueError("No pmultiqc plugin selected; skipping.")

def get_module(module_name, class_name):
    module = import_module(f"..{module_name}", __package__)
    return getattr(module, class_name)