from importlib import import_module

from matplotlib.colors import LinearSegmentedColormap, to_hex
from multiqc import BaseMultiqcModule, config


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
                and <a href='https://aptila.bio'>DIA-NN</a>
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
        heatmap_color_list = [
            [s, to_hex(color_map(s))] for s in [round(i * 0.1, 1) for i in range(11)]
        ]

        # section groups
        self.section_group_dict = dict()
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
        }

        # Parse ProteoBench results
        if config.kwargs.get("proteobench_plugin", False):

            ProteoBenchModule = get_module("proteobench", "ProteoBenchModule")
            pb = ProteoBenchModule(self.find_log_files, None, None)

            if pb.get_data():
                pb.draw_plots()

        # Parse MaxQuant results
        elif config.kwargs.get("maxquant_plugin", False):

            MaxQuantModule = get_module("maxquant", "MaxQuantModule")
            mq = MaxQuantModule(self.find_log_files, self.sub_sections, heatmap_color_list)

            if mq.get_data():
                mq.draw_plots()

        # Parse mzIdentML results
        elif config.kwargs.get("mzid_plugin", False):

            MzIdentMLModule = get_module("mzidentml", "MzIdentMLModule")
            mzid = MzIdentMLModule(self.find_log_files, self.sub_sections, heatmap_color_list)
            if mzid.get_data():
                mzid.draw_plots()

        # Parse DIA-NN results
        elif config.kwargs.get("diann_plugin", False):
            DiannModule = get_module("diann", "DiannModule")
            DiannModule(self.find_log_files, self.sub_sections, heatmap_color_list)

        # quantms, DIA-NN results
        elif config.kwargs.get("quantms_plugin", False):
            QuantMSModule = get_module("quantms", "QuantMSModule")
            QuantMSModule(self.find_log_files, self.sub_sections, heatmap_color_list)
        else:
            raise ValueError("No plugin defined")


def get_module(module_name, class_name):
    module = import_module(f"..{module_name}", __package__)
    return getattr(module, class_name)