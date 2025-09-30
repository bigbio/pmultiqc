from multiqc import BaseMultiqcModule, config
from importlib import import_module

from .section_groups import SUB_SECTIONS

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

        # Halt execution if we've disabled the plugin
        if config.kwargs.get("disable_plugin", True):
            return None

        # Parse ProteoBench results
        if config.kwargs.get("parse_proteobench", False):

            ProteoBenchModule = get_module("proteobench", "ProteoBenchModule")
            pb = ProteoBenchModule(self.find_log_files)

            if pb.get_proteobench():
                pb.draw_report_plots()

            return None

        # section groups
        self.sub_sections = SUB_SECTIONS

        # Parse MaxQuant results
        if config.kwargs.get("parse_maxquant", False):

            MaxQuantModule = get_module("maxquant", "MaxQuantModule")
            mq = MaxQuantModule(
                self.find_log_files,
                self.sub_sections
            )

            if mq.get_data():
                mq.draw_report_plots()

            return None

        # MzIdentML
        if config.kwargs.get("mzid_plugin", False):
            # Use MzIdentMLModule for mzid plugin
            MzIdentMLModule = get_module("mzidentml", "MzIdentMLModule")
            MzIdentMLModule(
                self.find_log_files,
                self.sub_sections
            )

            return None
        
        # DIA-NN
        if config.kwargs.get("diann_plugin", False):
            # Use DiaNNModule for regular quantms processing
            DiaNNModule = get_module("diann", "DiaNNModule")
            DiaNNModule(
                self.find_log_files,
                self.sub_sections
            )

            return None

        # quantms (include quantms Dia)
        else:
            # Use QuantMSModule for regular quantms processing
            QuantMSModule = get_module("quantms", "QuantMSModule")
            QuantMSModule(
                self.find_log_files,
                self.sub_sections
            )

def get_module(module_name, class_name):
    module = import_module(f"..{module_name}", __package__)
    return getattr(module, class_name)
