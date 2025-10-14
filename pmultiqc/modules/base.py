
import logging
from abc import abstractmethod, ABC


class BasePMultiqcModule(ABC):

    def __init__(self, find_log_files_func, sub_sections, heatmap_colors):
        self.find_log_files = find_log_files_func
        self.sub_sections = sub_sections
        self.heatmap_color_list = heatmap_colors

        # Initialize logging for this module via centralized logger
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)

    @abstractmethod
    def get_data(self) -> bool | None:
        pass

    @abstractmethod
    def draw_plots(self) -> None:
        pass