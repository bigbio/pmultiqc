from abc import abstractmethod
from abc import ABC
from pathlib import Path


class BaseParser(ABC):
    def __init__(self, file_paths: list[str, Path]):
        self.file_paths = file_paths
        

    @abstractmethod
    def parse(self, **kwargs) -> bool | None:
        pass

