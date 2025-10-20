from abc import abstractmethod
from abc import ABC
from pathlib import Path


class BaseParser(ABC):
    """
    Abstract base class for MS data format parsers.
    
    Subclasses must implement the parse() method to process their specific file formats.
    
    Attributes:
        file_paths: List of file paths to parse (strings or Path objects)
    """

    def __init__(self, file_paths: list[str | Path]):
        self.file_paths = file_paths

    @abstractmethod
    def parse(self, **kwargs) -> bool | None:
        return None