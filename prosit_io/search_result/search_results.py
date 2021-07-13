from abc import abstractmethod
import pandas as pd
from prosit_io.file.csv import write_file


class SearchResults:
    """
       Handle search results from different software.
    """
    path: str
    orig_res: pd.DataFrame
    fake_msms: pd.DataFrame

    def __init__(self, path):
        self.path = path

    @abstractmethod
    def read_result(self):
        raise NotImplementedError

    def handle_result(self):
        self.read_result()
        write_file(self.result, self.path)
