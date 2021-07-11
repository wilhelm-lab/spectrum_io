from abc import abstractmethod
import pandas as pd
from io.file.csv import write_file


class SearchResults:
    """
       Handle search results from different software.
    """
    path: str
    orig_res: pd.DataFrame
    fake_msms: pd.DataFrame

    def __init__(self, path):
        self.path = path

    def read_result(self):
        pass

    @abstractmethod
    def gen_fake(self):
        # set fake_msms
        pass

    def handle_result(self):
        self.read_result()
        self.gen_fake()
        write_file(self.fake_msms, self.path)
