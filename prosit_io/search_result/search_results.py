from abc import abstractmethod
from typing import Optional
import os
import logging

import pandas as pd

from prosit_io.file import csv

logger = logging.getLogger(__name__)


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

    def generate_internal(self, out_path: Optional[str] = None):
        if out_path is None:
            out_path = f"{os.path.splitext(self.path)[0]}.prosit"
        
        if os.path.isfile(out_path):
            logger.info(f"Found search results in internal format at {out_path}, skipping conversion")
            return out_path
        
        df = self.read_result(self.path)          
        csv.write_file(df, out_path)
          
        return out_path

    def read_internal(path):
        return csv.read_file(path)
