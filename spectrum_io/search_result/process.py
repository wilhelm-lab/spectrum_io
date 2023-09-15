import logging
import re
from abc import abstractmethod
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from spectrum_io.file import csv

from .mascot import read_mascot
from .maxquant import read_maxquant
from .msamanda import read_msamanda
from .msfragger import read_msfragger

logger = logging.getLogger(__name__)


def read_search_results(search_results: Union[str, Path], search_type: str, tmt_labeled: str):
    """Read seach results."""
    if search_type.lower() == "maxquant":
        read_maxquant(search_results, tmt_labeled)
    elif search_type.lower() == "mascot":
        read_mascot(search_results, tmt_labeled)
    elif search_type.lower() == "msfragger":
        read_msfragger(search_results, tmt_labeled)
    elif search_type.lower() == "msamanda":
        read_msamanda(search_results)
    else:
        raise ValueError(f"Unknown search_type provided: {search_type}")


def generate_internal(self, tmt_labeled: str, out_path: Optional[Union[str, Path]] = None) -> Path:
    """
    Generate df and save to out_path.

    :param out_path: path to output
    :param tmt_labeled: tmt label as str
    :return: path to output file
    """
    if out_path is None:
        out_path = self.path.with_suffix(".prosit")
    if isinstance(out_path, str):
        out_path = Path(out_path)

    if out_path.is_file():
        logger.info(f"Found search results in internal format at {out_path}, skipping conversion")
        return out_path

    df = self.read_result(self.path, tmt_labeled)
    csv.write_file(df, out_path)

    return out_path


def read_internal(self, path: Union[str, Path]) -> pd.DataFrame:
    """
    Read file from path.

    :param path: path to file
    :return: dataframe after reading the file
    """
    return csv.read_file(path)
