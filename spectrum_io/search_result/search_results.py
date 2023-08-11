import logging
import re
from abc import abstractmethod
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from spectrum_io.file import csv

logger = logging.getLogger(__name__)


def filter_valid_prosit_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter valid Prosit sequences.

    :param df: df to filter
    :return: df after filtering out unsupported peptides
    """
    logger.info(f"#sequences before filtering for valid prosit sequences: {len(df.index)}")
    # retain only peptides that fall within [7, 30] length supported by Prosit
    df = df[(df["PEPTIDE_LENGTH"] <= 30) & (df["PEPTIDE_LENGTH"] >= 7)]
    # remove unsupported mods to exclude
    unsupported_mods = [r"Acetyl \(Protein N\-term\)", "ac", r"\[[0-9]+\]"]
    exclude_mods_pattern = re.compile("|".join(unsupported_mods))
    df = df[~df["MODIFIED_SEQUENCE"].str.contains(exclude_mods_pattern)]
    # remove non-canonical aas
    df = df[(~df["SEQUENCE"].str.contains("U|O"))]
    # remove precursor charges greater than 6
    df = df[df["PRECURSOR_CHARGE"] <= 6]
    logger.info(f"#sequences after filtering for valid prosit sequences: {len(df.index)}")

    return df


class SearchResults:
    """Handle search results from different software."""

    orig_res: pd.DataFrame
    fake_msms: pd.DataFrame

    def __init__(self, path: Union[str, Path]):
        """
        Init Searchresults object.

        :param path: path to file
        """
        if isinstance(path, str):
            path = Path(path)
        self.path = path

    @abstractmethod
    def read_result(self, path: Union[str, Path], tmt_labeled: str):
        """Read result."""
        raise NotImplementedError

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
