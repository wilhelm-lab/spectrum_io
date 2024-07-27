import logging
import re
from abc import abstractmethod
from pathlib import Path
from typing import Optional, Union, Dict, Tuple
from spectrum_fundamentals.constants import update_custom_mods

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
    unsupported_mods = [r"Acetyl \(Protein N\-term\)", "ac", r"\[[0-9]+\]", r"\+"]
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
    def read_result(self, tmt_labeled: str, stat_mods: Optional[Dict[str, str]] = None, 
                    var_mods: Optional[Dict[str, str]] = None):
        """Read result.

        :param tmt_labeled: tmt label as str
        :param var_mods: variable modifications with custom identifier and respecitve internal equivalent 
        :param stat_mods: static modifications with custom identifier and respecitve internal equivalent

        """
        raise NotImplementedError

    def generate_internal(self, tmt_labeled: str, out_path: Optional[Union[str, Path]] = None, custom_mods: Optional[Dict[str, Dict[str, Tuple[str, float]]]] = None) -> pd.DataFrame:
        """
        Generate df and save to out_path if provided.

        :param out_path: path to output
        :param tmt_labeled: tmt label as str
        :param custom_mods: dict with static and variable custom modifications, their internal identifier and mass
        :raises AssertionError: if custom modification with illegal mass was provided
        :return: path to output file
        """
        stat_mods: Dict[str, str] = {key: value[0] for key, value in (custom_mods.get("stat_mods") or {}).items()}
        var_mods: Dict[str, str] = {key: value[0] for key, value in (custom_mods.get("var_mods") or {}).items()}

        if out_path is None:
            # convert and return
            return self.read_result(tmt_labeled, stat_mods=stat_mods, var_mods=var_mods)

        if isinstance(out_path, str):
            out_path = Path(out_path)

        if out_path.is_file():
            # only read converted and return
            logger.info(f"Found search results in internal format at {out_path}, skipping conversion")
            #TODO: internal_to_unimod
            return csv.read_file(out_path)

        # convert, save and return
        df = self.read_result(tmt_labeled, stat_mods=stat_mods, var_mods=var_mods)
        csv.write_file(df, out_path)
        return df

    def read_internal(self) -> pd.DataFrame:
        """
        Read file from path.

        :return: dataframe after reading the file
        """
        return csv.read_file(self.path)
