import logging
from pathlib import Path
from typing import Union

import pandas as pd
import spectrum_fundamentals.constants as c
from pyteomics import pepxml
from spectrum_fundamentals.mod_string import internal_without_mods, msfragger_to_internal
from tqdm import tqdm

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    @staticmethod
    def read_result(path: Union[str, Path], tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param path: path to pepXML folder or single pepXML file to read
        :param tmt_labeled: tmt label as str
        :raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :return: pd.DataFrame with the formatted data
        """
        if isinstance(path, str):
            path = Path(path)

        if path.is_file():
            file_list = [path]
        elif path.is_dir():
            file_list = list(path.rglob("*.pepXML"))
        else:
            raise FileNotFoundError(f"{path} could not be found.")

        ms_frag_results = []
        for pep_xml_file in tqdm(file_list):
            ms_frag_results.append(pepxml.DataFrame(str(pep_xml_file)))

        df = pd.concat(ms_frag_results)

        df = update_columns_for_prosit(df, tmt_labeled)
        return filter_valid_prosit_sequences(df)


def update_columns_for_prosit(df, tmt_labeled: str) -> pd.DataFrame:
    """
    Update columns of df to work with Prosit.

    :param df: df to modify
    :param tmt_labeled: True if tmt labeled
    :return: modified df as pd.DataFrame
    """
    df["REVERSE"] = df["protein"].apply(lambda x: "rev" in str(x))
    df["RAW_FILE"] = df["spectrum"].apply(lambda x: x.split(".")[0])
    df["MASS"] = df["precursor_neutral_mass"]
    df["PEPTIDE_LENGTH"] = df["peptide"].apply(lambda x: len(x))

    if tmt_labeled != "":
        unimod_tag = c.TMT_MODS[tmt_labeled]
        logger.info("Adding TMT fixed modifications")
        df["MODIFIED_SEQUENCE"] = msfragger_to_internal(
            df["modified_peptide"].to_list(),
            fixed_mods={"C": "C[UNIMOD:4]", r"n[\d+]": f"{unimod_tag}-", "K": f"K{unimod_tag}"},
        )
    else:
        df["MODIFIED_SEQUENCE"] = msfragger_to_internal(df["modified_peptide"].to_list())

    df.rename(
        columns={
            "assumed_charge": "PRECURSOR_CHARGE",
            "index": "SCAN_EVENT_NUMBER",
            "peptide": "SEQUENCE",
            "start_scan": "SCAN_NUMBER",
            "hyperscore": "SCORE",
        },
        inplace=True,
    )
    df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
    return df[
        [
            "RAW_FILE",
            "SCAN_NUMBER",
            "MODIFIED_SEQUENCE",
            "PRECURSOR_CHARGE",
            "SCAN_EVENT_NUMBER",
            "MASS",
            "SCORE",
            "REVERSE",
            "SEQUENCE",
            "PEPTIDE_LENGTH",
        ]
    ]
