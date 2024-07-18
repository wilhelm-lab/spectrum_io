import logging
from pathlib import Path
from typing import Optional, Union, Dict, Tuple

import pandas as pd
import spectrum_fundamentals.constants as c
from pyteomics import pepxml
from spectrum_fundamentals.mod_string import internal_without_mods, msfragger_or_custom_to_internal
from spectrum_fundamentals.constants import MSFRAGGER_VAR_MODS
from tqdm import tqdm

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    def read_result(self, tmt_labeled: str, stat_mods: Optional[Dict[str, str]] = None, 
                    var_mods: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :param var_mods: dict with custom variable identifier and respecitve internal equivalent 
        :param stat_mods: dict with custom static identifier and respecitve internal equivalent:raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :return: pd.DataFrame with the formatted data
        """
        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = list(self.path.rglob("*.pepXML"))
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")

        ms_frag_results = []
        for pep_xml_file in tqdm(file_list):
            ms_frag_results.append(pepxml.DataFrame(str(pep_xml_file)))

        df = pd.concat(ms_frag_results)

        df = update_columns_for_prosit(df, tmt_labeled, stat_mods=stat_mods, var_mods=var_mods)
        return filter_valid_prosit_sequences(df)


def update_columns_for_prosit(df, tmt_labeled: str, stat_mods: Optional[Dict[str, str]] = None, 
                              var_mods: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """
    Update columns of df to work with Prosit.

    :param df: df to modify
    :param tmt_labeled: True if tmt labeled
    :param var_mods: dict with custom variable identifier and respecitve internal equivalent 
    :param stat_mods: dict with custom static identifier and respecitve internal equivalent
    :return: modified df as pd.DataFrame
    """
    df["PROTEINS"] = df["protein"]
    df["PROTEINS"].fillna("UNKNOWN", inplace=True)
    df["REVERSE"] = df["protein"].apply(lambda x: "rev" in str(x))
    df["RAW_FILE"] = df["spectrum"].apply(lambda x: x.split(".")[0])
    df["MASS"] = df["precursor_neutral_mass"]
    df["PEPTIDE_LENGTH"] = df["peptide"].apply(lambda x: len(x))

    mods = {**(MSFRAGGER_VAR_MODS), **(stat_mods or {}), **(var_mods or {})}


    if tmt_labeled != "":
        unimod_tag = c.TMT_MODS[tmt_labeled]
        logger.info("Adding TMT fixed modifications")
        mods = {**{"C": "C[UNIMOD:4]", r"n[\d+]": f"{unimod_tag}-", "K": f"K{unimod_tag}"}, **mods}
        df["MODIFIED_SEQUENCE"] = msfragger_or_custom_to_internal(
            df["modified_peptide"].to_list(), mods=mods)
    else:
        #By default, i.e. if nothing is supplied to fixed_mods, carbamidomethylation on cystein will be included
        # in the fixed modifications. If you want to have no fixed modifictions at all, supply fixed_mods={}
        mods = {**{"C": "C[UNIMOD:4]"}, **mods}
        df["MODIFIED_SEQUENCE"] = msfragger_or_custom_to_internal(df["modified_peptide"].to_list(), mods=mods)

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
    df["PROTEINS"] = df["PROTEINS"].apply(lambda x: ";".join(x))

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
            "PROTEINS",
        ]
    ]
