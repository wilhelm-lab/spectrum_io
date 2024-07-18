import logging
from pathlib import Path
from typing import Optional, Union, Dict, Tuple

import pandas as pd
from spectrum_fundamentals.constants import MOD_MASSES_SAGE
from spectrum_fundamentals.mod_string import sage_to_internal

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class Sage(SearchResults):
    """Handle search results from Sage."""

    def read_result(self, tmt_labeled: str = "", stat_mods: Optional[Dict[str, str]] = None, 
                    var_mods: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Function to read a msms tsv and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :param var_mods: Variable modifications with custom identifiers and their respective internal equivalents
        :param stat_mods: Static modifications with custom identifiers and their respective internal equivalents
        :return: pd.DataFrame with the formatted data
        """
        logger.info(f"Reading {self.path}")
        df = pd.read_csv(
            self.path,
            usecols=["filename", "scannr", "peptide", "charge", "hyperscore", "calcmass", "label", "proteins"],
            sep="\t",
        )
        logger.info(f"Finished reading {self.path}")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = Sage.update_columns_for_prosit(df, tmt_labeled, stat_mods=stat_mods, var_mods=var_mods)
        return filter_valid_prosit_sequences(df)

    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame, tmt_labeled: str, stat_mods: Optional[Dict[str, str]] = None, 
                                  var_mods: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :param tmt_labeled: True if tmt labeled, ignored
        :param var_mods: Variable modifications with custom identifiers and their respective internal equivalents
        :param stat_mods: Static modifications with custom identifiers and their respective internal equivalents
        :return: modified df as pd.DataFrame
        """
        df = df.rename(
            columns={
                "FILENAME": "RAW_FILE",
                "SCANNR": "SCAN_NUMBER",
                "PEPTIDE": "MODIFIED_SEQUENCE",
                "CHARGE": "PRECURSOR_CHARGE",
                "CALCMASS": "MASS",
                "HYPERSCORE": "SCORE",
                "LABEL": "REVERSE",
            }
        )
        mods = {**(MOD_MASSES_SAGE), **(stat_mods or {}), **(var_mods or {})}

        # removing .mzML
        df["RAW_FILE"] = df["RAW_FILE"].str.replace(r"\.mz[M|m][l|L]", "", regex=True)
        # extracting only the scan number
        df["SCAN_NUMBER"] = [int(x.rsplit("=", 1)[-1]) for x in df["SCAN_NUMBER"]]
        # creating a column of decoys and targets
        df["REVERSE"] = df["REVERSE"] < 0
        # removing modification to create the unmodified sequences
        df["SEQUENCE"] = df["MODIFIED_SEQUENCE"].str.replace(r"\-|\[.*?\]", "", regex=True)
        # length of the peptide
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].str.len()
        # converting sage to unimod
        df["MODIFIED_SEQUENCE"] = sage_to_internal(df["MODIFIED_SEQUENCE"], mods=mods)
        df["PROTEINS"].fillna("UNKNOWN", inplace=True)

        return df
