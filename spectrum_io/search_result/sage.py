import logging
from pathlib import Path
from typing import Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import sage_to_internal

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class Sage(SearchResults):
    """Handle search results from Sage."""

    @staticmethod
    def read_result(path: Union[str, Path], tmt_labeled: str = "") -> pd.DataFrame:
        """
        Function to read a msms tsv and perform some basic formatting.

        :param path: path to msms.tsv to read
        :param tmt_labeled: tmt label as str
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading msms.tsv file")
        df = pd.read_csv(
            path,
            usecols=["filename", "scannr", "peptide", "charge", "hyperscore", "calcmass", "proteins"],
            sep="\t",
        )
        logger.info("Finished reading msms.tsv file")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = Sage.update_columns_for_prosit(df, tmt_labeled)
        return filter_valid_prosit_sequences(df)

    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame, tmt_labeled: str) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :param tmt_labeled: True if tmt labeled, ignored
        :return: modified df as pd.DataFrame
        """
        df = df.rename(
            columns={
                "FILENAME": "RAW_FILE",
                "SCANNR": "SCAN_NUMBER",
                "PEPTIDE": "MODIFIED_SEQUENCE",
                "CHARGE": "PRECURSOR_CHARGE",
            }
        )

        # removing .mzML
        df["RAW_FILE"] = df["RAW_FILE"].str.replace(".mzML", "", regex=True)
        # extracting only the scan number
        df["SCAN_NUMBER"] = [int(x.rsplit("=", 1)[-1]) for x in df["SCAN_NUMBER"]]
        # creating a column of decoys and targets
        df["REVERSE"] = df["PROTEINS"].str.startswith("rev_")
        # removing modification to create the unmodified sequences
        df["SEQUENCE"] = df["MODIFIED_SEQUENCE"].str.replace(r"\[.*?\]", "", regex=True)
        # length of the peptide
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].str.len()
        # mass of the peptide
        df["MASS"] = df["CALCMASS"]
        # score of the peptide
        df["SCORE"] = df["HYPERSCORE"]
        # converting proforma to unimode
        print(df)
        df["MODIFIED_SEQUENCE"] = sage_to_internal(df["MODIFIED_SEQUENCE"])

        print(df.columns)
        return df
