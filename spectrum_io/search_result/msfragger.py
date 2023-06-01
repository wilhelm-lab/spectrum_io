import logging
from pathlib import Path
from typing import Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    @staticmethod
    def read_result(path: Union[str, Path], tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param path: path to msms.txt to read
        :param tmt_labeled: tmt label as str
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading msfragger xlsx file")
        df = pd.read_excel(
            path,
            usecols=lambda x: x.upper()
            in [
                "SCANID",
                "PEPTIDE SEQUENCE",
                "PRECURSOR CHARGE",
                "PRECURSOR NEUTRAL MASS (DA)",
                "HYPERSCORE",
                "PROTEIN",
                "VARIABLE MODIFICATIONS DETECTED (STARTS WITH M, SEPARATED BY |, FORMATED AS POSITION,MASS)",
            ],
        )
        logger.info("Finished reading msfragger xlsx file")

        df.rename(
            columns={
                "ScanID": "SCAN NUMBER",
                "Peptide Sequence": "MODIFIED SEQUENCE",
                "Precursor neutral mass (Da)": "MASS",
                "Hyperscore": "SCORE",
                "Variable modifications detected (starts with M, separated by |, formated as position,mass)": "MODIFICATIONS",
            },
            inplace=True,
        )

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = MSFragger.update_columns_for_prosit(df, tmt_labeled)
        return filter_valid_prosit_sequences(df)

    @staticmethod
    def update_columns_for_prosit(df, tmt_labeled: str) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :param tmt_labeled: True if tmt labeled
        :return: modified df as pd.DataFrame
        """
        df.rename(columns={"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)

        df["REVERSE"] = df["PROTEIN"].str.contains("Reverse")
        # df["RAW_FILE"] = df.iloc[0]["PROTEIN"]
        df["RAW_FILE"] = "01625b_GA6-TUM_first_pool_41_01_01-DDA-1h-R2"
        logger.info("Converting MSFragger  peptide sequence to internal format")

        mod_masses_reverse = {round(float(v), 3): k for k, v in c.MOD_MASSES.items()}
        sequences = []
        for _, row in df.iterrows():
            modifications = row["MODIFICATIONS"].split("|")[1:]
            if len(modifications) == 0:
                sequences.append(row["MODIFIED_SEQUENCE"])
            else:
                sequence = row["MODIFIED_SEQUENCE"]
                skip = 0
                for mod in modifications:
                    pos, mass = mod.split("$")
                    sequence = (
                        sequence[: int(pos) + 1 + skip]
                        + mod_masses_reverse[round(float(mass), 3)]
                        + sequence[int(pos) + 1 + skip :]
                    )
                    skip = skip + len(mod_masses_reverse[round(float(mass), 3)])
                sequences.append(sequence)

        df["MODIFIED_SEQUENCE"] = sequences

        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].apply(lambda x: len(x))

        return df
