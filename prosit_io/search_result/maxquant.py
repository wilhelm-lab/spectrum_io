import pandas as pd
import numpy as np
from typing import List
from .search_results import SearchResults
from fundamentals import ALPHABET_MODS, MAXQUANT_VAR_MODS


class MaxQuant(SearchResults):

    @staticmethod
    def read_result(path):
        """
        Function to read a msms txt and perform some basic formatting
        :prarm path: Path to msms.txt to read
        :return: DataFrame
        """
        df = pd.read_csv(path,
                         usecols=lambda x: x.upper() in['RAW FILE',
                                  'SCAN NUMBER',
                                  'MODIFIED SEQUENCE',
                                  'CHARGE',
                                  'FRAGMENTATION',
                                  'MASS ANALYZER',
                                  'MASS',
                                  'SCORE',
                                  'REVERSE'],
                         sep = "\t")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")


        df["REVERSE"].fillna(False,inplace=True)
        df["REVERSE"].replace("+", True, inplace=True)
        df["MODIFIED_SEQUENCE"] = MaxQuant.transform_maxquant_mod_string(df["MODIFIED_SEQUENCE"])
        return df

    @staticmethod
    def transform_maxquant_mod_string(
        sequences: np.ndarray,
        fixed_mods: List[str] = [],
    ):
        """
        Function to translate a MaxQuant modstring to the Prosit format
        :param str: List, np.array or pd.Series of String sequences

        :return: pd.Series
        """
        assert all(x in ALPHABET_MODS.keys() for x in fixed_mods), f"Provided illegal fixed mod, supported modifications are {str(ALPHABET_MODS.keys())}."
        modified_sequences = []
        for seq in sequences:
            seq = seq.replace("_", "")
            for mod, unimode in MAXQUANT_VAR_MODS.items():
                seq = seq.replace(mod, unimode)
            for mod in fixed_mods:
                aa, remainder = mod.split("(")
                unimode = remainder.split(')')[0]
                seq = seq.replace(aa, mod)
            modified_sequences.append(seq)
        return modified_sequences
