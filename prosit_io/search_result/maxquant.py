import pandas as pd
import numpy as np
from .search_results import SearchResults
from fundamentals.mod_string import maxquant_to_internal


class MaxQuant(SearchResults):

    @staticmethod
    def read_result(path: str):
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
        df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy())
        return df
