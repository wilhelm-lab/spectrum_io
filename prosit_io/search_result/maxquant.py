import pandas as pd
import numpy as np
from .search_results import SearchResults
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import fundamentals.constants as C

class MaxQuant(SearchResults):

    @staticmethod
    def add_tmt_mod(mass, seq):
        num_of_tmt = seq.count('UNIMOD:737')
        mass += (num_of_tmt * C.MOD_MASSES['[UNIMOD:737]'])
        return mass

    @staticmethod
    def read_result(path: str):
        """
        Function to read a msms txt and perform some basic formatting
        :prarm path: Path to msms.txt to read
        :return: DataFrame
        """
        df = pd.read_csv(path,
                         usecols=lambda x: x.upper() in ['RAW FILE',
                                                         'SCAN NUMBER',
                                                         'MODIFIED SEQUENCE',
                                                         'CHARGE',
                                                         'FRAGMENTATION',
                                                         'MASS ANALYZER',
                                                         'MASS', # Experimental Precursor mass # TODO actually get column with experimental Precursor mass instead
                                                         'SCORE',
                                                         'REVERSE',
                                                         'RETENTION TIME'],
                         sep="\t")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)

        if "MASS_ANALYZER" not in df.columns:
            df['MASS_ANALYZER'] = 'FTMS'
        if "FRAGMENTATION" not in df.columns:
            df['FRAGMENTATION'] = 'HCD'

        #df['RETENTION_TIME'] = [x for x in range(len(df))]
        df["REVERSE"].fillna(False, inplace=True)
        df["REVERSE"].replace("+", True, inplace=True)
        df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy())
        df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE), axis=1)
        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        #Filter sequences remove sequences with length bigger than 30
        df['PEPTIDE_LENGTH'] = df["SEQUENCE"].apply(lambda x: len(x))
        df = df[df['PEPTIDE_LENGTH']<=30]
        return df

