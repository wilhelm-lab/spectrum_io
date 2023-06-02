import logging
import re
import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Plink2(SearchResults):
    """Handle search results from Plink2."""

    @staticmethod
    def read_result(path: str) -> pd.DataFrame:
        """
        Function to read a csv_file and perform some basic formatting.

        :param path: path to csv_file to read
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading csv_file file")
        df = pd.read_csv(
            path,
            usecols=lambda x: x.upper()
            in [
                "Title",
                "Charge",
                "Precursor_Mass",
                "Peptide",
                "Linker"
                "Modifications"  
                "Score",
                "LabelID",
            ],
            sep="\t",
        )
        logger.info("Finished reading csv_file file")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = Plink2.update_columns_for_prosit(df)
        df = Plink2.filter_valid_prosit_sequences(df)

        return df

    def add_mod_sequence(seq_a: str,
                         seq_b: str,
                         mod: str,
                         peptide_a_length: int,
                         crosslinker_position_a: int,
                         crosslinker_position_b: int):
        """
        Function adds modification in peptide sequence for xl-prosit 
        
        :seq_a: unmodified peptide a
        :seq_b: unmodified peptide b
        :mod: all modifications
        :crosslinker_position_a: crosslinker position of peptide a
        :crosslinker_position_b: crosslinker position of peptide b
        :return: modified sequence a and b
        """
        split_seq_a = [x for x in seq_a]
        split_seq_b = [x for x in seq_b]

        for mod in mod.split(";"):
            if mod.startswith("Oxidation"):
                modification = "M[UNIMOD:35]"
                pos_mod = int(mod.split("(")[-1][:-1])
                if peptide_a_length >= pos_mod :
                    split_seq_a[pos_mod-1] = modification
                else:
                    pos_mod = pos_mod - peptide_a_length - 3
                    split_seq_b[pos_mod-1] = modification
                
            elif mod.startswith("Carbamidomethyl"):
                modification = "C[UNIMOD:4]"
                pos_mod = int(mod.split("(")[-1][:-1])
                if pos_mod <= peptide_a_length :
                    split_seq_a[pos_mod-1] = modification
                else:
                    pos_mod = pos_mod - peptide_a_length - 3
                    split_seq_b[pos_mod-1] = modification
            elif mod == "null":
                break
            else:
                raise AssertionError(f"unknown modification provided:{mod}")

        split_seq_a[crosslinker_position_a-1] = "K[UNIMOD:1896]"
        split_seq_b[crosslinker_position_b-1] = "K[UNIMOD:1896]"
     
        seq_mod_a = ''.join(split_seq_a)    
        seq_mod_b = ''.join(split_seq_b)      

        return seq_mod_a, seq_mod_a


    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :return: modified df as pd.DataFrame
        """
        df.rename(columns={"Title": "RAW_FILE", 
                           "Precursor_Mass": "MASS" #Experimental Mass of crosslinked peptides
                           "Charge": "PRECURSOR_CHARGE",
                           "Peptide": "MODIFIED_SEQUENCE",
                           "Linker": "crosslinker_type",
                           "Score":"SCORE",
                           "LabelID": "REVERSE"},
                             inplace=True)

        if "MASS_ANALYZER" not in df.columns:
            df["MASS_ANALYZER"] = "FTMS"

        if "FRAGMENTATION" not in df.columns:
            df["FRAGMENTATION"] = "HCD"

        df['RAW_FILE'] = df['RAW_FILE'].str.split('.').str[0]
        df['SCAN_NUMBER'] = df['RAW_FILE'].str.split('.').str[1]

        logger.info("Converting Plink2 peptide sequence to internal format")
        df['SEQUENCE_A'] = df['crosslink'].apply(lambda x: re.split(r'[\(\)-]', x)[0])
        df['SEQUENCE_B'] = df['crosslink'].apply(lambda x: re.split(r'[\(\)-]', x)[3])
        df['crosslinker_position_A'] = df['crosslink'].apply(lambda x: re.split(r'[\(\)-]', x)[1])
        df['crosslinker_position_B'] = df['crosslink'].apply(lambda x: re.split(r'[\(\)-]', x)[4])

        df["PEPTIDE_LENGTH_A"] = df['SEQUENCE_A'].apply(len) 
        df["PEPTIDE_LENGTH_B"] = df['SEQUENCE_B'].apply(len)

        logger.info("Converting MaxQuant peptide sequence to internal format")
        df[['MODIFIED_SEQUENCE_A','MODIFIED_SEQUENCE_B']] = df.apply(lambda row: Plink2.add_mod_sequence(row['SEQUENCE_A'], 
                                                                                 row['SEQUENCE_B'],
                                                                                 row['PEPTIDE_LENGTH_A'],
                                                                                 row['crosslinker_position_A'],
                                                                                 row['crosslinker_position_B']), axis=1)
        
        #df["REVERSE"].fillna(False, inplace=True)
        df["REVERSE"].replace("1", True, inplace=True)
        return df

    @staticmethod
    def filter_valid_prosit_sequences(df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter valid Prosit sequences.

        :param df: df to filter
        :return: df after filtration
        """
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(df.index)}")

        df = df[(df["PEPTIDE_LENGTH_A"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_A"] >= 6]
        df = df[(df["PEPTIDE_LENGTH_B"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_B"] >= 6]
        df = df[(~df["SEQUENCE_A"].str.contains("U"))]
        df = df[(~df["SEQUENCE_B"].str.contains("U"))]
        df = df[df["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(df.index)}")

        return df



