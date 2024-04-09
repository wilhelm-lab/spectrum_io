import logging
import re
import pandas as pd
import spectrum_fundamentals.constants as c
from .search_results import SearchResults
import os
import glob
import numpy as np


logger = logging.getLogger(__name__)

class Xisearch(SearchResults):
    """Handle search results from xisearch."""

    def read_result(self, tmt_labeled: str = "") -> pd.DataFrame:
        """
        Function to read a csv of CSMs and perform some basic formatting.

        :param path: path to msms.csv to read
        :return: pd.DataFrame with the formatted data
        """
        
        logger.info("Reading msms.csv file")
        columns_to_read = ["run_name",
                        "scan_number",
                        "precursor_mass",
                        "precursor_charge",
                        "crosslinker_name",
                        "decoy_p1",
                        "base_sequence_p1",
                        "aa_len_p1",
                        "link_pos_p1",
                        "linked_aa_p1",
                        "mods_p1",
                        "mod_pos_p1",
                        "decoy_p2",
                        "base_sequence_p2",
                        "aa_len_p2",
                        "link_pos_p2",
                        "linked_aa_p2",
                        "mods_p2",
                        "mod_pos_p2",
                        "linear",
                        "match_score"] 
        
        # Initialize path variable
        path = self.path

        if str(self.path).endswith(".txt"):
            path = self.path.with_suffix('.tsv')
        df = pd.read_csv(path,  sep='\t', usecols= columns_to_read)
        logger.info("Finished reading msms.tsv file")
        # Standardize column names
        df = Xisearch.filter_xisearch_result(df)
        df = Xisearch.update_columns_for_prosit(df)
        df = Xisearch.filter_valid_prosit_sequences(df)
        df.to_csv("/cmnfs/home/m.kalhor/wilhelmlab/spectrum_io/tests/unit_tests/data/xisearch_output_internal.tsv",  index=False)
        return df

    def filter_xisearch_result (df: pd.DataFrame) -> pd.DataFrame:
        """
        remove unsupported modifications and keep only k-k as linked amino acid .

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """ 
        df = df[df['linear'] != True]
        df = df[df['linked_aa_p1'].notna() & df['linked_aa_p1'].str.contains('K')]
        df = df[df['linked_aa_p2'].notna() & df['linked_aa_p2'].str.contains('K')]
        df = df[~df['mods_p1'].str.contains('dsso-hyd', na=False)]
        df = df[~df['mods_p2'].str.contains('dsso-hyd', na=False)]
        valid_modifications = ['cm', 'ox', pd.NA]
        df = df[df['mods_p1'].apply(lambda x: any(mod in str(x).split(';') if pd.notnull(x) else mod is pd.NA for mod in valid_modifications))]
        df = df[df['mods_p2'].apply(lambda x: any(mod in str(x).split(';') if pd.notnull(x) else mod is pd.NA for mod in valid_modifications))]
        
        return df
    
    def add_mod_sequence(seq_a: str,
                         seq_b: str,
                         mod_a: str,
                         mod_b: str,
                         crosslinker_position_a: int,
                         crosslinker_position_b: int,
                         mod_a_positions: str,
                         mod_b_positions: str
                         ):
        
        """
        Function adds modification in peptide sequence for xl-prosit 
        
        :seq_a: unmodified peptide a
        :seq_b: unmodified peptide b
        :mod_a: all modifications of pep a
        :mod_b: all modifications of pep b
        :crosslinker_position_a: crosslinker position of peptide a
        :crosslinker_position_b: crosslinker position of peptide b
        :crosslinker_type: crosslinker tpe eg. DSSO, DSBU
        :mod_a_positions: position of all modifications of peptide a
        :mod_b_positions: position of all modifications of peptide b
        :return: modified sequence a and b
        """
        
        split_seq_a = [x for x in seq_a]
        split_seq_b = [x for x in seq_b]
        mod_a_positions = str(mod_a_positions)
        mod_b_positions = str(mod_b_positions)
        

        if mod_a_positions not in ["nan", "null"]:
            if ";" in mod_a_positions:
                split_pos_mod_a = [int(num) for num in mod_a_positions.split(";")]
                split_mod_a = [str(mod) for mod in mod_a.split(";")]
                for index, pos_a in enumerate(split_pos_mod_a):
                    if split_mod_a[index] == "ox":
                        modification = "M[UNIMOD:35]"
                    if split_mod_a[index] == "cm":
                        modification = "C[UNIMOD:4]"
                    pos_mod_a = int(pos_a)
                    split_seq_a[pos_mod_a-1] = modification
            else:
                if mod_a == "ox" :
                    modification = "M[UNIMOD:35]"
                if mod_a == "cm" :
                    modification = "C[UNIMOD:4]"
                try:
                    mod_a_positions_float = float(mod_a_positions)
                    split_seq_a[int(mod_a_positions_float)-1] = modification
                except ValueError:
                    print(f"Error occurred with mod_a_positions value: {mod_a_positions}")

        if mod_b_positions not in ["nan", "null"]:
            if ";" in mod_b_positions:
                split_pos_mod_b = [int(num) for num in mod_b_positions.split(";")]
                split_mod_b = [str(mod) for mod in mod_b.split(";")]
                for index, pos_b in enumerate(split_pos_mod_b):
                    if split_mod_b[index] == "ox":
                        modification = "M[UNIMOD:35]"
                    if split_mod_b[index] == "cm":
                        modification = "C[UNIMOD:4]"
                    pos_mod_b = int(pos_b)
                    split_seq_b[pos_mod_b-1] = modification
            else:
                if mod_b == "ox" :
                    modification = "M[UNIMOD:35]"
                if mod_b == "cm" :
                    modification = "C[UNIMOD:4]"
                try:
                    mod_b_positions_float = float(mod_b_positions)
                    split_seq_b[int(mod_b_positions_float)-1] = modification
                except ValueError:
                    print(f"Error occurred with mod_a_positions value: {mod_b_positions}")
                  
        split_seq_a[int(crosslinker_position_a)-1] = "K[UNIMOD:1896]"
        split_seq_b[int(crosslinker_position_b)-1] = "K[UNIMOD:1896]"

        seq_mod_a = ''.join(split_seq_a)    
        seq_mod_b = ''.join(split_seq_b)      

        return seq_mod_a, seq_mod_b

    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame) -> pd.DataFrame:
        """
        Update columns of df to work with xl-prosit.

        :param df: df to modify
        :return: modified df as pd.DataFrame
        """

        df['decoy'] = df['decoy_p1'] | df['decoy_p2']
        df["RAW_FILE"] = df["run_name" ]
        df["MASS"] = df["precursor_mass"]
        df["PRECURSOR_CHARGE"] = df["precursor_charge"]
        df["CROSSLINKER_TYPE"] = df["crosslinker_name"]
        df["SCORE"] = df["match_score"]
        df["REVERSE"] = df["decoy"]
        df["SCAN_NUMBER"] = df["scan_number"]
        df["SEQUENCE_A"] = df["base_sequence_p1"]
        df["SEQUENCE_B"] = df["base_sequence_p2"]
        df["Modifications_A"] = df["mods_p1"]
        df["Modifications_B"] = df["mods_p2"]
        df["CROSSLINKER_POSITION_A"] = df["link_pos_p1"]
        df["CROSSLINKER_POSITION_B"] = df["link_pos_p2"]
        df["ModificationPositions1"] = df["mod_pos_p1"]
        df["ModificationPositions2"] = df["mod_pos_p2"]
        df["PEPTIDE_LENGTH_A"] = df["aa_len_p1"]
        df["PEPTIDE_LENGTH_B"] = df["aa_len_p2"]
        logger.info("Converting xisearch peptide sequence to internal format")

        df['RAW_FILE'] = df['RAW_FILE'].str.replace('.raw', '')                          
        df['Modifications_A'] = df['Modifications_A'].astype('str') 
        df['Modifications_B'] = df['Modifications_B'].astype('str')
        
        df['CROSSLINKER_POSITION_A'] = df['CROSSLINKER_POSITION_A'].astype('int')
        df['CROSSLINKER_POSITION_B'] = df['CROSSLINKER_POSITION_B'].astype('int')
        
          
        df[['MODIFIED_SEQUENCE_A','MODIFIED_SEQUENCE_B']] = df.apply(lambda row: Xisearch.add_mod_sequence(row['SEQUENCE_A'], 
                                                                                 row['SEQUENCE_B'],
                                                                                 row['Modifications_A'],
                                                                                 row['Modifications_B'],
                                                                                 row['CROSSLINKER_POSITION_A'],
                                                                                 row['CROSSLINKER_POSITION_B'],
                                                                                 row['ModificationPositions1'],
                                                                                 row['ModificationPositions2']), axis=1, result_type='expand')
        
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
