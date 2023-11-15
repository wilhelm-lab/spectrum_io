
import logging
import re
import pandas as pd
import spectrum_fundamentals.constants as c
from .search_results import SearchResults
import os
from spectrum_io.spectral_library import digest 
import glob
import numpy as np

#from search_results import SearchResults

logger = logging.getLogger(__name__)

class Xisearch(SearchResults):
    """Handle search results from xisearch."""

    @staticmethod
    def read_result(path: str, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a csv of CSMs and perform some basic formatting.

        :param path: path to msms.csv to read
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading msms.csv file")
        columns_to_read = ["Run",
                        "Scan",
                        "PrecursorMass",
                        "PrecoursorCharge",
                        "decoy",
                        "Crosslinker",
                        "decoy",
                        "Protein1decoy",
                        "BasePeptide1",
                        "LengthPeptide1",
                        "Link1",
                        "Linked AminoAcid 1",
                        "Modifications1",
                        "ModificationPositions1",
                        "Protein2decoy",
                        "BasePeptide2",
                        "LengthPeptide2",
                        "Link2",
                        "Linked AminoAcid 2",
                        "Modifications2",
                        "ModificationPositions2",
                        "match score"] 

        # I must update codes and remmove these three lines
        path = str(path)  
        if path.endswith(".txt"):
            path = path[:-4] + ".csv"
        
        df = pd.read_csv(path, usecols=columns_to_read)
        logger.info("Finished reading msms.csv file")
        # Standardize column names
        df = Xisearch.filter_xisearch_result(df)
        df = Xisearch.update_columns_for_prosit(df)
        df = Xisearch.filter_valid_prosit_sequences(df)
        df.to_csv('/cmnfs/home/m.kalhor/wilhelmlab/notebooks/notebooks/xiserach_fdr/df.csv', index=False)
        return df

    def filter_xisearch_result (df: pd.DataFrame) -> pd.DataFrame:
        """
        remove unsupported modifications and keep only k-k as linked amino acid .

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """
        df = df[df['Linked AminoAcid 1'].notna() & df['Linked AminoAcid 1'].str.contains('K')]
        df = df[df['Linked AminoAcid 2'].notna() & df['Linked AminoAcid 2'].str.contains('K')]
        df = df[df['Modifications1'].astype(str).str.strip().isin(['Mox', 'Mox;Mox', 'Mox;Mox;Mox']) | pd.isna(df['Modifications1'])]
        df = df[df['Modifications2'].astype(str).str.strip().isin(['Mox', 'Mox;Mox', 'Mox;Mox;Mox']) | pd.isna(df['Modifications2'])]

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
        #print(mod_a_positions)
        #print(mod_a)

        if mod_a_positions not in ["nan", "null"]:
            if ";" in mod_a_positions:
                split_pos_mod_a = [int(num) for num in mod_a_positions.split(";")]
                split_mod_a = [str(mod) for mod in mod_a.split(";")]
                for index, pos_a in enumerate(split_pos_mod_a):
                    if split_mod_a[index] == "Mox":
                        modification = "M[UNIMOD:35]"
                        pos_mod_a = int(pos_a)
                        split_seq_a[pos_mod_a-1] = modification
            else:
                split_seq_a[int(mod_a_positions)-1] = "M[UNIMOD:35]"

        if mod_b_positions not in ["nan", "null"]:
            if ";" in mod_b_positions:
                split_pos_mod_b = [int(num) for num in mod_b_positions.split(";")]
                split_mod_b = [str(mod) for mod in mod_b.split(";")]
                for index, pos_b in enumerate(split_pos_mod_b):
                    if split_mod_b[index] == "Mox":
                        modification = "M[UNIMOD:35]"
                        pos_mod_b = int(pos_b)
                        split_seq_b[pos_mod_b-1] = modification
            else:
                split_seq_b[int(mod_b_positions)-1] = "M[UNIMOD:35]"
                
    
        if "C" in split_seq_a:
            c_index_pep_a = split_seq_a.index('C')
            split_seq_a[c_index_pep_a] = 'C[UNIMOD:4]'

        if "C" in split_seq_b:
            c_index_pep_b = split_seq_b.index('C')
            split_seq_b[c_index_pep_b] = 'C[UNIMOD:4]'
        
           
        split_seq_a[int(crosslinker_position_a)-1] = "K[UNIMOD:1896]"
        split_seq_b[int(crosslinker_position_b)-1] = "K[UNIMOD:1896]"

        seq_mod_a = ''.join(split_seq_a)    
        seq_mod_b = ''.join(split_seq_b)      

        return seq_mod_a, seq_mod_b


    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :return: modified df as pd.DataFrame
        """
        df.rename(columns={"Run": "RAW_FILE", 
                           "PrecursorMass": "MASS", #Experimental Mass of crosslinked peptides
                           "PrecoursorCharge": "PRECURSOR_CHARGE",
                           "Crosslinker": "CROSSLINKER_TYPE",
                           "match score":"SCORE",
                           "decoy": "REVERSE",
                           "Scan": "SCAN_NUMBER",
                           "BasePeptide1": "SEQUENCE_A",
                           "BasePeptide2": "SEQUENCE_B",
                           "Modifications1": "Modifications_A",
                           "Modifications2": "Modifications_B",
                           "Link1": "CROSSLINKER_POSITION_A",
                           "Link2": "CROSSLINKER_POSITION_B",
                           "LengthPeptide1": "PEPTIDE_LENGTH_A",
                           "LengthPeptide2": "PEPTIDE_LENGTH_B"},
                             inplace=True)
        
        logger.info("Converting xisearch peptide sequence to internal format")

                                    
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





