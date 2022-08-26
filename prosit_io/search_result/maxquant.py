import pandas as pd
import numpy as np
import logging

from .search_results import SearchResults
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import fundamentals.constants as C

logger = logging.getLogger(__name__)


class MaxQuant(SearchResults):

    @staticmethod
    def add_tmt_mod(mass, seq, unimod_tag):
        num_of_tmt = seq.count(unimod_tag)
        mass += (num_of_tmt * C.MOD_MASSES[f'{unimod_tag}'])
        return mass

    @staticmethod
    def read_result(path: str, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting
        :prarm path: Path to msms.txt to read
        :return: DataFrame
        """
        logger.info("Reading msms.txt file")
        df = pd.read_csv(path,
                         usecols=lambda x: x.upper() in ['RAW FILE',
                                                         'SCAN NUMBER',
                                                         'MODIFIED SEQUENCE',
                                                         'CHARGE',
                                                         'FRAGMENTATION',
                                                         'MASS ANALYZER',
                                                         'SCAN EVENT NUMBER',
                                                         'LABELING STATE',
                                                         'MASS', # = Calculated Precursor mass; TODO get column with experimental Precursor mass instead
                                                         'SCORE',
                                                         'REVERSE',
                                                         'RETENTION TIME'],
                         sep="\t")
        logger.info("Finished reading msms.txt file")
        
        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = MaxQuant.update_columns_for_prosit(df, tmt_labeled)
        df = MaxQuant.filter_prosit_valid_sequences(df)
        
        return df
    
    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame, tmt_labeled: str) -> pd.DataFrame:    
        df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)

        if "MASS_ANALYZER" not in df.columns:
            df['MASS_ANALYZER'] = 'FTMS'
        if "FRAGMENTATION" not in df.columns:
            df['FRAGMENTATION'] = 'HCD'

        df["REVERSE"].fillna(False, inplace=True)
        df["REVERSE"].replace("+", True, inplace=True)
        logger.info("Converting MaxQuant peptide sequence to internal format")
        if tmt_labeled != "":
            unimod_tag = C.TMT_MODS[tmt_labeled]
            logger.info("Adding TMT fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                           '^_':f'_{unimod_tag}', 
                                                                                                           'K': f'K{unimod_tag}'})
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, unimod_tag), axis=1)
            if 'msa' in tmt_labeled:
                logger.info("Replacing phospho by dehydration for Phospho-MSA")
                df["MODIFIED_SEQUENCE_MSA"] = df["MODIFIED_SEQUENCE"].str.replace("[UNIMOD:21]", "[UNIMOD:23]", regex=False)
        elif "LABELING_STATE" in df.columns:
            logger.info("Adding SILAC fixed modifications")
            df.loc[df['LABELING_STATE'] == 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(df[df['LABELING_STATE'] == 1]["MODIFIED_SEQUENCE"].to_numpy(), 
                                                                                          fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                      'K': 'K[UNIMOD:259]', 
                                                                                                      'R': 'R[UNIMOD:267]'})
            df.loc[df['LABELING_STATE'] != 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(df[df['LABELING_STATE'] != 1]["MODIFIED_SEQUENCE"].to_numpy())
            df.drop(columns=['LABELING_STATE'], inplace=True)
        else:
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy())
        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        df['PEPTIDE_LENGTH'] = df["SEQUENCE"].apply(lambda x: len(x))
        
        return df
    
    @staticmethod
    def filter_prosit_valid_sequences(df: pd.DataFrame) -> pd.DataFrame:
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(df.index)}")
        df = df[(df['PEPTIDE_LENGTH'] <= 30)]
        df = df[(~df['MODIFIED_SEQUENCE'].str.contains('\(ac\)'))]
        df = df[
            (~df['MODIFIED_SEQUENCE'].str.contains('\(Acetyl \(Protein N-term\)\)'))]
        df = df[(~df['SEQUENCE'].str.contains('U'))]
        df = df[df['PRECURSOR_CHARGE'] <= 6]
        df = df[df['PEPTIDE_LENGTH'] >= 7]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(df.index)}")

        return df

