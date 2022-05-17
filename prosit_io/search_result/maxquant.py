import pandas as pd
import numpy as np
import logging

from .search_results import SearchResults
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods, parse_modstrings
import fundamentals.constants as C

logger = logging.getLogger(__name__)


class MaxQuant(SearchResults):

    @staticmethod
    def add_tmt_mod(mass, seq, tag):
        if tag == "tmt":
            num_of_tmt = seq.count('UNIMOD:737')
            mass += (num_of_tmt * C.MOD_MASSES['[UNIMOD:737]'])
        elif tag == "tmtpro":
            num_of_tmt = seq.count('UNIMOD:2016')
            mass += (num_of_tmt * C.MOD_MASSES['[UNIMOD:2016]'])
        elif tag == "itraq4":
            num_of_tmt = seq.count('UNIMOD:214')
            mass += (num_of_tmt * C.MOD_MASSES['[UNIMOD:214]'])
        elif tag == "itraq8":
            num_of_tmt = seq.count('UNIMOD:730')
            mass += (num_of_tmt * C.MOD_MASSES['[UNIMOD:730]'])
        return mass

    @staticmethod
    def read_result(path: str, tmt_labeled=''):
        """
        Function to read a msms txt and perform some basic formatting
        :prarm path: Path to msms.txt to read
        :return: DataFrame
        """
        logger.info("Reading msms.txt file " + path)
        df = pd.read_csv(path,
                         usecols=lambda x: x.upper() in ['RAW FILE',
                                                         'SCAN NUMBER',
                                                         'ALL MODIFIED SEQUENCES',
                                                         'MODIFIED SEQUENCE',
                                                         'CHARGE',
                                                         'FRAGMENTATION',
                                                         'MASS ANALYZER',
                                                         'SCAN EVENT NUMBER',
                                                         'LABELING STATE',
                                                         'MASS', # = Calculated Precursor mass; TODO get column with experimental Precursor mass instead
                                                         'SCORE',
                                                         'REVERSE',
                                                         'RETENTION TIME',
                                                         'COLLISION ENERGY',
                                                         'LOCALIZATION PROB'],
                         sep="\t")
        logger.info("Finished reading msms.txt file")
        
        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df.rename(columns = {"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)
        df.rename(columns = {"COLIISION_ENERGY": "ORIGINAL_COLLISION_ENERGY"}, inplace=True)
        #df['MODIFIED_SEQUENCE'] = df['ALL_MODIFIED_SEQUENCES'].apply(lambda x: x.split(';'))
        #df = df.explode('MODIFIED_SEQUENCE')
        if "MASS_ANALYZER" not in df.columns:
            df['MASS_ANALYZER'] = 'FTMS'
        if "FRAGMENTATION" not in df.columns:
            df['FRAGMENTATION'] = 'HCD'

        df["REVERSE"].fillna(False, inplace=True)
        df["REVERSE"].replace("+", True, inplace=True)
        #df = df[(~df['MODIFIED_SEQUENCE'].str.contains('(tm)'))]

        logger.info("Converting MaxQuant peptide sequence to internal format")
        if tmt_labeled == "tmt":
            logger.info("Adding TMT fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                           '^_(tm)':'_[UNIMOD:737]',
                                                                                                           'K(tm)': 'K[UNIMOD:737]'})
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, tmt_labeled), axis=1)
        elif tmt_labeled == "tmtpro":
            logger.info("Adding TMTpro fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                           '^_':'_[UNIMOD:2016]', 
                                                                                                           'K': 'K[UNIMOD:2016]'})
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, tmt_labeled), axis=1)
        elif tmt_labeled == "itraq4":
            logger.info("Adding iTRAQ4 fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                           '^_':'_[UNIMOD:214]', 
                                                                                                           'K': 'K[UNIMOD:214]'})
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, tmt_labeled), axis=1)
        elif tmt_labeled == "itraq8":
            logger.info("Adding iTRAQ8 fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                           '^_':'_[UNIMOD:730]', 
                                                                                                           'K': 'K[UNIMOD:730]'})
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, tmt_labeled), axis=1)
        elif "LABELING_STATE" in df.columns:
            logger.info("Adding SILAC fixed modifications")
            df.loc[df['LABELING_STATE'] == 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(df[df['LABELING_STATE'] == 1]["MODIFIED_SEQUENCE"].to_numpy(), 
                                                                                          fixed_mods={'C': 'C[UNIMOD:4]',
                                                                                                      'K': 'K[UNIMOD:259]', 
                                                                                                      'R': 'R[UNIMOD:267]'})
            df.loc[df['LABELING_STATE'] != 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(df[df['LABELING_STATE'] != 1]["MODIFIED_SEQUENCE"].to_numpy())
            df.drop(columns=['LABELING_STATE'], inplace=True)
        else:
            df["MODIFIED_SEQUENCE_INT"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy(), {}, False)
            #enum_char = enumerate(
             #   parse_modstrings(df["MODIFIED_SEQUENCE_INT"].to_numpy(), alphabet=C.ALPHABET, translate=True, filter=False))
            #sequences = []
            #for _, seq in enum_char:
             #   sequences.append(seq)
            #df["MODIFIED_SEQUENCE_INT"] = sequences
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy())
        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        df['PEPTIDE_LENGTH'] = df["SEQUENCE"].apply(lambda x: len(x))

        logger.info(f"No of sequences before Filtering is {len(df['PEPTIDE_LENGTH'])}")
        #df = df[(df['PEPTIDE_LENGTH'] <= 30)]

        #df = df[(~df['SEQUENCE'].str.contains('U'))]
        #df = df[df['PRECURSOR_CHARGE'] <= 6]
        #df = df[df['PEPTIDE_LENGTH'] >= 7]
        #df = df[df['SCORE'] > 70]
        logger.info(f"No of sequences after Filtering is {len(df['PEPTIDE_LENGTH'])}")

        return df

