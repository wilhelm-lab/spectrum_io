import pandas as pd
import numpy as np
import logging
import sqlite3

from .search_results import SearchResults
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import fundamentals.constants as C

logger = logging.getLogger(__name__)


class Mascot(SearchResults):

    @staticmethod
    def read_result(path: str, tmt_labeled):
        """
        Function to read a msms txt and perform some basic formatting
        :prarm path: Path to msms.txt to read
        :return: DataFrame
        """
        logger.info("Reading mascot msf file")
        connection = sqlite3.connect(path)
        cursor = connection.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        df = pd.read_sql(f"SELECT * FROM MSnSpectrumInfo", connection)[["SpectrumID", "SpectrumFileName", "RetentionTime", "Mass", "Charge"]]
        df_id_map = pd.read_sql(f"SELECT * FROM TargetPsmsMSnSpectrumInfo", connection)[["MSnSpectrumInfoSpectrumID", "TargetPsmsPeptideID"]]
        df = df.merge(df_id_map, left_on = "SpectrumID", right_on = "MSnSpectrumInfoSpectrumID")
        df_target_psms = pd.read_sql(f"SELECT * FROM TargetPsms", connection)[["PeptideID", "Sequence", "ModifiedSequence", "Modifications", "XCorr"]]
        df = df.merge(df_target_psms, left_on="TargetPsmsPeptideID", right_on="PeptideID")
        df_modif = pd.read_sql(f"SELECT * FROM TargetPsmsFoundModifications", connection)[["TargetPsmsPeptideID", "FoundModificationsModificationID", "Position"]]
        df_modif_mass = pd.read_sql(f"SELECT * FROM FoundModifications", connection)[["ModificationID", "DeltaMonoisotopicMass"]]
        df_modif = df_modif.merge(df_modif_mass, left_on="FoundModificationsModificationID", right_on="ModificationID")
        df = df.merge(df_modif, on="TargetPsmsPeptideID")

        logger.info("Finished reading mascot msf file")

        df.rename(columns = {"SpectrumID": "SCAN NUMBER", "ModifiedSequence": "MODIFIED SEQUENCE", "Charge": "PRECURSOR CHARGE",
                             "XCorr": "SCORE", "RetentionTime": "RETENTION TIME", "SpectrumFileName": "RAW FILE"}, inplace=True)
        
        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")
        #TODO reverse
        df["REVERSE"] = df["SEQUENCE"].str.contains("Reverse") 
        logger.info("Converting MSFragger  peptide sequence to internal format")
        df["RAW_FILE"] = df["RAW_FILE"].str.replace(".raw", "")
        df["MODIFICATIONS"] = (df["POSITION"].astype(int) - 1).astype(str) + "$" + df["DELTAMONOISOTOPICMASS"].astype(str)
        df = df.groupby('SCAN_NUMBER', as_index=False).apply(lambda x: x.sort_values('POSITION'))
        df = df.groupby(["SCAN_NUMBER", "PRECURSOR_CHARGE", "SCORE", "RETENTION_TIME", "RAW_FILE", "SEQUENCE", "REVERSE"], as_index=False).agg({"MODIFICATIONS": "|".join})
        mod_masses_reverse = {round(float(v), 3): k for k, v in C.MOD_MASSES.items()}

        sequences = []
        for index, row in df.iterrows():
            modifications = row["MODIFICATIONS"].split("|")
            if len(modifications) == 0:
                sequences.append(row["SEQUENCE"])
            else:
                sequence = row["SEQUENCE"]
                skip = 0
                for mod in modifications:
                    pos, mass = mod.split("$")
                    sequence = sequence[:int(pos)+1+skip] + mod_masses_reverse[round(float(mass), 3)] + sequence[int(pos)+1+skip:]
                    skip = skip + len(mod_masses_reverse[round(float(mass), 3)])
                sequences.append(sequence)

        df["MODIFIED_SEQUENCE"] = sequences

        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        df['PEPTIDE_LENGTH'] = df["SEQUENCE"].apply(lambda x: len(x))
        logger.info(f"No of sequences before Filtering is {len(df['PEPTIDE_LENGTH'])}")
        df = df[(df['PEPTIDE_LENGTH'] <= 30)]
        df = df[(~df['MODIFIED_SEQUENCE'].str.contains('\(ac\)'))]
        df = df[
            (~df['MODIFIED_SEQUENCE'].str.contains('\(Acetyl \(Protein N-term\)\)'))]
        df = df[(~df['SEQUENCE'].str.contains('U'))]
        df = df[df['PRECURSOR_CHARGE'] <= 6]
        df = df[df['PEPTIDE_LENGTH'] >= 7]
        logger.info(f"No of sequences after Filtering is {len(df['PEPTIDE_LENGTH'])}")
        return df
