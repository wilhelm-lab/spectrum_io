import logging
import re
import sqlite3
from pathlib import Path
from typing import Optional, Union, Dict, Tuple

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class Mascot(SearchResults):
    """Handle search results from Mascot."""

    def read_result(self, tmt_labeled: str, custom_stat_mods: Dict[str, Tuple[str, float]] = None, 
                    custom_var_mods: Dict[str, Tuple[str, float]] = None) -> pd.DataFrame:
        """
        Function to read a mascot msf file and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading mascot msf file")
        connection = sqlite3.connect(self.path)
        # cursor = connection.cursor()
        # cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        df = pd.read_sql("SELECT * FROM MSnSpectrumInfo", connection)[
            ["SpectrumID", "SpectrumFileName", "RetentionTime", "Mass", "Charge"]
        ]
        df_id_map = pd.read_sql("SELECT * FROM TargetPsmsMSnSpectrumInfo", connection)[
            ["MSnSpectrumInfoSpectrumID", "TargetPsmsPeptideID"]
        ]
        df = df.merge(df_id_map, left_on="SpectrumID", right_on="MSnSpectrumInfoSpectrumID")
        df_target_psms = pd.read_sql("SELECT * FROM TargetPsms", connection)[
            ["PeptideID", "Sequence", "ModifiedSequence", "Modifications", "XCorr"]
        ]
        df = df.merge(df_target_psms, left_on="TargetPsmsPeptideID", right_on="PeptideID")
        df_modif = pd.read_sql("SELECT * FROM TargetPsmsFoundModifications", connection)[
            ["TargetPsmsPeptideID", "FoundModificationsModificationID", "Position"]
        ]
        df_modif_mass = pd.read_sql("SELECT * FROM FoundModifications", connection)[
            ["ModificationID", "DeltaMonoisotopicMass"]
        ]
        df_modif = df_modif.merge(df_modif_mass, left_on="FoundModificationsModificationID", right_on="ModificationID")
        df = df.merge(df_modif, on="TargetPsmsPeptideID")

        logger.info("Finished reading mascot msf file.")

        df.rename(
            columns={
                "SpectrumID": "SCAN NUMBER",
                "ModifiedSequence": "MODIFIED SEQUENCE",
                "Charge": "PRECURSOR CHARGE",
                "XCorr": "SCORE",
                "SpectrumFileName": "RAW FILE",
            },
            inplace=True,
        )

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")
        # TODO reverse
        df["REVERSE"] = df["SEQUENCE"].str.contains("Reverse")
        logger.info("Converting MSFragger  peptide sequence to internal format")
        df["RAW_FILE"] = df["RAW_FILE"].str.replace(".raw", "")
        df["MODIFICATIONS"] = (
            (df["POSITION"].astype(int) - 1).astype(str) + "$" + df["DELTAMONOISOTOPICMASS"].astype(str)
        )
        df = df.groupby("SCAN_NUMBER", as_index=False).apply(lambda x: x.sort_values("POSITION"))
        df = df.groupby(
            ["SCAN_NUMBER", "PRECURSOR_CHARGE", "SCORE", "RAW_FILE", "SEQUENCE", "REVERSE"],
            as_index=False,
        ).agg({"MODIFICATIONS": "|".join})
        mod_masses_reverse = {round(float(v), 3): k for k, v in c.MOD_MASSES.items()}


        def custom_regex_escape(key: str) -> str:
            """
            Subfunction to escape only normal brackets in the modstring.

            :param key: The match to escape
            :return: match with escaped special characters
            """
            for k, v in {"[": r"\[", "]": r"\]", "(": r"\(", ")": r"\)"}.items():
                key = key.replace(k, v)
            return key
        
        def find_replacement(match: re.Match, seq: str) -> str:
            """
        Subfunction to find the corresponding substitution for a match.

        :param match: an re.Match object found by re.sub
        :return: substitution string for the given match
        """
            key = match.string[match.start() : match.end()]
            if custom_var_mods is not None and key in custom_var_mods.keys():
                assert isinstance(custom_mods[key][0], str), f"Provided illegal custom mod format, expected dict-values are (str, float), recieved {(type(custom_mods[key][0]).__name__), (type(custom_mods[key][1]).__name__)}."
                end = match.span()[1]
                if end < len(seq) and (seq[end] == "[" or seq[end]== "("):
                    return key
                if not custom_mods[key][0].startswith(key):
                        return key + custom_mods[key][0]
                return custom_mods[key][0]
            elif custom_stat_mods is not None and key in custom_stat_mods.keys():
                assert isinstance(custom_mods[key][0], str), f"Provided illegal custom mod format, expected dict-values are (str, float), recieved {(type(replacements[key][0]).__name__), (type(replacements[key][1]).__name__)}."
                return custom_mods[key][0]
            return custom_mods[key]
        

        custom_mods = {}

        if custom_var_mods is not None:
            custom_mods.update(custom_var_mods)
        if custom_stat_mods is not None:
            custom_mods.update(custom_stat_mods)

        if custom_mods:
            regex = re.compile("|".join(map(custom_regex_escape, custom_mods.keys())))

        sequences = []
        for _, row in df.iterrows():
            modifications = row["MODIFICATIONS"].split("|")
            sequence = row["SEQUENCE"]
            if custom_mods:
                sequence = regex.sub(lambda match: find_replacement(match, sequence), sequence) 

            if len(modifications) == 0:
                sequences.append(sequence)
            else:                 
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

        return filter_valid_prosit_sequences(df)
