import logging
import sqlite3
from typing import Dict, Optional

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Mascot(SearchResults):
    """Handle search results from Mascot."""

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {}

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: Optional[Dict[str, int]] = None,
        ptm_unimod_id: Optional[int] = 0,
        ptm_sites: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Function to read a mascot msf file and perform some basic formatting.

        :param tmt_label: tmt label as str
        :param custom_mods: dict with custom variable and static identifier and respecitve internal equivalent and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :raises NotImplementedError: always
        :return: pd.DataFrame with the formatted data
        """
        raise NotImplementedError

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
        mod_masses = c.update_mod_masses()
        mod_masses_reverse = {round(float(v), 3): k for k, v in mod_masses.items()}

        parsed_mods = self.standard_mods | (custom_mods or {})

        sequences = []
        for _, row in df.iterrows():
            modifications = row["MODIFICATIONS"].split("|")
            sequence = row["SEQUENCE"]
            if parsed_mods:
                sequence = sequence.replace(parsed_mods, sequence, regex=True)

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

        return self.filter_valid_prosit_sequences()
