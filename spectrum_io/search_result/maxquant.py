import logging
from pathlib import Path
from typing import Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal

from .search_results import SearchResults, filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


class MaxQuant(SearchResults):
    """Handle search results from MaxQuant."""

    def __init__(self, path: Union[str, Path]):
        """
        Init Searchresults object.

        :param path: path to file
        """
        if isinstance(path, str):
            path = Path(path)

        if path.is_file() and path.name == "msms.txt":
            path = path.parent
        self.path = path

    @staticmethod
    def add_tmt_mod(mass: float, seq: str, unimod_tag: str) -> float:
        """
        Add tmt modification.

        :param mass: mass without tmt modification
        :param seq: sequence of the peptide
        :param unimod_tag: UNIMOD tag for the modification
        :return: mass as float
        """
        num_of_tmt = seq.count(unimod_tag)
        mass += num_of_tmt * c.MOD_MASSES[f"{unimod_tag}"]
        return mass

    def read_result(self, tmt_labeled: str) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_labeled: tmt label as str
        :return: pd.DataFrame with the formatted data
        """
        logger.info("Reading msms.txt file")
        df = pd.read_csv(
            self.path / "msms.txt",
            usecols=lambda x: x.upper()
            in [
                "RAW FILE",
                "SCAN NUMBER",
                "MODIFIED SEQUENCE",
                "CHARGE",
                "SCAN EVENT NUMBER",
                "LABELING STATE",
                "MASS",  # = Calculated Precursor mass; TODO get column with experimental Precursor mass instead
                "SCORE",
                "REVERSE",
                "PROTEINS",
            ],
            sep="\t",
        )
        logger.info("Finished reading msms.txt file")

        # Standardize column names
        df.columns = df.columns.str.upper()
        df.columns = df.columns.str.replace(" ", "_")

        df = MaxQuant.update_columns_for_prosit(df, tmt_labeled)
        return filter_valid_prosit_sequences(df)

    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame, tmt_labeled: str) -> pd.DataFrame:
        """
        Update columns of df to work with Prosit.

        :param df: df to modify
        :param tmt_labeled: True if tmt labeled
        :return: modified df as pd.DataFrame
        """
        df.rename(columns={"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)

        df["REVERSE"].fillna(False, inplace=True)
        df["REVERSE"].replace("+", True, inplace=True)
        logger.info("Converting MaxQuant peptide sequence to internal format")
        if tmt_labeled != "":
            unimod_tag = c.TMT_MODS[tmt_labeled]
            logger.info("Adding TMT fixed modifications")
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(
                df["MODIFIED_SEQUENCE"].to_numpy(),
                fixed_mods={"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}-", "K": f"K{unimod_tag}"},
            )
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, unimod_tag), axis=1)
            if "msa" in tmt_labeled:
                logger.info("Replacing phospho by dehydration for Phospho-MSA")
                df["MODIFIED_SEQUENCE_MSA"] = df["MODIFIED_SEQUENCE"].str.replace(
                    "[UNIMOD:21]", "[UNIMOD:23]", regex=False
                )
        elif "LABELING_STATE" in df.columns:
            logger.info("Adding SILAC fixed modifications")
            df.loc[df["LABELING_STATE"] == 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(
                df[df["LABELING_STATE"] == 1]["MODIFIED_SEQUENCE"].to_numpy(),
                fixed_mods={"C": "C[UNIMOD:4]", "K": "K[UNIMOD:259]", "R": "R[UNIMOD:267]"},
            )
            df.loc[df["LABELING_STATE"] != 1, "MODIFIED_SEQUENCE"] = maxquant_to_internal(
                df[df["LABELING_STATE"] != 1]["MODIFIED_SEQUENCE"].to_numpy()
            )
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, "[UNIMOD:259]"), axis=1)
            df["MASS"] = df.apply(lambda x: MaxQuant.add_tmt_mod(x.MASS, x.MODIFIED_SEQUENCE, "[UNIMOD:267]"), axis=1)
            df.drop(columns=["LABELING_STATE"], inplace=True)
        else:
            df["MODIFIED_SEQUENCE"] = maxquant_to_internal(df["MODIFIED_SEQUENCE"].to_numpy())
        df["SEQUENCE"] = internal_without_mods(df["MODIFIED_SEQUENCE"])
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].apply(lambda x: len(x))
        df["PROTEINS"].fillna("UNKNOWN", inplace=True)

        return df

    def generate_internal_timstof_metadata(self):
        """
        Load information files required for correct aggregation of spectra in timsTOF experiments.

        :return: dataframe containing the columns RAW_FILE, SCANNUMBER, PRECURSOR, FRAME, SCANNUMBEGIN, SCANNUMEND, CollisionEnergy
        """
        df_msms = pd.read_csv(self.path / "msms.txt", sep="\t", usecols=["Raw file", "Scan number"])
        df_msms.columns = ["RAW_FILE", "SCAN_NUMBER"]

        df_precursors = pd.read_csv(
            self.path / "accumulatedMsmsScans.txt", sep="\t", usecols=["Raw file", "Scan number", "PASEF precursor IDs"]
        )
        df_precursors.columns = ["RAW_FILE", "SCAN_NUMBER", "PRECURSOR"]
        df_precursors.query("SCAN_NUMBER in @df_msms.SCAN_NUMBER", inplace=True)
        df_precursors["PRECURSOR"] = df_precursors["PRECURSOR"].str.split(";")
        df_precursors = df_precursors.explode("PRECURSOR")
        df_precursors["PRECURSOR"] = df_precursors["PRECURSOR"].astype("int")

        df_pasef = pd.read_csv(
            self.path / "pasefMsmsScans.txt",
            sep="\t",
            usecols=["Raw file", "Frame", "Precursor", "ScanNumBegin", "ScanNumEnd", "CollisionEnergy"],
        )
        df_pasef.columns = ["RAW_FILE", "FRAME", "PRECURSOR", "SCAN_NUM_BEGIN", "SCAN_NUM_END", "COLLISION_ENERGY"]

        return df_pasef.merge(df_precursors).sort_values(["FRAME", "PRECURSOR"])
