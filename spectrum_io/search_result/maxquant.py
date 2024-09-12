import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods

from .search_results import SearchResults, filter_valid_prosit_sequences, parse_mods

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

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {
            "C": 4,
            "M(ox)": 35,
            "M(Oxidation (M))": 35,
        }

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

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: Optional[Dict[str, int]] = None,
    ) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_label: optional tmt label as str
        :param custom_mods: optional dictionary mapping MaxQuant-specific mod pattern to UNIMOD IDs.
            If None, static carbamidomethylation of cytein and variable oxidation of methionine
            are mapped automatically. To avoid this, explicitely provide an empty dictionary.
        :return: pd.DataFrame with the formatted data
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))
        if tmt_label:
            unimod_tag = c.TMT_MODS[tmt_label]
            parsed_mods["K"] = f"K{unimod_tag}"
            parsed_mods["^_"] = f"_{unimod_tag}-"

        logger.info("Reading msms.txt file")
        self.results = pd.read_csv(
            self.path / "msms.txt",
            usecols=[
                "Raw file",
                "Scan number",
                "Modified sequence",
                "Charge",
                "Scan event number",
                "Mass",  # = Calculated Precursor mass; TODO get column with experimental Precursor mass instead
                "Score",
                "Reverse",
                "Proteins",
            ],
            sep="\t",
        )

        logger.info("Finished reading msms.txt file")

        self.convert_to_internal(mods=parsed_mods)
        return filter_valid_prosit_sequences(self.results)

    def convert_to_internal(self, mods: Dict[str, str]):
        """
        Convert all columns in the MaxQuant output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping MaxQuant-specific mod patterns (keys) to ProForma standard (values)
        """
        df = self.results
        # Standardize column names
        # df.columns = df.columns.str.upper()
        # df.columns = df.columns.str.replace(" ", "_")
        # df.rename(columns={"CHARGE": "PRECURSOR_CHARGE"}, inplace=True)

        mods["_"] = ""

        df.fillna({"Reverse": "", "Proteins": "UNKNOWN"}, inplace=True)
        df["Reverse"] = df["Reverse"].astype(bool)
        df.replace({"Modified sequence": mods}, regex=True, inplace=True)

        df["Sequence"] = internal_without_mods(df["Modified sequence"])
        df["PEPTIDE_LENGTH"] = df["Sequence"].str.len()

        df.rename(
            columns={
                "Reverse": "REVERSE",
                "Sequence": "SEQUENCE",
                "Modified sequence": "MODIFIED_SEQUENCE",
                "Proteins": "PROTEINS",
                "Charge": "PRECURSOR_CHARGE",
                "Raw file": "RAW_FILE",
                "Scan number": "SCAN_NUMBER",
                "Scan event number": "SCAN_EVENT_NUMBER",
                "Mass": "MASS",
                "Score": "SCORE",
            },
            inplace=True,
        )

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
