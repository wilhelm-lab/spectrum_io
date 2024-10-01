from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import add_permutations, internal_without_mods

from .search_results import SearchResults, parse_mods

logger = logging.getLogger(__name__)


class MaxQuant(SearchResults):
    """Handle search results from MaxQuant."""

    def __init__(self, path: str | Path):
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
            "R(Citrullination)": 7,
            "Q(Deamidation (NQ))": 7,
            "N(Deamidation (NQ))": 7,
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

    def filter_valid_prosit_sequences(self):
        """Filter valid Prosit sequences."""
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results.index)}")
        # retain only peptides that fall within [7, 30] length supported by Prosit
        self.results = self.results[(self.results["PEPTIDE_LENGTH"] <= 30) & (self.results["PEPTIDE_LENGTH"] >= 7)]
        # remove unsupported mods to exclude
        self.results = self.results[~self.results["MODIFIED_SEQUENCE"].str.contains(r"\(", regex=True)]
        # remove precursor charges greater than 6
        self.results = self.results[self.results["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(self.results.index)}")

        return self.results

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Function to read a msms txt and perform some basic formatting.

        :param tmt_label: optional tmt label as str
        :param custom_mods: optional dictionary mapping MaxQuant-specific mod pattern to UNIMOD IDs.
            If None, static carbamidomethylation of cytein and variable oxidation of methionine
            are mapped automatically. To avoid this, explicitely provide an empty dictionary.
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
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

        self.convert_to_internal(mods=parsed_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites)
        return self.filter_valid_prosit_sequences()

    def convert_to_internal(self, mods: dict[str, str], ptm_unimod_id: int | None, ptm_sites: list[str] | None):
        """
        Convert all columns in the MaxQuant output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping MaxQuant-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
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
        if ptm_unimod_id != 0:

            # PTM permutation generation
            if ptm_unimod_id == 7:
                allow_one_less_modification = True
            else:
                allow_one_less_modification = False

            df["Modified sequence"] = df["Modified sequence"].apply(
                add_permutations,
                unimod_id=ptm_unimod_id,
                residues=ptm_sites,
                allow_one_less_modification=allow_one_less_modification,
            )
            df = df.explode("Modified sequence", ignore_index=True)

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
        self.results = df

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
