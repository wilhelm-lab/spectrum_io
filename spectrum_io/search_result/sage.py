from __future__ import annotations

import logging

import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import internal_without_mods

from .search_results import SearchResults, parse_mods

logger = logging.getLogger(__name__)


class Sage(SearchResults):
    """Handle search results from Sage."""

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {
            "C[+57.0215]": 4,
            "M[+15.9949]": 35,
            "M[+15.994]": 35,
            "R[+0.98402]": 7,
            "Q[+0.98402]": 7,
            "N[+0.98402]": 7,
        }

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Function to read a msms tsv and perform some basic formatting.

        :param tmt_label: optional tmt label as str
        :param custom_mods: optional dictionary mapping Sage-specific mod pattern to UNIMOD IDs.
            If None, static carbamidomethylation of cytein and variable oxidation of methionine
            are mapped automatically. To avoid this, explicitely provide an empty dictionary.
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :return: pd.DataFrame with the formatted data
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))
        if tmt_label:
            unimod_tag = c.TMT_MODS[tmt_label]
            parsed_mods[r"K\[\+\d+\.\d+\]"] = f"K{unimod_tag}"
            parsed_mods[r"^\[\+\d+\.\d+\]"] = f"{unimod_tag}"

        logger.info(f"Reading {self.path}...")
        self.results = pd.read_csv(
            self.path,
            usecols=["filename", "scannr", "peptide", "charge", "hyperscore", "calcmass", "label", "proteins"],
            sep="\t",
        )
        logger.info(f"Finished reading {self.path}.")

        self.convert_to_internal(mods=parsed_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites)
        return self.filter_valid_prosit_sequences()

    def filter_valid_prosit_sequences(self):
        """Filter valid Prosit sequences."""
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results.index)}")
        # retain only peptides that fall within [7, 30] length supported by Prosit
        self.results = self.results[(self.results["PEPTIDE_LENGTH"] <= 30) & (self.results["PEPTIDE_LENGTH"] >= 7)]
        # remove unsupported mods to exclude
        self.results = self.results[~self.results["MODIFIED_SEQUENCE"].str.contains(r"\[\d+\]", regex=True)]
        # remove precursor charges greater than 6
        self.results = self.results[self.results["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(self.results.index)}")

        return self.results

    def convert_to_internal(self, mods: dict[str, str], ptm_unimod_id: int | None, ptm_sites: list[str] | None):
        """
        Convert all columns in the Sage output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping Sage-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        df = self.results

        df.fillna({"proteins": "UNKNOWN"}, inplace=True)
        df.replace({"filename": {r"\.mz[M|m][l|L]": ""}, "peptide": mods}, regex=True, inplace=True)
        if not df["scannr"].dtype == int:
            df["scannr"] = df["scannr"].str.rsplit(pat="=", n=1).str[1].astype("int64")
        df["label"] = df["label"] < 0
        df["SEQUENCE"] = internal_without_mods(df["peptide"])
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].str.len()

        df.rename(
            columns={
                "filename": "RAW_FILE",
                "scannr": "SCAN_NUMBER",
                "peptide": "MODIFIED_SEQUENCE",
                "charge": "PRECURSOR_CHARGE",
                "calcmass": "MASS",
                "hyperscore": "SCORE",
                "label": "REVERSE",
                "proteins": "PROTEINS",
            },
            inplace=True,
        )
