from __future__ import annotations

import logging
from typing import Dict, Optional

import pandas as pd
import spectrum_fundamentals.constants as c
from pyteomics import pepxml
from spectrum_fundamentals.constants import MSFRAGGER_VAR_MODS
from spectrum_fundamentals.mod_string import add_permutations, internal_without_mods
from tqdm import tqdm

from .search_results import SearchResults, parse_mods

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {"C[160]": 4, "M[147]": 35, "R[157]": 7, "Q[129]": 7, "N[115]": 7}

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
        :param custom_mods: optional dictionary mapping MSFragger-specific mod pattern to UNIMOD IDs.
            If None, static carbamidomethylation of cytein and variable oxidation of methionine
            are mapped automatically. To avoid this, explicitely provide an empty dictionary.
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :raises FileNotFoundError: in case the given path is neither a file, nor a directory.
        :return: pd.DataFrame with the formatted data
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))
        if tmt_label:
            unimod_tag = c.TMT_MODS[tmt_label]
            parsed_mods["K"] = f"K{unimod_tag}"
            parsed_mods[r"^n\[\d+\]"] = f"{unimod_tag}-"
        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = list(self.path.rglob("*.pepXML"))
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")

        ms_frag_results = []
        for pep_xml_file in tqdm(file_list):
            ms_frag_results.append(pepxml.DataFrame(str(pep_xml_file)))

        self.results = pd.concat(ms_frag_results)

        self.convert_to_internal(mods=parsed_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites)
        return self.filter_valid_prosit_sequences()

    @staticmethod
    def check_decoys(protein_names: str):
        """
        Check if all protein names in a given string correspond to decoy proteins.

        :param protein_names: A string containing one or more protein names separated by semicolons (';').
                            Each protein name is checked for the presence of the substring 'rev'.
        :return: `True` if all proteins are decoy proteins (i.e., if all protein names contain 'rev'),
                otherwise `False`.
        """
        all_proteins = protein_names.split(";")
        reverse = True
        for protein in all_proteins:
            if "rev" not in protein:
                reverse = False
                break
        return reverse

    def convert_to_internal(self, mods: dict[str, str], ptm_unimod_id: int | None, ptm_sites: list[str] | None):
        """
        Convert all columns in the MSFragger output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping MSFragger-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        df = self.results
        df["protein"] = df["protein"].fillna("UNKNOWN").apply(lambda x: ";".join(x))

        df["REVERSE"] = df["protein"].apply(lambda x: MSFragger.check_decoys(x))

        df["spectrum"] = df["spectrum"].str.split(pat=".", n=1).str[0]
        df["PEPTIDE_LENGTH"] = df["peptide"].str.len()
        df.replace({"modified_peptide": mods}, regex=True, inplace=True)
        df["peptide"] = internal_without_mods(df["modified_peptide"])
        if ptm_unimod_id != 0:

            # PTM permutation generation
            if ptm_unimod_id == 7:
                allow_one_less_modification = True
            else:
                allow_one_less_modification = False

            df["modified_peptide"] = df["modified_peptide"].apply(
                add_permutations,
                unimod_id=ptm_unimod_id,
                residues=ptm_sites,
                allow_one_less_modification=allow_one_less_modification,
            )
            df = df.explode("modified_peptide", ignore_index=True)

        df.rename(
            columns={
                "assumed_charge": "PRECURSOR_CHARGE",
                "index": "SCAN_EVENT_NUMBER",
                "peptide": "SEQUENCE",
                "start_scan": "SCAN_NUMBER",
                "hyperscore": "SCORE",
                "modified_peptide": "MODIFIED_SEQUENCE",
                "protein": "PROTEINS",
                "precursor_neutral_mass": "MASS",
                "spectrum": "RAW_FILE",
            },
            inplace=True,
        )
        self.results = df
        """
        return df[
            [
                "RAW_FILE",
                "SCAN_NUMBER",
                "MODIFIED_SEQUENCE",
                "PRECURSOR_CHARGE",
                "SCAN_EVENT_NUMBER",
                "MASS",
                "SCORE",
                "REVERSE",
                "SEQUENCE",
                "PEPTIDE_LENGTH",
                "PROTEINS",
            ]
        ]
        """
