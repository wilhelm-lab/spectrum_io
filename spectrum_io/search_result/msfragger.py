from __future__ import annotations

import logging

import pandas as pd
import spectrum_fundamentals.constants as c
from pyteomics import pepxml
from spectrum_fundamentals.mod_string import internal_without_mods
from tqdm import tqdm

from .search_results import SearchResults, parse_mods

logger = logging.getLogger(__name__)


class MSFragger(SearchResults):
    """Handle search results from MSFragger."""

    @property
    def standard_mods(self):
        """Standard modifications that are always applied if not otherwise specified."""
        return {
            "C[160]": 4,
            "M[147]": 35,
            "R[157]": 7,
            "Q[129]": 7,
            "N[115]": 7,
        }

    @staticmethod
    def apply_modifications(peptide: str, modifications: list[str]) -> str:
        """
        Applies modifications to a peptide sequence with specific formatting rules:
        - Masses with whole number 170 always show 2 decimal places
        - Other masses show decimals only if non-zero

        Args:
            peptide: Unmodified peptide string (e.g., "KAAGASKR")
            modifications: List of modification strings in "mass@position" format

        Returns:
            Modified peptide string with formatted modifications
        """

        def format_mass(mass_str: str) -> str:
            """Formats mass strings according to specified rules"""
            if "." in mass_str:
                integer_part, decimal_part = mass_str.split(".", 1)
            else:
                integer_part = mass_str

            # if integer_part == "170":
            #    decimal = (decimal_part.ljust(2, "0"))[:2]  # Pad/truncate to 2 decimals
            #    return f"170.{decimal}"
            # else:
            return integer_part

        # Create dictionary of {position: formatted_mass} (0-based)
        mod_dict = {}
        for mod in modifications:
            mass, pos = mod.split("@")
            mod_dict[int(pos) - 1] = format_mass(mass)  # Convert to 0-based indexing

        # Build modified sequence
        modified_sequence = []
        if -1 in mod_dict:
            modified_sequence.append(f"n[{mod_dict[-1]}]")
        for idx, aa in enumerate(peptide):
            if idx in mod_dict:
                modified_sequence.append(f"{aa}[{mod_dict[idx]}]")
            else:
                modified_sequence.append(aa)

        return "".join(modified_sequence)

    def filter_valid_prosit_sequences(self):
        """Filter valid Prosit sequences."""
        # logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results.index)}")
        # retain only peptides that fall within [7, 30] length supported by Prosit
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results.index)}")

        self.results = self.results[(self.results["PEPTIDE_LENGTH"] <= 30) & (self.results["PEPTIDE_LENGTH"] >= 6)]
        # remove unsupported mods to exclude
        # self.results = self.results[~self.results["MODIFIED_SEQUENCE"].isna()]

        self.results = self.results[~self.results["MODIFIED_SEQUENCE"].str.contains(r"\[\d+\]", regex=True)]

        # self.results = self.results[~self.results['MODIFIED_SEQUENCE'].str.contains(r'C(?!\[)', regex=True)]
        # remove precursor charges greater than 6
        self.results = self.results[self.results["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(self.results.index)}")

        return self.results

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | dict[int, list[str]] | None = 0,
        ptm_sites: list[str] | None = None,
        ptm_model: bool = False,
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

        if ptm_model:
            parsed_mods = dict(c.MSFRAGGER_VAR_MODS_IO)
        if tmt_label:
            unimod_tag = c.TMT_MODS[tmt_label]
            # apply_modifications writes the labelled lysine as K[<mass>]; also map a bare K
            k_tmt_mass = int(c.AA_MASSES["K"] + c.MOD_MASSES[unimod_tag])
            parsed_mods[rf"K(\[{k_tmt_mass}\])?"] = f"K{unimod_tag}"
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
        self.results["modified_peptide"] = self.results[["peptide", "modifications"]].apply(
            lambda x: MSFragger.apply_modifications(x["peptide"], x["modifications"]), axis=1
        )

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

    def convert_to_internal(
        self, mods: dict[str, str], ptm_unimod_id: int | dict[int, list[str]] | None, ptm_sites: list[str] | None
    ):
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
        df = df[~df["modified_peptide"].str.contains(r"U", regex=True)]
        df["PEPTIDE_LENGTH"] = df["peptide"].str.len()

        df.replace({"modified_peptide": mods}, regex=True, inplace=True)
        df["peptide"] = internal_without_mods(df["modified_peptide"])
        # PTM permutation generation (only if localization was requested)
        df = self.add_ptm_permutations(df, ptm_unimod_id, ptm_sites, sequence_column="modified_peptide")

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
                "num_matched_ions": "Annotated_Ions_MSF",
                "tot_num_ions": "Total_Ions_MSF",
                "massdiff": "MZ_diff_MSF",
                "nextscore": "NEXT_SCORE",
                "expect": "EXPECT",
            },
            inplace=True,
        )
        self.results = df[
            [
                "RAW_FILE",
                "SCAN_NUMBER",
                "MODIFIED_SEQUENCE",
                "PRECURSOR_CHARGE",
                "SCAN_EVENT_NUMBER",
                "MASS",
                "SCORE",
                "Annotated_Ions_MSF",
                "Total_Ions_MSF",
                "MZ_diff_MSF",
                "EXPECT",
                "NEXT_SCORE",
                "REVERSE",
                "SEQUENCE",
                "PEPTIDE_LENGTH",
                "PROTEINS",
            ]
            # localization of the search engine, if PTM site permutations were generated
            + [column for column in ["SEARCH_MODIFIED_SEQUENCE"] if column in df.columns]
        ]
