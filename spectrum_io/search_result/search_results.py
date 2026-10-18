from __future__ import annotations

import logging
import re
from abc import abstractmethod
from pathlib import Path

import pandas as pd
from spectrum_fundamentals.mod_string import add_permutations

from spectrum_io.file import csv

logger = logging.getLogger(__name__)


def parse_mods(mods: dict[str, int]) -> dict[str, str]:
    """
    Parse provided mapping of custom modification pattern to ProForma standard.

    This function takes a dictionary mapping custom modification pattern for specific aminoacids (keys) to a
    UNIMOD ID (values). The pattern is translated to ProForma standard and a new dictionary mapping the custom
    modification patterns to the ProForma standard is returned.
    The pattern for the custom modifications must start with the one-letter code for an aminoacid or '^' / '$',
    to describe n- / c-terminal modifications, respectively, followed by an optional pattern (which can be
    empty).
    This means that 'X'  or 'X(custom_pattern)', is both mapped to 'X[UNIMOD:#]'.
    For the n-terminus, an additional dash will be added automatically, which maps 'X(custom_pattern)' to
    'X[UNIMOD:#]-'. If the sequence to apply the transformation on already contains the dash, it needs to be part
    of the custom_pattern (i.e. 'X(custom_pattern)-'), to avoid adding an additional dash.

    :param mods: Dictionary mapping custom modification patterns (keys) to UNIMOD IDs (values)
    :raises TypeError: if keys are not strings or values are not integers
    :raises ValueError: if the keys do not start with [A-Z,a-z,^,$]
    :return: A dictionary mapping custom modification patterns (keys) to the ProForma standard (values)
    """
    key_pattern = (
        "'X' or 'X<mod_pattern>' where X is either the one-letter code of an aminoacid or '^' / '$' defining the"
        " n- or c-terminus, respectively, followed by an optional pattern identifying a specific modification."
    )
    unimod_regex_map = {}
    for k, v in mods.items():
        if not isinstance(v, int):
            raise TypeError(f"UNIMOD id {v} for replacement {k} not understood. UNIMOD IDs must be integers.")
        if not isinstance(k, str):
            raise TypeError(
                f"Replacement {k} not understood. Replacements must be strings and follow the pattern {key_pattern}"
            )
        if k[0].isalpha():
            unimod_regex_map[re.escape(k)] = f"{k[0].upper()}[UNIMOD:{v}]"
            continue
        if k[0] == "^":
            unimod_regex_map[f"^{re.escape(k[1:])}"] = f"[UNIMOD:{v}]-"
            continue
        raise ValueError(
            f"Replacement {k} not understood. {k[0]} is not a valid aminoacid. "
            f"Replacements most follow the pattern {key_pattern}"
        )
    # sort descending by length to avoid replacement of single-aa keys within larger keys
    return dict(sorted(unimod_regex_map.items(), key=lambda x: -len(x[0])))


class SearchResults:
    """Handle search results from different software."""

    orig_res: pd.DataFrame
    fake_msms: pd.DataFrame

    def __init__(self, path: str | Path):
        """
        Init Searchresults object.

        :param path: path to file
        """
        if isinstance(path, str):
            path = Path(path)
        self.path = path

    @abstractmethod
    def filter_valid_prosit_sequences(self):
        """Filter valid Prosit sequences."""
        raise NotImplementedError

    @staticmethod
    def get_ptm_sites(
        ptm_unimod_id: int | dict[int, list[str]] | None = 0, ptm_sites: list[str] | None = None
    ) -> dict[int, list[str]]:
        """
        Normalize the PTM localization settings to a mapping of unimod id to possible residues.

        :param ptm_unimod_id: unimod id of a single PTM, or a mapping of unimod id to possible residues
        :param ptm_sites: possible residues if a single unimod id is given
        :return: mapping of unimod id to possible residues, empty if no localization is requested
        """
        if isinstance(ptm_unimod_id, dict):
            return {int(unimod_id): list(sites) for unimod_id, sites in ptm_unimod_id.items() if sites}
        if ptm_unimod_id and ptm_sites:
            return {int(ptm_unimod_id): list(ptm_sites)}
        return {}

    @staticmethod
    def add_ptm_permutations(
        df: pd.DataFrame,
        ptm_unimod_id: int | dict[int, list[str]] | None = 0,
        ptm_sites: list[str] | None = None,
        sequence_column: str = "MODIFIED_SEQUENCE",
    ) -> pd.DataFrame:
        """
        Generate all possible localizations of the given PTM(s) for each PSM.

        Each PSM carrying a PTM is duplicated once per possible combination of modified sites, so that all
        localization candidates can be scored. PSMs without the PTM(s) are left unchanged. If several PTMs are
        given, the permutations of all PTMs are combined.

        :param df: search results containing the modified sequences in internal format
        :param ptm_unimod_id: unimod id of a single PTM, or a mapping of unimod id to possible residues,
            0 or None to disable
        :param ptm_sites: amino acids the PTM can be located on, if a single unimod id is given
        :param sequence_column: name of the column containing the modified sequences
        :return: search results with one row per possible PTM localization, and the localization of the search engine
            in the column SEARCH_MODIFIED_SEQUENCE
        """
        ptm_site_map = SearchResults.get_ptm_sites(ptm_unimod_id, ptm_sites)
        if not ptm_site_map:
            return df
        df = df.copy()
        df["SEARCH_MODIFIED_SEQUENCE"] = df[sequence_column]
        for unimod_id, residues in ptm_site_map.items():
            # for citrullination / deamidation also consider that the mod mass was mistakenly picked as monoisotopic peak
            df[sequence_column] = df[sequence_column].apply(
                add_permutations,
                unimod_id=unimod_id,
                residues=residues,
                allow_one_less_modification=unimod_id == 7,
            )
            df = df.explode(sequence_column, ignore_index=True)
        return df

    @abstractmethod
    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ):
        """Read result.

        :param tmt_label: tmt label as str
        :param custom_mods: dict with custom variable and static identifier and respecitve internal equivalent and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        raise NotImplementedError

    def generate_internal(
        self,
        tmt_label: str = "",
        out_path: str | Path | None = None,
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | dict[int, list[str]] | None = 0,
        ptm_sites: list[str] | None = None,
        ptm_model: bool = False,
        xl: bool = False,
    ) -> pd.DataFrame:
        """
        Generate df and save to out_path if provided.

        :param out_path: path to output
        :param tmt_label: tmt label as str
        :param custom_mods: dict with static and variable custom modifications, their internal identifier and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :param xl: set to True for crosslinking data
        :return: path to output file
        """
        if out_path is None:
            # convert and return
            filtered_df = self.read_result(
                tmt_label, custom_mods=custom_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites
            )
            return filtered_df

        if isinstance(out_path, str):
            out_path = Path(out_path)

        if out_path.is_file():
            # only read converted and return
            logger.info(f"Found search results in internal format at {out_path}, skipping conversion")
            # TODO: internal_to_unimod
            return csv.read_file(out_path)

        # convert, save and return
        df = self.read_result(
            tmt_label, custom_mods=custom_mods, ptm_unimod_id=ptm_unimod_id, ptm_sites=ptm_sites, ptm_model=ptm_model
        )
        csv.write_file(df, out_path)
        return df

    def read_internal(self) -> pd.DataFrame:
        """
        Read file from path.

        :return: dataframe after reading the file
        """
        return csv.read_file(self.path)

    @abstractmethod
    def convert_to_internal(self, mods: dict[str, str], ptm_unimod_id: int | None, ptm_sites: list[str] | None):
        """
        Convert all columns in the search engine-specific output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping search engine-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        raise NotImplementedError
