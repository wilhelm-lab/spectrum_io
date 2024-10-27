from __future__ import annotations

import glob
import logging
import os
import re
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
import spectrum_fundamentals.constants as c
from spectrum_fundamentals.mod_string import xisearch_to_internal

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Scout(SearchResults):
    """Handle search results from xisearch."""

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ) -> pd.DataFrame:
        """
        Function to read a csv of CSMs and perform some basic formatting.

        :param tmt_label: tmt label as str
        :param custom_mods: dict with custom variable and static identifier and respecitve internal equivalent and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :raises NotImplementedError: if TMT label is provided
        :return: pd.DataFrame with the formatted data
        """
        if tmt_label != "":
            raise NotImplementedError("TMT is not supported for Scout")

        logger.info("Reading search results file...")
        columns_to_read = [
            "ScanNumber",
            "Charge",
            "ExperimentalMZ",
            "AlphaPeptide",
            "BetaPeptide",
            "AlphaPos",
            "BetaPos",
            "AlphaMappings",
            "BetaMappings",
            "ClassificationScore",
            "Peptide Position 1",
            "Peptide Position 2",
            "Protein 1",
            "Protein 2",
            "FileName",
        ]

        self.results = pd.read_csv(self.path, usecols=columns_to_read)
        logger.info("Finished reading search results file.")
        # Standardize column names
        self.convert_to_internal(mods={})
        return self.filter_valid_prosit_sequences()
        # df = Scout._filter_duplicates(df)

    @staticmethod
    def _filter_duplicates(df: pd.DataFrame) -> pd.DataFrame:
        """
        Keep csm with higher score and remove duplicate (only top ranks).

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """
        repetitive_combinations = df[df.duplicated(subset=["ScanNumber", "RAW_FILE"], keep=False)]
        filtered_df = repetitive_combinations.groupby(["ScanNumber", "RAW_FILE"]).apply(
            lambda x: x.loc[x["ClassificationScore"].idxmax()]
        )
        filtered_df.reset_index(drop=True, inplace=True)
        final_df = pd.concat([df.drop_duplicates(subset=["ScanNumber", "RAW_FILE"], keep=False), filtered_df])
        final_df.reset_index(drop=True, inplace=True)
        df = final_df
        return df

    @staticmethod
    def _extract_modifications(peptide_seq: str):
        modifications = []
        # Find all matches of modifications
        matches = re.findall(r"([CM])\(\+([\d.]+)\)", peptide_seq)
        for match in matches:
            mod, _ = match
            # Add modification to the list
            if mod == "C":
                modifications.append("cm")
            elif mod == "M":
                modifications.append("ox")
        return ";".join(modifications)

    @staticmethod
    def _extract_modification_positions(peptide_seq: str):
        pattern = r"([A-Z])(\(\+\d+\.\d+\))?"
        matches = re.findall(pattern, peptide_seq)
        split_peptide = []
        for match in matches:
            amino_acid = match[0]
            modification = match[1] if match[1] else ""
            split_peptide.append(amino_acid + modification)
        positions = [str(i + 1) for i, component in enumerate(split_peptide) if "+" in component]
        return ";".join(positions)

    @staticmethod
    def _self_or_between_mp(df: pd.DataFrame) -> pd.DataFrame:
        df["tmp_id"] = df.index
        df_expl = df.copy()
        df_expl.loc[:, "AlphaMappings"] = df_expl["AlphaMappings"].str.split(";")
        df_expl.loc[:, "BetaMappings"] = df_expl["BetaMappings"].str.split(";")
        df_expl = df_expl.explode("AlphaMappings")
        df_expl = df_expl.explode("BetaMappings")
        df_expl.loc[:, "self"] = False
        df_expl.loc[df_expl["AlphaMappings"] == df_expl["BetaMappings"], "self"] = True
        id_to_self = df_expl.groupby("tmp_id", dropna=False).agg({"self": "max"}).reset_index()
        df = df.drop(["self"], axis=1, errors="ignore").merge(id_to_self, on=["tmp_id"], validate="1:1")
        df.loc[:, "fdr_group"] = df["self"].apply(lambda x: "self" if x else "between")
        return df

    def convert_to_internal(
        self, mods: dict[str, str], ptm_unimod_id: int | None = None, ptm_sites: list[str] | None = None
    ):
        """
        Convert all columns in the search engine-specific output to the internal format used by Oktoberfest.

        :param mods: dictionary mapping search engine-specific mod patterns (keys) to ProForma standard (values)
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        """
        # Filter csms that does not contain any "k"
        df = self.results
        df = df[(df["AlphaPeptide"].str.contains("K")) & (df["BetaPeptide"].str.contains("K"))]
        df["decoy_p1"] = df["AlphaMappings"].str.contains("Reverse").astype(bool)
        df["decoy_p2"] = df["BetaMappings"].str.contains("Reverse").astype(bool)
        df["protein_p1"] = df["AlphaMappings"]
        df["protein_p2"] = df["BetaMappings"]
        df["decoy"] = df["decoy_p1"] | df["decoy_p2"]
        df["REVERSE"] = df["decoy"]
        df["RAW_FILE"] = df["FileName"].apply(lambda x: x.split("\\")[-1])
        df["MASS"] = df["ExperimentalMZ"]
        df["PRECURSOR_CHARGE"] = df["Charge"]
        df["CROSSLINKER_TYPE"] = "DSSO"
        df["crosslinker_name"] = "DSSO"
        df["linked_aa_p1"] = "K"
        df["linked_aa_p2"] = "K"
        df["linear"] = "False"
        df["match_score"] = "ClassificationScore"
        df["SCORE"] = df["ClassificationScore"]
        df["SCAN_NUMBER"] = df["ScanNumber"]
        df["SEQUENCE_A"] = df["AlphaPeptide"].apply(lambda x: re.sub(r"\([^)]*\)", "", x))
        df["SEQUENCE_B"] = df["BetaPeptide"].apply(lambda x: re.sub(r"\([^)]*\)", "", x))
        df["base_sequence_p1"] = df["SEQUENCE_A"]
        df["base_sequence_p2"] = df["SEQUENCE_B"]
        df = df[df.apply(lambda row: row["SEQUENCE_A"][row["AlphaPos"]] == "K", axis=1)]
        df = df[df.apply(lambda row: row["SEQUENCE_B"][row["BetaPos"]] == "K", axis=1)]
        df["Modifications_A"] = df["AlphaPeptide"].apply(Scout._extract_modifications)
        df["Modifications_B"] = df["BetaPeptide"].apply(Scout._extract_modifications)
        df["mods_p1"] = df["AlphaPeptide"].apply(Scout._extract_modifications)
        df["mods_p2"] = df["BetaPeptide"].apply(Scout._extract_modifications)
        df["ModificationPositions1"] = df["AlphaPeptide"].apply(Scout._extract_modification_positions)
        df["ModificationPositions2"] = df["BetaPeptide"].apply(Scout._extract_modification_positions)
        df["CROSSLINKER_POSITION_A"] = df["AlphaPos"] + 1
        df["CROSSLINKER_POSITION_B"] = df["BetaPos"] + 1
        df["mod_pos_p1"] = df["AlphaPos"] + 1
        df["mod_pos_p2"] = df["BetaPos"] + 1
        df["link_pos_p1"] = df["AlphaPos"] + 1
        df["link_pos_p2"] = df["BetaPos"] + 1
        df["PEPTIDE_LENGTH_A"] = df["SEQUENCE_A"].apply(len)
        df["PEPTIDE_LENGTH_B"] = df["SEQUENCE_B"].apply(len)
        df["aa_len_p1"] = df["SEQUENCE_A"].apply(len)
        df["aa_len_p2"] = df["SEQUENCE_B"].apply(len)
        df = Scout._self_or_between_mp(df)
        df["fdr_group"] = np.where(
            df["AlphaMappings"].str.replace("Reverse_", "") == df["BetaMappings"].str.replace("Reverse_", ""),
            "self",
            "between",
        )
        df.drop(columns=["self"], inplace=True)
        df.drop(columns=["tmp_id"], inplace=True)
        logger.info("Converting Scout peptide sequence to internal format...")
        df["RAW_FILE"] = df["RAW_FILE"].str.replace(".raw", "")
        df["MODIFIED_SEQUENCE_A"] = df.apply(
            lambda row: xisearch_to_internal(
                xl=row["CROSSLINKER_TYPE"],
                seq=row["SEQUENCE_A"],
                mod=row["Modifications_A"],
                crosslinker_position=row["CROSSLINKER_POSITION_A"],
                mod_positions=row["ModificationPositions1"],
            ),
            axis=1,
            result_type="expand",
        )
        df["MODIFIED_SEQUENCE_B"] = df.apply(
            lambda row: xisearch_to_internal(
                xl=row["CROSSLINKER_TYPE"],
                seq=row["SEQUENCE_B"],
                mod=row["Modifications_B"],
                crosslinker_position=row["CROSSLINKER_POSITION_B"],
                mod_positions=row["ModificationPositions2"],
            ),
            axis=1,
            result_type="expand",
        )
        new_column_names = {
            "FileName": "run_name",
            "ScanNumber": "scan_number",
            "ExperimentalMZ": "precursor_mass",
            "Charge": "precursor_charge",
            "scan_number": "ScanNumber",
        }
        self.results = df.rename(columns=new_column_names)

    def filter_valid_prosit_sequences(self) -> pd.DataFrame:
        """
        Filter valid Prosit sequences.

        :return: df after filtration
        """
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results)}")
        self.results = self.results[(self.results["PEPTIDE_LENGTH_A"] <= 30)]
        self.results = self.results[self.results["PEPTIDE_LENGTH_A"] >= 6]
        self.results = self.results[(self.results["PEPTIDE_LENGTH_B"] <= 30)]
        self.results = self.results[self.results["PEPTIDE_LENGTH_B"] >= 6]
        self.results = self.results[(~self.results["SEQUENCE_A"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        self.results = self.results[(~self.results["SEQUENCE_B"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        self.results = self.results[self.results["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(self.results)}")

        return self.results
