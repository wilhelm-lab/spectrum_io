from __future__ import annotations

import logging
import os
import re
from math import ceil

# import multiprocess as mp
from multiprocessing import pool
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from spectrum_fundamentals.mod_string import xisearch_to_internal

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Xisearch(SearchResults):
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
            raise NotImplementedError("TMT is not supported for XIsearch")

        logger.info("Reading search results file...")
        # assume xiSEARCH2 style
        try:
            columns_to_read = [
                "run_name",
                "scan_number",
                "precursor_mass",
                "precursor_charge",
                "crosslinker_name",
                "decoy_p1",
                "base_sequence_p1",
                "sequence_p1",
                "sequence_p2",
                "start_pos_p1",
                "start_pos_p2",
                "aa_len_p1",
                "link_pos_p1",
                "linked_aa_p1",
                "mods_p1",
                "mod_pos_p1",
                "protein_p1",
                "decoy_p2",
                "base_sequence_p2",
                "aa_len_p2",
                "link_pos_p2",
                "linked_aa_p2",
                "mods_p2",
                "mod_pos_p2",
                "protein_p2",
                "linear",
                "match_score",
            ]

            converters = {
                "mods_p1": str,
                "mods_p2": str,
                "mod_pos_p1": str,
                "mod_pos_p2": str,
                "start_pos_p1": str,
                "start_pos_p2": str,
            }

            self.results = pd.read_csv(self.path, sep="\t", usecols=columns_to_read, converters=converters)
        except ValueError:
            # assume xiSEARCH1 style
            column_mapping = {
                "Run": "run_name",
                "Scan": "scan_number",
                "PrecursorMass": "precursor_mass",
                "PrecoursorCharge": "precursor_charge",
                "Crosslinker": "crosslinker_name",
                "BasePeptide1": "base_sequence_p1",
                "Peptide1": "sequence_p1",
                "Peptide2": "sequence_p2",
                "Start1": "start_pos_p1",
                "Start2": "start_pos_p2",
                "LengthPeptide1": "aa_len_p1",
                "Link1": "link_pos_p1",
                "Linked AminoAcid 1": "linked_aa_p1",
                "Modifications1": "mods_p1",
                "ModificationPositions1": "mod_pos_p1",
                "Protein1": "protein_p1",
                "Protein2decoy": "decoy_p2",
                "BasePeptide2": "base_sequence_p2",
                "LengthPeptide2": "aa_len_p2",
                "Link2": "link_pos_p2",
                "Linked AminoAcid 2": "linked_aa_p2",
                "Modifications2": "mods_p2",
                "ModificationPositions2": "mod_pos_p2",
                "Protein2": "protein_p2",
                "match score": "match_score",
            }

            converters = {
                "Modifications1": str,
                "Modifications2": str,
                "ModificationPositions1": str,
                "ModificationPositions2": str,
                "Start1": str,
                "Start2": str,
            }
            # read in the xi1 columns
            # could be csv
            try:
                self.results = pd.read_csv(
                    self.path, sep=",", usecols=[*column_mapping.keys(), "decoy"], converters=converters
                )
            except ValueError:
                # or tsv
                self.results = pd.read_csv(
                    self.path, sep="\t", usecols=[*column_mapping.keys(), "decoy"], converters=converters
                )

            # convert ot xi2 column names
            self.results.rename(columns=column_mapping, inplace=True)
            # add the decoy columns
            self.results["decoy_p1"] = self.results["decoy"] & (
                self.results["protein_p1"].str.contains("REV_") | self.results["protein_p1"].str.contains("RAN_")
            )
            self.results["decoy_p2"] = self.results["decoy"] & (
                self.results["protein_p2"].str.contains("REV_") | self.results["protein_p2"].str.contains("RAN_")
            )
            # flag linears
            self.results["Linear"] = self.results["protein_p2"].isna()

        logger.info("Finished reading search results file.")
        # Standardize column names
        self.filter_xisearch_result()
        self.convert_to_internal(mods={})
        self.filter_valid_prosit_sequences()
        # df = Xisearch._filter_duplicates(df)
        self.results = Xisearch._fdr_group(self.results, fdr_group_col=None)

        return self.results

    def filter_xisearch_result(self):
        """Remove unsupported modifications and keep only k-k as linked amino acid."""
        df = self.results
        df["linear"] = df["linear"].fillna(True)
        df["linear"] = df["linear"].astype(bool)
        df = df[~df["linear"]]
        df = df[df["linked_aa_p1"].notna() & df["linked_aa_p1"].str.contains("K")]
        df = df[df["linked_aa_p2"].notna() & df["linked_aa_p2"].str.contains("K")]
        df = df[~df["mods_p1"].str.contains("dsso-hyd", na=False)]
        df = df[~df["mods_p2"].str.contains("dsso-hyd", na=False)]

        self.results = df

    @staticmethod
    def _self_or_between(df):
        df.loc[:, "protein_p1_arr"] = (
            df.loc[:, "protein_p1"].astype(str).str.replace("REV_", "").str.split(";").map(np.unique).map(list)
        )
        df.loc[:, "protein_p2_arr"] = (
            df.loc[:, "protein_p2"].astype(str).str.replace("REV_", "").str.split(";").map(np.unique).map(list)
        )
        df.loc[:, "proteins_arr"] = df.loc[:, "protein_p1_arr"] + df.loc[:, "protein_p2_arr"]
        df.loc[:, "proteins_arr_unique"] = df.loc[:, "proteins_arr"].map(np.unique)
        df.loc[:, "is_between"] = ~(
            df.loc[:, "proteins_arr"].map(len) - df.loc[:, "proteins_arr_unique"].map(len)
        ).astype(bool)
        df.drop(["protein_p1_arr", "protein_p2_arr", "proteins_arr", "proteins_arr_unique"], axis=1, inplace=True)
        return df.loc[:, "is_between"].map({True: "between", False: "self"})

    @staticmethod
    def _self_or_between_mp(df):
        pool_size = min([10, os.cpu_count()])
        slice_size = ceil(len(df) / pool_size)
        cols = ["protein_p1", "protein_p2"]
        df = df[cols]
        df_slices = [df.iloc[i * slice_size : (i + 1) * slice_size][cols] for i in range(pool_size)]
        print("slicing done")
        print(f"Pool size: {pool_size}")
        with pool.Pool(processes=pool_size) as self_or_between_pool:
            map_res = self_or_between_pool.map(Xisearch._self_or_between, df_slices)
        return pd.concat(map_res).copy()

    @staticmethod
    def _fdr_group(df, fdr_group_col=None):
        if fdr_group_col is None:
            df.loc[:, "fdr_group"] = Xisearch._self_or_between_mp(df)
        return df

    @staticmethod
    def _filter_duplicates(df: pd.DataFrame) -> pd.DataFrame:
        """
        Keep csm with higher score and remove duplicate (only top ranks).

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """
        repetitive_combinations = df[df.duplicated(subset=["scan_number", "run_name"], keep=False)]
        filtered_df = repetitive_combinations.groupby(["scan_number", "run_name"]).apply(
            lambda x: x.loc[x["match_score"].idxmax()]
        )
        filtered_df.reset_index(drop=True, inplace=True)
        final_df = pd.concat([df.drop_duplicates(subset=["scan_number", "run_name"], keep=False), filtered_df])
        final_df.reset_index(drop=True, inplace=True)
        df = final_df
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
        df = self.results
        df["crosslinker_name"] = df["crosslinker_name"].replace(to_replace="*", value="DSSO")
        df["decoy"] = df["decoy_p1"] | df["decoy_p2"]
        df["run_name"] = df["run_name"].str.replace("-", "_")
        df["RAW_FILE"] = df["run_name"]
        df["MASS"] = df["precursor_mass"]
        df["PRECURSOR_CHARGE"] = df["precursor_charge"]
        df["CROSSLINKER_TYPE"] = df["crosslinker_name"]
        df["SCORE"] = df["match_score"]
        df["REVERSE"] = df["decoy"]
        df["SCAN_NUMBER"] = df["scan_number"]
        df["SEQUENCE_A"] = df["base_sequence_p1"]
        df["SEQUENCE_B"] = df["base_sequence_p2"]
        df["Modifications_A"] = df["mods_p1"]
        df["Modifications_B"] = df["mods_p2"]
        df["CROSSLINKER_POSITION_A"] = df["link_pos_p1"]
        df["CROSSLINKER_POSITION_B"] = df["link_pos_p2"]
        df["ModificationPositions1"] = df["mod_pos_p1"]
        df["ModificationPositions2"] = df["mod_pos_p2"]
        df["PEPTIDE_LENGTH_A"] = df["aa_len_p1"]
        df["PEPTIDE_LENGTH_B"] = df["aa_len_p2"]
        logger.info("Converting Xisearch peptide sequence to internal format...")

        df["RAW_FILE"] = df["RAW_FILE"].str.replace(".raw", "", regex=False)

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

        self.results = df

    def filter_valid_prosit_sequences(self) -> pd.DataFrame:
        """
        Filter valid Prosit sequences.

        :return: df after filtration
        """
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(self.results)}")
        df = self.results
        df = df[(df["PEPTIDE_LENGTH_A"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_A"] >= 6]
        df = df[(df["PEPTIDE_LENGTH_B"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_B"] >= 6]
        df = df[(~df["SEQUENCE_A"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        df = df[(~df["SEQUENCE_B"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        df = df[df["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(df)}")
        self.results = df
        return self.results
