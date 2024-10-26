import logging

import os
import re
from pathlib import Path
from typing import Union
from math import ceil
import numpy as np

from typing import Dict, Optional


import pandas as pd
from spectrum_fundamentals.mod_string import xisearch_to_internal
import multiprocess as mp

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Xisearch(SearchResults):
    """Handle search results from xisearch."""

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: Optional[Dict[str, int]] = None,
        ptm_unimod_id: Optional[int] = 0,
        ptm_sites: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Function to read a csv of CSMs and perform some basic formatting.

        :param tmt_label: tmt label as str
        :param custom_mods: dict with custom variable and static identifier and respecitve internal equivalent and mass
        :param ptm_unimod_id: unimod id used for site localization
        :param ptm_sites: possible sites that the ptm can exist on
        :raises NotImplementedError: if a tmt label is provided
        :return: pd.DataFrame with the formatted data
        """
        if tmt_label != "":
            raise NotImplementedError("TMT is not supported for XIsearch")

        logger.info("Reading search results file...")
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

        converters = {"mods_p1": str, "mods_p2": str, "mod_pos_p1": str, "mod_pos_p2": str, "start_pos_p1": str, "start_pos_p2": str}

        df = pd.read_csv(self.path, sep="\t", usecols=columns_to_read, converters=converters)

        logger.info("Finished reading search results file.")
        # Standardize column names
        df = Xisearch.filter_xisearch_result(df)
        df = Xisearch.update_columns_for_prosit(df)
        df = Xisearch.filter_valid_prosit_sequences(df)
        #df = Xisearch.filter_duplicates(df)
        df = Xisearch.fdr_group(df, fdr_group_col=None, decoy_class=None)

        return df

    @staticmethod
    def filter_xisearch_result(df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove unsupported modifications and keep only k-k as linked amino acid .

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """
        df['linear'] = df['linear'].fillna(True)
        df['linear'] = df['linear'].astype(bool)
        df = df[~df["linear"]]
        df = df[df["linked_aa_p1"].notna() & df["linked_aa_p1"].str.contains("K")]
        df = df[df["linked_aa_p2"].notna() & df["linked_aa_p2"].str.contains("K")]
        df = df[~df["mods_p1"].str.contains("dsso-hyd", na=False)]
        df = df[~df["mods_p2"].str.contains("dsso-hyd", na=False)]

        return df

    @staticmethod
    def self_or_between(df):
        df.loc[:,'protein_p1_arr'] = df.loc[:,'protein_p1'].astype(str).str.replace('REV_','').str.split(';').map(np.unique).map(list)
        df.loc[:,'protein_p2_arr'] = df.loc[:,'protein_p2'].astype(str).str.replace('REV_','').str.split(';').map(np.unique).map(list)
        df.loc[:,'proteins_arr'] = (df.loc[:,'protein_p1_arr'] + df.loc[:,'protein_p2_arr'])
        df.loc[:,'proteins_arr_unique'] = df.loc[:,'proteins_arr'].map(np.unique)
        df.loc[:,'is_between'] = ~(df.loc[:,'proteins_arr'].map(len) - df.loc[:,'proteins_arr_unique'].map(len)).astype(bool)
        df.drop(
            [
                'protein_p1_arr',
                'protein_p2_arr',
                'proteins_arr',
                'proteins_arr_unique'
            ],
            axis=1,
            inplace=True
    )
        return df.loc[:,'is_between'].map({True: 'between', False: 'self'})

    @staticmethod
    def self_or_between_mp(df):
        pool_size = min([10,os.cpu_count()])
        slice_size = ceil(len(df)/pool_size)
        cols = ['protein_p1', 'protein_p2']
        df = df[cols]
        df_slices = [
            df.iloc[i*slice_size:(i+1)*slice_size][cols]
            for i in range(pool_size)
        ]
        print('slicing done')
        print(f"Pool size: {pool_size}")
        with mp.Pool(processes=pool_size) as pool:
            map_res = pool.map(Xisearch.self_or_between, df_slices)
        return pd.concat(map_res).copy()


    @staticmethod
    def fdr_group(df, fdr_group_col=None, decoy_class=None):
        original_order = df.index
        if fdr_group_col is None:
            df.loc[:,'fdr_group'] = Xisearch.self_or_between_mp(df)
        return df



    @staticmethod
    def filter_duplicates(df: pd.DataFrame) -> pd.DataFrame:
        """
        keep csm with higher score and remove duplicate (only top ranks) .

        :param df: df to filter
        :return: filtered df as pd.DataFrame
        """
        repetitive_combinations = df[df.duplicated(subset=['scan_number', 'run_name'], keep=False)]
        filtered_df = repetitive_combinations.groupby(['scan_number', 'run_name']).apply(lambda x: x.loc[x['match_score'].idxmax()])
        filtered_df.reset_index(drop=True, inplace=True)
        final_df = pd.concat([df.drop_duplicates(subset=['scan_number', 'run_name'], keep=False), filtered_df])
        final_df.reset_index(drop=True, inplace=True)
        df = final_df
        return df

    @staticmethod
    def update_columns_for_prosit(df: pd.DataFrame) -> pd.DataFrame:
        """
        Update columns of df to work with xl-prosit.

        :param df: df to modify
        :return: modified df as pd.DataFrame
        """
        df['crosslinker_name'] = df['crosslinker_name'].replace(to_replace='*', value='DSSO')
        df["decoy"] = df["decoy_p1"] | df["decoy_p2"]
        df['run_name'] = df['run_name'].str.replace('-', '_')
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

        return df

    @staticmethod
    def filter_valid_prosit_sequences(df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter valid Prosit sequences.

        :param df: df to filter
        :return: df after filtration
        """
        logger.info(f"#sequences before filtering for valid prosit sequences: {len(df.index)}")

        df = df[(df["PEPTIDE_LENGTH_A"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_A"] >= 6]
        df = df[(df["PEPTIDE_LENGTH_B"] <= 30)]
        df = df[df["PEPTIDE_LENGTH_B"] >= 6]
        df = df[(~df["SEQUENCE_A"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        df = df[(~df["SEQUENCE_B"].str.contains(r"B|\*|\.|U|O|X|Z|\(|\)"))]
        df = df[df["PRECURSOR_CHARGE"] <= 6]
        logger.info(f"#sequences after filtering for valid prosit sequences: {len(df.index)}")

        return df
