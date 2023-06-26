import logging
from pathlib import Path
from typing import Union

import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES

from .search_results import filter_valid_prosit_sequences

logger = logging.getLogger(__name__)


def _update_columns_for_prosit(df: pd.DataFrame):
    """
    Update columns of df to work with Prosit.

    :param df: df to modify
    :return: modified df as pd.DataFrame
    """

    def _check_unsupported_mods(modstring: str):
        if modstring:
            for mod in modstring.split(";"):
                if not mod.split("(", 1)[1].startswith(("Oxidation", "Carbamidomethyl")):
                    raise AssertionError(f"Unsupported modification: {mod}")

    def _make_mod_seq(seq: str):
        seq = seq.replace("m", "M[UNIMOD:35]")
        seq = seq.replace("c", "C[UNIMOD:4]")
        return seq

    df["REVERSE"] = df["Protein Accessions"].str.startswith("REV_")
    df["MASS"] = df["m/z"] * df["PRECURSOR_CHARGE"] - (df["PRECURSOR_CHARGE"] * PARTICLE_MASSES["PROTON"])
    df["Modifications"].fillna("").apply(_check_unsupported_mods)
    df["MODIFIED_SEQUENCE"] = df["SEQUENCE"].apply(_make_mod_seq)
    df["SEQUENCE"] = df["SEQUENCE"].str.upper()
    df["PEPTIDE_LENGTH"] = df["SEQUENCE"].str.len()
    df["RAW_FILE"] = df["RAW_FILE"].str.replace(r"\.\w+$", "", regex=True)

    return df[
        [
            "SCORE",
            "MASS",
            "PRECURSOR_CHARGE",
            "RAW_FILE",
            "SCAN_NUMBER",
            "REVERSE",
            "MODIFIED_SEQUENCE",
            "SEQUENCE",
            "PEPTIDE_LENGTH",
        ]
    ]


def _remove_decoys_in_targets(full_df):
    duplicated_psms = full_df[["SCAN_NUMBER", "RAW_FILE", "MODIFIED_SEQUENCE"]].duplicated(keep=False)
    logger.info(f"Removing {sum(duplicated_psms)} duplicated PSMs...")
    full_df = full_df[~(full_df["REVERSE"] & duplicated_psms)]
    if any(full_df[["SCAN_NUMBER", "RAW_FILE", "MODIFIED_SEQUENCE"]].duplicated(keep=False)):
        raise AssertionError("There are duplicate target PSMs. This is a bug of msamanda!")
    logger.info(f"{len(full_df)} remaining PSMs")

    return full_df


def read_result(path: Union[str, Path], suffix: str = "output.csv") -> pd.DataFrame:
    """
    Function to read a msms txt and perform some basic formatting.

    :param path: path to msms.txt to read
    :param suffix: Optional suffix to determine which fileresult files should be taken from the supplied path
    :raises FileNotFoundError: If the supplied path is not found
    :raises AssertionError: If the supplied path does not contain any files matching the provided suffix.
    :return: pd.DataFrame with the formatted data
    """
    if isinstance(path, str):
        path = Path(path)
    if path.is_file():
        pathlist = [path]
    elif path.is_dir():
        pathlist = list(path.glob(f"*{suffix}"))
        if not pathlist:
            raise AssertionError(f"The directory does not contain any files that match the pattern *{suffix}")
    else:
        raise FileNotFoundError(f"{path} does not exist.")

    df_list = []
    for output_file in pathlist:
        logger.info(f"Reading {output_file}...")
        df = pd.read_csv(
            output_file,
            sep="\t",
            skiprows=1,
            usecols=[
                "Scan Number",
                "Sequence",
                "Modifications",
                "Filename",
                "Amanda Score",
                "m/z",
                "Charge",
                "Protein Accessions",
            ],
        )
        df.columns = [
            "SCAN_NUMBER",
            "SEQUENCE",
            "Modifications",
            "Protein Accessions",
            "SCORE",
            "m/z",
            "PRECURSOR_CHARGE",
            "RAW_FILE",
        ]
        logger.info(f"Finished reading {output_file}")

        df = _update_columns_for_prosit(df)
        df = filter_valid_prosit_sequences(df)
        df_list.append(df)

    if len(df_list) == 1:
        full_df = df_list[0]
    full_df = pd.concat(df_list, axis=0)
    full_df = _remove_decoys_in_targets(full_df)
    return full_df
