import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def filter_valid_prosit_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter valid Prosit sequences.

    :param df: df to filter
    :return: df after filtering out unsupported peptides
    """
    logger.info(f"#sequences before filtering for valid prosit sequences: {len(df.index)}")
    # retain only peptides that fall within [7, 30] length supported by Prosit
    df = df[(df["PEPTIDE_LENGTH"] <= 30) & (df["PEPTIDE_LENGTH"] >= 7)]
    # remove unsupported mods to exclude
    unsupported_mods = [r"Acetyl \(Protein N\-term\)", "ac", r"\[[0-9]+\]"]
    exclude_mods_pattern = re.compile("|".join(unsupported_mods))
    df = df[~df["MODIFIED_SEQUENCE"].str.contains(exclude_mods_pattern)]
    # remove non-canonical aas
    df = df[(~df["SEQUENCE"].str.contains("U|O"))]
    # remove precursor charges greater than 6
    df = df[df["PRECURSOR_CHARGE"] <= 6]
    logger.info(f"#sequences after filtering for valid prosit sequences: {len(df.index)}")

    return df
