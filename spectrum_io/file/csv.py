from pathlib import Path
from typing import Union

import pandas as pd


def read_file(path: Union[str, Path]) -> pd.DataFrame:
    """
    Read csv file and return df with contents.

    :param path: path to file to read
    :return: df with contents as pd.DataFrame
    """
    df = pd.read_csv(path, sep=",")
    return df


def write_file(df: pd.DataFrame, path: Union[str, Path]):
    """
    Write dataframe to csv file.

    :param df: df with contents as pd.DataFrame
    :param path: path to file to write
    """
    df.to_csv(path, index=False)
