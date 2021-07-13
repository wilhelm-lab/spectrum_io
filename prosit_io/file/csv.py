import pandas as pd


def read_file(path: str):
    """
    Read csv file and return df with contents
    :param path:
    :return
    """
    df = pd.read_csv(path, sep=",")
    return df


def write_file(df, path):
    """
    Write dataframe to csv file
    :param df:
    :param path:
    """
