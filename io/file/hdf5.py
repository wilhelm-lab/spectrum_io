import pandas as pd
import h5py


def read(path, start_index=0, end_index=-1):
    """
    Read hdf5 file and return df with contents.
    With possibility to partial load for memory issues.
    :param start_index:
    :param end_index:
    :return:
    :param path:
    :return:
    """
    df = pd.DataFrame()
    return df


def write(df, path):
    """
    Write dataframe to hdf5 file
    :param df:
    :param path:
    """
