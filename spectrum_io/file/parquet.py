import logging
from pathlib import Path
from typing import List, Union

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy

# TODO add sparse matrix / anndata support
# TODO add speed benchmarks
# TODO add support for HuggingFace datasets API

Pathlike = Union[Path, str]
Dataset = Union[pd.DataFrame, scipy.sparse.spmatrix]

logger = logging.getLogger(__name__)


def read_file(path: Pathlike) -> pd.DataFrame:
    """
    Read a Parquet file and return a Pandas DataFrame with its contents.

    :param path: Path to the Parquet file to read
    :return: a Pandas DataFrame with the contents of the file
    """
    return pd.read_parquet(path)


def read_partition(path: Pathlike, dataset_name: str) -> pd.DataFrame:
    """
    Read a single table from a partitioned dataset.

    :param path: Root path of the partitioned dataset
    :param dataset_name: Name of the dataset to extract

    :return: a Pandas DataFrame of the specified table from the partitioned dataset
    """
    try:
        dataset = pq.ParquetDataset(path, filters=[("dataset", "=", dataset_name)])
        df = dataset.read().to_pandas().drop("dataset", axis=1)
    except pa.lib.ArrowNotImplementedError as e1:
        if dataset_name.isdigit():
            try:
                logger.warning("Failed to read from partition using string of integer as key, trying with int...")
                dataset = pq.ParquetDataset(path, filters=[("dataset", "=", int(dataset_name))])
                df = dataset.read().to_pandas().drop("dataset", axis=1)
            except pa.lib.ArrowNotImplementedError as e2:
                logger.exception(e2)
        else:
            logger.exception(e1)

    return df


def write_file(data: Dataset, path: Pathlike) -> None:
    """Writes a single DataFrame or matrix to a Parquet file.

    :param data: Data to store
    :param path: Path to write the Parquet file to

    :raises NotImplementedError: if anything else but a Pandas DataFrame is used as the dataset
    """
    if isinstance(data, pd.DataFrame):
        data.to_parquet(path)
    else:
        raise NotImplementedError


def write_partition(datasets: List[Dataset], path: Pathlike, dataset_names: List[str]) -> None:
    """
    Write several datasets to a Parquet dataset as a directory containing subdirectories partinioned by dataset name.

    :param datasets: Datasets to write
    :param path: Root path to write the partitioned dataset to
    :param dataset_names: Names to assign to the datasets for retrieval. Careful: If all of these are strings of ints,
        Parquet will convert them to raw integers!

    :raises NotImplementedError: if anything else but a Pandas DataFrame is used as the dataset
    """
    if all(isinstance(x, pd.DataFrame) for x in datasets):
        df = pd.concat([dataset.assign(dataset=name) for dataset, name in zip(datasets, dataset_names)])
        table = pa.Table.from_pandas(df)
    else:
        raise NotImplementedError

    if isinstance(path, str):
        path = Path(path)
    path.mkdir(exist_ok=True)

    pq.write_to_dataset(
        table,
        root_path=path,
        partition_cols=["dataset"],
        existing_data_behavior="delete_matching",
    )
