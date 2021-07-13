import pandas as pd
import h5py
from scipy.sparse import coo_matrix
import pandas as pd

def read_file(path, start_index=0, end_index=-1):
    """
    Read hdf5 file and return df with contents.
    With possibility to partial load for memory issues.
    :param path:
    :param start_index:
    :param end_index:
    :return:
    """
    with h5py.File(path, 'r') as f:
        df = None
        for index in range(start_index, end_index):
            ## TODO indexing
            ## TODO type checks and numpy convert
            ## categorical data ?
            if f[index].startswith("sparse"):
                values = f[index]["values"]
                i = f[index]["i"]
                j = f[index]["j"]
                shape = f[index]["shape"]
                sparse_data = coo_matrix((values, (i, j)), shape=shape)
                current_df = pd.DataFrame.sparse.from_spmatrix(sparse_data)
            else:
                current_df = pd.DataFrame(f[index])
            df = current_df if df is None else pd.concat(df, current_df, axis=1)

    return df


def write_file(
    data:Union[pd.DataFrame, scipy.sparse.coo_matrix],
    path: str,
    dataset_name: str,
    mode: str = 'w',
):
    ## TODO categorical data?
    with h5py.File(path, mode) as f:
        if isinstance(data, pd.DataFrame):
            f.create_dataset(dataset_name, data=np.asarray(data, compression="gzip", chunks=True, maxshape=(None,))
        elif isinstance(data, scipy.sparse.coo_matrix):
            group_name = f"sparse_{dataset_name}"
            f.create_group(group_name)
            i,j, val = scipy.sparse.find(data)
            shape = data.shape
            f.create_dataset(f"{group_name}/i", i)
            f.create_dataset(f"{group_name}/j", j)
            f.create_dataset(f"{group_name}/values", values)
            f.create_dataset(f"{group_name}/shape", shape)
        else:
            assert False, "Only pd.DataFrame and scipy.sparse.coo_matrix are supported."
