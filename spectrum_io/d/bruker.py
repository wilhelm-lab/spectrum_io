import logging
import os
import pickle
from pathlib import Path
from typing import Optional, Tuple, Union

import alphatims
import alphatims.bruker
import alphatims.utils
import numpy as np
import pandas as pd

from .masterSpectrum import MasterSpectrum

logger = logging.getLogger(__name__)


def _sanitize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sanitize DataFrame column names.

    This function replaces spaces with underscores, converts column names to
    uppercase, and removes specific suffixes ('_INDICES' and '_VALUES') from
    column names.

    :param df: pd.DataFrame whose column names need to the sanitized.
    :return: pd.DataFrame with sanitized column names.
    """
    df.columns = [c.replace(" ", "_") for c in df.columns]
    df.columns = [c.upper() for c in df.columns]
    df.columns = [c.replace("_INDICES", "") for c in df.columns]
    df.columns = [c.replace("_VALUES", "") for c in df.columns]
    return df


def binning(inp: pd.DataFrame, ignore_charges: bool, rescoring_path: Path) -> pd.DataFrame:
    """
    Perform binning on the input MasterSpectrum.

    This function loads a MasterSpectrum from the provided input, performs binning, and exports the results
    to a temporary text file in the given rescoring path. It then reads the temporary file as a DataFrame,
    modifies the DataFrame by adding SCAN_NUMBER and dropping specified columns before returning it.

    :param inp: Input data used to perform binning.
    :param ignore_charges: indicating whether charges should be ignored during binning.
    :param rescoring_path: Path where the temporary file will be exported.
    :return: Modified DataFrame after binning with added SCAN_NUMBER and specified columns dropped.
    """
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignore_charges)
    ms.export_to_csv(rescoring_path / "temp.txt")
    comb_ms = pd.read_csv(rescoring_path / "temp.txt")
    scan = inp["SCAN_NUMBER"]
    comb_ms["SCAN_NUMBER"] = scan
    comb_ms = comb_ms.drop(
        columns=[
            "counts",
            "left border",
            "right border",
            "start_mz",
            "ms1_charge",
            "rel_intensity_ratio",
            "counts_ratio",
        ]
    )
    return comb_ms


def aggregate_timstof(raw_spectra: pd.DataFrame, temp_path: Path, chunk_size: int = 1000) -> pd.DataFrame:
    """
    Combine spectra from the provided pd.DataFrame and perform binning on chunks.

    This function splits the input pd.DataFrame into chunks, performs binning on each chunk of data,
    merges the binning results with the original data, processes the combined spectra, and returns
    the combined and processed spectra as a pd.DataFrame.

    :param raw_spectra: pd.DataFrame containing spectra information.
    :param temp_path: Path used for intermediate results during binning.
    :param chunk_size: Size of chunks used for splitting the DataFrame.
    :return: pd.DataFrame containing combined and processed spectra.
    """
    chunk_list = []
    # Split both DataFrames into chunks
    chunks = [raw_spectra[i : i + chunk_size] for i in range(0, len(raw_spectra), chunk_size)]
    for chunk in chunks:
        bin_result_list = []
        for _, line in chunk.iterrows():
            bin_result = binning(line, True, temp_path)
            bin_result_list.append(bin_result)
        bin_result_df = pd.concat(bin_result_list)
        bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
        scans_combined = pd.merge(chunk, bin_result_df_collapsed, on="SCAN_NUMBER")
        scans_comb = scans_combined.drop(columns=["combined_INTENSITIES", "combined_MZ"]).rename(
            columns={"intensity": "INTENSITIES", "mz": "MZ", "median_CE": "COLLISION_ENERGY"}
        )
        # Convert lists into arrays
        scans_comb["INTENSITIES"] = scans_comb["INTENSITIES"].apply(lambda x: np.array(x))
        scans_comb["MZ"] = scans_comb["MZ"].apply(lambda x: np.array(x))
        # Sort the MZ (and linked INTENSITIES)
        for i in range(len(scans_comb)):
            zipped_list = zip(scans_comb.iloc[i]["MZ"], scans_comb.iloc[i]["INTENSITIES"])
            sorted_pair = sorted(zipped_list)
            tuples = zip(*sorted_pair)
            list_1, list_2 = (np.array(tuple) for tuple in tuples)
            scans_comb.at[i, "MZ"] = list_1
            scans_comb.at[i, "INTENSITIES"] = list_2
        scans_comb["MZ_RANGE"] = "0"
        chunk_list.append(scans_comb)
    chunk_comb = pd.concat(chunk_list, ignore_index=True)
    return chunk_comb


def read_timstof(d_path, scan_to_precursor_map):
    # preparation of filter

    df_frame_group = (
        scan_to_precursor_map[["FRAME", "PRECURSOR"]]
        .drop_duplicates()
        .groupby("FRAME", as_index=False)
        .agg(
            {
                "PRECURSOR": tuple,
            }
        )
        .groupby("PRECURSOR", as_index=False)
        .agg({"FRAME": tuple})
    )

    # load filtered stuff
    data = alphatims.bruker.TimsTOF(str(d_path), slice_as_dataframe=False)

    raw_idx = []
    for frames, precursors in zip(df_frame_group["FRAME"], df_frame_group["PRECURSOR"]):
        raw_idx.extend(data[frames, :, precursors])

    df = data.as_dataframe(
        raw_idx,
        raw_indices=False,
        frame_indices=True,
        scan_indices=True,
        quad_indices=False,
        tof_indices=False,
        precursor_indices=True,
        rt_values=True,
        rt_values_min=False,
        mobility_values=True,
        quad_mz_values=False,
        push_indices=False,
        mz_values=True,
        intensity_values=True,
        corrected_intensity_values=False,
        raw_indices_sorted=False,
    )

    df.columns = ["FRAME", "SCAN", "PRECURSOR", "RETENTION_TIME", "INV_ION_MOBILITY", "MZ", "INTENSITY"]

    # aggregation
    df_combined_grouped = (
        df.merge(
            scan_to_precursor_map[
                ["SCAN_NUM_BEGIN", "SCAN_NUM_END", "PRECURSOR", "FRAME", "COLLISION_ENERGY"]
            ].drop_duplicates()
        )
        .query("SCAN_NUM_BEGIN <= SCAN <= SCAN_NUM_END")  # can probably be skipped
        .groupby(["PRECURSOR", "FRAME"], as_index=False)  # aggregate fragments per precursor in FRAME
        .agg(
            {
                "INTENSITY": list,
                "MZ": list,
                "RETENTION_TIME": "first",
                "COLLISION_ENERGY": "first",
                "INV_ION_MOBILITY": "first",
            }
        )
        .merge(scan_to_precursor_map.reset_index())
        .groupby("SCAN_NUMBER", as_index=False)  # aggregate PRECURSORS for same SCAN_NUMBER
        .agg(
            median_CE=("COLLISION_ENERGY", "median"),
            combined_INTENSITIES=("INTENSITY", lambda x: [item for sublist in x for item in sublist]),
            combined_MZ=("MZ", lambda x: [item for sublist in x for item in sublist]),
            median_RETENTION_TIME=("RETENTION_TIME", "median"),
            median_INV_ION_MOBILITY=("INV_ION_MOBILITY", "median"),
        )
    )

    return df_combined_grouped


def convert_d_hdf(
    input_path: Union[Path, str],
    output_path: Optional[Union[Path, str]] = None,
):
    data = alphatims.bruker.TimsTOF(str(input_path))
    data.save_as_hdf(directory=str(output_path.parent), filename=str(output_path.name))


def read_and_aggregate_timstof(source: Path, scan_to_precursor_map: Path):
    """
    Read raw spectra from timstof hdf spectra file and aggregate to MS2 spectra.

    :param source: Path to the hdf file
    :param scan_to_precursor_map: Dataframe mapping scan numbers to precursors
    """
    raw_spectra = read_timstof(source, scan_to_precursor_map)
    df_combined = aggregate_timstof(raw_spectra, temp_path=Path("/tmp"))
    return df_combined
