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

# from .masterSpectrum import MasterSpectrum

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


def _chunk_merge(df1: pd.DataFrame, df2: pd.DataFrame, common_column: str, chunk_size: int = 100000) -> pd.DataFrame:
    """
    Merge two DataFrames in chunks based on a common column.

    Splits both input DataFrames into chunks of specified size and merges corresponding chunks
    based on the provided common column. It then filters the merged chunks based on specific conditions
    before concatenating them to obtain the final merged DataFrame.

    :param df1: First pd.DataFrame to be merged.
    :param df2: Second pd.DataFrame to be merged.
    :param common_column: Name of the column on which the DataFrames will be merged.
    :param chunk_size: Size of chunks used for splitting the DataFrames.
    :return: Merged pd.DataFrame containing the results of merging the input DataFrames in chunks.
    """
    merged_chunks = []
    # Split both DataFrames into chunks
    chunks_df1 = [df1[i : i + chunk_size] for i in range(0, len(df1), chunk_size)]
    chunks_df2 = [df2[i : i + chunk_size] for i in range(0, len(df2), chunk_size)]
    # Iterate through the chunks and merge them
    for chunk1 in chunks_df1:
        for chunk2 in chunks_df2:
            merged_chunk = pd.merge(chunk1, chunk2, on=common_column)
            merged_chunk = merged_chunk[
                (merged_chunk["SCANNUMBEGIN"] <= merged_chunk["SCAN"])
                & (merged_chunk["SCAN"] <= merged_chunk["SCANNUMEND"])
            ]
            merged_chunks.append(merged_chunk)
    # Concatenate the merged chunks to get the final result
    merged_df = pd.concat(merged_chunks, ignore_index=True)
    return merged_df


def _load_maxquant_txt(maxquant_out_path: Path, file_name: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load information files from MaxQuant output.

    This function searches for specific text files in the provided MaxQuant output folder,
    filters entries based on the given file name, and returns the content of these files as
    pd.DataFrames.

    :param maxquant_out_path: Path to the MaxQuant output folder.
    :param file_name: Raw file name used to filter entries in the loaded txt files.
    :raises AssertionError: If any of the requested text files does not does not contain any entries
        for the provided file name.
    :return: Tuple containing pd.DataFrames containing the data from the loaded txt files
        (msms.txt, accumulatedMsmsScans.txt, pasefMsmsScans.txt).
    """
    df_msms = pd.read_csv(maxquant_out_path / "msms.txt", sep="\t")
    df_precursors = pd.read_csv(maxquant_out_path / "accumulatedMsmsScans.txt", sep="\t")
    df_pasef = pd.read_csv(maxquant_out_path / "pasefMsmsScans.txt", sep="\t")
    for df in [df_msms, df_precursors, df_pasef]:
        df = _sanitize_columns(df)
    df_msms = df_msms[df_msms["RAW_FILE"] == file_name]
    df_precursors = df_precursors[df_precursors["RAW_FILE"] == file_name]
    df_pasef = df_pasef[df_pasef["RAW_FILE"] == file_name]
    if df_msms.empty:
        raise AssertionError(f"No entries for {file_name} found in msms.txt. " "Check the rawfile name.")
    if df_precursors.empty:
        raise AssertionError(
            f"No entries for {file_name} found in accumulatedMsmsScans.txt. " "Check the rawfile name."
        )
    if df_pasef.empty:
        raise AssertionError(f"No entries for {file_name} found in pasefMsmsScans.txt. " "Check the rawfile name.")
    return df_msms, df_precursors, df_pasef


def _load_metadata(
    out_path: Path, file_name: str, search_engine: str
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load metadata based on the search engine used in MaxQuant output.

    Depending on the specified search engine, this function loads metadata from the MaxQuant
    output or raises an error if the search engine is not supported.

    :param out_path: Path to the search engine output folder.
    :param file_name: Raw file name used to filter entries in the loaded txt files.
    :param search_engine: Search engine used for processing the data (e.g., 'maxquant', 'sage', 'msfragger').
    :raises NotImplementedError: If the search engine is 'sage' or 'msfragger'.
    :raises ValueError: If the search engine is not recognized.
    :return: Tuple containing pd.DataFrames containing metadata based on the search engine used.
    """
    if search_engine == "maxquant":
        return _load_maxquant_txt(out_path, file_name)
    if search_engine == "sage":
        raise NotImplementedError("Oktoberfest does not yet supper results from SAGE.")
    if search_engine == "msfragger":
        raise NotImplementedError("Oktoberfest does not yet supper results from MSFragger.")
    raise ValueError(f"{search_engine} is not recognized.")


def load_timstof(d_path: Path, out_path: Path, search_engine: str = "maxquant") -> pd.DataFrame:
    """
    Load timsTOF data and merge with metadata from search engine output.

    This function loads raw timsTOF data and merges it with metadata (precursor, frame, scan information)
    obtained from the search engine output.

    :param d_path: Path to the raw timsTOF data.
    :param out_path: Path to the search engine output folder.
    :param search_engine: Search engine used for processing the data ('maxquant', 'sage', 'msfragger').
    :return: DataFrame containing combined information from the timsTOF data and search engine metadata.
    """
    # Load the raw bruker data
    logger.info("Loading data...")
    data = alphatims.bruker.TimsTOF(str(d_path))
    # Get the filename to ensure the correct data is mapped
    file_name = d_path.stem
    # Load the PSMs, precursor and frame information
    df_msms, df_precursors, df_pasef = _load_metadata(out_path, file_name, search_engine)
    df_precursors = df_precursors[df_precursors["SCAN_NUMBER"].isin(df_msms.SCAN_NUMBER)]
    df_precursors["PRECURSOR"] = df_precursors["PASEF_PRECURSOR_IDS"].str.split(";")
    df_precursors = df_precursors.explode("PRECURSOR", ignore_index=True).drop_duplicates()
    df_precursors["PRECURSOR"] = df_precursors["PRECURSOR"].astype(int)
    # Get which precursors are in which scan
    scan_precursor_map = df_precursors[["SCAN_NUMBER", "PRECURSOR"]].drop_duplicates()
    df_pasef = df_pasef[df_pasef["PRECURSOR"].isin(df_precursors.PRECURSOR)]
    df_pasef = df_pasef.rename(columns={"COLLISIONENERGY": "COLLISION_ENERGY"})
    # Get where each frame starts and ends
    df_pasef = df_pasef[["PRECURSOR", "FRAME", "SCANNUMBEGIN", "SCANNUMEND", "COLLISION_ENERGY"]].drop_duplicates()
    # Get the frames from the raw bruker data
    df_raw = data[: df_pasef.FRAME]
    df_raw = data[
        {
            "frame_indices": list(df_pasef.FRAME),
        }
    ]
    df_raw = _sanitize_columns(df_raw)
    df_raw = pd.DataFrame(df_raw)
    df_raw = df_raw.rename(columns={"MOBILITY": "INV_ION_MOBILITY", "RT": "RETENTION_TIME"})
    df_raw = df_raw[["FRAME", "SCAN", "TOF", "INTENSITY", "MZ", "INV_ION_MOBILITY", "RETENTION_TIME"]]
    # Map scan information to the raw bruker data
    df_raw_mapped = _chunk_merge(df1=df_raw, df2=df_pasef, common_column="FRAME")
    df_raw_mapped = df_raw_mapped.drop_duplicates()
    # Combine the MZ and INTENSITY information on frame-level
    df_raw_mapped_test = (
        df_raw_mapped.groupby(["PRECURSOR", "FRAME"])
        .agg(
            {
                "INTENSITY": list,
                "MZ": list,
                "RETENTION_TIME": "first",
                "COLLISION_ENERGY": "first",
                "INV_ION_MOBILITY": "first",
            }
        )
        .reset_index()
    )
    # Combine the MZ and INTENSITY information on the summed scan-level
    df_scans = (
        df_raw_mapped_test.merge(scan_precursor_map)
        .groupby("SCAN_NUMBER")
        .agg(
            median_CE=("COLLISION_ENERGY", "median"),
            combined_INTENSITIES=("INTENSITY", lambda x: [item for sublist in x for item in sublist]),
            combined_MZ=("MZ", lambda x: [item for sublist in x for item in sublist]),
            median_RETENTION_TIME=("RETENTION_TIME", "median"),
            median_INV_ION_MOBILITY=("INV_ION_MOBILITY", "median"),
        )
        .reset_index()
    )
    # Get the CHARGE, MASS_ANALYZER, RAW_FILE, and FRAGMENTATION from the msms.txt file
    df_msms_scans = pd.merge(df_scans, df_msms, on="SCAN_NUMBER")
    df_msms_scans = df_msms_scans[
        [
            "RAW_FILE",
            "SCAN_NUMBER",
            "combined_INTENSITIES",
            "combined_MZ",
            "MASS_ANALYZER",
            "FRAGMENTATION",
            "median_RETENTION_TIME",
            "median_INV_ION_MOBILITY",
            "median_CE",
            "CHARGE",
        ]
    ]
    logger.info("Done loading data.")
    return df_msms_scans


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


def combine_spectra(df_msms_scans: pd.DataFrame, temp_path: Path, chunk_size: int = 1000) -> pd.DataFrame:
    """
    Combine spectra from the provided pd.DataFrame and perform binning on chunks.

    This function splits the input pd.DataFrame into chunks, performs binning on each chunk of data,
    merges the binning results with the original data, processes the combined spectra, and returns
    the combined and processed spectra as a pd.DataFrame.

    :param df_msms_scans: pd.DataFrame containing spectra information.
    :param temp_path: Path used for intermediate results during binning.
    :param chunk_size: Size of chunks used for splitting the DataFrame.
    :return: pd.DataFrame containing combined and processed spectra.
    """
    chunk_list = []
    # Split both DataFrames into chunks
    chunks = [df_msms_scans[i : i + chunk_size] for i in range(0, len(df_msms_scans), chunk_size)]
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


def convert_d_pkl(d_path: Path, search_output_path: Path, output_path: Path):
    """
    Converts a .d folder to pkl.

    :param d_path: Path to the .d folder
    :param search_output_path: Path to the output folder from the search engine
    :param output_path: file path of the pkl file
    """
    df_msms_scans = load_timstof(d_path, search_output_path)
    df_combined = combine_spectra(df_msms_scans, output_path)
    # Write to pickle
    file_name = df_combined["RAW_FILE"][0]
    df_combined.to_pickle(output_path / f"{file_name}.pkl")
