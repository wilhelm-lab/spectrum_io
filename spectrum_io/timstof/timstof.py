import numpy as np
import pandas as pd
import os
from typing import Optional, Union
from pathlib import Path
from mgf_filter.masterSpectrum import MasterSpectrum
import pickle
import alphatims
import alphatims.utils
import alphatims.bruker

logger = logging.getLogger(__name__)

def rename_columns(df):
    df.columns = [c.replace(' ', '_') for c in df.columns]
    df.columns = [c.upper() for c in df.columns]
    df.columns = [c.replace('_INDICES', '') for c in df.columns]
    df.columns = [c.replace('_VALUES', '') for c in df.columns]
    return df

def chunk_merge(df1, df2, common_column, chunk_size = 100000):
    merged_chunks = []
    # Split both DataFrames into chunks
    chunks_df1 = [df1[i:i + chunk_size] for i in range(0, len(df1), chunk_size)]
    chunks_df2 = [df2[i:i + chunk_size] for i in range(0, len(df2), chunk_size)]
    # Iterate through the chunks and merge them
    for chunk1 in chunks_df1:
        for chunk2 in chunks_df2:
            merged_chunk = pd.merge(chunk1, chunk2, on=common_column)
            merged_chunk = merged_chunk[(merged_chunk['SCANNUMBEGIN'] <= merged_chunk['SCAN']) & (merged_chunk['SCAN'] <= merged_chunk['SCANNUMEND'])]
            merged_chunks.append(merged_chunk)
    # Concatenate the merged chunks to get the final result
    merged_df = pd.concat(merged_chunks, ignore_index=True)
    return merged_df

def load_maxquant_txt(txt_path, txt_file, raw_file_name):
    df = pd.read_csv(os.path.join(txt_path, txt_file), sep="\t")
    df = rename_columns(df)
    df = df[df["RAW_FILE"] == file_name]
    return df

def load_timstof(d_path, txt_path):
    # Load the raw bruker data
    data = alphatims.bruker.TimsTOF(d_path)
    # Get the filename to ensure the correct data is mapped
    raw_file_name = os.path.splitext(os.path.basename(d_path))[0]
    # Load the PSMs
    df_msms = load_maxquant_txt(txt_path, "msms.txt", raw_file_name)
    # Load the precursor information
    df_precursors = load_maxquant_txt(txt_path, "accumulatedMsmsScans.txt", raw_file_name)
    df_precursors = df_precursors[df_precursors["SCAN_NUMBER"].isin(df_msms.SCAN_NUMBER)]
    df_precursors["PRECURSOR"] = df_precursors["PASEF_PRECURSOR_IDS"].str.split(";")
    df_precursors = df_precursors.explode('PRECURSOR', ignore_index=True).drop_duplicates()
    df_precursors["PRECURSOR"] = df_precursors["PRECURSOR"].astype(int)
    # Get which precursors are in which scan
    scan_precursor_map = df_precursors[["SCAN_NUMBER", "PRECURSOR"]].drop_duplicates()
    # Load the frame information
    df_pasef = load_maxquant_txt(txt_path, "pasefMsmsScans.txt", raw_file_name)
    df_pasef = df_pasef[df_pasef["PRECURSOR"].isin(df_precursors.PRECURSOR)]
    df_pasef = df_pasef.rename(columns={"COLLISIONENERGY":"COLLISION_ENERGY"})
    # Get where each frame starts and ends
    df_pasef = df_pasef[["PRECURSOR", "FRAME", "SCANNUMBEGIN", "SCANNUMEND", "COLLISION_ENERGY"]].drop_duplicates()
    # Get the frames from the raw bruker data
    df_raw = data[:df_pasef.FRAME]
    df_raw = data[{"frame_indices": list(df_pasef.FRAME),}]
    df_raw = rename_columns(df_raw)
    df_raw = pd.DataFrame(df_raw)
    df_raw = df_raw.rename(columns={"MOBILITY":"INV_ION_MOBILITY", "RT":"RETENTION_TIME"})
    df_raw = df_raw[['FRAME', 'SCAN', 'TOF', 'INTENSITY', 'MZ', 'INV_ION_MOBILITY','RETENTION_TIME']]
    # Map scan information to the raw bruker data
    df_raw_mapped = chunk_merge(df1=df_raw, df2=df_pasef, common_column="FRAME")
    df_raw_mapped = df_raw_mapped.drop_duplicates()
    # Combine the MZ and INTENSITY information on frame-level
    df_raw_mapped_test = df_raw_mapped.groupby(['PRECURSOR', 'FRAME']).agg({
            'INTENSITY': list,
            'MZ': list,
            'RETENTION_TIME': 'first',
            'COLLISION_ENERGY': 'first',
            'INV_ION_MOBILITY': 'first'
        }).reset_index()
    # Combine the MZ and INTENSITY information on the summed scan-level
    df_scans = df_raw_mapped_test.merge(scan_precursor_map).groupby('SCAN_NUMBER').agg(
            median_CE=('COLLISION_ENERGY', 'median'),
            combined_INTENSITIES=('INTENSITY', lambda x: [item for sublist in x for item in sublist]),
            combined_MZ=('MZ', lambda x: [item for sublist in x for item in sublist]),
            median_RETENTION_TIME=('RETENTION_TIME', 'median'),
            median_INV_ION_MOBILITY=('INV_ION_MOBILITY', 'median')).reset_index()
    # Get the CHARGE, MASS_ANALYZER, RAW_FILE, and FRAGMENTATION from the msms.txt file
    df_msms_scans = pd.merge(df_scans, df_msms, on="SCAN_NUMBER")
    df_msms_scans = df_msms_scans[["RAW_FILE", "SCAN_NUMBER", "combined_INTENSITIES",
            "combined_MZ", "MASS_ANALYZER", "FRAGMENTATION",
            "median_RETENTION_TIME", "median_INV_ION_MOBILITY",
            "median_CE", "CHARGE"]]
    return df_msms_scans

def binning(inp, ignoreCharges, rescoring_path):
    ms = MasterSpectrum()
    ms.load_from_tims(inp, ignoreCharges)
    ms.export_to_csv(rescoring_path + "/temp.txt")
    comb_ms = pd.read_csv(rescoring_path + "/temp.txt")
    scan = inp["SCAN_NUMBER"]
    comb_ms["SCAN_NUMBER"] = scan
    comb_ms = comb_ms.drop(columns=["counts", "left border", "right border", "start_mz", "ms1_charge", "rel_intensity_ratio", "counts_ratio"])
    return comb_ms

def combine_spectra(df_msms_scans, rescoring_path, chunk_size = 1000):
    chunk_list = []
    # Split both DataFrames into chunks
    chunks = [df_msms_scans[i:i + chunk_size] for i in range(0, len(df_msms_scans), chunk_size)]
    for chunk in chunks:
        bin_result_list = []
        for index, line in chunk.iterrows():
            bin_result = binning(line, True, rescoring_path)
            bin_result_list.append(bin_result)
        bin_result_df = pd.concat(bin_result_list)
        bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
        scans_combined = pd.merge(chunk, bin_result_df_collapsed, on="SCAN_NUMBER")
        scans_comb = scans_combined.drop(columns = ["combined_INTENSITIES", "combined_MZ"]).rename(columns={"intensity":"INTENSITIES", "mz":"MZ", "median_CE":"COLLISION_ENERGY"})
        # Convert lists into arrays
        scans_comb['INTENSITIES'] = scans_comb['INTENSITIES'].apply(lambda x: np.array(x))
        scans_comb['MZ'] = scans_comb['MZ'].apply(lambda x: np.array(x))
        # Sort the MZ (and linked INTENSITIES)
        for i in range (len(scans_comb)):
            zipped_list = zip(scans_comb.iloc[i]["MZ"], scans_comb.iloc[i]["INTENSITIES"])
            sorted_pair = sorted(zipped_list)
            tuples = zip(*sorted_pair)
            list_1, list_2 = [np.array(tuple) for tuple in tuples]
            scans_comb.at[i, "MZ"] = list_1
            scans_comb.at[i, "INTENSITIES"] = list_2
        scans_comb["MZ_RANGE"] = "0"
        chunk_list.append(scans_comb)
    chunk_comb = pd.concat(chunk_list, ignore_index=True)
    return chunk_comb

def convert_d_pkl(
    input_path: Union[Path, str],
    txt_path: Union[Path, str],
    output_path: Optional[Union[Path, str]] = None):
    """Converts a .d folder to pkl.
    :param input_path: path of the .d folder
    :param txt_path: path of the txt folder from MaxQuant,
    containing msms.txt, accumulatedMsmsScans.txt, and pasefMsmsScans.txt
    :param output_path: file path of the pkl file
    :return: path to converted file as string
    """
    df_msms_scans = load_timstof(input_path, txt_path)
    df_combined = combine_spectra(df_msms_scans, output_path)
    # Write to pickle
    file_name = df_combined["RAW_FILE"][0]
    df_combined.to_pickle(output_path + "/" + file_name + ".pkl")
