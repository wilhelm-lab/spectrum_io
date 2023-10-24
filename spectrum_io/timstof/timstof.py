import numpy as np
import pandas as pd
import os

import rpy2
import rpy2.robjects as robjects #ro
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import pandas2ri

from fundamentals import constants
from fundamentals.fragments import initialize_peaks
from fundamentals.annotation.annotation import annotate_spectra
from fundamentals.mod_string import maxquant_to_internal, internal_without_mods
import argparse, pathlib
from fundamentals.mod_string import parse_modstrings, maxquant_to_internal
from fundamentals.constants import ALPHABET
from mgf_filter.util import timeStamped
from mgf_filter.masterSpectrum import MasterSpectrum
import pickle

logger = logging.getLogger(__name__)

def load_timstof(txt_path, path_to_bruker_dll, d_path):
    script_directory = os.path.dirname(os.path.abspath(__file__))
    r_script = f'source("{script_directory}/extract_d.R")'
    base = importr("base")
    utils = importr("utils")
    dplyr = importr("dplyr")
    data_table = importr("data.table")
    opentimsr = importr("opentimsr")
    tidyverse = importr("tidyverse")

    pandas2ri.activate()
    robjects.r(r_script)

    extract_d = robjects.globalenv["extract_d"] # da = name of function
    res = extract_d(d_path, txt_path, path_to_bruker_dll, "AspN_F1_E1_1_4468")

    df_raw = robjects.conversion.rpy2py(res[0])
    df_pasef = robjects.conversion.rpy2py(res[1])
    scan_precursor_map = robjects.conversion.rpy2py(res[2])
    df_msms = robjects.conversion.rpy2py(res[3])

    df_raw_pasef = pd.merge(df_raw, df_pasef, on="FRAME")
    df_raw_mapped = df_raw_pasef[(df_raw_pasef['SCANNUMBEGIN'] <= df_raw_pasef['SCAN']) &
        (df_raw_pasef['SCAN'] <= df_raw_pasef['SCANNUMEND'])]

    df_raw_mapped = df_raw_mapped.sort_values(by='MZ')

    df_raw_scans = df_raw_mapped.groupby(['PRECURSOR', 'FRAME']).agg({
        'INTENSITY': list,
        'MZ': list,
        'RETENTION_TIME': 'first',
        'COLLISION_ENERGY': 'first',
        'INV_ION_MOBILITY': 'first'
    }).reset_index()

    df_scan_group = pd.merge(df_raw_scans, scan_precursor_map)

    df_scans = df_scan_group.groupby('SCAN_NUMBER').agg(
        median_CE=('COLLISION_ENERGY', 'median'),
        combined_INTENSITIES=('INTENSITY', lambda x: [item for sublist in x for item in sublist]),
        combined_MZ=('MZ', lambda x: [item for sublist in x for item in sublist]),
        median_RETENTION_TIME=('RETENTION_TIME', 'median'),
        median_INV_ION_MOBILITY=('INV_ION_MOBILITY', 'median')
    ).reset_index()

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

def combine_spectra(df_msms_scans):
    bin_result_list = []
    for index, line in df_msms_scans.iterrows():
        bin_result = binning(line, True, rescoring_path)
        bin_result_list.append(bin_result)

    bin_result_df = pd.concat(bin_result_list)
    bin_result_df_collapsed = bin_result_df.groupby("SCAN_NUMBER").agg(list)
    scans_combined = pd.merge(df_msms_scans, bin_result_df_collapsed, on="SCAN_NUMBER")
    scans_comb = scans_combined.drop(columns = ["combined_INTENSITIES", "combined_MZ"]).rename(columns={"intensity":"INTENSITIES", "mz":"MZ", "median_CE":"COLLISION_ENERGY"})

    scans_comb['INTENSITIES'] = scans_comb['INTENSITIES'].apply(lambda x: np.array(x))
    scans_comb['MZ'] = scans_comb['MZ'].apply(lambda x: np.array(x))

    total_list = []
    for i in range (len(scans_comb)):
        zipped_list = zip(scans_comb.iloc[i]["MZ"], scans_comb.iloc[i]["INTENSITIES"])
        sorted_pair = sorted(zipped_list)
        tuples = zip(*sorted_pair)
        list_1, list_2 = [np.array(tuple) for tuple in tuples]
        scans_comb.at[i, "MZ"] = list_1
        scans_comb.at[i, "INTENSITIES"] = list_2
        total_list.append(scans_comb)

    total_df = pd.concat(total_list, ignore_index=True)
    total_df["MZ_RANGE"] = "0"
    file_name = total_df["RAW_FILE"][0]
    # Write to pickle
    total_df.to_pickle(rescoring_path + "/" + file_name + ".pkl")
