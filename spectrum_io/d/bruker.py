import logging
import os
from pathlib import Path
from typing import List, Optional, Tuple, Union

import alphatims
import alphatims.bruker
import alphatims.utils
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .masterSpectrum import MasterSpectrum

logger = logging.getLogger(__name__)


def binning(mzs: List[float], intensities: List[int], ignore_charges: bool) -> Tuple[List[float], List[float]]:
    """
    Perform binning on the input MasterSpectrum.

    This function loads a MasterSpectrum from the provided input, performs binning, and exports the results
    to a temporary text file in the given rescoring path. It then reads the temporary file as a DataFrame,
    modifies the DataFrame by adding SCAN_NUMBER and dropping specified columns before returning it.

    :param intensities: Input data used to perform binning.
    :param mzs: Path where the temporary file will be exported.
    :param ignore_charges: indicating whether charges should be ignored during binning.
    :return: Tuple containing the sorted list of fragment mzs and associated intensities
    """
    ms = MasterSpectrum()
    ms.load_from_tims(intensities, mzs, ignore_charges)

    mzs_out = [mp.mz for key in ms.spectrum[0].keys() for mp in ms.spectrum[0][key]]
    intensities_out = [mp.intensity for key in ms.spectrum[0].keys() for mp in ms.spectrum[0][key]]

    return mzs_out, intensities_out


def aggregate_timstof(raw_spectra: pd.DataFrame) -> pd.DataFrame:
    """
    Combine spectra from the provided pd.DataFrame and perform binning on chunks.

    This function splits the input pd.DataFrame into chunks, performs binning on each chunk of data,
    merges the binning results with the original data, processes the combined spectra, and returns
    the combined and processed spectra as a pd.DataFrame.

    :param raw_spectra: pd.DataFrame containing spectra information.
    :return: pd.DataFrame containing combined and processed spectra.
    """
    for i, (combined_intensities, combined_mzs) in tqdm(
        enumerate(zip(raw_spectra["INTENSITIES"], raw_spectra["MZ"])),
        total=len(raw_spectra),
        desc="Aggregating spectra",
    ):
        mz, intensity = binning(combined_mzs, combined_intensities, True)
        raw_spectra.at[i, "INTENSITIES"] = intensity
        raw_spectra.at[i, "MZ"] = mz

    return raw_spectra


def read_timstof(hdf_file: Path, scan_to_precursor_map: pd.DataFrame) -> pd.DataFrame:
    """
    Read selected spectra from a given timstof hdf file.

    This function queries a given hdf file for spectra that are provided within a scan to precursor map.
    #TODO elaborate

    :param hdf_file: Path to hdf file containing spectra
    :param scan_to_precursor_map: Dataframe containing metadata to select spectra

    :return: Dataframe containing the relevant spectra read from the hdf file
    """
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
    data = alphatims.bruker.TimsTOF(str(hdf_file), slice_as_dataframe=False)

    raw_idx = []
    for frames, precursors in zip(df_frame_group["FRAME"], df_frame_group["PRECURSOR"]):
        raw_idx.extend(data[frames, :, precursors])

    # read spectra
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
    df.columns = ["FRAME", "SCAN", "PRECURSOR", "RETENTION_TIME", "INV_ION_MOBILITY", "MZ", "INTENSITIES"]

    # converting RETENTION TIME from seconds to minutes
    df["RETENTION_TIME"] = df["RETENTION_TIME"].div(60)

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
                "INTENSITIES": list,
                "MZ": list,
                "RETENTION_TIME": "first",
                "COLLISION_ENERGY": "first",
                "INV_ION_MOBILITY": "first",
            }
        )
        .merge(scan_to_precursor_map.reset_index())
        .groupby("SCAN_NUMBER", as_index=False)  # aggregate PRECURSORS for same SCAN_NUMBER
        .agg(
            COLLISION_ENERGY=("COLLISION_ENERGY", "median"),
            INTENSITIES=("INTENSITIES", lambda x: [item for sublist in x for item in sublist]),
            MZ=("MZ", lambda x: [item for sublist in x for item in sublist]),
            RETENTION_TIME=("RETENTION_TIME", "median"),
            median_INV_ION_MOBILITY=("INV_ION_MOBILITY", "median"),
        )
    )

    return df_combined_grouped


def convert_d_hdf(
    input_path: Union[Path, str],
    output_path: Union[Path, str],
):
    """
    Convert a bruker d folder to hdf format.

    # TODO long description
    :param input_path: Path to the d folder to be converted
    :param output_path: Path to the desired output location of the converted hdf file
    """
    if isinstance(output_path, str):
        output_path = Path(input_path)
    if output_path.is_file():
        logger.info(f"Found converted file at {output_path}, skipping conversion")
        return
    logger.info("Converting bruker d to hdf using alphatims...")
    data = alphatims.bruker.TimsTOF(str(input_path))
    data.save_as_hdf(directory=str(output_path.parent), file_name=str(output_path.name))


def read_and_aggregate_timstof(source: Path, tims_meta_file: Path) -> pd.DataFrame:
    """
    Read raw spectra from timstof hdf spectra file and aggregate to MS2 spectra.

    :param source: Path to the hdf file
    :param tims_meta_file: Path to metadata mapping scan numbers to precursors / frames
    :return: Dataframe containing the MS2 spectra
    """
    scan_to_precursor_map = pd.read_csv(tims_meta_file)
    raw_spectra = read_timstof(source, scan_to_precursor_map)
    df_combined = aggregate_timstof(raw_spectra)
    df_combined["RAW_FILE"] = source.stem
    df_combined["MASS_ANALYZER"] = "TOF"
    df_combined["FRAGMENTATION"] = "HCD"
    df_combined["INSTRUMENT_TYPES"] = "TIMSTOF"

    return df_combined
