import logging
import os
import warnings
from abc import abstractmethod
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import pymzml
from pyteomics import mzml
from spectrum_fundamentals.constants import MZML_DATA_COLUMNS

logger = logging.getLogger(__name__)


class MSRaw:
    """Main to read mzml file and generate dataframe containing intensities and m/z values."""

    path: Optional[str]
    output_path: Optional[str]

    def __init__(self, path: Optional[str] = None, output_path: Optional[str] = None):
        """
        Initialize a MSRaw object.

        :param path: path to mzml file
        :param output_path: path to save the output file
        """
        self.path = path
        self.output_path = output_path

    @staticmethod
    def read_mzml(
        source: Union[str, List[str]],
        ext: str = "mzml",
        package: str = "pyteomics",
        search_type: str = "Maxquant",
        scanidx: Optional[List] = None,
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param source: a directory containing mzml files, a list of files or a single file
        :param ext: file extension for searching a specified directory
        :param package: package for parsing the mzml file. Can eiter be "pymzml" or "pyteomics"
        :param scanidx: optional list of scan numbers to extract. if not specified, all scans will be extracted
        :param search_type: type of the search (Maxquant, Mascot, Msfragger)
        :param args: additional positional arguments
        :param kwargs: additional keyword arguments
        :raises AssertionError: if package has an unexpected type
        :return: pd.DataFrame with intensities and m/z values
        """
        if isinstance(source, str):
            file_list = MSRaw.get_file_list(source, ext)
            source = file_list
        data = {}  # type: Dict[str, Any]
        if package == "pymzml":
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ImportWarning)
                for file_path in source:
                    logger.info(f"Reading mzML file: {file_path}")
                    MSRaw._get_scans_pymzml(file_path, data, scanidx, *args, **kwargs)
        elif package == "pyteomics":
            for file_path in source:
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = mzml.read(source=file_path, *args, **kwargs)
                file_name = os.path.splitext(os.path.basename(file_path))[0]
                for spec in data_iter:
                    id = spec["id"].split("scan=")[-1]
                    mass_analyzer = spec["scanList"]["scan"][0]["filter string"].split()[0]
                    fragmentation = spec["scanList"]["scan"][0]["filter string"].split("@")[1][:3]
                    mz_range = spec["scanList"]["scan"][0]["filter string"].split("[")[1][:-1]
                    key = f"{file_name}_{id}"
                    if search_type == "maxquant":
                        data[key] = [file_name, id, spec["intensity array"], spec["m/z array"], mz_range]
                    else:
                        data[key] = [
                            file_name,
                            id,
                            spec["intensity array"],
                            spec["m/z array"],
                            mz_range,
                            mass_analyzer,
                            fragmentation,
                        ]
                data_iter.close()
        else:
            raise AssertionError("Choose either 'pymzml' or 'pyteomics'")
        if search_type == "maxquant":
            data = pd.DataFrame.from_dict(data, orient="index", columns=MZML_DATA_COLUMNS)
        else:
            data = pd.DataFrame.from_dict(
                data, orient="index", columns=MZML_DATA_COLUMNS + ["MASS_ANALYZER", "FRAGMENTATION"]
            )
        data["SCAN_NUMBER"] = pd.to_numeric(data["SCAN_NUMBER"])
        return data

    @staticmethod
    def get_file_list(source: Union[str, List[str]], ext: str = "mzml"):
        """
        Get list of files from source.

        :param source: a directory containing mzml files, a list of files or a single file
        :param ext: file extension for searching a specified directory
        :return: list of files
        """
        file_list = []
        if isinstance(source, str) and os.path.isdir(source):
            # if string is provided and is a directory, search all mzml files with provided extension
            for file in os.listdir(source):
                if file.lower().endswith(ext.lower()):
                    file_list.append(file)

        else:
            if isinstance(source, list):
                file_list.extend(source)
            else:
                file_list.append(source)
        return file_list

    @staticmethod
    def _get_scans_pymzml(file_path: str, data: Dict, scanidx: Optional[List] = None, *args, **kwargs) -> None:
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param file_path: path to a single mzml file.
        :param data: dictionary to be added to by this function
        :param scanidx: optional list of scan numbers to extract. if not specified, all scans will be extracted
        :param args: additional positional arguments
        :param kwargs: additional keyword arguments
        """
        data_iter = pymzml.run.Reader(file_path, args=args, kwargs=kwargs)
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        if scanidx is None:
            for spec in data_iter:
                key = f"{file_name}_{spec.ID}"
                data[key] = [file_name, spec.ID, spec.i, spec.mz]
        else:
            for idx in scanidx:
                spec = data_iter[idx]
                # this does not work if some spectra are filtered out, e.g. mzML files with only MS2 spectra, see:
                # https://github.com/pymzml/pymzML/blob/a883ff0e61fd97465b0a74667233ff594238e335/pymzml/file_classes
                # /standardMzml.py#L81-L84
                key = f"{file_name}_{spec.ID}"
                data[key] = [file_name, spec.ID, spec.i, spec.mz]
        data_iter.close()
