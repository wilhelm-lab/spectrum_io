import re
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import pandas as pd
import numpy as np
from pyteomics import mgf


class MGFparser:

    @staticmethod
    def get_file_list(source: Union[str, Path, List[Union[str, Path]]]) -> List[Path]:
        """
        Get list of files from source.

        :param source: a directory containing mzml files, a list of files or a single file
        :param ext: file extension for searching a specified directory
        :raises FileNotFoundError: if one of the files given by source does not exist
        :raises TypeError: if source is not given as a str, Path or list object
        :return: list of files
        """
        file_list = []
        if isinstance(source, str):
            source = Path(source)
        if isinstance(source, Path):
            if source.is_file():
                file_list = [source]
            elif source.is_dir():
                file_list = list(source.glob("*[mM][gG][fF]"))
            else:
                raise FileNotFoundError(f"{source} does not exist.")
        elif isinstance(source, list):
            for elem in source:
                if isinstance(elem, str):
                    elem = Path(elem)
                if elem.is_file():
                    file_list.append(elem)
                else:
                    raise FileNotFoundError(f"{elem} does not exist,")
        else:
            raise TypeError("source can only be a single str or Path or a list of files.")
        return file_list

    @staticmethod
    def parse(mgf_file: Path):
        """Read MGF file generated from MSConvert"""

        scan = lambda s: int(s[s.find("scan=") + 5 :].split('"')[0].split()[0]) if "scan=" in s else None

        with mgf.read(str(mgf_file)) as reader:
            file_name = mgf_file.name.removesuffix(".mgf")
            spectra = []
            for i, spec in enumerate(reader):
                idx = i
                title = spec["params"]["title"]
                scan_number = scan(title)

                parsed = {
                    "SCAN_IDX": np.int64(idx),
                    "SCAN_NUMBER": np.int64(scan_number) if scan_number is not None else None,
                    "RAW_FILE": str(file_name),
                    "RETENTION_TIME": np.float64(spec["params"]["rtinseconds"]),
                    "MZ": np.array(spec["m/z array"]),
                    "INTENSITIES": np.array(spec["intensity array"]),
                    "PRECURSOR_MZ": np.float64(spec["params"]["pepmass"][0]),
                    "PRECURSOR_CHARGE": np.int64(spec["params"]["charge"][0]),
                    "PRECURSOR_INTENSITY": np.float64(spec["params"]["pepmass"][1]),
                }

                spectra.append(parsed)

        spectra_df = pd.DataFrame(
            spectra,
            columns=[
                "SCAN_IDX",
                "SCAN_NUMBER",
                "RAW_FILE",
                "RETENTION_TIME",
                "MZ",
                "INTENSITIES",
                "PRECURSOR_CHARGE",
                "PRECURSOR_MZ",
                "PRECURSOR_INTENSITY",
            ],
        )

        return spectra_df

    @staticmethod
    def read_mgf(source: Union[str, Path, List[Union[str, Path]]]) -> pd.DataFrame:
        file_list = MGFparser.get_file_list(source)
        return pd.concat([MGFparser.parse(file) for file in file_list])

    def write_file():
        # mgf.write(spectra=spectra, header=header)
        pass
