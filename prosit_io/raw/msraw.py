import logging
import os
from abc import abstractmethod
from typing import Union, List, Optional
import pandas as pd
from fundamentals.constants import MZML_DATA_COLUMNS

logger = logging.getLogger(__name__)

class MSRaw:
    path: str
    output_path: str

    def __init__(self, path: Optional[str] = None, output_path: Optional[str] = None):
        self.path = path
        self.output_path = output_path

    @abstractmethod
    def convert_raw_mzml(self, input_path, output_path):
        """
        Use https://github.com/compomics/ThermoRawFileParser for conversion
        """
        raise NotImplementedError

    @staticmethod
    def read_mzml(
        source: Union[str, List[str]],
        ext: str = 'mzml',
        package: str = 'pymzml',
        scanidx: Optional[List] = None,
        *args,
        **kwargs
    ):
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param source: A directory containing mzml files, a list of files or a single file.
        :param ext: File extension for searching a specified directory.
        :param package: Package for parsing the mzml file. Can eiter be "pymzml" or "pyteomics"

        :return: Pandas DataFrame
        """

        if isinstance(source, str):
            file_list = []
            if os.path.isdir(source):
                # if string is provided and is a directory, search all mzml files with provided extension
                for file in os.listdir(source):
                    if file.lower().endswith(ext.lower()):
                        file_list.append(file)

            else:
                file_list = [source]
            source = file_list
        data = {}
        if package == 'pymzml':
            import pymzml
            for file_path in source:
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = pymzml.run.Reader(file_path, args=args, kwargs=kwargs)
                file_name = os.path.splitext(os.path.basename(file_path))[0]
                if scanidx is None:
                    for spec in data_iter:
                        key = f"{file_name}_{spec.ID}"
                        data[key] = [file_name, spec.ID, spec.i, spec.mz]
                else:
                    for idx in scanidx:
                        spec = data_iter[idx]
                        key = f"{file_name}_{spec.ID}"
                        data[key] = [file_name, spec.ID, spec.i, spec.mz]


        elif package == 'pyteomics':
            from pyteomics import mzml
            for file_path in source:
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = mzml.read(source=file_path, *args, **kwargs)
                file_name = os.path.splitext(os.path.basename(file_path))[0]
                for spec in data_iter:
                    id = spec['id'].split('scan=')[-1]
                    key = f"{file_name}_{id}"
                    data[key] = [file_name, id, spec['intensity array'], spec['m/z array']]
        else:
            assert False, "Choose either 'pymzml' or 'pyteomics'"

        data = pd.DataFrame.from_dict(data, orient='index', columns=MZML_DATA_COLUMNS)
        data['SCAN_NUMBER'] = pd.to_numeric(data['SCAN_NUMBER'])
        return data
