import logging
import warnings
from abc import abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from xml.etree import ElementTree

import pandas as pd
import pymzml
from pyteomics import mzml
from spectrum_fundamentals.constants import MZML_DATA_COLUMNS

logger = logging.getLogger(__name__)


def check_analyzer(mass_analyzers: Dict[str, str]) -> Dict[str, str]:
    """
    Convert mass analyzer accession ids to internal format.

    :param mass_analyzers: dictionary with instrumentConfigurationRef, analyzer accession
    :raises AssertionError: if the mass analyzer metadata cannot be found in the file or the search
        was conducted with an unsupported mass analyzer.
    :return: dictionary with instrumentConfigurationRef, one of (ITMS, FTMS, TOF)
    """
    for elem in mass_analyzers.keys():
        accession = mass_analyzers[elem]
        if accession in ["MS:1000079", "MS:1000484"]:  # fourier transform ion cyclotron, orbitrap
            mass_analyzers[elem] = "FTMS"
        elif accession in ["MS:1000082", "MS:1000264"]:  # quadrupole ion-trap, io-trap
            mass_analyzers[elem] = "ITMS"
        elif accession in ["MS:1000084"]:  # TOF
            mass_analyzers[elem] = "TOF"
        else:
            raise AssertionError(f"The mass analyzer with accession {accession} is not supported.")
    return mass_analyzers


def get_mass_analyzer(file_path: Path) -> Dict[str, str]:
    """
    Retrieve mass analyzer information from mzml file.

    This is using the description of the mzml format to check for specific accessions in the mzml file
    that are not covered by pyteomics or pymzml. The documentation can be found here:
    https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo

    :param file_path: The path to the mzml file to parse
    :return: A dictionary with instrumentConfigurationId, mass analyzer (ITMS, FTMS or TOF) to represent
        the respective mass analyzer category for each MS level that is present in the mzml file, i.e. MS1/MS2/MS3.
    """
    context = ElementTree.iterparse(file_path, events=("start", "end"))
    _, root = next(context)  # Get the root element

    namespace = "{http://psi.hupo.org/ms/mzml}"
    mass_analyzers = {}

    config_id = None
    accession = None
    within_instrument_configuration_list = False
    within_instrument_configuration = False

    for event, element in context:
        if event == "start" and element.tag == f"{namespace}instrumentConfigurationList":
            within_instrument_configuration_list = True
        elif event == "end" and element.tag == f"{namespace}instrumentConfigurationList":
            within_instrument_configuration_list = False
            break
        elif (
            within_instrument_configuration_list
            and event == "start"
            and element.tag == f"{namespace}instrumentConfiguration"
        ):
            within_instrument_configuration = True
            config_id = element.attrib.get("id")
        elif (
            within_instrument_configuration
            and element.tag == f"{namespace}analyzer"
            and element.attrib.get("order") == "2"
            and event == "start"
        ):
            accession = None  # Reset accession
        elif (
            within_instrument_configuration
            and element.tag == f"{namespace}cvParam"
            and element.get("accession") is not None
        ):
            if accession is None:
                accession = element.attrib.get("accession")
        elif (
            within_instrument_configuration and event == "end" and element.tag == f"{namespace}instrumentConfiguration"
        ):
            within_instrument_configuration = False
            if config_id is not None and accession is not None:
                mass_analyzers[config_id] = accession

    return check_analyzer(mass_analyzers)


class MSRaw:
    """Main to read mzml file and generate dataframe containing intensities and m/z values."""

    def __init__(self, path: Optional[Union[str, Path]] = None, output_path: Optional[Union[str, Path]] = None):
        """
        Initialize a MSRaw object.

        :param path: path to mzml file
        :param output_path: path to save the output file
        """
        if isinstance(path, str):
            path = Path(path)
        if isinstance(output_path, str):
            output_path = Path(output_path)
        self.path = path
        self.output_path = output_path

    @staticmethod
    def read_mzml(
        source: Union[str, Path, List[Union[str, Path]]],
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
        file_list = MSRaw.get_file_list(source, ext)
        data = {}  # type: Dict[str, Any]
        if package == "pymzml":
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ImportWarning)
                for file_path in file_list:
                    logger.info(f"Reading mzML file: {file_path}")
                    MSRaw._get_scans_pymzml(file_path, data, scanidx, *args, **kwargs)
        elif package == "pyteomics":
            for file_path in file_list:
                mass_analyzer = get_mass_analyzer(file_path)
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = mzml.read(source=str(file_path), *args, **kwargs)
                file_name = file_path.stem
                for spec in data_iter:
                    if spec["ms level"] != 1:  # filter out ms1 spectra if there are any
                        spec_id = spec["id"].split("scan=")[-1]
                        instrument_configuration_ref = spec["scanList"]["scan"][0]["instrumentConfigurationRef"]
                        fragmentation = spec["scanList"]["scan"][0]["filter string"].split("@")[1][:3].upper()
                        mz_range = spec["scanList"]["scan"][0]["filter string"].split("[")[1][:-1]
                        rt = spec["scanList"]["scan"][0]["scan start time"]
                        key = f"{file_name}_{spec_id}"
                        data[key] = [
                            file_name,
                            spec_id,
                            spec["intensity array"],
                            spec["m/z array"],
                            mz_range,
                            rt,
                            mass_analyzer[instrument_configuration_ref],
                            fragmentation,
                        ]
                data_iter.close()
        else:
            raise AssertionError("Choose either 'pymzml' or 'pyteomics'")

        data = pd.DataFrame.from_dict(data, orient="index", columns=MZML_DATA_COLUMNS)
        data["SCAN_NUMBER"] = pd.to_numeric(data["SCAN_NUMBER"])
        return data

    @staticmethod
    def get_file_list(source: Union[str, Path, List[Union[str, Path]]], ext: str = "mzml") -> List[Path]:
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
                file_list = list(source.glob("*[mM][zZ][mM][lL]"))
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
    def _get_scans_pymzml(
        file_path: Union[str, Path], data: Dict, scanidx: Optional[List] = None, *args, **kwargs
    ) -> None:
        """
        Reads mzml and generates a dataframe containing intensities and m/z values.

        :param file_path: path to a single mzml file.
        :param data: dictionary to be added to by this function
        :param scanidx: optional list of scan numbers to extract. if not specified, all scans will be extracted
        :param args: additional positional arguments
        :param kwargs: additional keyword arguments
        """
        if isinstance(file_path, str):
            file_path = Path(file_path)
        data_iter = pymzml.run.Reader(file_path, args=args, kwargs=kwargs)
        file_name = file_path.stem
        mass_analyzer = get_mass_analyzer(file_path)
        if scanidx is None:
            for spec in data_iter:
                if spec.ms_level != 1:  # filter out ms1 spectra if there are any
                    key = f"{file_name}_{spec.ID}"
                    instrument_configuration_ref = spec["scanList"]["scan"][0]["instrumentConfigurationRef"]
                    filter_string = str(spec.element.find(".//*[@accession='MS:1000512']").get("value"))
                    fragmentation = filter_string.split("@")[1][:3].upper()
                    mz_range = filter_string.split("[")[1][:-1]
                    data[key] = [
                        file_name,
                        spec.ID,
                        spec.i,
                        spec.mz,
                        mz_range,
                        spec.scan_time_in_minutes(),
                        mass_analyzer[instrument_configuration_ref],
                        fragmentation,
                    ]
        else:
            for idx in scanidx:
                spec = data_iter[idx]
                # this does not work if some spectra are filtered out, e.g. mzML files with only MS2 spectra, see:
                # https://github.com/pymzml/pymzML/blob/a883ff0e61fd97465b0a74667233ff594238e335/pymzml/file_classes
                # /standardMzml.py#L81-L84
                key = f"{file_name}_{spec.ID}"
                instrument_configuration_ref = spec["scanList"]["scan"][0]["instrumentConfigurationRef"]
                filter_string = str(spec.element.find(".//*[@accession='MS:1000512']").get("value"))
                fragmentation = filter_string.split("@")[1][:3].upper()
                mz_range = filter_string.split("[")[1][:-1]
                data[key] = [
                    file_name,
                    spec.ID,
                    spec.i,
                    spec.mz,
                    mz_range,
                    spec.scan_time_in_minutes(),
                    mass_analyzer[instrument_configuration_ref],
                    fragmentation,
                ]
        data_iter.close()
