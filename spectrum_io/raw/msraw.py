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
    :return: dictionary with instrumentConfigurationRef, one of "ITMS", "FTMS", "TOF" or "unknown",
        in case the analyzer accession cannot be mapped to one of the three groups.
    """
    for elem in mass_analyzers.keys():
        accession = mass_analyzers[elem]
        if accession in ["MS:1000079", "MS:1000484"]:  # fourier transform ion cyclotron, orbitrap
            mass_analyzers[elem] = "FTMS"
        elif accession in ["MS:1000082", "MS:1000264", "MS:1000078"]:  # quadrupole ion-trap, ion-trap, linear ion-trap
            mass_analyzers[elem] = "ITMS"
        elif accession in ["MS:1000084"]:  # TOF
            mass_analyzers[elem] = "TOF"
        else:
            mass_analyzers[elem] = "unsupported"
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
        :param args: additional positional arguments
        :param kwargs: additional keyword arguments
        :raises AssertionError: if package has an unexpected type
        :return: pd.DataFrame with intensities and m/z values
        """
        file_list = MSRaw.get_file_list(source, ext)

        if package == "pymzml":
            data = MSRaw._read_mzml_pymzml(file_list, scanidx, *args, **kwargs)
        elif package == "pyteomics":
            data = MSRaw._read_mzml_pyteomics(file_list, *args, **kwargs)
        else:
            raise AssertionError("Choose either 'pymzml' or 'pyteomics'")

        data["SCAN_NUMBER"] = pd.to_numeric(data["SCAN_NUMBER"])
        return data

    @staticmethod
    def _read_mzml_pymzml(file_list: List[Path], scanidx: Optional[List] = None, *args, **kwargs) -> pd.DataFrame:
        data_dict = {}
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ImportWarning)
            for file_path in file_list:
                logger.info(f"Reading mzML file: {file_path}")
                data_iter = pymzml.run.Reader(file_path, args=args, kwargs=kwargs)
                file_name = file_path.stem
                mass_analyzer = get_mass_analyzer(file_path)
                namespace = "{http://psi.hupo.org/ms/mzml}"
                instrument_name = data_iter.info["referenceable_param_group_list_element"][0][0].get("name")

                if scanidx is None:
                    spectra = data_iter
                else:
                    # this does not work if some spectra are filtered out, e.g. mzML files with only MS2 spectra, see:
                    # https://github.com/pymzml/pymzML/blob/a883ff0e61fd97465b0a74667233ff594238e335/pymzml/file_classes
                    # /standardMzml.py#L81-L84
                    spectra = (data_iter[idx] for idx in scanidx)

                for spec in spectra:
                    if spec.ms_level != 2:
                        continue  # filter out ms1 spectra if there are any
                    key = f"{file_name}_{spec.ID}"
                    scan = spec.get_element_by_path(["scanList", "scan"])[0]
                    instrument_configuration_ref = scan.get("instrumentConfigurationRef", "")
                    activation = spec.get_element_by_path(["precursorList", "precursor", "activation"])[0]
                    fragmentation = "unknown"
                    collision_energy = 0.0
                    for cv_param in activation:
                        name = cv_param.get("name")
                        if name == "collision energy":
                            collision_energy = float(cv_param.get("value"))
                            continue
                        if "beam-type" in name:
                            fragmentation = "HCD"
                        elif "collision-induced dissociation" in name:
                            fragmentation = "CID"
                        else:
                            fragmentation = name
                    scan_window = scan.find(f".//{namespace}scanWindow")
                    scan_lower_limit = float(
                        scan_window.find(f'./{namespace}cvParam[@accession="MS:1000501"]').get("value")
                    )
                    scan_upper_limit = float(
                        scan_window.find(f'./{namespace}cvParam[@accession="MS:1000500"]').get("value")
                    )
                    mz_range = f"{scan_lower_limit}-{scan_upper_limit}"
                    data_dict[key] = [
                        file_name,
                        spec.ID,
                        spec.i,
                        spec.mz,
                        mz_range,
                        spec.scan_time_in_minutes(),
                        mass_analyzer.get(instrument_configuration_ref, "unknown"),
                        fragmentation,
                        collision_energy,
                        instrument_name,
                    ]
                data_iter.close()
        data = pd.DataFrame.from_dict(data_dict, orient="index", columns=MZML_DATA_COLUMNS)
        return data

    @staticmethod
    def _read_mzml_pyteomics(file_list: List[Path], *args, **kwargs) -> pd.DataFrame:
        data_dict = {}
        for file_path in file_list:
            mass_analyzer = get_mass_analyzer(file_path)
            logger.info(f"Reading mzML file: {file_path}")
            data_iter = mzml.read(str(file_path), *args, **kwargs)
            file_name = file_path.stem
            try:
                instrument_params = data_iter.get_by_id("commonInstrumentParams")
            except KeyError:
                instrument_params = data_iter.get_by_id("CommonInstrumentParams")
            instrument_name = list(instrument_params.keys())[1]
            for spec in data_iter:
                if spec["ms level"] != 2:
                    continue  # filter out ms1 spectra if there are any
                spec_id = spec["id"].split("scan=")[-1]
                scan = spec["scanList"]["scan"][0]
                instrument_configuration_ref = scan.get("instrumentConfigurationRef", "")
                activation = spec["precursorList"]["precursor"][0]["activation"]
                fragmentation = "unknown"
                collision_energy = 0.0
                for key, value in activation.items():
                    if key == "collision energy":
                        collision_energy = value
                    elif "beam-type" in key:
                        fragmentation = "HCD"
                    elif "collision-induced dissociation" in key:
                        fragmentation = "CID"
                    else:
                        fragmentation = key
                scan_lower_limit = scan["scanWindowList"]["scanWindow"][0]["scan window lower limit"]
                scan_upper_limit = scan["scanWindowList"]["scanWindow"][0]["scan window upper limit"]
                mz_range = f"{scan_lower_limit}-{scan_upper_limit}"
                rt = spec["scanList"]["scan"][0]["scan start time"]
                key = f"{file_name}_{spec_id}"
                data_dict[key] = [
                    file_name,
                    spec_id,
                    spec["intensity array"],
                    spec["m/z array"],
                    mz_range,
                    rt,
                    mass_analyzer.get(instrument_configuration_ref, "unknown"),
                    fragmentation,
                    collision_energy,
                    instrument_name,
                ]
            data_iter.close()
        data = pd.DataFrame.from_dict(data_dict, orient="index", columns=MZML_DATA_COLUMNS)
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
