import logging
import subprocess
from pathlib import Path
from sys import platform
from typing import Any, List, Optional, Tuple, Union

from .msraw import MSRaw

logger = logging.getLogger(__name__)


def _type_check(var: Any, varname: str, types: Union[type, Tuple[type, ...]]):
    if isinstance(var, types):
        return
    if isinstance(types, type):
        possible_types_str = types.__name__
    else:
        type_names = [t.__name__ for t in types]
        possible_types_str = "{} or {}".format(", ".join(type_names[:-1]), type_names[-1])
    raise TypeError(f"{varname} must be of type {possible_types_str}")


def _assemble_arg_list(
    input_path: Path, output_path: Path, ms_level: List[int], gzip: bool, thermo_exe: Path
) -> List[Union[str, Path]]:
    exec_arg_list: List[Union[str, Path]] = [
        thermo_exe,
        f"--msLevel={','.join([str(level) for level in ms_level])}",
        "-i",
        input_path.resolve(),
        "-b",
        output_path,
    ]
    if gzip:
        exec_arg_list.append("-g")
    if "linux" in platform or platform == "darwin":
        exec_arg_list.insert(0, "mono")

    return exec_arg_list


class ThermoRaw(MSRaw):
    """Main to convert a ThermoRaw file into mzml file."""

    @staticmethod
    def convert_raw_mzml(
        input_path: Union[Path, str],
        gzip: bool = False,
        ms_level: Union[int, List[int]] = 2,
        output_path: Optional[Union[Path, str]] = None,
        thermo_exe: Union[Path, str] = "ThermoRawFileParser.exe",
    ) -> Path:
        """Converts a ThermoRaw file to mzML.

        Use https://github.com/compomics/ThermoRawFileParser for conversion.

        :param input_path: file path of the Thermo Rawfile
        :param gzip: whether to gzip the file
        :param ms_level: level of MS, can be a single integer (1, 2, 3) or any combination of that provided as a list
        :param output_path: file path of the mzML path
        :param thermo_exe: path to the executable of ThermoRawFileParser. Default: ThermoRawFileParser.exe
        :raises subprocess.CalledProcessError: if the subprocess for conversion failed
        :raises ValueError: if ms_level(s) provided are other than 1, 2 or 3.
        :return: path to converted file as string
        """
        _type_check(input_path, "input_path", (Path, str))
        input_path = Path(input_path)
        if output_path is None:
            output_path = input_path.with_suffix(".mzML")
        _type_check(output_path, "output_path", (Path, str))
        output_path = Path(output_path)

        _type_check(thermo_exe, "thermo_exe", (Path, str))
        thermo_exe = Path(thermo_exe)

        _type_check(ms_level, "ms_level", (int, list))
        if isinstance(ms_level, int):
            ms_level = [ms_level]
        for level in ms_level:
            _type_check(level, "all ms_levels in list", int)
            if not 1 <= level <= 3:
                raise ValueError(f"Value of all ms_levels must be within [1,3]. Got {level}")

        if output_path.is_file():
            logger.info(f"Found converted file at {output_path}, skipping conversion")
            return output_path

        exec_arg_list = _assemble_arg_list(input_path, output_path, ms_level, gzip, thermo_exe)

        logger.info(
            f"Converting thermo rawfile to mzml with the command: {' '.join([str(arg) for arg in exec_arg_list])}"
        )

        try:
            subprocess.run(exec_arg_list, shell=False, check=True)
        except subprocess.CalledProcessError:
            if output_path.is_file():
                output_path.unlink()
            raise  # reraise only after removing a corrupted file

        return output_path


if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2:
        converter = ThermoRaw().convert_raw_mzml(argv[1], ms_level=[1, 2])
    else:
        print("Please specify a rawfile")
