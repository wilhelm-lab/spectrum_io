import logging
import os
import pathlib
import subprocess
from sys import platform
from typing import Optional, Union

from .msraw import MSRaw

logger = logging.getLogger(__name__)


class ThermoRaw(MSRaw):
    """Main to convert a ThermoRaw file into mzml file."""

    @staticmethod
    def convert_raw_mzml(
        input_path: Union[pathlib.Path, str],
        gzip: bool = False,
        ms_level: str = "2",
        output_path: Optional[Union[pathlib.Path, str]] = None,
    ) -> str:
        """
        Converts a ThermoRaw file to mzML.

        :param input_path: file path of the Thermo Rawfile
        :param gzip: whether to gzip the file
        :param ms_level: level of MS
        :param output_path: file path of the mzML path
        :raises subprocess.CalledProcessError: if the subprocess for conversion failed
        :return: path to converted file as string
        """
        if isinstance(input_path, str):
            input_path = pathlib.Path(input_path)
        if output_path is None:
            output_path = input_path.with_suffix(".mzML")
        if isinstance(output_path, str):
            output_path = pathlib.Path(output_path)

        if os.path.isfile(output_path):
            logger.info(f"Found converted file at {output_path}, skipping conversion")
            return output_path

        exec_path = pathlib.Path(__file__).parent.absolute()  # get path of parent directory of this file
        exec_path /= "/utils/ThermoRawFileParser/ThermoRawFileParser.exe"

        exec_arg_list = [exec_path, f"{'-g ' if gzip else ''}--msLevel {ms_level} -i", input_path, "-b", output_path]
        if "linux" in platform:
            exec_arg_list = ["mono"] + exec_arg_list

        logger.info(
            f"Converting thermo rawfile to mzml with the command: {' '.join([str(arg) for arg in exec_arg_list])}"
        )

        try:
            subprocess.run(exec_arg_list, shell=True, check=True)
        except subprocess.CalledProcessError:
            if os.path.isfile(output_path):
                os.remove(output_path)
            raise  # reraise only after removing a corrupted file

        return output_path


if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2:
        converter = ThermoRaw().convert_raw_mzml(argv[1], ms_level="1,2")
    else:
        print("Please specify a rawfile")
