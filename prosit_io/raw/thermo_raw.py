import logging
import os
import pathlib
import subprocess
from sys import platform
from typing import Optional

from .msraw import MSRaw

logger = logging.getLogger(__name__)


class ThermoRaw(MSRaw):
    """Main to convert a ThermoRaw file into mzml file."""

    @staticmethod
    def convert_raw_mzml(
        input_path: str, output_path: Optional[str] = None, gzip: bool = False, ms_level: str = "2"
    ) -> str:
        """
        Converts a ThermoRaw file to mzML.

        :param input_path: file path of the Thermo Rawfile
        :param output_path: file path of the mzML path
        :param gzip: whether to gzip the file
        :param ms_level: level of MS
        :return: path to converted file as string
        """
        if output_path is None:
            output_path = os.path.splitext(input_path)[0] + ".mzML"

        if os.path.isfile(output_path):
            logger.info(f"Found converted file at {output_path}, skipping conversion")
            return output_path

        if gzip:
            gzip_opt = "-g"
        else:
            gzip_opt = ""

        if "linux" in platform:
            mono = "mono"
        elif "win" in platform:
            mono = ""

        exec_path = pathlib.Path(__file__).parent.absolute()  # get path of parent directory of current file
        exec_command = f"{mono} {exec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe \
                        {gzip_opt} --msLevel {ms_level} -i {input_path} -b {output_path}.tmp"
        logger.info(f"Converting thermo rawfile to mzml with the command: '{exec_command}'")
        subprocess.run(exec_command, shell=True, check=True)

        # only rename the file now, so that we don't have a partially converted file if something fails
        os.rename(f"{output_path}.tmp", output_path)

        return output_path


if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2:
        converter = ThermoRaw().convert_raw_mzml(argv[1], ms_level="1,2")
    else:
        print("Please specify a rawfile")
