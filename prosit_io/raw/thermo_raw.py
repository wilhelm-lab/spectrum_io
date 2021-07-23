from typing import Optional
from .msraw import MSRaw
import os
import pathlib
from sys import platform

import logging
logger = logging.getLogger(__name__)

class ThermoRaw(MSRaw):

    @staticmethod
    def convert_raw_mzml(input_path: str, output_path: Optional[str] = None, gzip = False):
        """
        Converts a ThermoRaw file to mzML

        :param input_path: File path of the Thermo Rawfile
        :param output_path: File path of the mzML path
        """
        if output_path is not None:
            output_path = f"-b {output_path}"
        else:
            output_path = ""

        if gzip:
            gzip = "-g"
        else:
            gzip = ""

        if "linux" in platform:
            mono = "mono"
        elif "win" in platform:
            mono = ""

        exec_path = pathlib.Path(__file__).parent.absolute() # get path of parent directory of current file
        exec_command = f"{mono} {exec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe {gzip} -i {input_path} {output_path}"
        logger.info(f"Converting thermo rawfile to mzml with the command: '{exec_command}'")
        os.system(exec_command)
        return output_path
