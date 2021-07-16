from typing import Optional
from .msraw import MSRaw
import os
import pathlib

class ThermoRaw(MSRaw):

    @staticmethod
    def convert_raw_mzml(input_path: str, output_path: Optional[str] = None):
        """
        Converts a ThermoRaw file to mzML

        :param input_path: File path of the Thermo Rawfile
        :param output_path: File path of the mzML path
        """
        if output_path is None:
            output_path = f"{os.path.splitext(input_path)[0]}.mzml"

        exec_path = pathlib.Path(__file__).parent.absolute() # get path of parent directory of current file
        exec_command = f"{exec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe -i {input_path} -b {output_path}"
        os.system(exec_command)
        return output_path
