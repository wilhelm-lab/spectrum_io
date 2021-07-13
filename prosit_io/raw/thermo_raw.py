from .msraw import MSRaw
import os
import pathlib

class ThermoRaw(MSRaw):

    @staticmethod
    def raw_mzml(input_path: str, output_path: str):
        """
        Converts a ThermoRaw file to mzML

        :param input_path: File path of the Thermo Rawfile
        :param output_path: File path of the mzML path
        """
        exec_path = pathlib.Path(__file__).parent.absolute() # get path of parent directory of current file
        exec_command = f"mono {excec_path}/utils/ThermoRawFileParser/ThermoRawFileParser.exe -i {input_path} -b {output_path}"
        os.system(exec_command)
