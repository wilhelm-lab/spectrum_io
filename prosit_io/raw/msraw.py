from abc import abstractmethod


class MSRaw:
    path: str
    output_path: str

    def __init__(self, path=None, output_path=None):
        self.path = path
        self.output_path = output_path

    @abstractmethod
    def raw_mzml(self, input_path, output_path):
        """
        Use https://github.com/compomics/ThermoRawFileParser for conversion
        """
        raise NotImplementedError


    def read_mzml(self):
        """
        read mzml and generate peaks for each mzml file.
        columns for peaks csv [scan_number,scan_type,masses,intensities]

        writes peaks file

        """
        pass

    def read_peaks(self):
        """
        read peaks generated for each raw file.

        Read a csv for peaks

        """
        pass
