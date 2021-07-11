from abc import abstractmethod


class MSRaw:
    path: str
    output_path: str

    def __init__(self, path, output_path):
        self.path = path
        self.output_path = output_path

    @abstractmethod
    def raw_mzml(self):
        """
        Use https://github.com/compomics/ThermoRawFileParser for conversion
        """
        pass

    def read_mzml(self):
        """
        read mzml and generate peaks for each mzml file.
        """
        pass

    def read_peaks(self):
        """
        read peaks generated for each raw file.
        """
        pass