from abc import abstractmethod

import pandas as pd
import numpy as np


class SpectralLibrary:
    # Check https://gitlab.lrz.de/proteomics/prosit_tools/converter for old code
    path: str
    out_path: str
    spectra_output: pd.DataFrame

    def __init__(self, path, output_path):
        self.path = path
        self.out_path = output_path

    def load(self):
        """
        Load predictions from hdf5 file
        """

    @abstractmethod
    def write(self):
        pass

    @abstractmethod
    def prepare_spectrum(self):
        pass
