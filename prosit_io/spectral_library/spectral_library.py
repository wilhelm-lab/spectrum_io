from abc import abstractmethod

import pandas as pd
import numpy as np


class SpectralLibrary:
    # Check https://gitlab.lrz.de/proteomics/prosit_tools/converter for old code
    spectra_input: pd.DataFrame
    grpc_output: dict
    spectra_output: pd.DataFrame
    out_path: str

    def __init__(self, input_dataframe, grpc_dict, output_path):
        """
        @param input_dataframe: Dataframe of sequences, charges, and masses of all library peptides
        @param grpc_dict: GRPC client output dictionary with spectrum, irt, and proteotypicity prediction
        @param output_path: Path to output file including file name
        """
        self.spectra_input = input_dataframe
        self.grpc_output = grpc_dict
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
