from abc import abstractmethod
from typing import Optional, Union

import numpy as np
import pandas as pd


class SpectralLibrary:
    """Main to initialze a SpectralLibrary obj."""

    # Check https://gitlab.lrz.de/proteomics/prosit_tools/converter for old code
    spectra_input: pd.DataFrame
    grpc_output: dict
    spectra_output: pd.DataFrame
    out_path: str

    def __init__(self, input_dataframe: pd.DataFrame, grpc_dict: dict, output_path: str):
        """
        Initialize a SpectralLibrary obj.

        :param input_dataframe: dataframe of sequences, charges, and masses of all library peptides
        :param grpc_dict: GRPC client output dictionary with spectrum, irt, and proteotypicity prediction
        :param output_path: path to output file including file name
        """
        self.spectra_input = input_dataframe
        self.grpc_output = grpc_dict
        self.out_path = output_path

    def load(self):
        """Load predictions from hdf5 file."""

    @abstractmethod
    def write(self, chunksize: Optional[Union[None, int]]):
        """Write predictions."""
        pass

    @abstractmethod
    def prepare_spectrum(self):
        """Prepare spectrum."""
        pass
