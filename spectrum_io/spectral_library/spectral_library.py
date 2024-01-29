from abc import abstractmethod
from multiprocessing import Queue
from multiprocessing.managers import ValueProxy
from pathlib import Path
from typing import IO, Dict, Optional, Union

import numpy as np
import pandas as pd


class SpectralLibrary:
    """Main to initialze a SpectralLibrary obj."""

    def __init__(
        self,
        output_path: Union[str, Path],
        mode: str = "w",
        min_intensity_threshold: float = 5e-4,
        chunksize: Optional[int] = None,
    ):
        """
        Initialize a SpectralLibrary obj.

        :param output_path: path to output file including file name
        :param mode: Whether to append ('a') to or overwrite ('w') an extisting file
            at the provided output path (if it is present).
        :param min_intensity_threshold: optional filter for low intensity peaks
        :param chunksize: optional chunksize for dlib
        """
        if isinstance(output_path, str):
            output_path = Path(output_path)
        self.out_path = output_path
        self.mode = mode
        self.min_intensity_threshold = min_intensity_threshold
        self.chunksize = chunksize

    def load(self):
        """Load predictions from hdf5 file."""

    def write(self, *args, **kwargs):
        """
        Write content to the output file.

        :param args: Positional arguments to be passed to the internal _write method.
        :param kwargs: Keyword arguments to be passed to the internal _write method.
        """
        with open(self.out_path, self.mode) as out:
            self._write_header(out)
            self._write(out, *args, **kwargs)

    def async_write(self, queue: Queue, progress: ValueProxy):
        """
        Asynchronously write content to the output file from a queue.

        :param queue: A queue from which content will be retrieved for writing.
        :param progress: An integer value representing the progress of the writing process.
        """
        with open(self.out_path, self.mode) as out:
            self._write_header(out)
            while True:
                content = queue.get()
                if content is None:
                    break
                self._write(out, *content)
                progress.value += 1

    def _fragment_filter_passed(
        self, f_mz: Union[np.ndarray, float], f_int: Union[np.ndarray, float]
    ) -> Union[np.ndarray, bool]:
        """
        Return if a fragment passes a common filter.

        This function returns a boolean to determine if a fragment can exist (as defined by its mz != -1)
        and if its (predicted) intensity is at least the minimal threshold given during instantiation of
        the writer class.

        :param f_mz: The fragment's mass to charge ratio
        :param f_int: The fragment's (predicted) intensity

        :return: boolean (array) that determines if the fragment passes the filter
        """
        return (f_mz != -1) & (f_int >= self.min_intensity_threshold)

    @abstractmethod
    def _write(self, out: IO, data: Dict[str, np.ndarray], metadata: pd.DataFrame):
        """
        Internal writer function.

        This function takes a batch of data and corresponding metadata and writes it to the provided
        output file handle. The file handle logic is provided using the public write or async_write
        functions.

        :param out: file handle accepting the data to be written to disk
        :param data: Dictionary containing TODO keys and corresponding values as numpy array
        :param metadata: a dataframe that contains the columns TODO
        """
        pass

    @abstractmethod
    def _write_header(self, out: IO):
        pass

    @abstractmethod
    def prepare_spectrum(self):
        """Prepare spectrum."""
        pass
