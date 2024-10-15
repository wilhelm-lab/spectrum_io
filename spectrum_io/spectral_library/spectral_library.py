import re
from abc import abstractmethod
from multiprocessing import Queue
from multiprocessing.managers import ValueProxy
from pathlib import Path
from sqlite3 import Connection
from typing import IO, Dict, Optional, Union

import numpy as np
import pandas as pd


def parse_mods(mods: Dict[str, int]) -> Dict[str, str]:
    """
    Parse provided mapping of custom modification pattern to ProForma standard.

    This function takes a dictionary mapping custom modification pattern for specific aminoacids (keys) to a
    UNIMOD ID (values). The pattern is translated to ProForma standard and a new dictionary mapping the custom
    modification patterns to the ProForma standard is returned.
    The pattern for the custom modifications must start with the one-letter code for an aminoacid or '^' / '$',
    to describe n- / c-terminal modifications, respectively, followed by an optional pattern (which can be
    empty).
    This means that 'X'  or 'X(custom_pattern)', is both mapped to 'X[UNIMOD:#]'.
    For the n-terminus, an additional dash will be added automatically, which maps 'X(custom_pattern)' to
    'X[UNIMOD:#]-'. If the sequence to apply the transformation on already contains the dash, it needs to be part
    of the custom_pattern (i.e. 'X(custom_pattern)-'), to avoid adding an additional dash.

    :param mods: Dictionary mapping custom modification patterns (keys) to UNIMOD IDs (values)
    :raises TypeError: if keys are not strings or values are not integers
    :raises ValueError: if the keys do not start with [A-Z,a-z,^,$]
    :return: A dictionary mapping custom modification patterns (keys) to the ProForma standard (values)
    """
    key_pattern = (
        "'X' or 'X<mod_pattern>' where X is either the one-letter code of an aminoacid or '^' / '$' defining the"
        " n- or c-terminus, respectively, followed by an optional pattern identifying a specific modification."
    )
    unimod_regex_map = {}
    for k, v in mods.items():
        if not isinstance(v, int):
            raise TypeError(f"UNIMOD id {v} for replacement {k} not understood. UNIMOD IDs must be integers.")
        if not isinstance(k, str):
            raise TypeError(
                f"Replacement {k} not understood. Replacements must be strings and follow " f"the pattern {key_pattern}"
            )
        if k[0].isalpha():
            unimod_regex_map[re.escape(f"{k[0]}[UNIMOD:{v}]")] = k
            continue
        if k[0] == "^":
            to_escape = re.escape(f"[UNIMOD:{v}]-")
            unimod_regex_map[f"^{to_escape}"] = k[1:]
            continue
        raise ValueError(
            f"Replacement {k} not understood. {k[0]} is not a valid aminoacid. "
            f"Replacements most follow the pattern {key_pattern}"
        )
    return unimod_regex_map


class SpectralLibrary:
    """Main to initialze a SpectralLibrary obj."""

    @property
    @abstractmethod
    def standard_mods(self) -> Dict[str, int]:
        """Standard modifications that are always applied if not otherwise specified."""
        pass

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

    def write(self, *, custom_mods: Optional[Dict[str, int]] = None, **kwargs):
        """
        Write content to the output file.

        :param custom_mods: optional dictionary mapping libary format-specific modification patterns (keys)
            to UNIMOD IDs (values)
        :param kwargs: Keyword arguments to be passed to the internal _write method.
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))
        with self._get_handle() as out:
            self._initialize(out)
            self._write(out, mods=parsed_mods, **kwargs)

    def _get_handle(self):
        return open(self.out_path, self.mode)

    def async_write(self, queue: Queue, progress: ValueProxy, custom_mods: Optional[Dict[str, int]] = None):
        """
        Asynchronously write content to the output file from a queue.

        :param queue: A queue from which content will be retrieved for writing.
        :param progress: An integer value representing the progress of the writing process.
        :param custom_mods: dict with custom variable and static identifier and respecitve internal equivalent and mass
        """
        parsed_mods = parse_mods(self.standard_mods | (custom_mods or {}))

        with self._get_handle() as out:
            self._initialize(out)
            while True:
                content = queue.get()
                if content is None:
                    break
                data, metadata = content
                self._write(out, data=data, metadata=metadata, mods=parsed_mods)
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
    def _write(
        self,
        out: Union[IO, Connection],
        data: Dict[str, np.ndarray],
        metadata: pd.DataFrame,
        mods: Dict[str, str],
    ):
        """
        Internal writer function.

        This function takes a batch of data and corresponding metadata and writes it to the provided
        output file handle. The file handle logic is provided using the public write or async_write
        functions.

        :param out: file handle accepting the data to be written to disk
        :param data: Dictionary containing TODO keys and corresponding values as numpy array
        :param metadata: a dataframe that contains the columns TODO
        :param mods: optional dictionary mapping libary format-specific modification patterns (keys)
            to UNIMOD IDs (values)

        """
        pass

    @abstractmethod
    def _initialize(self, out: Union[IO, Connection]):
        pass
