import pickle
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import MZML_DATA_COLUMNS

import spectrum_io.raw.msraw as msraw


def _test_read_mzml(package: str):
    source = Path(__file__).parent / "data/test.mzml"

    msraw_obj = msraw.MSRaw(source, source)
    df = msraw_obj.read_mzml(source, package=package)
    # pickle.dump(df, open(Path(__file__).parent / "data/testdf.pkl", "wb"))

    target_df = pickle.load(open(Path(__file__).parent / "data/testdf.pkl", "rb"))
    pd.testing.assert_frame_equal(df, target_df)


class TestMsraw(unittest.TestCase):
    """Class to test msraw."""

    def test_read_mzml_with_pyteomics(self):
        """Test read_mzml."""
        _test_read_mzml(package="pyteomics")

    def test_read_mzml_with_pymzml(self):
        """Test read_mzml."""
        _test_read_mzml(package="pymzml")
