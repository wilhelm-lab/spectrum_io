import pickle
import unittest
from pathlib import Path

import pandas as pd

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


class TestCollisionEnergyFromFilterString(unittest.TestCase):
    """Test the collision energy fallback for mzML files without a collision energy cvParam."""

    def test_collision_energy_from_filter_string(self):
        """Test the collision energy is read from Thermo filter strings."""
        from spectrum_io.raw.msraw import _collision_energy_from_filter_string

        self.assertEqual(
            _collision_energy_from_filter_string("FTMS + p NSI d Full ms2 419.3156@hcd32.00 [110.0000-849.0000]"), 32.0
        )
        self.assertEqual(_collision_energy_from_filter_string("ITMS + c NSI d Full ms2 500.20@cid35.00 [50-1000]"), 35.0)
        self.assertEqual(_collision_energy_from_filter_string("FTMS + p NSI Full ms [350.0000-1400.0000]"), 0.0)
        self.assertEqual(_collision_energy_from_filter_string(None), 0.0)
