import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import MZML_DATA_COLUMNS

import spectrum_io.raw.msraw as msraw


class TestMsraw(unittest.TestCase):
    """Class to test msraw."""

    def test_read_mzml(self):
        """Test read_mzml."""
        source = Path(__file__).parent / "data/test.mzml"

        msraw_obj = msraw.MSRaw(source, source)
        df = msraw_obj.read_mzml(source)

        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(df.shape[1], 8)
        self.assertListEqual(list(df.columns), MZML_DATA_COLUMNS)
