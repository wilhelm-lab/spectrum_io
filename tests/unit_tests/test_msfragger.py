import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result.msfragger import MSFragger


class TestMSFragger(unittest.TestCase):
    """Class to test MSFragger."""

    def test_read_result(self):
        """Test read_result for MSFragger."""
        msfragger = MSFragger(str(Path(__file__).parent / "data/"))
        df = msfragger.read_result(str(Path(__file__).parent / "data/psm.tsv"), "")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertTrue("SEQUENCE" in df.columns)
        self.assertTrue("MODIFIED_SEQUENCE" in df.columns)
        self.assertTrue("MASS" in df.columns)
        self.assertTrue("SCORE" in df.columns)
        self.assertTrue("PRECURSOR_CHARGE" in df.columns)
        self.assertTrue("RAW_FILE" in df.columns)
        self.assertTrue("SCAN_NUMBER" in df.columns)
        self.assertTrue("REVERSE" in df.columns)
        self.assertTrue("PEPTIDE_LENGTH" in df.columns)
        self.assertTrue("SPECTRUM" in df.columns)
        self.assertTrue("RETENTION_TIME" in df.columns)
