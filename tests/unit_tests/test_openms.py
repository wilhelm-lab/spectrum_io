import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result.openms import OpenMS

COLUMNS = [
    "SCAN_NUMBER",
    "MODIFIED_SEQUENCE",
    "PRECURSOR_CHARGE",
    "MASS",
    "SCORE",
    "REVERSE",
    "SEQUENCE",
    "PEPTIDE_LENGTH",
    "PROTEINS",
]


class TestOpenMS(unittest.TestCase):
    """Class to test OpenMS."""

    def test_read_result(self):
        """Test read_result for MSFragger."""
        openms = OpenMS(Path(__file__).parent / "data/openms.idXML")
        df = openms.read_result("")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertTrue("RAW_FILE" in df.columns)
        self.assertTrue("SCAN_NUMBER" in df.columns)
        self.assertTrue("PRECURSOR_CHARGE" in df.columns)
        self.assertTrue("SCAN_EVENT_NUMBER" in df.columns)
        self.assertTrue("MODIFIED_SEQUENCE" in df.columns)
        self.assertTrue("MASS" in df.columns)
        self.assertTrue("SCORE" in df.columns)
        self.assertTrue("REVERSE" in df.columns)
        self.assertTrue("SEQUENCE" in df.columns)
        self.assertTrue("PEPTIDE_LENGTH" in df.columns)
        self.assertTrue("PROTEINS" in df.columns)

    def test_read_openms(self):
        #Test function for reading OpenMS results and transforming to Prosit format.
        expected_openms_internal_path = Path(__file__).parent / "data" / "openms.csv"
        OpenMS(Path(__file__).parent / "data" / "openms.idXML").generate_internal(out_path=Path(__file__).parent / "data" / "openms.csv")

        internal_search_results_df = OpenMS(Path(__file__).parent / "data" / "openms.idXML").read_result().reset_index(drop=True)
        expected_df = pd.read_csv(expected_openms_internal_path)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

