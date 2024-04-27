import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result.msfragger import MSFragger


class TestMSFragger(unittest.TestCase):
    """Class to test MSFragger."""

    def test_read_result(self):
        """Test read_result for MSFragger."""
        msfragger = MSFragger(Path(__file__).parent / "data/psm.pepXML")
        df = msfragger.read_result("")
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

    def test_read_msfragger(self):
        """Test function for reading sage results and transforming to Prosit format."""
        expected_msfragger_internal_path = Path(__file__).parent / "data" / "psm_tmt_internal.csv"

        internal_search_results_df = MSFragger(Path(__file__).parent / "data" / "psm_tmt.pepXML").read_result(
            tmt_labeled="tmtpro"
        )
        expected_df = pd.read_csv(expected_msfragger_internal_path, index_col=0)

        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
