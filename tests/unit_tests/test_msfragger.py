import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result.msfragger import MSFragger

COLUMNS = [
    "RAW_FILE",
    "SCAN_NUMBER",
    "MODIFIED_SEQUENCE",
    "PRECURSOR_CHARGE",
    "SCAN_EVENT_NUMBER",
    "MASS",
    "SCORE",
    "REVERSE",
    "SEQUENCE",
    "PEPTIDE_LENGTH",
    "PROTEINS",
]


class TestMSFragger(unittest.TestCase):
    """Class to test MSFragger."""

    def test_check_decoys(self):
        """Test only proteins with the decoy prefix rev are decoys, not proteins with rev elsewhere in the name."""
        self.assertTrue(MSFragger.check_decoys("rev_sp|P12345|A_HUMAN"))
        self.assertTrue(MSFragger.check_decoys("rev_sp|P12345|A_HUMAN;rev_sp|P67890|B_HUMAN"))
        self.assertFalse(MSFragger.check_decoys("rev_sp|P12345|A_HUMAN;sp|P67890|B_HUMAN"))
        self.assertFalse(MSFragger.check_decoys("sp|Q9BREV|PREV_HUMAN"))
        self.assertFalse(MSFragger.check_decoys("sp|P12345|SREV1_HUMAN Protein sorev"))

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
        self.assertTrue("PROTEINS" in df.columns)

    def test_read_msfragger(self):
        """Test function for reading msfragger results and transforming to Prosit format."""
        expected_msfragger_internal_path = Path(__file__).parent / "data" / "psm_tmt_internal.csv"

        internal_search_results_df = MSFragger(Path(__file__).parent / "data" / "psm_tmt.pepXML").read_result(
            tmt_label="tmtpro"
        )
        expected_df = pd.read_csv(expected_msfragger_internal_path, index_col=0)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

    def test_read_msfragger_with_custom_mods(self):
        """Test function for reading msfragger results and transforming to Prosit format with custom mods."""
        expected_msfragger_internal_path = Path(__file__).parent / "data" / "psm_tmt_internal.csv"
        custom_mods = {
            "M[147]": 35,
            "C": 4,
        }

        internal_search_results_df = MSFragger(Path(__file__).parent / "data" / "psm_tmt.pepXML").read_result(
            tmt_label="tmtpro", custom_mods=custom_mods
        )
        expected_df = pd.read_csv(expected_msfragger_internal_path, index_col=0)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

    def test_read_msfragger_with_custom_mods_with_tmt(self):
        """Test function for reading msfragger results and transforming to Prosit format with custom mods and explicit TMT."""
        expected_msfragger_internal_path = Path(__file__).parent / "data" / "psm_tmt_internal.csv"
        custom_mods = {"^n[305]": 2016, "K[432]": 2016}

        internal_search_results_df = MSFragger(Path(__file__).parent / "data" / "psm_tmt.pepXML").read_result(
            custom_mods=custom_mods
        )
        expected_df = pd.read_csv(expected_msfragger_internal_path, index_col=0)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])
