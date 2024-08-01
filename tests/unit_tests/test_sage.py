import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result import Sage

COLUMNS = [
    "RAW_FILE",
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


class TestSage(unittest.TestCase):
    """Test class to check Sage search result processing."""

    def test_read_sage(self):
        """Test function for reading sage results and transforming to Prosit format."""
        expected_sage_internal_path = Path(__file__).parent / "data" / "sage_output_internal.csv"
        internal_search_results_df = Sage(Path(__file__).parent / "data" / "sage_output.tsv").read_result(
            tmt_label="tmt"
        )
        expected_df = pd.read_csv(expected_sage_internal_path)
        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

    def test_read_sage_with_custom_mods(self):
        """Test function for reading sage results with custom mods and transforming to Prosit format ."""
        custom_mods = {
            "M[+15.9948]": 35,
        }
        expected_sage_internal_path = Path(__file__).parent / "data" / "sage_output_internal_mods.csv"
        internal_search_results_df = Sage(Path(__file__).parent / "data" / "sage_output_mods.tsv").read_result(
            custom_mods=custom_mods, tmt_label="tmt"
        )
        expected_df = pd.read_csv(expected_sage_internal_path)
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)

    def test_read_sage_with_custom_mods_with_tmt(self):
        """Test function for reading sage results with custom mods and transforming to Prosit format ."""
        custom_mods = {"M[+15.9948]": 35, "K[+229.1629]": 737, "^[+229.1629]-": 737}
        expected_sage_internal_path = Path(__file__).parent / "data" / "sage_output_internal_mods.csv"
        internal_search_results_df = Sage(Path(__file__).parent / "data" / "sage_output_mods.tsv").read_result(
            custom_mods=custom_mods
        )
        expected_df = pd.read_csv(expected_sage_internal_path)
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
