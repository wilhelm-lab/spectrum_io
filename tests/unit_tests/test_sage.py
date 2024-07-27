import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result import Sage


class TestSage(unittest.TestCase):
    """Test vlass to check Sage search result processing."""

    def test_read_sage(self):
        """Test function for reading sage results and transforming to Prosit format."""
        expected_sage_internal_path = Path(__file__).parent / "data" / "sage_output_internal.csv"
        internal_search_results_df = (
            Sage(Path(__file__).parent / "data" / "sage_output.tsv").read_result().reset_index(drop=True)
        )
        expected_df = pd.read_csv(expected_sage_internal_path)
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)

    def test_read_sage_custom(self):
        """Test function for reading sage results with custom mods and transforming to Prosit format ."""
        stat_mods = {"15.9948": "[UNIMOD:35]"}
        expected_sage_internal_path = Path.cwd() / "data" / "sage_output_internal_mods.csv"
        internal_search_results_df = (
            Sage(Path.cwd() / "data" / "sage_output_mods.tsv").read_result(stat_mods=stat_mods).reset_index(drop=True)
        )
        expected_df = pd.read_csv(expected_sage_internal_path)
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
