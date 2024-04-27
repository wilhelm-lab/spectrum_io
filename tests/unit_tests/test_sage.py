import unittest
from pathlib import Path

import pandas as pd

from spectrum_io.search_result import Sage


class TestSage(unittest.TestCase):
    """Test vlass to check Sage search result processing."""

    def test_read_sage(self):
        """Test function for reading sage results and transforming to Prosit format."""
        expected_sage_internal_path = Path(__file__).parent / "data" / "sage_output_internal.csv"
        internal_search_results_df = Sage(Path(__file__).parent / "data" / "sage_output.tsv").read_result()
        print(internal_search_results_df.columns)
        expected_df = pd.read_csv(expected_sage_internal_path, index_col=0)
        print(expected_df.columns)

        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
