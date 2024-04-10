import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from spectrum_io.search_result import Xisearch


class TestXisearch(unittest.TestCase):
    """Test vlass to check Xisearch search result processing."""

    def test_read_xisearch(self):
        """Test function for reading Xisearch results and transforming to Prosit format."""
        expected_xisearch_internal_path = Path(__file__).parent / "data" / "xisearch_output_internal.tsv"

        internal_search_results_df = Xisearch(Path(__file__).parent / "data" / "xisearch_output.tsv").read_result()
        internal_search_results_df.reset_index(drop=True, inplace=True)
        expected_df = pd.read_csv(expected_xisearch_internal_path)
        expected_df["Modifications_A"] = expected_df["Modifications_A"].fillna("nan")
        expected_df["Modifications_B"] = expected_df["Modifications_B"].fillna("nan")
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
