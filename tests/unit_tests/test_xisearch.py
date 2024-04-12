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

        internal_search_results_df = Xisearch(path=Path(__file__).parent / "data" / "xisearch_output.tsv").read_result()
        internal_search_results_df.reset_index(drop=True, inplace=True)

        converters = {
            "Modifications_A": str,
            "Modifications_B": str,
            "ModificationPositions1": str,
            "ModificationPositions2": str,
            "mods_p1": str,
            "mods_p2": str,
            "mod_pos_p1": str,
            "mod_pos_p2": str,
        }

        expected_df = pd.read_csv(expected_xisearch_internal_path, converters=converters)
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
