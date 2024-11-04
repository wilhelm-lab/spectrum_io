import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from spectrum_io.search_result import Scout


class TestXisearch(unittest.TestCase):
    """Test vlass to check Xisearch search result processing."""

    def test_read_scout(self):
        """Test function for reading Scout results and transforming to Prosit format."""
        expected_scout_internal_path = Path(__file__).parent / "data" / "scout_output_internal.csv"
        internal_search_results_df = (
            Scout(path=Path(__file__).parent / "data" / "scout_output.csv").read_result().reset_index(drop=True)
        )
        internal_search_results_df["linear"] = internal_search_results_df.linear.astype(bool)
        expected_df = pd.read_csv(
            expected_scout_internal_path,
            converters={
                "linear": bool,
                "Modifications_A": str,
                "Modifications_B": str,
                "mods_p1": str,
                "mods_p2": str,
                "ModificationPositions1": str,
                "ModificationPositions2": str,
            },
        )
        pd.testing.assert_frame_equal(internal_search_results_df, expected_df)
