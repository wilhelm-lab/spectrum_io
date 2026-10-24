import unittest

import pandas as pd

from spectrum_io.search_result.search_results import SearchResults


class TestPTMPermutations(unittest.TestCase):
    """Test the generation of PTM localization candidates."""

    def test_permutations_keep_search_engine_localization(self):
        """All localizations of the phosphorylation are generated, the search engine localization is kept."""
        df = pd.DataFrame({"MODIFIED_SEQUENCE": ["PES[UNIMOD:21]TK", "PEPK"], "SCAN_NUMBER": [1, 2]})
        result = SearchResults.add_ptm_permutations(df, {21: ["S", "T"]})
        self.assertEqual(
            sorted(result.loc[result["SCAN_NUMBER"] == 1, "MODIFIED_SEQUENCE"]),
            ["PEST[UNIMOD:21]K", "PES[UNIMOD:21]TK"],
        )
        self.assertTrue(
            (result.loc[result["SCAN_NUMBER"] == 1, "SEARCH_MODIFIED_SEQUENCE"] == "PES[UNIMOD:21]TK").all()
        )
        self.assertEqual(result.loc[result["SCAN_NUMBER"] == 2, "MODIFIED_SEQUENCE"].tolist(), ["PEPK"])

    def test_no_permutations_without_ptms(self):
        """Without PTMs to localize, the search results are unchanged."""
        df = pd.DataFrame({"MODIFIED_SEQUENCE": ["PES[UNIMOD:21]TK"]})
        pd.testing.assert_frame_equal(SearchResults.add_ptm_permutations(df, {}), df)
