import io
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from spectrum_io.search_result.maxquant import MaxQuant

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


class TestAddTMTMod:
    """Class to test tmt modification addition."""

    def test_add_tmt_mod(self):
        """Test addition of tmt modification."""
        assert (
            MaxQuant.add_tmt_mod(1.0, "[UNIMOD:2016]ABC[UNIMOD:4]K[UNIMOD:2016]", "[UNIMOD:2016]")
            == 1.0 + 2 * 304.207146
        )


class TestMaxQuant(unittest.TestCase):
    """Class to test MSFragger."""

    def test_read_maxquant(self):
        """Test function for reading msfragger results and transforming to Prosit format."""
        expected_df_path = Path(__file__).parent / "data" / "msms_internal.csv"

        internal_search_results_df = MaxQuant(Path(__file__).parent / "data" / "msms.txt").read_result()
        expected_df = pd.read_csv(expected_df_path)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

    def test_read_maxquant_with_custom_mods(self):
        """Test function for reading msfragger results and transforming to Prosit format with custom mods."""
        expected_df_path = Path(__file__).parent / "data" / "msms_internal.csv"
        custom_mods = {"M(Oxidation (M)": 35, "C": 4, "S(Phospho (STY))": 21}
        internal_search_results_df = MaxQuant(Path(__file__).parent / "data" / "msms.txt").read_result(
            custom_mods=custom_mods
        )
        expected_df = pd.read_csv(expected_df_path)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])

    def test_read_maxquant_with_custom_mods_with_tmt(self):
        """Test function for reading msfragger results and transforming to Prosit format with custom mods."""
        expected_df_path = Path(__file__).parent / "data" / "msms_internal_tmt.csv"
        internal_search_results_df = MaxQuant(Path(__file__).parent / "data" / "msms.txt").read_result(tmt_label="")
        custom_mods = {
            "K": 2016,
            "^": 2016,
        }

        internal_search_results_df = MaxQuant(Path(__file__).parent / "data" / "msms.txt").read_result(
            custom_mods=custom_mods
        )
        expected_df = pd.read_csv(expected_df_path)

        pd.testing.assert_frame_equal(internal_search_results_df[COLUMNS], expected_df[COLUMNS])
