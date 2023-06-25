import io

import numpy as np
import pandas as pd
import pytest

import spectrum_io.search_result.maxquant as mq
from spectrum_io.search_result.search_results import filter_valid_prosit_sequences


class TestAddTMTMod:
    """Class to test tmt modification addition."""

    def test_add_tmt_mod(self):
        """Test addition of tmt modification."""
        assert (
            mq.MaxQuant.add_tmt_mod(1.0, "[UNIMOD:2016]ABC[UNIMOD:4]K[UNIMOD:2016]", "[UNIMOD:2016]")
            == 1.0 + 2 * 304.207146
        )


class TestUpdateColumns:
    """Class to test update columns."""

    def test_update_columns(self, maxquant_df: pd.DataFrame):
        """
        Test column update.

        :param maxquant_df: maxquant df as pd.DataFrame
        """
        prosit_df = mq.MaxQuant.update_columns_for_prosit(maxquant_df, tmt_labeled="")
        assert not prosit_df["REVERSE"][0]
        assert prosit_df["REVERSE"][3]

        assert prosit_df["MODIFIED_SEQUENCE"][0] == "DS[UNIMOD:21]DS[UNIMOD:21]WDADAFSVEDPVRK"
        assert prosit_df["MODIFIED_SEQUENCE"][3] == "SS[UNIMOD:21]PTPES[UNIMOD:21]PTMLTK"

        assert prosit_df["SEQUENCE"][0] == "DSDSWDADAFSVEDPVRK"
        assert prosit_df["SEQUENCE"][3] == "SSPTPESPTMLTK"

        assert prosit_df["PEPTIDE_LENGTH"][0] == 18
        assert prosit_df["PEPTIDE_LENGTH"][3] == 13

    def test_update_columns_silac(self, maxquant_df: pd.DataFrame):
        """
        Test column update silac.

        :param maxquant_df: maxquant df as pd.DataFrame
        """
        maxquant_df["LABELING_STATE"] = [1, 1, 1, 2, 2]
        prosit_df = mq.MaxQuant.update_columns_for_prosit(maxquant_df, tmt_labeled="")
        assert prosit_df["MODIFIED_SEQUENCE"][0] == "DS[UNIMOD:21]DS[UNIMOD:21]WDADAFSVEDPVR[UNIMOD:267]K[UNIMOD:259]"
        assert prosit_df["MODIFIED_SEQUENCE"][3] == "SS[UNIMOD:21]PTPES[UNIMOD:21]PTMLTK"

        assert prosit_df["MASS"][0] == 1.0 + 8.014199 + 10.008269
        assert prosit_df["MASS"][3] == 2.0

    def test_update_columns_tmt(self, maxquant_df: pd.DataFrame):
        """
        Test column update tmt.

        :param maxquant_df: maxquant df as pd.DataFrame
        """
        prosit_df = mq.MaxQuant.update_columns_for_prosit(maxquant_df, tmt_labeled="tmt")
        assert prosit_df["MODIFIED_SEQUENCE"][0] == "[UNIMOD:737]DS[UNIMOD:21]DS[UNIMOD:21]WDADAFSVEDPVRK[UNIMOD:737]"
        assert prosit_df["MODIFIED_SEQUENCE"][3] == "[UNIMOD:737]SS[UNIMOD:21]PTPES[UNIMOD:21]PTMLTK[UNIMOD:737]"

        assert prosit_df["MASS"][0] == 1.0 + 2 * 229.162932
        assert prosit_df["MASS"][3] == 2.0 + 2 * 229.162932

    def test_update_columns_tmt_msa(self, maxquant_df: pd.DataFrame):
        """
        Test column update tmt msa.

        :param maxquant_df: maxquant df as pd.DataFrame
        """
        prosit_df = mq.MaxQuant.update_columns_for_prosit(maxquant_df, tmt_labeled="tmt_msa")
        assert (
            prosit_df["MODIFIED_SEQUENCE_MSA"][0] == "[UNIMOD:737]DS[UNIMOD:23]DS[UNIMOD:23]WDADAFSVEDPVRK[UNIMOD:737]"
        )
        assert prosit_df["MODIFIED_SEQUENCE_MSA"][3] == "[UNIMOD:737]SS[UNIMOD:23]PTPES[UNIMOD:23]PTMLTK[UNIMOD:737]"

    def test_filter_valid_prosit_sequences(self, invalid_df: pd.DataFrame):
        """Test filter_valid_prosit_sequences."""
        filtered_df = filter_valid_prosit_sequences(invalid_df)
        assert filtered_df["MODIFIED_SEQUENCE"][0] == "ABCDEFG"
        assert len(filtered_df) == 1
        assert "(ac)" not in filtered_df["MODIFIED_SEQUENCE"]
        assert "(Acetyl (Protein N-term))" not in filtered_df["MODIFIED_SEQUENCE"]
        assert "U" not in filtered_df["SEQUENCE"]
        assert filtered_df["PEPTIDE_LENGTH"].min() >= 7
        assert filtered_df["PRECURSOR_CHARGE"].max() <= 6


@pytest.fixture
def maxquant_df():
    """Create dataframes from strings: https://towardsdatascience.com/67b0c2b71e6a."""
    df_string = """              MODIFIED_SEQUENCE; REVERSE; MASS;
_DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK_;        ;  1.0;
_DS(Phospho (STY))DS(Phospho (STY))WDADAFSVEDPVRK_;        ;  1.0;
_DS(Phospho (STY))DSWDADAFS(Phospho (STY))VEDPVRK_;        ;  1.0;
     _SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK_;       +;  2.0;
     _SS(Phospho (STY))PTPES(Phospho (STY))PTMLTK_;       +;  2.0;"""
    df = pd.read_csv(io.StringIO(df_string), delimiter=";", skipinitialspace=True)
    df["Charge"] = 2
    return df


@pytest.fixture
def invalid_df():
    """Create invalid df."""
    df = pd.DataFrame(
        {
            "PEPTIDE_LENGTH": [7, 7, 6, 32],
            "MODIFIED_SEQUENCE": [
                "ABCDEFG",
                "GHD(ac)IJKL",
                "MN(Acetyl (Protein N-term))OPQR",
                "STUVWDEFSTUVWDEFSTUVWDEFSTUVWDEF",
            ],
            "SEQUENCE": ["ABCDEFG", "GHDIJKL", "MNOPQR", "STUVWDEFSTUVWDEFSTUVWDEFSTUVWDEF"],
            "PRECURSOR_CHARGE": [2, 5, 7, 6],
        }
    )
    return df
