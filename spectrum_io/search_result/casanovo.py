from __future__ import annotations

import logging
import re

import pandas as pd
from pyteomics.mztab import MzTab
from tqdm import tqdm

from .search_results import SearchResults

logger = logging.getLogger(__name__)


class Casanovo(SearchResults):
    """Class to handle search results from Casanovo."""

    MOD_PATTERNS = {
        r"C\[Carbamidomethyl\]": "C[UNIMOD:4]",
        r"M\[Oxidation\]": "M[UNIMOD:35]",
    }

    EXCLUDE_MODS = ["N[Deamidated]", "Q[Deamidated]", "[Acetyl]", "[Carbamyl]", "[Ammonia-loss]"]

    def parse_modifications(self, sequence: str) -> str:
        """
        Parsing using regex to avoid accidental replacements.
        """

        if any(ex in sequence for ex in self.EXCLUDE_MODS):
            return float("nan")

        for pattern, replacement in self.MOD_PATTERNS.items():
            sequence = re.sub(pattern, replacement, sequence)

        return sequence

    def map_metadata(self, metadata):
        """
        Build mapping of ms_run keys to raw filenames.
        Example:
            'ms_run[1]-location': 'file:///path/to/raw1.mgf'
            → {'ms_run[1]': 'raw1'}
        """
        filename_mapping = {}
        for k, v in metadata.items():
            if "ms_run" in k and "-location" in k:
                ms_run_key = k.split("-")[0]
                filename = v.split("/")[-1].split(".")[0]
                filename_mapping[ms_run_key] = filename
        self.filename_mapping = filename_mapping
        return filename_mapping

    def read_result(
        self,
        tmt_label: str = "",
        custom_mods: dict[str, int] | None = None,
        ptm_unimod_id: int | None = 0,
        ptm_sites: list[str] | None = None,
    ):
        """
        Read Casanovo results from a file or directory of .mztab files.
        Returns a concatenated DataFrame of results.
        """
        # TODO: write custom_mods aware functions
        if self.path.is_file():
            file_list = [self.path]
        elif self.path.is_dir():
            file_list = [f for f in self.path.rglob("*") if f.suffix.lower() == ".mztab"]
        else:
            raise FileNotFoundError(f"{self.path} could not be found.")

        casanovo_results = []
        for mztab_file in tqdm(file_list, desc="Reading Casanovo results"):
            try:
                tables = MzTab(str(mztab_file))
                result_df = self.convert_to_internal(tables)
                casanovo_results.append(result_df)
            except Exception as e:
                logger.error(f"Error processing {mztab_file}: {e}")
                continue

        return pd.concat(casanovo_results, ignore_index=True)

    def convert_to_internal(self, tables: MzTab):
        """
        Convert MzTab tables into Casanovo's internal DataFrame format.
        """
        # Build mapping from metadata
        table_metadata = self.map_metadata(tables.metadata)

        # Spectrum match table from pyteomics
        df = tables.spectrum_match_table
        df.loc[df["search_engine_score[1]"] < 0.0, "search_engine_score[1]"] += 1.0
        df["opt_ms_run[1]_aa_scores"] = (
            df["opt_ms_run[1]_aa_scores"].str.split(",").apply(lambda x: "|".join(f"{float(i):.2f}" for i in x))
        )
        spectra_ref_parts = df["spectra_ref"].str.split(":", n=1, expand=True)
        df["RAW_FILE"] = spectra_ref_parts[0].map(table_metadata)
        df["SCAN_NUMBER"] = (
            spectra_ref_parts[1].str.split("scan=").str[-1] if len(spectra_ref_parts.columns) > 1 else ""
        )
        df["MODIFIED_SEQUENCE"] = df["opt_ms_run[1]_proforma"].map(self.parse_modifications)
        df.dropna(subset=["MODIFIED_SEQUENCE"], inplace=True)
        df["SEQUENCE"] = df["MODIFIED_SEQUENCE"].str.replace(r"\[[^\]]+\]", "", regex=True)
        df["PEPTIDE_LENGTH"] = df["SEQUENCE"].str.len()
        df = df.query("PEPTIDE_LENGTH > 6 and PEPTIDE_LENGTH < 31").copy()
        df["REVERSE"] = False
        df["PROTEINS"] = "UNKNOWN"

        df.rename(
            columns={
                "exp_mass_to_charge": "MASS",
                "calc_mass_to_charge": "CALC_MASS",
                "search_engine_score[1]": "SCORE",
                "charge": "PRECURSOR_CHARGE",
                "opt_ms_run[1]_aa_scores": "AA_SCORE",
            },
            inplace=True,
        )

        result_columns = [
            "RAW_FILE",
            "SCAN_NUMBER",
            "MODIFIED_SEQUENCE",
            "PRECURSOR_CHARGE",
            "MASS",
            "CALC_MASS",
            "SCORE",
            "AA_SCORE",
            "REVERSE",
            "SEQUENCE",
            "PEPTIDE_LENGTH",
            "PROTEINS",
        ]

        return df[result_columns]
