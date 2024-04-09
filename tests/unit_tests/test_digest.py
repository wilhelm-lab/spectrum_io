import unittest
from pathlib import Path

import numpy as np

import spectrum_io.spectral_library.digest as digest


class TestDigest(unittest.TestCase):
    """Class to test digest."""

    def test_main(self):
        """Test digest."""
        args = [
            "--enzyme",
            "trypsinp",
            "--cleavages",
            "0",
            "--fragmentation",
            "CID",
            "--prosit_input",
            str(Path(__file__).parent / "data/prosit_input.csv"),
            "--fasta",
            str(Path(__file__).parent / "data/fasta.fasta"),
            "--peptide_protein_map",
            str(Path(__file__).parent / "data/peptide_protein_map.tsv"),
        ]

        digest.main(args)
