import unittest
from pathlib import Path

import numpy as np

import spectrum_io.spectral_library.digest as digest


class TestCharge(unittest.TestCase):
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
        ]

        digest.main(args)
