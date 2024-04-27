import unittest
from pathlib import Path

import numpy as np

import spectrum_io.spectral_library.digest as digest


class TestDigest(unittest.TestCase):
    """Class to test digest."""

    def test_main(self):
        """Test digest."""
        prosit_input = Path(__file__).parent / "data/prosit_input.csv"
        prosit_input_with_proteins = Path(__file__).parent / "data/prosit_input_with_proteins.csv"

        fasta = Path(__file__).parent / "data/fasta.fasta"
        pep_prot_map = Path(__file__).parent / "data/peptide_protein_map.tsv"
        pep_prot_params = Path(__file__).parent / "data/peptide_protein_map.tsv.params.txt"
        args = [
            "--enzyme",
            "trypsinp",
            "--cleavages",
            "0",
            "--fragmentation",
            "CID",
            "--prosit_input",
            str(prosit_input),
            "--fasta",
            str(fasta),
            "--peptide_protein_map",
            str(pep_prot_map),
        ]

        digest.main(args)

        prosit_input.unlink()
        pep_prot_map.unlink()
        pep_prot_params.unlink()
        prosit_input_with_proteins.unlink()
