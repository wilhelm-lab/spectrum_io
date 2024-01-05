import os
import re
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES
from spectrum_fundamentals.mod_string import internal_to_spectronaut, internal_without_mods

from .spectral_library import SpectralLibrary


def split_annotations(annotations: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Split annotation strings into components and return separate arrays.

    :param annotations: A 2D numpy array containing annotation strings.
    :return: Four numpy arrays representing fragment numbers, ion types, charges, and losses.
    :raises ValueError: If an annotation has an invalid format.
    """
    f_ns = np.zeros_like(annotations, dtype=int)
    f_types = np.empty_like(annotations, dtype="U1")
    f_charges = np.zeros_like(annotations, dtype=int)
    f_losses = np.empty_like(annotations, dtype="U10")

    with np.nditer(
        [annotations.astype("S", copy=False), f_ns, f_types, f_charges, f_losses],
        flags=["multi_index"],
        op_flags=[["readonly"], ["writeonly"], ["writeonly"], ["writeonly"], ["writeonly"]],
    ) as it:

        for annot, f_n, f_type, f_charge, f_loss in it:
            match = re.match(rb"([by])(\d+)\+(\d)(?:-(\w+))?", annot)
            if match:
                f_type[...] = match.group(1)
                f_n[...] = int(match.group(2))
                f_charge[...] = int(match.group(3))
                f_loss[...] = match.group(4) if match.group(4) else "noloss"
            else:
                raise ValueError(f"Invalid annotation format: {annot}")

    return f_ns, f_types, f_charges, f_losses


class Spectronaut(SpectralLibrary):
    """Main to initialze a Spectronaut obj."""

    # Check spectronaut folder for output format.

    def _write(self, out: str, data: Dict[str, np.ndarray], metadata: pd.DataFrame):
        # prepare metadata
        seqs = metadata["SEQUENCE"]
        modseqs = internal_to_spectronaut(metadata["peptide_sequences"].apply(lambda x: "_" + x + "_"))
        p_charges = metadata["precursor_charges"]
        p_mzs = (metadata["MASS"] + (p_charges * PARTICLE_MASSES["PROTON"])) / p_charges
        ces = metadata["collision_energies"]

        # prepare spectra
        irts = data["irt"][:, 0]  # should create a 1D view of the (n_peptides, 1) shaped array
        f_mzss = data["mz"]
        f_intss = data["intensities"]
        f_annotss = data["annotation"]

        f_nss, f_typess, f_chargess, f_lossess = split_annotations(f_annotss)

        for f_ints, f_mzs, modseq, seq, p_charge, p_mz, irt, ce, f_ns, f_types, f_charges, f_losses in zip(
            f_intss, f_mzss, modseqs, seqs, p_charges, p_mzs, irts, ces, f_nss, f_typess, f_chargess, f_lossess
        ):
            line_start = f"{modseq},{seq},{seq},{p_charge},{p_mz},{irt},{ce},"
            for f_int, f_mz, f_n, f_type, f_charge, f_loss in zip(f_ints, f_mzs, f_ns, f_types, f_charges, f_losses):
                if f_mz != -1:
                    out.write(line_start)
                    out.write(f"{f_int},{f_mz},{f_n},{f_type},{f_charge},{f_loss}\n")

    def _write_header(self, out: str):
        if self.mode == "w":
            out.write(
                "ModifiedPeptide,LabeledPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,CollisionEnergy,"
                "RelativeFragmentIntensity,FragmentMz,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType\n"
            )
