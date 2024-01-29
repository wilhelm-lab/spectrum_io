import os
import re
from itertools import chain, cycle
from typing import IO, Dict, Tuple

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES
from spectrum_fundamentals.mod_string import internal_to_spectronaut, internal_without_mods

from .spectral_library import SpectralLibrary


class Spectronaut(SpectralLibrary):
    """Main to initialze a Spectronaut obj."""

    # Check spectronaut folder for output format.

    @staticmethod
    def _assemble_fragment_string(f_int: float, f_mz: float, f_annot: bytes):
        m = re.match(r"([by])(\d+)\+(\d)(?:-(\w+))?", f_annot.decode())
        if m is None:
            raise ValueError(f"Malformed annotation string encountered: {f_annot.decode()}")
        return (
            f"{f_int:.4f},{f_mz:.8f},{m.group(2)},{m.group(1)},{m.group(3)},{m.group(4) if m.group(4) else 'noloss'}\n"
        )

    def _write(self, out: IO, data: Dict[str, np.ndarray], metadata: pd.DataFrame):
        # prepare metadata
        seqs = metadata["SEQUENCE"]
        modseqs = internal_to_spectronaut(metadata["MODIFIED_SEQUENCE"].apply(lambda x: "_" + x + "_"))
        p_charges = metadata["PRECURSOR_CHARGE"]
        p_mzs = (metadata["MASS"] + (p_charges * PARTICLE_MASSES["PROTON"])) / p_charges
        ces = metadata["COLLISION_ENERGY"]

        # prepare spectra
        irts = data["irt"][:, 0]  # should create a 1D view of the (n_peptides, 1) shaped array
        f_mzss = data["mz"]
        f_intss = data["intensities"]
        f_annotss = data["annotation"].astype("S", copy=False)

        vec_assemble = np.vectorize(Spectronaut._assemble_fragment_string)

        for f_ints, f_mzs, modseq, seq, p_charge, p_mz, irt, ce, f_annots in zip(
            f_intss, f_mzss, modseqs, seqs, p_charges, p_mzs, irts, ces, f_annotss
        ):
            cond = self._fragment_filter_passed(f_mzs, f_ints)
            line_start = [f"{modseq},{seq},{seq},{p_charge},{p_mz:.8f},{irt:.2f},{ce},"]
            fragment_list = vec_assemble(f_ints[cond], f_mzs[cond], f_annots[cond])
            out.writelines(chain.from_iterable(zip(cycle(line_start), fragment_list)))

    def _write_header(self, out: IO):
        if self.mode == "w":
            out.write(
                "ModifiedPeptide,LabeledPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,CollisionEnergy,"
                "RelativeFragmentIntensity,FragmentMz,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType\n"
            )
