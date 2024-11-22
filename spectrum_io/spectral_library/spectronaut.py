import re
from itertools import chain, cycle
from sqlite3 import Connection
from typing import IO, Dict, Union

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES

from .spectral_library import SpectralLibrary


class Spectronaut(SpectralLibrary):
    """Main to initialze a Spectronaut obj."""

    # Check spectronaut folder for output format.

    @property
    def standard_mods(self) -> Dict[str, int]:
        """Standard modifications that are always applied if not otherwise specified."""
        return {
            "C[Carbamidomethyl (C)]": 4,
            "M[Oxidation (O)]": 35,
            "^[TMT_6]": 737,
            "K[TMT_6]": 737,
            "^[TMT_Pro]": 2016,
            "K[TMT_Pro]": 2016,
        }

    @staticmethod
    def _assemble_fragment_string(f_int: float, f_mz: float, f_annot: bytes):
        m = re.match(r"([by])(\d+)\+(\d)(?:-(\w+))?", f_annot.decode())
        if m is None:
            raise ValueError(f"Malformed annotation string encountered: {f_annot.decode()}")
        return (
            f"{f_int:.4f},{f_mz:.8f},{m.group(2)},{m.group(1)},{m.group(3)},{m.group(4) if m.group(4) else 'noloss'}\n"
        )

    def _write(
        self,
        out: Union[IO, Connection],
        data: Dict[str, np.ndarray],
        metadata: pd.DataFrame,
        mods: Dict[str, str],
    ):
        # prepare metadata
        if isinstance(out, Connection):
            raise TypeError("Not supported. Use DLib if you want to write a database file.")
        seqs = metadata["SEQUENCE"]

        modseqs = metadata["MODIFIED_SEQUENCE"].replace(mods, regex=True).apply(lambda x: "_" + x + "_")
        # modseqs = internal_to_spectronaut(metadata["MODIFIED_SEQUENCE"].apply(lambda x: "_" + x + "_"))
        p_charges = metadata["PRECURSOR_CHARGE"]
        p_mzs = (metadata["MASS"] + (p_charges * PARTICLE_MASSES["PROTON"])) / p_charges
        ces = metadata["COLLISION_ENERGY"]
        pr_ids = metadata["PROTEINS"]

        # prepare spectra
        irts = data["irt"][:, 0]  # should create a 1D view of the (n_peptides, 1) shaped array
        f_mzss = data["mz"]
        f_intss = data["intensities"]
        f_annotss = data["annotation"].astype("S", copy=False)

        vec_assemble = np.vectorize(Spectronaut._assemble_fragment_string)

        for f_ints, f_mzs, modseq, seq, p_charge, p_mz, irt, ce, pr_id, f_annots in zip(
            f_intss, f_mzss, modseqs, seqs, p_charges, p_mzs, irts, ces, pr_ids, f_annotss
        ):
            cond = self._fragment_filter_passed(f_mzs, f_ints)
            line_start = [f"{modseq},{seq},{seq},{p_charge},{p_mz:.8f},{irt:.2f},{ce},{pr_id},"]
            fragment_list = vec_assemble(f_ints[cond], f_mzs[cond], f_annots[cond])
            out.writelines(chain.from_iterable(zip(cycle(line_start), fragment_list)))

    def _initialize(self, out: Union[IO, Connection]):
        if isinstance(out, Connection):
            raise TypeError("Not supported. Use DLib if you want to write a database file.")
        if self.mode == "w":
            out.write(
                "ModifiedPeptide,LabeledPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,CollisionEnergy,ProteinIds,"
                "RelativeFragmentIntensity,FragmentMz,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType\n"
            )
