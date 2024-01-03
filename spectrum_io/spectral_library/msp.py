from typing import Dict

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES
from spectrum_fundamentals.mod_string import internal_to_mod_names, internal_without_mods

from .spectral_library import SpectralLibrary


class MSP(SpectralLibrary):
    """Main to initialze a MSP obj."""

    def _write(self, out: str, data: Dict[str, np.ndarray], metadata: pd.DataFrame):
        # prepare metadata
        stripped_peptides = metadata["SEQUENCE"]
        modss = internal_to_mod_names(metadata["peptide_sequences"])
        p_charges = metadata["precursor_charges"]
        p_mzs = (metadata["MASS"] + (p_charges * PARTICLE_MASSES["PROTON"])) / p_charges
        ces = metadata["collision_energies"]

        # prepare spectra
        irts = data["irt"][:, 0]  # should create a 1D view of the (n_peptides, 1) shaped array
        f_mzss = data["mz"]
        f_intss = data["intensities"]
        f_annotss = data["annotation"]

        for stripped_peptide, p_charge, p_mz, ce, mods, irt, f_mzs, f_ints, f_annots in zip(
            stripped_peptides, p_charges, p_mzs, ces, modss, irts, f_mzss, f_intss, f_annotss
        ):
            out.write(f"Name: {stripped_peptide}/{p_charge}\nMW: {p_mz}\n")
            out.write(
                f"Comment: Parent={p_mz} Collision_energy={ce} Mods={mods[0]} "
                f"ModString={mods[1]}/{p_charge} iRT={irt}\n"
            )
            fragment_list = []
            for f_mz, f_int, f_annot in zip(
                f_mzs,
                f_ints,
                f_annots,
            ):
                if f_mz != -1:
                    if f_annot.endswith(b"1"):
                        annot = f_annot[:-2]
                    else:
                        annot = f_annot.replace(b"+", b"^")
                    fragment_list.append(f'{f_mz}\t{f_int}\t"{annot.decode()}/0.0ppm"\n')

            out.write(f"Num peaks: {len(fragment_list)}\n")
            for line in fragment_list:
                out.write(line)

    def _write_header(self, out: str):
        pass
