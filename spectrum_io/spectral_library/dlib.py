import sqlite3
import zlib
from typing import IO, Dict, Union

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES

from .spectral_library import SpectralLibrary

DLIB_COL_NAMES = [
    "MassArray",
    "IntensityArray",
    "MassEncodedLength",
    "IntensityEncodedLength",
    "PrecursorCharge",
    "PeptideModSeq",
    "PeptideSeq",
    "RTInSeconds",
    "PrecursorMz",
]


class DLib(SpectralLibrary):
    """Main to init a DLib obj."""

    @property
    def standard_mods(self) -> Dict[str, int]:
        """Standard modifications that are always applied if not otherwise specified."""
        return {
            "M[+15.994915]": 35,
            "C[+57.021464]": 4,
            "^[+304.207146]": 2016,
            "K[+304.207146]": 2016,
            "^[+229.162932]": 737,
            "K[+229.162932]": 737,
        }

    def _initialize(self, out: Union[IO, sqlite3.Connection]):
        if isinstance(out, IO):
            raise TypeError("Not supported. Use msp/spectronaut if you want to write a text file.")
        if self.mode == "w":
            DLib._create_database(out)

    def _get_handle(self):
        return sqlite3.connect(self.out_path)

    def _calculate_masked_values(self, fragmentmz: np.ndarray, intensities: np.ndarray):
        """
        Internal function that masks, filters, byte encodes, swaps and compresses fragmentmz \
        and intensities.

        This will produce the data for the following columns in this order:
            - 'MassArray'
            - 'IntensityArray',
            - 'MassEncodedLength',
            - 'IntensityEncodedLength'.
        :param fragmentmz: fragmentmz provided in __init__
        :param intensities: intensities provided in __init__
        :return: 4 lists as described above
        """
        mz_bytes_list = []
        i_bytes_list = []
        mz_lengths = []
        i_lengths = []

        full_mask = self._fragment_filter_passed(fragmentmz, intensities)
        for mz, i, mask in zip(fragmentmz, intensities, np.array(full_mask)):
            # mask to only existing peaks
            sort_index = np.argsort(mz[mask])
            masked_mz_ordered = mz[mask][sort_index]
            masked_i_ordered = i[mask][sort_index]

            # convert arrays to big-endian byte blops
            masked_mz_ordered = np.array(masked_mz_ordered, dtype=np.double)
            masked_i_ordered = np.array(masked_i_ordered * 100, dtype=np.float32)
            masked_mz_ordered.byteswap(inplace=True)
            masked_i_ordered.byteswap(inplace=True)
            bytes_mz = bytes(masked_mz_ordered)
            bytes_i = bytes(masked_i_ordered)
            mz_bytes_list.append(zlib.compress(bytes_mz))
            i_bytes_list.append(zlib.compress(bytes_i))
            mz_lengths.append(len(bytes_mz))
            i_lengths.append(len(bytes_i))
        return mz_bytes_list, i_bytes_list, mz_lengths, i_lengths

    @staticmethod
    def _create_database(conn: sqlite3.Connection):
        """
        Creates the database file with prefab tables entries, peptidetoprotein (p2p) and metadata, according to the \
        dlib specification.

        :param conn: specifies the path of the created database file
        """
        sql_create_entries = """
            CREATE TABLE IF NOT EXISTS entries
            (
                PrecursorMz REAL NOT NULL,
                PrecursorCharge INTEGER NOT NULL,
                PeptideModSeq TEXT NOT NULL,
                PeptideSeq TEXT NOT NULL,
                Copies INTEGER NOT NULL DEFAULT 1,
                RTInSeconds REAL NUT NULL,
                Score REAL NOT NULL DEFAULT 0,
                MassEncodedLength INTEGER NOT NULL,
                MassArray BLOB NOT NULL,
                IntensityEncodedLength INTEGER NOT NULL,
                IntensityArray BLOB NOT NULL,
                CorrelationEncodedLength INTEGER,
                CorrelationArray BLOB,
                RTInSecondsStart REAL,
                RTInSecondsStop REAL,
                MedianChromatogramEncodedLength INTEGER,
                MedianChromatogramArray BLOB,
                SourceFile TEXT NOT NULL DEFAULT 'Oktoberfest'
            )
        """
        sql_create_p2p = """
            CREATE TABLE IF NOT EXISTS peptidetoprotein
            (
                PeptideSeq TEXT NOT NULL,
                isDecoy BOOL DEFAULT FALSE,
                ProteinAccession TEXT NOT NULL DEFAULT 'UNKNOWN'
            )
        """
        sql_create_meta = """
            CREATE TABLE IF NOT EXISTS metadata (Key string not null, Value string not null)
        """
        sql_insert_meta = "INSERT INTO metadata VALUES (?,?)"
        c = conn.cursor()
        c.execute(sql_create_entries)
        c.execute(sql_create_p2p)
        c.execute(sql_create_meta)
        c.execute(sql_insert_meta, ["version", "0.1.14"])
        c.execute(sql_insert_meta, ["staleProteinMapping", "true"])
        conn.commit()

    def _write(
        self,
        out: Union[IO, sqlite3.Connection],
        data: Dict[str, np.ndarray],
        metadata: pd.DataFrame,
        mods: Dict[str, str],
    ):
        if isinstance(out, IO):
            raise TypeError("Not supported. Use msp/spectronaut if you want to write a text file.")
        seqs = metadata["SEQUENCE"]
        modseqs = metadata["MODIFIED_SEQUENCE"].replace(mods, regex=True)
        # mass_mod_sequences = internal_to_mod_mass(modseqs)#, custom_mods)
        # print(mass_mod_sequences, "masmodseq")

        p_charges = metadata["PRECURSOR_CHARGE"]
        p_mzs = (metadata["MASS"] + (p_charges * PARTICLE_MASSES["PROTON"])) / p_charges
        # ces = metadata["COLLISION_ENERGY"]

        pr_ids = metadata["PROTEINS"]

        # prepare spectra
        irts = data["irt"][:, 0]  # should create a 1D view of the (n_peptides, 1) shaped array
        f_mzss = data["mz"]
        f_intss = data["intensities"]
        # f_annotss = data["annotation"].astype("S", copy=False)

        masked_values = self._calculate_masked_values(f_mzss, f_intss)

        data_list = [*masked_values, p_charges, modseqs, seqs, irts, p_mzs]
        entries = pd.DataFrame(dict(zip(DLIB_COL_NAMES, data_list)))
        p2p = pd.DataFrame({"PeptideSeq": seqs, "ProteinAccession": pr_ids})

        out.execute("BEGIN")

        entries.to_sql(index=False, name="entries", con=out, if_exists="append", method="multi")
        p2p.to_sql(index=False, name="peptidetoprotein", con=out, if_exists="append", method="multi")

        out.commit()
        # conn.close()

        # def prepare_spectrum(self):
        """Converts grpc output and metadata dataframe into dlib format."""
        # precursor_mz: Union[List[float], np.ndarray],
        # precursor_charges: Union[List[int], np.ndarray],
        # modified_sequences: List[str],
        # retention_times: Union[List[float], np.ndarray],
        # fragmentmz: List[np.ndarray],
        # intensities: List[np.ndarray],

        # intensities = self.grpc_output[list(self.grpc_output)[0]]["intensity"]
        # fragment_mz = self.grpc_output[list(self.grpc_output)[0]]["fragmentmz"]
        # annotation = self.grpc_output[list(self.grpc_output)[0]]["annotation"]
        # irt = self.grpc_output[list(self.grpc_output)[1]]
        # retention_times = irt.flatten()
        # modified_sequences = self.spectra_input["MODIFIED_SEQUENCE"]

        # precursor_charges = self.spectra_input["PRECURSOR_CHARGE"]
        # precursor_masses = self.spectra_input["MASS"]
        # precursor_mz = (precursor_masses + (precursor_charges * PARTICLE_MASSES["PROTON"])) / precursor_charges

        # self.create_database(self.out_path)

        # gather all values for the entries table and create pandas DataFrame
        # masked_values = self._calculate_masked_values(fragment_mz, intensities)
        # mass_mod_sequences = internal_to_mod_mass(modified_sequences)
        # sequences = internal_without_mods(modified_sequences)
        # data_list = [*masked_values, precursor_charges, mass_mod_sequences, sequences, retention_times, precursor_mz]
        # self.entries = pd.DataFrame(dict(zip(DLIB_COL_NAMES, data_list)))

        # hardcoded entries that we currently not use.
        # Visit https://bitbucket.org/searleb/encyclopedia/wiki/EncyclopeDIA%20File%20Formats for dlib specs
        # self.entries["Copies"] = 1  # this is hardcorded for now and unused
        # self.entries["Score"] = 0
        # self.entries["CorrelationEncodedLength"] = None
        # self.entries["CorrelationArray"] = None
        # self.entries["RTInSecondsStart"] = None
        # self.entries["RTInSecondsStop"] = None
        # self.entries["MedianChromatogramEncodedLength"] = None
        # self.entries["MedianChromatogramArray"] = None
        # self.entries["SourceFile"] = "Prosit"

        # gather all values for the p2p table and create pandas DataFrame
        # self.p2p = pd.DataFrame({"PeptideSeq": sequences, "isDecoy": False, "ProteinAccession": "unknown"})
