import sqlite3
import zlib
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from spectrum_fundamentals.constants import PARTICLE_MASSES
from spectrum_fundamentals.mod_string import internal_to_mod_mass, internal_without_mods

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

    def _calculate_masked_values(self, fragmentmz: List[np.ndarray], intensities: List[np.ndarray]):
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
        for mz, i in zip(fragmentmz, intensities):
            # mask to only existing peaks
            mask = i >= self.min_intensity_threshold
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
            i_bytes_list.append(zlib.compress(bytes(bytes_i)))
            mz_lengths.append(len(bytes_mz))
            i_lengths.append(len(bytes_i))
        return mz_bytes_list, i_bytes_list, mz_lengths, i_lengths

    @staticmethod
    def create_database(path: Union[str, Path]):
        """
        Creates the database file with prefab tables entries, peptidetoprotein (p2p) and metadata, according to the \
        dlib specification.

        :param path: specifies the path of the created database file
        """
        sql_create_entries = """
            CREATE TABLE entries
            (   PrecursorMz double not null,
                PrecursorCharge int not null,
                PeptideModSeq string not null,
                PeptideSeq string not null,
                Copies int not null,
                RTInSeconds double not null,
                Score double not null,
                MassEncodedLength int not null,
                MassArray blob not null,
                IntensityEncodedLength int not null,
                IntensityArray blob not null,
                CorrelationEncodedLength int,
                CorrelationArray blob,
                RTInSecondsStart double,
                RTInSecondsStop double,
                MedianChromatogramEncodedLength int,
                MedianChromatogramArray blob,
                SourceFile string not null
            )
        """
        sql_create_p2p = """
            CREATE TABLE peptidetoprotein
            (PeptideSeq string not null, isDecoy boolean, ProteinAccession string not null)
        """
        sql_create_meta = """
            CREATE TABLE metadata (Key string not null, Value string not null)
        """
        sql_insert_meta = "INSERT INTO metadata VALUES (?,?)"
        conn = sqlite3.connect(path)
        c = conn.cursor()
        c.execute(sql_create_entries)
        c.execute(sql_create_p2p)
        c.execute(sql_create_meta)
        c.execute(sql_insert_meta, ["version", "0.1.14"])
        c.execute(sql_insert_meta, ["staleProteinMapping", "true"])
        conn.commit()

    def write(self):
        """Writes the entries ad p2p table to file."""
        self._write_entries(index=False, if_exists="append", method="multi", chunksize=self.chunksize)
        self._write_p2p(index=False, if_exists="append", method="multi", chunksize=self.chunksize)

    def _write_entries(self, *args, **kwargs):
        """
        Internal function to write the entries table.

        :param args: forwarded to pandas.to_sql
        :param kwargs: forwarded to pandas.to_sql
        """
        conn = sqlite3.connect(self.out_path)
        self.entries.to_sql(name="entries", con=conn, *args, **kwargs)
        conn.commit()

    def _write_p2p(self, *args, **kwargs):
        """
        Internal function to write the p2p table.

        :param args: forwarded to pandas.to_sql
        :param kwargs: forwarded to pandas.to_sql
        """
        conn = sqlite3.connect(self.out_path)
        self.p2p.to_sql(name="peptidetoprotein", con=conn, *args, **kwargs)
        conn.commit()

    def prepare_spectrum(self):
        """Converts grpc output and metadata dataframe into dlib format."""
        # precursor_mz: Union[List[float], np.ndarray],
        # precursor_charges: Union[List[int], np.ndarray],
        # modified_sequences: List[str],
        # retention_times: Union[List[float], np.ndarray],
        # fragmentmz: List[np.ndarray],
        # intensities: List[np.ndarray],

        intensities = self.grpc_output[list(self.grpc_output)[0]]["intensity"]
        fragment_mz = self.grpc_output[list(self.grpc_output)[0]]["fragmentmz"]
        # annotation = self.grpc_output[list(self.grpc_output)[0]]["annotation"]
        irt = self.grpc_output[list(self.grpc_output)[1]]
        retention_times = irt.flatten()
        modified_sequences = self.spectra_input["MODIFIED_SEQUENCE"]

        precursor_charges = self.spectra_input["PRECURSOR_CHARGE"]
        precursor_masses = self.spectra_input["MASS"]
        precursor_mz = (precursor_masses + (precursor_charges * PARTICLE_MASSES["PROTON"])) / precursor_charges

        self.create_database(self.out_path)

        # gather all values for the entries table and create pandas DataFrame
        masked_values = self._calculate_masked_values(fragment_mz, intensities)
        mass_mod_sequences = internal_to_mod_mass(modified_sequences)
        sequences = internal_without_mods(modified_sequences)
        data_list = [*masked_values, precursor_charges, mass_mod_sequences, sequences, retention_times, precursor_mz]
        self.entries = pd.DataFrame(dict(zip(DLIB_COL_NAMES, data_list)))

        # hardcoded entries that we currently not use.
        # Visit https://bitbucket.org/searleb/encyclopedia/wiki/EncyclopeDIA%20File%20Formats for dlib specs
        self.entries["Copies"] = 1  # this is hardcorded for now and unused
        self.entries["Score"] = 0
        self.entries["CorrelationEncodedLength"] = None
        self.entries["CorrelationArray"] = None
        self.entries["RTInSecondsStart"] = None
        self.entries["RTInSecondsStop"] = None
        self.entries["MedianChromatogramEncodedLength"] = None
        self.entries["MedianChromatogramArray"] = None
        self.entries["SourceFile"] = "Prosit"

        # gather all values for the p2p table and create pandas DataFrame
        self.p2p = pd.DataFrame({"PeptideSeq": sequences, "isDecoy": False, "ProteinAccession": "unknown"})
