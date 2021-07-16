from .spectral_library import SpectralLibrary
import sqlite3
import pandas as pd
from typing import Union, List, Optional
import numpy as np
import zlib
from fundamentals.mod_string import internal_to_mod_mass, internal_without_mods

DLIB_COL_NAMES = [
    'MassArray',
    'IntensityArray',
    'MassEncodedLength',
    'IntensityEncodedLength',
    'PrecursorCharge',
    'PeptideModSeq',
    'PeptideSeq',
    'RTInSeconds',
    'PrecursorMz'
]


class DLib(SpectralLibrary):

    def __init__(
        self,
        precursor_mz: Union[List[float], np.ndarray],
        precursor_charges: Union[List[int], np.ndarray],
        modified_sequences: Union[List[str], np.ndarray],
        retention_times: Union[List[float], np.ndarray],
        fragmentmz: List[np.ndarray],
        intensities: List[np.ndarray],
        path: str,
        min_intensity_threshold: Optional[float] = 0.05,
    ):
        self.path = path
        self.create_database(self.path)

        # gather all values for the entries table and create pandas DataFrame
        masked_values = self._calculate_masked_values(fragmentmz, intensities, min_intensity_threshold)
        mass_mod_sequences = internal_to_mod_mass(modified_sequences)
        sequences = internal_without_mods(modified_sequences)
        data_list = [*masked_values, precursor_charges, mass_mod_sequences, sequences, retention_times, precursor_mz]
        self.entries = pd.DataFrame(dict(zip(DLIB_COL_NAMES, data_list)))

        # hardcoded entries that we currently not use.
        # Visit https://bitbucket.org/searleb/encyclopedia/wiki/EncyclopeDIA%20File%20Formats for dlib specs
        self.entries['Copies'] = 1 # this is hardcorded for now and unused
        self.entries['Score'] = 0
        self.entries['CorrelationEncodedLength'] = None
        self.entries['CorrelationArray'] = None
        self.entries['RTInSecondsStart'] = None
        self.entries['RTInSecondsStop'] = None
        self.entries['MedianChromatogramEncodedLength'] = None
        self.entries['MedianChromatogramArray'] = None
        self.entries['SourceFile'] = 'Prosit'

        # gather all values for the p2p table and create pandas DataFrame
        self.p2p = pd.DataFrame({
            'PeptideSeq': sequences,
            'isDecoy': False,
            'ProteinAccession': 'unknown'
        })

    @staticmethod
    def _calculate_masked_values(
        fragmentmz: Union[List[np.ndarray]],
        intensities: Union[List[np.ndarray]],
        intensity_min_threshold: Optional[float] = 0.05
    ):
        mz_bytes_list = []
        i_bytes_list = []
        mz_lengths = []
        i_lengths = []
        for mz, i in zip(fragmentmz, intensities):
            # mask to only existing peaks
            mask = i >= intensity_min_threshold
            print(mask)
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
    def create_database(path: str):
        SQL_CREATE_ENTRIES = """
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
        SQL_CREATE_P2P = """
            CREATE TABLE peptidetoprotein
            (PeptideSeq string not null, isDecoy boolean, ProteinAccession string not null)
        """
        SQL_CREATE_META = """
            CREATE TABLE metadata (Key string not null, Value string not null)
        """
        SQL_INSERT_META = "INSERT INTO metadata VALUES (?,?)"
        conn = sqlite3.connect(path)
        c = conn.cursor()
        c.execute(SQL_CREATE_ENTRIES)
        c.execute(SQL_CREATE_P2P)
        c.execute(SQL_CREATE_META)
        c.execute(SQL_INSERT_META, ["version", "0.1.14"])
        c.execute(SQL_INSERT_META, ["staleProteinMapping", "true"])
        conn.commit()
        return conn, c

    def write(
        self,
        chunksize: Optional[Union[None, int]]
    ):
        self._write_entries(
            index=False,
            if_exists='append',
            method='multi',
            chunksize=chunksize)
        self._write_p2p(
            index=False,
            if_exists='append',
            method='multi',
            chunksize=chunksize)

    def _write_entries(self, *args, **kwargs):
        conn = sqlite3.connect(self.path)
        self.entries.to_sql(name='entries', con=conn, *args, **kwargs)
        conn.commit()

    def _write_p2p(self, *args, **kwargs):
        conn = sqlite3.connect(self.path)
        self.p2p.to_sql(name='peptidetoprotein', con=conn, *args, **kwargs)
        conn.commit()

    def prepare_spectrum(self):
        pass
