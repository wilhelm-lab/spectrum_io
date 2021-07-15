from spectral_library import SpectralLibrary


class DLib(SpectralLibrary):


    def __init__(
        self,
        grpc_out: dict,
        intensity_model: str,
        iRT_model: str
        sequence: str,
        charge: str,
        collision_energy: str
        path: str,
    ):
        self.intensities = grpc_out[intensity_model]['intensity']
        self.masses = grpc_out[intensity_model]['fragmentmz']
        self.iRT = grpc_out[iRT_model]
        self.sequence = sequence
        self.charge = charge
        self.collision_energy = collision_energy
        self.path = path

    def insert_entries(self, entries):
        SQL_ENTRY = 'INSERT INTO entries VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        conn = sqlite3.connect(self.path)
        c = conn.cursor()
        c.executemany(SQL_ENTRY, [tuple(e) for e in entries])
        return conn.commit()

    def insert_p2p(self, p2p):
        SQL_P2P = 'INSERT INTO peptidetoprotein VALUES (?,?,?)'
        conn = sqlite3.connect(self.path)
        c = conn.cursor()
        c.executemany(SQL_P2P, [tuple(e) for e in p2p])
        return conn.commit()

    def create_database(self):
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
        #SQL_CREATE_META = """
        #    CREATE TABLE metadata (Key string not null, Value string not null)
        #"""
        #SQL_INSERT_META = "INSERT INTO metadata VALUES (?,?)"
        conn = sqlite3.connect(self.path)
        c = conn.cursor()
        c.execute(SQL_CREATE_ENTRIES)
        c.execute(SQL_CREATE_P2P)
        #c.execute(SQL_CREATE_META)
        #c.execute(SQL_INSERT_META, ["version", "0.1.14"])
        #c.execute(SQL_INSERT_META, ["staleProteinMapping", "true"])
        conn.commit()
        return conn, c


    # Check dlib folder for output format.
    def write(self):
        pass

    def prepare_spectrum(self):
        pass
