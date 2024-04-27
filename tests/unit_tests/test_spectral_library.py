from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from spectrum_io.spectral_library import MSP, Spectronaut


class TestMSP:
    """Class to test msp."""

    def test_write(self, spectra_input, data, metadata):
        """Test write to file."""
        out_file = Path(__file__).parent / "test.msp"
        msp_lib = MSP(out_file)
        msp_lib.write(data, metadata)
        with open(out_file) as test_file:
            file_content = test_file.read()

        file_content = file_content.replace("\r", "")  # explicitly remove to work for windows
        anticipated_content = (
            "Name: AAACCCCKR/1\n"
            "MW: 124.407276467\n"
            "Comment: Parent=124.40727647 Collision_energy=10.0 Mods=2/3,C,Carbamidomethyl/5,C,Carbamidomethyl "
            "ModString=AAACCCCKR//Carbamidomethyl@C3; Carbamidomethyl@C5/1 iRT=982.12\n"
            "Num peaks: 2\n"
            '0.80000000	0.2000	"b1/0.0ppm"\n'
            '0.30000000	0.8000	"b2^2/0.0ppm"\n'
            "Name: AAACILKKR/2\n"
            "MW: 1617.057276467\n"
            "Comment: Parent=1617.05727647 Collision_energy=20.0 Mods=0 ModString=AAACILKKR///2 iRT=382.12\n"
            "Num peaks: 3\n"
            '0.50000000	0.5000	"b1/0.0ppm"\n'
            '0.40000000	0.6000	"y2^2/0.0ppm"\n'
            '0.30000000	0.0010	"b2^2/0.0ppm"\n'
        )
        assert file_content == anticipated_content

        out_file.unlink()


class TestSpectronaut:
    """Class to test msp."""

    def test_write(self, spectra_input, data, metadata):
        """Test write to file."""
        out_file = Path(__file__).parent / "test.csv"
        msp_lib = Spectronaut(out_file)
        msp_lib.write(data, metadata)
        with open(out_file) as test_file:
            file_content = test_file.read()

        file_content = file_content.replace("\r", "")  # explicitly remove to work for windows
        anticipated_content = (
            "ModifiedPeptide,LabeledPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,CollisionEnergy,"
            "RelativeFragmentIntensity,FragmentMz,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType\n"
            "_AAAC[Carbamidomethyl (C)]CC[Carbamidomethyl (C)]CKR_,AAACCCCKR,AAACCCCKR,1,124.40727647,982.12,"
            "10.0,0.2000,0.80000000,1,b,1,noloss\n"
            "_AAAC[Carbamidomethyl (C)]CC[Carbamidomethyl (C)]CKR_,AAACCCCKR,AAACCCCKR,1,124.40727647,982.12,"
            "10.0,0.8000,0.30000000,2,b,2,noloss\n"
            "_AAACILKKR_,AAACILKKR,AAACILKKR,2,1617.05727647,382.12,20.0,0.5000,0.50000000,1,b,1,noloss\n"
            "_AAACILKKR_,AAACILKKR,AAACILKKR,2,1617.05727647,382.12,20.0,0.6000,0.40000000,2,y,2,noloss\n"
            "_AAACILKKR_,AAACILKKR,AAACILKKR,2,1617.05727647,382.12,20.0,0.0010,0.30000000,2,b,2,noloss\n"
        )
        assert file_content == anticipated_content

        out_file.unlink()


@pytest.fixture
def spectra_input():
    """Test spectra input."""
    spectra_input = pd.DataFrame()
    spectra_input["MODIFIED_SEQUENCE_SPEC"] = ["AAACCCC", "AAACILK"]
    spectra_input["MODIFIED_SEQUENCE"] = ["AAAC[UNIMOD:4]CC[UNIMOD:4]CKR", "AAACILKKR"]
    spectra_input["MASS"] = [123.4, 3232.1]
    spectra_input["COLLISION_ENERGY"] = [10.0, 20.0]
    spectra_input["PRECURSOR_CHARGE"] = [1, 2]
    return spectra_input


@pytest.fixture
def data():
    """Creates data dictionary."""
    return {
        "intensities": np.array([[1e-5, 0.2, 0.3, 0.8], [-1, 0.5, 0.6, 0.001]]),
        "mz": np.array([[0.9, 0.8, -1, 0.3], [0.6, 0.5, 0.4, 0.3]]),
        "annotation": np.array([[b"y1+1", b"b1+1", b"y2+2", b"b2+2"], [b"y1+1", b"b1+1", b"y2+2", b"b2+2"]]),
        "irt": np.array([[982.12], [382.12]]),
    }


@pytest.fixture
def metadata():
    """Creates metadata df."""
    return pd.DataFrame(
        {
            "SEQUENCE": ["AAACCCCKR", "AAACILKKR"],
            "MODIFIED_SEQUENCE": ["AAAC[UNIMOD:4]CC[UNIMOD:4]CKR", "AAACILKKR"],
            "PRECURSOR_CHARGE": [1, 2],
            "MASS": [123.4, 3232.1],
            "COLLISION_ENERGY": [10.0, 20.0],
        }
    )
