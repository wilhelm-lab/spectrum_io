import tempfile

import pytest
import numpy as np
import pandas as pd

import prosit_io.spectral_library.msp as msp


class TestMspPrepareSpectrum:
    def test_prepare_spectrum(self, spectra_input, grpc_dict):
        output_path = ""
        msp_lib = msp.MSP(spectra_input, grpc_dict, output_path)
        msp_lib.prepare_spectrum()
    
    def test_write(self, spectra_input, grpc_dict):
        out_file = tempfile.NamedTemporaryFile(delete=False)
        print(out_file.name)
        msp_lib = msp.MSP(spectra_input, grpc_dict, out_file.name)
        msp_lib.prepare_spectrum()
        msp_lib.write()
        assert out_file.read() == b"""Name: AAACCCCKR/1
MW: 124.4
Comment: Parent=124.4 Collision_energy=10.0 Mods=2/3,C,Carbamidomethyl/5,C,Carbamidomethyl ModString=AAACCCCKR//Carbamidomethyl@C3; Carbamidomethyl@C5/1 iRT=982.12 proteotypicity=123.1
Num peaks: 3
0.9	0.1	"b1/0.0ppm"
0.8	0.2	"y1^2/0.0ppm"
0.7	0.3	"N2^3/0.0ppm"
Name: AAACILKKR/2
MW: 1617.05
Comment: Parent=1617.05 Collision_energy=20.0 Mods=0 ModString=AAACILKKR///2 iRT=382.12 proteotypicity=234.2
Num peaks: 3
0.6	0.4	"b1^2/0.0ppm"
0.5	0.5	"y3^3/0.0ppm"
0.4	0.6	"N5/0.0ppm"
"""
        
@pytest.fixture
def spectra_input():
    spectra_input = pd.DataFrame()
    spectra_input['MODIFIED_SEQUENCE_SPEC'] = ['AAACCCC', 'AAACILK']
    spectra_input['MODIFIED_SEQUENCE'] = ['AAAC(U:4)CC(U:4)CKR', 'AAACILKKR']
    spectra_input['MASS'] = [123.4, 3232.1]
    spectra_input['COLLISION_ENERGY'] = [10.0, 20.0]
    spectra_input['PRECURSOR_CHARGE'] = [1, 2]
    return spectra_input

@pytest.fixture
def grpc_dict():
    grpc_dict = { 
        "model" : {
            "intensity" : np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
            "fragmentmz" : np.array([[0.9, 0.8, 0.7], [0.6, 0.5, 0.4]]),
            "annotation": {
                "charge" : np.array([[1, 2, 3], [2, 3, 1]]),
                "number" : np.array([[1, 1, 2], [1, 3, 5]]),
                "type": np.array([['b', 'y', 'N'], ['b', 'y', 'N']])
            }
        },
        "model_irt" : np.array([[982.12], [382.12]]),
        "model_proteotypicity" : np.array([[123.1], [234.2]])
    }
    return grpc_dict
