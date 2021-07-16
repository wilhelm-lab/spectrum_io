import pytest
import numpy as np
import pandas as pd

import prosit_io.msp as msp


class TestMspPrepareSpectrum:
    def test_prepare_spectrum(self):        
        spectra_input = pd.DataFrame()
        spectra_input['MODIFIED_SEQUENCE_SPEC'] = ['AAACCCC', 'AAACILK']
        spectra_input['MODIFIED_SEQUENCE'] = ['AAACCCCKR', 'AAACILKKR']
        spectra_input['MASS'] = [123.4, 3232.1]
        spectra_input['COLLISION_ENERGY'] = [10.0, 20.0]
        spectra_input['PRECURSOR_CHARGE'] = [1, 2]
        
        grpc_dict = { 
            "model" : {
                "intensity" : [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
                "fragmentmz" : [[0.9, 0.8, 0.7], [0.6, 0.5, 0.4]],
                "annotation": {
                    "charge" : [[1, 2, 3], [2, 3, 1]],
                    "number" : [[1, 1, 2], [1, 3, 5]],
                    "type": [['b', 'y', 'N'], ['b', 'y', 'N']]
                }
            },
            "model_irt" : [[982.12], [382.12]],
            "model_proteotypicity" : [[123.1], [234.2]]
        }
        output_path = ""
        msp = MSP(spectra_input, grpc_dict, output_path)
        msp.prepare_spectrum()
