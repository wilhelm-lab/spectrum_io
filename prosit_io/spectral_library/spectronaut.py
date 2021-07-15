import pandas as pd
import numpy as np
from .spectral_library import SpectralLibrary


class Spectronaut(SpectralLibrary):
    # Check spectronaut folder for output format.

    def write(self):

        out = open(self.out_path, "w")
        out.truncate()
        out.close()

        n = 10000  # split df into chunks of size n
        initial = True
        for group, segment in self.spectra_output.groupby(np.arange(len(self.spectra_output)) // n):
            segment = segment.explode(['intensities', 'fragment_mz', 'fragment_types', 'fragment_numbers', 'fragment_charges'])
            segment = segment[segment['intensities'] > 0]  # set to >= if 0 should be kept
            segment.to_csv(self.out_path, mode='a', header=initial)
            if initial:
                initial = False

    def prepare_spectrum(self):
        intensities = self.grpc_output[list(self.grpc_output)[0]]['intensity']
        fragment_mz = self.grpc_output[list(self.grpc_output)[0]]['fragmentmz']
        annotation = self.grpc_output[list(self.grpc_output)[0]]['annotation']
        fragment_types = annotation['type']
        fragment_numbers = annotation['number']
        fragment_charges = annotation['charge']
        irt = self.grpc_output[list(self.grpc_output)[1]]
        proteotypicity = self.grpc_output[list(self.grpc_output)[2]]

        modified_sequences = self.spectra_input['MODIFIED_SEQUENCE']

        def label_sequences(sequences):
            return sequences

        def strip_sequences(sequences):
            return sequences

        labelled_sequences = label_sequences(modified_sequences)
        stripped_peptide = strip_sequences(modified_sequences)
        charges = self.spectra_input['PRECURSOR_CHARGE']
        precursor_masses = self.spectra_input['MASS']
        precursor_mz = (precursor_masses + charges) / charges

        inter_df = pd.DataFrame(data={'ModifiedPeptide': modified_sequences, 'LabeledPeptide': labelled_sequences,
                                      'StrippedPeptide': stripped_peptide, 'PrecursorCharge': charges,
                                      'PrecursorMz': precursor_mz})
        inter_df['iRT'], inter_df['proteotypicity'] = irt.tolist(), proteotypicity.tolist()
        inter_df['intensities'], inter_df['fragment_mz'] = intensities.tolist(), fragment_mz.tolist()
        inter_df['fragment_types'] = fragment_types.tolist()
        inter_df['fragment_numbers'] = fragment_numbers.tolist()
        inter_df['fragment_charges'] = fragment_charges.tolist()

        self.spectra_output = inter_df
