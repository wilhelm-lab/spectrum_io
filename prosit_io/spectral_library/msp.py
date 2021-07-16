import re

from spectral_library import SpectralLibrary


class MSP(SpectralLibrary):
    # Check msp folder for output format.
    def write(self):
        """
        Writing method; writes intermediate dataframe as msp format spectra
        :return: None
        """
        out = open(self.out_path, "w")
        
        modification_regexp = re.compile(
        for spectrum in self.spectra_output.iterrows():
            out.write(f"Name: {spectrum['StrippedPeptide']}/{spectrum['PrecursorCharge']}\n")
            out.write(f"MW: {spectrum['PrecursorMz']}\n")
            modifications = MSP.get_modifications(spectrum['ModifiedPeptide'])
            out.write(f"Comment: Parent={spectrum['PrecursorMz']} "
                      f"Collision_energy={spectrum['CollisionEnergy']} "
                      f"Mods={len(modifications)} "
                      f"ModString={}/{spectrum['PrecursorCharge']} "
                      f"iRT={spectrum['iRT']} "
                      f"proteotypicity={spectrum['proteotypicity']}\n") 
            for fmz, fintensity, ftype, fcharge, fnumber in zip(spectrum['fragment_mz'], spectrum['intensities'], spectrum['fragment_types'], spectrum['fragment_charges'], spectrum['fragment_numbers']):
                out.write(f'{fmz}\t{fintensity}\t"{ftype}{fnumber}"' + (f'^{fcharge}' if fcharge != 1 else '') + '/0.0ppm')
        out.close()            

    def prepare_spectrum(self):
        """
        Converts grpc output and metadata dataframe into msp format
        :return: 
        """
        intensities = self.grpc_output[list(self.grpc_output)[0]]['intensity']
        fragment_mz = self.grpc_output[list(self.grpc_output)[0]]['fragmentmz']
        annotation = self.grpc_output[list(self.grpc_output)[0]]['annotation']
        fragment_types = annotation['type']
        fragment_numbers = annotation['number']
        fragment_charges = annotation['charge']
        irt = self.grpc_output[list(self.grpc_output)[1]]
        proteotypicity = self.grpc_output[list(self.grpc_output)[2]]
        
        modified_sequences = self.spectra_input['MODIFIED_SEQUENCE']
        collision_energies = self.spectra_input['COLLISION_ENERGY']

        stripped_peptide = internal_without_mods(modified_sequences, remove_underscores=True)
        charges = self.spectra_input['PRECURSOR_CHARGE']
        precursor_masses = self.spectra_input['MASS']
        precursor_mz = (precursor_masses + charges) / charges

        inter_df = pd.DataFrame(data={'ModifiedPeptide': modified_sequences,
                                      'StrippedPeptide': stripped_peptide, 'PrecursorCharge': charges,
                                      'PrecursorMz': precursor_mz, 'PrecursorMass': precursor_masses,
                                      'CollisionEnergy': collision_energies})
        inter_df['iRT'], inter_df['proteotypicity'] = irt.tolist(), proteotypicity.tolist()
        inter_df['intensities'], inter_df['fragment_mz'] = intensities.tolist(), fragment_mz.tolist()
        inter_df['fragment_types'] = fragment_types.tolist()
        inter_df['fragment_numbers'] = fragment_numbers.tolist()
        inter_df['fragment_charges'] = fragment_charges.tolist()
        
        self.spectra_output = inter_df
    
    @staticmethod
    def get_modifications(modified_sequence):
        re.
        
