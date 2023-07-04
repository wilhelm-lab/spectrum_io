# Spectrum IO

Spectrum IO is a package primarily developed for usage within oktoberfest (https://github.com/wilhelm-lab/oktoberfest). It handles file conversions and input / output operations for oktoberfest.

## Installation

### Prerequisites

If you want to convert raw files to mzml, make sure you have ThermoRawFileParser (https://github.com/compomics/ThermoRawFileParser) installed.

If you are on linux or MacOS, make sure mono (https://www.mono-project.com/) is installed (for ThermoRawFileParser).

### Using pip

```bash
pip install oktoberfest
```

## Features

-   Read search results from different search engines (Mascot, MaxQuant, MSFragger, MS Amanda) and transform them to the internal format used by oktoberfest
-   Read thermo raw files and convert them to mzml, required by oktoberfest
-   Read a fasta file and digest with various configurations (protease, missed cleavages, length of peptides, fragmentation, ...) for spectral library generation
-   Create spectral libraries from peptide lists and output as dlib, msp or spectronaut(csv) format
-   read and write data created as part of oktoberfest in hdf5 libraries

## Documentation

Please refer to https://spectrum-io.readthedocs.io for further documentation.

## License

The project is licensed under the [MIT license](https://github.com/wilhelm-lab/spectrum_io/blob/main/LICENSE).
