[![PyPI](https://img.shields.io/pypi/v/oktoberfest.svg)](https://pypi.org/project/oktoberfest/)
[![Python Version](https://img.shields.io/pypi/pyversions/oktoberfest)](https://pypi.org/project/oktoberfest)
[![License](https://img.shields.io/github/license/wilhelm-lab/oktoberfest)](https://opensource.org/licenses/MIT)
[![Read the Docs](https://img.shields.io/readthedocs/oktoberfest/latest.svg?label=Read%20the%20Docs)](https://oktoberfest.readthedocs.io/)
[![Build](https://github.com/wilhelm-lab/oktoberfest/workflows/Build%20oktoberfest%20Package/badge.svg)](https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Package)
[![Tests](https://github.com/wilhelm-lab/oktoberfest/workflows/Run%20oktoberfest%20Tests/badge.svg)](https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Tests)
[![Codecov](https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg)](https://codecov.io/gh/wilhelm-lab/oktoberfest)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

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
