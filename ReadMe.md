[![PyPI](https://img.shields.io/pypi/v/spectrum_io.svg)](https://pypi.org/project/spectrum_io/)
[![Python Version](https://img.shields.io/pypi/pyversions/spectrum_io)](https://pypi.org/project/spectrum_io)
[![License](https://img.shields.io/github/license/wilhelm-lab/spectrum_io)](https://opensource.org/licenses/MIT)
[![Read the Docs](https://img.shields.io/readthedocs/spectrum_io/latest.svg?label=Read%20the%20Docs)](https://spectrum-io.readthedocs.io/)
[![Build](https://github.com/wilhelm-lab/spectrum_io/workflows/Build%20spectrum_io%20Package/badge.svg)](https://github.com/wilhelm-lab/spectrum_io/actions?workflow=Package)
[![Tests](https://github.com/wilhelm-lab/spectrum_io/workflows/Run%20spectrum_io%20Tests/badge.svg)](https://github.com/wilhelm-lab/spectrum_io/actions?workflow=Tests)
[![Codecov](https://codecov.io/gh/wilhelm-lab/spectrum_io/branch/main/graph/badge.svg)](https://codecov.io/gh/wilhelm-lab/spectrum_io)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Spectrum IO: File / Data Conversion for Mass Spec data within the Oktoberfest ecosystem

spectrum_io is a package primarily developed for usage within the rescoring and spectral library generation pipeline oktoberfest (https://github.com/wilhelm-lab/oktoberfest).

It provides the following functionalities:

-   read search results from different search engines (MaxQuant, MSFragger, Sage, Xisearch) or a generic csv format and transform them to the internal format for rescoring with oktoberfest
-   extraction of MS2 level spectra from .RAW files and conversion to to mzML for rescoring with oktoberfest
-   spectra extraction from .d folders, conversion to .hdf5 format, and aggregation to MS2 level with metadata from a MaxQuant search for timsTOF rescoring with oktoberfest
-   in-silico digestion of a fasta file with various configuration options (protease, missed cleavages, length of peptides, fragmentation, ...) for spectral library generation with oktoberfest
-   write spectral libraries in dlib, msp, or spectronaut(csv) format
-   parquet file creation for peptide prediction model development and refinement within DLOmix

## Documentation

The official documentation can be found at https://spectrum-fundamentals.readthedocs.io

## How to cite

Please always cite the main publication:

[Oktoberfest] Picciani M, Gabriel W, Giurcoiu VG et al. (2023), _Oktoberfest: Open-source spectral library generation and rescoring pipeline based on Prosit_, [Proteomics](https://doi.org/10.1002/pmic.202300112)
