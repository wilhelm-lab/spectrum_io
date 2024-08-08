Spectrum IO: File / Data Conversion for Mass Spec data within the Oktoberfest ecosystem
=======================================================================================

|PyPI| |Python Version| |License| |Read the Docs| |Build| |Tests| |Codecov| |pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/spectrum_io.svg
   :target: https://pypi.org/project/spectrum_io/
   :alt: PyPI
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/spectrum_io
   :target: https://pypi.org/project/spectrum_io
   :alt: Python Version
.. |License| image:: https://img.shields.io/github/license/wilhelm-lab/spectrum_io
   :target: https://opensource.org/licenses/MIT
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/spectrum_io/latest.svg?label=Read%20the%20Docs
   :target: https://spectrum-io.readthedocs.io/
   :alt: Read the documentation at https://spectrum-io.readthedocs.io/
.. |Build| image:: https://github.com/wilhelm-lab/spectrum_io/workflows/Build%20spectrum_io%20Package/badge.svg
   :target: https://github.com/wilhelm-lab/spectrum_io/actions?workflow=Package
   :alt: Build Package Status
.. |Tests| image:: https://github.com/wilhelm-lab/spectrum_io/workflows/Run%20spectrum_io%20Tests/badge.svg
   :target: https://github.com/wilhelm-lab/spectrum_io/actions?workflow=Tests
   :alt: Run Tests Status
.. |Codecov| image:: https://codecov.io/gh/wilhelm-lab/spectrum_io/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/wilhelm-lab/spectrum_io
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

spectrum_io is a package primarily developed for usage within the rescoring and spectral library generation pipeline oktoberfest (https://github.com/wilhelm-lab/oktoberfest).

It provides the following functionalities:
 -   read search results from different search engines (MaxQuant, MSFragger, Sage, Xisearch) or a generic csv format and transform them to the internal format for rescoring with oktoberfest
 -   extraction of MS2 level spectra from .RAW files and conversion to to mzML for rescoring with oktoberfest
 -   spectra extraction from .d folders, conversion to .hdf5 format, and aggregation to MS2 level with metadata from a MaxQuant search for timsTOF rescoring with oktoberfest
 -   in-silico digestion of a fasta file with various configuration options (protease, missed cleavages, length of peptides, fragmentation, ...) for spectral library generation with oktoberfest
 -   write spectral libraries in dlib, msp, or spectronaut(csv) format
 -   parquet file creation for peptide prediction model development and refinement within DLOmix
