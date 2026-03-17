Usage
=====

Spectrum IO is a Python library for reading and processing mass spectrometry data in various formats.

Basic Example
~~~~~~~~~~~~~

.. code-block:: python

    from spectrum_io.raw import ThermoRaw
    from spectrum_io.search_result import MaxQuant

    # Read raw mass spectrometry data
    raw = ThermoRaw("path/to/file.raw")

    # Read search results
    results = MaxQuant("path/to/msms.txt")

For detailed API documentation, see the :doc:`reference`.
