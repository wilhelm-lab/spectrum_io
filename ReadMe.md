# Status

[![pipeline status](https://gitlab.lrz.de/compmass/prosit/tools/prosit_io/badges/develop/pipeline.svg)](https://gitlab.lrz.de/compmass/prosit/tools/prosit_io/-/commits/develop) [![coverage report](https://gitlab.lrz.de/compmass/prosit/tools/prosit_io/badges/develop/coverage.svg)](https://gitlab.lrz.de/compmass/prosit/tools/prosit_io/-/commits/develop)

# Updating fundamentals dependency

For some reason poetry does install updated versions of git repository dependencies even though they are written to the poetry.lock file (https://github.com/python-poetry/poetry/issues/2921).
To circumvent this, uninstall the package first before running `poetry install`:
```
poetry update fundamentals
pip uninstall fundamentals
poetry install
```
