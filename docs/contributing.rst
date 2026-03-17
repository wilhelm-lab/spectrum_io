Contributor Guide
=================

Thank you for your interest in improving this project.
This project is open-source under the `MIT license`_ and
highly welcomes contributions in the form of bug reports, feature requests, and pull requests.

Here is a list of important resources for contributors:

- `Source Code`_
- `Documentation`_
- `Issue Tracker`_
- `Code of Conduct`_

.. _MIT license: https://opensource.org/licenses/MIT
.. _Source Code: https://github.com/wilhelm-lab/spectrum_io
.. _Documentation: https://spectrum-io.readthedocs.io/
.. _Issue Tracker: https://github.com/wilhelm-lab/spectrum_io/issues

How to report a bug
-------------------

Report bugs on the `Issue Tracker`_.


How to request a feature
------------------------

Request features on the `Issue Tracker`_.


How to set up your development environment
------------------------------------------

You need Python 3.10+ and the following tools:

- Poetry_

You can install Poetry with:

.. code:: console

    $ pip install poetry

Install the package with development requirements:

.. code:: console

   $ make install

You can now run an interactive Python session,
or the command-line interface:

.. code:: console

   $ poetry run python
   $ poetry run spectrum_io

.. _Poetry: https://python-poetry.org/


How to test the project
-----------------------

Run the full test suite with all checks (formatting, type checking, and unit tests):

.. code:: console

   $ make check

Run only unit tests:

.. code:: console

   $ make test

Run type checking with mypy:

.. code:: console

   $ make typecheck

Run code formatting and linting checks:

.. code:: console

   $ make lint

Auto-format your code with ruff:

.. code:: console

   $ make format

Unit tests are located in the ``tests`` directory,
and are written using the pytest_ testing framework.

.. _pytest: https://pytest.readthedocs.io/

How to build and view the documentation
---------------------------------------

This project uses Sphinx_ together with several extensions to build the documentation.

To install all required dependencies including documentation tools, use `make install` which installs all dependencies from pyproject.toml:

.. code:: console

    $ make install

To build the documentation run:

.. code:: console

    $ make docs

This will build the documentation and automatically open it in your default browser.

Alternatively, to manually build the documentation:

.. code:: console

    $ cd docs
    $ poetry run make html

The generated static HTML files can be found in the `_build/html` folder.
Simply open them with your favorite browser.

.. _sphinx: https://www.sphinx-doc.org/en/master/

How to submit changes
---------------------

Open a `pull request`_ to submit changes to this project against the ``development`` branch.

Your pull request needs to meet the following guidelines for acceptance:

- The test suite must pass without errors (run `make check` locally).
- Include unit tests. This project maintains a high code coverage.
- If your changes add functionality, update the documentation accordingly.

To run linting and code formatting checks before committing your change, you can install pre-commit as a Git hook by running:

.. code:: console

   $ pre-commit install

It is recommended to open an issue before starting work on anything.
This will allow a chance to talk it over with the owners and validate your approach.

.. _pull request: https://github.com/wilhelm-lab/spectrum_io/pulls


How to make a release
---------------------

Releases are published to PyPI automatically when a GitHub Release is published.
The version string lives only in ``pyproject.toml`` — ``__version__`` is read from
the installed package metadata at runtime.

Release Drafter continuously updates a draft GitHub Release with an accumulated changelog
from merged PR labels and a suggested next version (e.g. ``0.9.1``). It is a changelog
generator — it never modifies any file in the repository.

**Branch model:** ``development`` is the integration branch; ``main`` mirrors exactly
what is published on PyPI. The release tag is created on ``development`` and subsequently
merged into ``main``.

1. **Check the draft release** on GitHub to see the suggested next version (e.g. ``0.10.0``).
   The version is inferred automatically from the labels on merged PRs since the last release.

2. **Bump the version on** ``development``:

   .. code:: console

      $ git checkout development && git pull
      $ poetry version <next-version>   # e.g. poetry version 0.10.0
      $ git add pyproject.toml
      $ git commit -m "bump version to $(poetry version -s)"
      $ git push origin development

3. **Publish the draft release** on GitHub.
   The draft already targets ``development`` (set via ``commitish: development`` in
   ``.github/release-drafter.yml``), so no target branch change is needed.
   Clicking **Publish release** triggers the CI workflow, which:

   - Re-runs the full CI suite as a hard gate.
   - Builds the wheel and sdist with ``poetry build``.
   - Publishes to PyPI via OIDC Trusted Publishing (no secrets required).
   - Creates the tag ``v<next-version>`` on ``development``.

4. **Merge the tagged commit into** ``main`` so that ``main`` reflects the release:

   .. code:: console

      $ git checkout main && git pull
      $ git merge v<next-version> --no-ff -m "release: v<next-version>"
      $ git push origin main

.. _Code of Conduct: CODE_OF_CONDUCT.rst
