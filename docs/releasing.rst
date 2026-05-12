Releasing spectrum_io
====================

This guide describes the process for releasing a new version of spectrum_io.

Prerequisites
-------------

- You have push access to the ``development`` and ``main`` branches
- **PyPI OIDC trusted publisher configured** (one-time setup):
  1. Go to https://pypi.org/manage/account/publishing/
  2. Click "Add a new pending publisher"
  3. Fill in:
     - **PyPI Project Name**: ``spectrum_io``
     - **Owner**: ``wilhelm-lab``
     - **Repository name**: ``spectrum_io``
     - **Workflow name**: ``publish_package.yml``
     - **Environment name**: ``release``
  4. Approve the pending publisher (you'll receive a notification on PyPI)
- All changes for the release are merged into the ``development`` branch

Release Process
---------------

1. **Update the version in pyproject.toml**

   Edit ``pyproject.toml`` and update the ``version`` field in the ``[tool.poetry]`` section:

   .. code:: toml

       [tool.poetry]
       name = "spectrum_io"
       version = "0.9.0"  # Update this to the new version

   The version is the single source of truth. It will be automatically read by:
   - ``spectrum_io.__init__.py`` (via importlib.metadata)
   - ``docs/conf.py`` (via importlib.metadata)
   - CI/CD workflows

2. **Commit the version change**

   .. code:: console

       $ git add pyproject.toml
       $ git commit -m "chore: bump version to 0.9.0"
       $ git push origin development

3. **Create a pull request to main**

   Create a PR from ``development`` to ``main``. Once approved and merged, this will trigger the release workflow.

4. **Release Drafter automatically creates a draft release**

   When you push to ``development``, the Release Drafter workflow automatically creates a draft release on GitHub with:
   - A dynamically determined version (via ``$RESOLVED_VERSION``)
   - Auto-generated changelog based on merged pull requests
   - Categorized changes (Features, Bug Fixes, Maintenance, Dependencies)

5. **Publish the release**

   - Go to the `Releases page <https://github.com/wilhelm-lab/spectrum_io/releases>`_ on GitHub
   - Edit the draft release
   - Verify the version, changelog, and notes
   - Click "Publish release"

   This will:
   - Create a git tag (e.g., ``v0.9.0``)
   - Trigger the ``Publish spectrum_io to PyPI`` workflow
   - Automatically build the package and upload to PyPI (if ``PYPI_TOKEN`` is configured)

6. **Verify the release**

   - Check that the new version is available on `PyPI <https://pypi.org/project/spectrum_io/>`_
   - Verify the documentation is updated on `ReadTheDocs <https://spectrum-io.readthedocs.io/>`_

Versioning
----------

Spectrum_io uses `Semantic Versioning <https://semver.org/>`_:

- **MAJOR.MINOR.PATCH** (e.g., 1.2.3)
- MAJOR: Incompatible API changes
- MINOR: New functionality, backwards compatible
- PATCH: Bug fixes, backwards compatible

Version Labeling
~~~~~~~~~~~~~~~~

Pull requests can be labeled to automatically determine version bumps:

- ``major`` label → Version bump: Increment MAJOR
- ``minor`` label → Version bump: Increment MINOR
- ``patch`` label (or no label) → Version bump: Increment PATCH

The Release Drafter uses these labels to automatically bump the version in draft releases.

Manual Version in pyproject.toml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you prefer manual versioning, simply edit ``pyproject.toml`` directly and the CI/CD pipelines will read the exact version you specify.

Automated Version Reading
~~~~~~~~~~~~~~~~~~~~~~~~~~

The version in ``pyproject.toml`` is dynamically read at runtime:

- **Python package**: ``spectrum_io.__version__`` uses ``importlib.metadata.version("spectrum_io")``
- **Documentation**: ``docs/conf.py`` reads the version dynamically
- **Workflows**: CI/CD workflows can read the version from the built package

This ensures a single source of truth and eliminates version sync issues.

Continuous Integration
----------------------

The following workflows run automatically:

- **run_tests.yml**: Tests on every push to ``development``, ``main``, and ``release/*`` branches (runs on Python 3.10, 3.11, 3.12 on Ubuntu)
- **build_package.yml**: Validates package building on every push and pull request
- **publish_docs.yml**: Builds and deploys docs to GitHub Pages on pushes to ``main`` and ``master`` branches
- **publish_package.yml**: Uploads to PyPI on release publication (uses OIDC for authentication - no token needed!)
- **release_drafter.yml**: Auto-generates draft releases with changelogs

Troubleshooting
---------------

**Release not appearing on PyPI after publishing?**

Ensure the OIDC trusted publisher is configured correctly in PyPI. See Prerequisites section. The ``publish_package.yml`` workflow will fail if the trusted publisher isn't set up properly.

If you previously used a token-based approach, you can verify OIDC is working by checking the workflow run details in the Actions tab.
