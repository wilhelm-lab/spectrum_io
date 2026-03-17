.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	if exist build rd /s /q build
	if exist build rd /s /q dist
	if exist build rd /s /q .eggs
	for /d /r . %%d in (*.egg-info) do @if exist "%%d" echo "%%d" && rd /s/q "%%d"
	del /q /s /f .\*.egg


clean-pyc: ## remove Python file artifacts
	del /s /f /q .\*.pyc
	del /s /f /q .\*.pyo
	del /s /f /q .\*~
	for /d /r . %%d in (*__pycache__) do @if exist "%%d" echo "%%d" && rd /s/q "%%d"

clean-test: ## remove test and coverage artifacts
	if exist .tox rd /s /q .tox
	if exist .coverage rd /s /q .coverage
	if exist htmlcov rd /s /q htmlcov
	if exist .pytest_cache rd /s /q .pytest_cache

lint: ## check style with ruff
	poetry run ruff check spectrum_io tests

format: ## auto-format code with ruff
	poetry run ruff format spectrum_io tests

typecheck: ## run type checking with mypy
	poetry run mypy spectrum_io

test: ## run tests quickly with the default Python
	poetry run pytest

check: format typecheck test ## run all checks (format, typecheck, test) - simulates CI
	@echo ✓ All checks passed!

coverage: ## check code coverage quickly with the default Python
	poetry run coverage run --source spectrum_io -m pytest
	poetry run coverage report -m
	poetry run coverage html
	$(BROWSER) htmlcov\index.html

docs: ## build HTML documentation with sphinx
	poetry run sphinx-build docs docs\_build\html

docs-serve: ## build docs and serve locally with live reload
	poetry run sphinx-autobuild docs docs\_build\html --open-browser

release: dist ## package and upload a release
	poetry release

dist: clean-build clean-pyc ## builds source and wheel package
	poetry build

install: ## install the package to the active Python's site-packages
	poetry install
