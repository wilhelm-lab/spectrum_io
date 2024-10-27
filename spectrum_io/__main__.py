#!/usr/bin/env python
"""Command-line interface."""
import click
from rich import traceback


@click.command()
@click.version_option(version="0.6.4", message=click.style("spectrum_io Version: 0.6.4"))
def main() -> None:
    """spectrum_io."""


if __name__ == "__main__":
    traceback.install()
    main(prog_name="spectrum_io")  # pragma: no cover
