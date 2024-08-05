"""SpectrumIO: File / Data Conversion for Mass Spec data within the Oktoberfest ecosystem."""

from datetime import datetime

__author__ = """The Oktoberfest development team (Wilhelmlab at Technical University of Munich)"""
__copyright__ = f"Copyright {datetime.now():%Y}, Wilhelmlab at Technical University of Munich"
__license__ = "MIT"
__version__ = "0.7.0"

import logging
import logging.handlers
import sys

from spectrum_io import d, file, raw

# from .search_result import MaxQuant
# from .spectral_library import DLib, Spectronaut

CONSOLE_LOG_LEVEL = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class _InfoWarningFilter(logging.Filter):
    def filter(self, record):
        return CONSOLE_LOG_LEVEL <= record.levelno <= logging.WARNING


if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    # add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.addFilter(_InfoWarningFilter())
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # add error handler
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)
    logger.addHandler(error_handler)
else:
    logger.info("Logger already initialized. Resuming normal operation.")
