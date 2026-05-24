"""
Centralized Logging and Transparency Module for Hydrosheaf.
Provides scientific audit logging and debugging capabilities.
"""

import logging
import sys
from pathlib import Path
from typing import Optional
import functools
import time


# Create a custom logger class to add scientific context
class ScienceLogger(logging.Logger):
    def science(self, msg, *args, **kwargs):
        """Log a scientific decision or significant model state change (Level 25)."""
        if self.isEnabledFor(25):
            self._log(25, msg, args, **kwargs)


logging.setLoggerClass(ScienceLogger)
logging.addLevelName(25, "SCIENCE")

_SETUP_DONE = False


def setup_logging(
    verbose: bool = False, log_file: Optional[str] = None, capture_warnings: bool = True
) -> None:
    """
    Configure the Hydrosheaf logging system.

    Parameters
    ----------
    verbose : bool
        If True, set console level to DEBUG. Default is INFO.
    log_file : str, optional
        Path to write a comprehensive log file.
    """
    global _SETUP_DONE
    if _SETUP_DONE:
        return

    logger = logging.getLogger("hydrosheaf")
    logger.setLevel(logging.DEBUG)  # capture everything, handlers filter it
    logger.propagate = False

    # 1. Console Handler (User facing)
    console_handler = logging.StreamHandler(sys.stdout)
    console_level = logging.DEBUG if verbose else logging.INFO
    console_handler.setLevel(console_level)

    # Format: [Component] Message
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # 2. File Handler (Audit facing)
    if log_file:
        path = Path(log_file)
        path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file

        # Format: Time | Level | Module | Message
        file_formatter = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(name)-20s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    if capture_warnings:
        logging.captureWarnings(True)
        py_warnings = logging.getLogger("py.warnings")
        py_warnings.addHandler(console_handler)
        if log_file:
            py_warnings.addHandler(file_handler)

    _SETUP_DONE = True

    logger.info("Hydrosheaf Logging Initialized")
    if log_file:
        logger.info(f"Audit log writing to: {log_file}")


def _normalize_logger_name(name: str) -> str:
    clean = name.strip(".")
    if clean == "hydrosheaf" or clean.startswith("hydrosheaf."):
        return clean
    return f"hydrosheaf.{clean}"


def get_logger(name: str) -> ScienceLogger:
    """Get a logger instance under the hydrosheaf namespace."""
    return logging.getLogger(_normalize_logger_name(name))  # type: ignore


# --- Decorators for Transparency ---


def log_step(description: str):
    """
    Decorator to log the start and end of a significant scientific step.
    Includes timing.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger = get_logger(func.__module__)
            logger.info(f"START: {description}")
            t0 = time.time()
            try:
                result = func(*args, **kwargs)
                dt = time.time() - t0
                logger.info(f"DONE:  {description} ({dt:.3f}s)")
                return result
            except Exception as e:
                dt = time.time() - t0
                logger.error(f"FAIL:  {description} ({dt:.3f}s) - {str(e)}")
                raise

        return wrapper

    return decorator
