"""
Centralized logging configuration for Hydrosheaf backend.
"""

import os
import sys
import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

# Log level from environment variable
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
LOG_FILE = os.getenv("LOG_FILE", "hydrosheaf.log")
LOG_MAX_BYTES = int(os.getenv("LOG_MAX_BYTES", str(10 * 1024 * 1024)))  # 10MB default
LOG_BACKUP_COUNT = int(os.getenv("LOG_BACKUP_COUNT", "5"))

# Create logs directory if it doesn't exist
LOGS_DIR = Path(__file__).parent.parent / "logs"
LOGS_DIR.mkdir(exist_ok=True)

# Formatter for consistent log output
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def get_logger(name: str) -> logging.Logger:
    """
    Get a configured logger instance.

    Args:
        name: Logger name (typically __name__ of the module)

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)

    # Only configure if not already configured
    if not logger.handlers:
        logger.setLevel(getattr(logging, LOG_LEVEL, logging.INFO))

        formatter = logging.Formatter(LOG_FORMAT, datefmt=DATE_FORMAT)

        # Console handler (stdout)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # File handler with rotation
        try:
            file_handler = RotatingFileHandler(
                LOGS_DIR / LOG_FILE,
                maxBytes=LOG_MAX_BYTES,
                backupCount=LOG_BACKUP_COUNT,
            )
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        except Exception as e:
            logger.warning(f"Could not create file handler: {e}")

        # Prevent propagation to root logger
        logger.propagate = False

    return logger


# Pre-configured loggers for common modules
analysis_logger = get_logger("hydrosheaf.analysis")
network_logger = get_logger("hydrosheaf.network")
samples_logger = get_logger("hydrosheaf.samples")
projects_logger = get_logger("hydrosheaf.projects")
