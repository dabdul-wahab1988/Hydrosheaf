"""Data parsing and transformation utilities."""
from datetime import datetime
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

import pandas as pd


def sample_list(
    samples: Union[Mapping[str, Any], Sequence[Any]]
) -> List[Mapping[str, Any]]:
    """Convert a sample input (dict or sequence) to a list of sample dicts."""
    if isinstance(samples, Mapping):
        return list(samples.values())
    if isinstance(samples, Sequence):
        return list(samples)
    raise TypeError("Unsupported samples input type.")


def _datetime_to_decimal_year(value: datetime) -> float:
    """Convert a datetime object to a decimal year float."""
    start = datetime(value.year, 1, 1, tzinfo=value.tzinfo)
    end = datetime(value.year + 1, 1, 1, tzinfo=value.tzinfo)
    span_seconds = (end - start).total_seconds()
    if span_seconds <= 0:
        return float(value.year)
    elapsed_seconds = (value - start).total_seconds()
    return float(value.year) + (elapsed_seconds / span_seconds)


def _numeric_to_decimal_year(value: float) -> Optional[float]:
    """Convert a numeric timestamp or year to decimal year."""
    if 1800.0 <= value <= 3000.0:
        return float(value)
    timestamp = value
    if timestamp > 1.0e12:  # milliseconds
        timestamp /= 1000.0
    if timestamp > 10000.0:  # seconds since epoch
        try:
            return _datetime_to_decimal_year(datetime.fromtimestamp(timestamp))
        except (OverflowError, OSError, ValueError):
            return None
    return None


def parse_decimal_year(value: object) -> Optional[float]:
    """Robustly parse a date/time value into a decimal year float.

    Uses pandas.to_datetime for flexible parsing of strings.
    """
    if value is None:
        return None
    if isinstance(value, datetime):
        return _datetime_to_decimal_year(value)
    if isinstance(value, (int, float)):
        return _numeric_to_decimal_year(float(value))

    text = str(value).strip()
    if not text:
        return None

    # Try numeric parsing first
    try:
        numeric = float(text)
        year = _numeric_to_decimal_year(numeric)
        if year is not None:
            return year
    except ValueError:
        pass

    # Use pandas for robust parsing
    try:
        dt = pd.to_datetime(text, errors='coerce')
        if not pd.isna(dt):
            # pandas Timestamp to datetime
            return _datetime_to_decimal_year(dt.to_pydatetime())
    except (ValueError, TypeError):
        pass

    return None


def extract_sample_decimal_year(sample: Mapping[str, object]) -> Optional[float]:
    """Extract decimal year from a sample dictionary using common date keys."""
    for key in (
        "sample_date",
        "collection_date",
        "date",
        "datetime",
        "timestamp",
        "sample_datetime",
        "collection_datetime",
        "year",
    ):
        if key in sample:
            parsed = parse_decimal_year(sample.get(key))
            if parsed is not None:
                return parsed
    return None
