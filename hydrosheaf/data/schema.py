"""Input schema helpers."""

from typing import Iterable, List, Mapping, Optional, Sequence, Tuple


def required_columns(ion_order: Iterable[str]) -> List[str]:
    return ["sample_id", "site_id", *ion_order, "EC", "TDS", "pH"]


def _columns_from_samples(samples: object) -> Iterable[str]:
    if hasattr(samples, "columns"):
        return samples.columns
    if isinstance(samples, Sequence) and samples:
        first = samples[0]
        if isinstance(first, Mapping):
            return first.keys()
    raise TypeError("Unsupported samples input type.")


def validate_required_columns(samples: object, ion_order: Iterable[str]) -> None:
    columns = set(_columns_from_samples(samples))
    missing = [col for col in required_columns(ion_order) if col not in columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def missing_required(sample: Mapping[str, object], ion_order: Iterable[str]) -> List[str]:
    required = required_columns(ion_order)
    missing = [col for col in required if col not in sample or sample[col] in ("", None)]
    return missing


def missing_ions(sample: Mapping[str, object], ion_order: Iterable[str]) -> List[str]:
    missing = [ion for ion in ion_order if ion not in sample or sample[ion] in ("", None)]
    return missing


def extract_vector(sample: Mapping[str, float], ion_order: Iterable[str]) -> list:
    return [float(sample[ion]) for ion in ion_order]


def safe_extract_vector(sample: Mapping[str, float], ion_order: Iterable[str]) -> list:
    if missing_required(sample, ion_order):
        raise ValueError("Sample is missing required ion values.")
    return extract_vector(sample, ion_order)


def parse_numeric(value: object, detection_policy: str) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    if text.startswith("<") or text.startswith(">"):
        prefix = text[0]
        number = text[1:].strip()
        try:
            limit = float(number)
        except ValueError:
            return None
        if detection_policy == "drop":
            return None
        if detection_policy == "zero":
            return 0.0 if prefix == "<" else limit
        if detection_policy == "half":
            return 0.5 * limit if prefix == "<" else limit
        return limit
    try:
        return float(text)
    except ValueError:
        return None


def normalize_sample(
    sample: Mapping[str, object],
    ion_order: Iterable[str],
    detection_policy: str,
) -> dict:
    numeric_keys = set(ion_order) | {"EC", "TDS", "pH", "K"}
    normalized: dict = {}
    for key, value in sample.items():
        if key in numeric_keys:
            normalized[key] = parse_numeric(value, detection_policy)
        else:
            normalized[key] = value
    return normalized


def normalize_samples(
    samples: Sequence[Mapping[str, object]],
    ion_order: Iterable[str],
    detection_policy: str,
) -> List[dict]:
    return [normalize_sample(sample, ion_order, detection_policy) for sample in samples]


def vector_from_sample(
    sample: Mapping[str, object],
    ion_order: Iterable[str],
    missing_policy: str,
    detection_policy: str,
) -> Tuple[Optional[List[float]], dict]:
    normalized = normalize_sample(sample, ion_order, detection_policy)
    values: List[float] = []
    for ion in ion_order:
        value = normalized.get(ion)
        if value is None:
            if missing_policy == "impute_zero":
                values.append(0.0)
            else:
                return None, normalized
        else:
            values.append(float(value))
    return values, normalized


def build_endmember_vectors(
    samples: Sequence[Mapping[str, float]],
    ids: Iterable[str],
    ion_order: Iterable[str],
) -> Mapping[str, List[float]]:
    lookup = {str(sample.get("site_id")): sample for sample in samples}
    endmembers = {}
    for end_id in ids:
        sample = lookup.get(str(end_id))
        if sample is None or missing_ions(sample, ion_order):
            continue
        endmembers[str(end_id)] = extract_vector(sample, ion_order)
    return endmembers
