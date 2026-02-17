"""Data validation and requirement checks."""
from dataclasses import replace
from typing import Dict, List, Mapping, Sequence

from ..config import Config
from .parsing import sample_list
from .schema import parse_numeric


def _any_numeric(
    samples: Sequence[Mapping[str, object]],
    key: str,
    detection_policy: str,
) -> bool:
    """Check if any sample has a valid numeric value for the given key."""
    for sample in samples:
        if parse_numeric(sample.get(key), detection_policy) is not None:
            return True
    return False


def _any_pair_numeric(
    samples: Sequence[Mapping[str, object]],
    key_a: str,
    key_b: str,
    detection_policy: str,
) -> bool:
    """Check if any sample has valid numeric values for both keys."""
    for sample in samples:
        if (
            parse_numeric(sample.get(key_a), detection_policy) is not None
            and parse_numeric(sample.get(key_b), detection_policy) is not None
        ):
            return True
    return False


def _missing_keys(
    samples: Sequence[Mapping[str, object]],
    keys: Sequence[str],
    detection_policy: str,
) -> List[str]:
    """Identify sample IDs missing required numeric keys."""
    missing_ids: List[str] = []
    for sample in samples:
        missing = False
        for key in keys:
            if parse_numeric(sample.get(key), detection_policy) is None:
                missing = True
                break
        if missing:
            sample_id = sample.get("site_id") or sample.get("sample_id") or "unknown"
            missing_ids.append(str(sample_id))
    return missing_ids


def validate_required_inputs(samples: object, config: Config) -> None:
    """Raise ValueError if required inputs for enabled modules are missing."""
    s_list = sample_list(samples)
    detection_policy = config.detection_limit_policy
    missing_reports: List[str] = []

    if config.phreeqc_enabled:
        required_phreeqc = [
            "pH",
            "Ca",
            "Mg",
            "Na",
            "K",
            "Cl",
            "SO4",
            "NO3",
            "F",
            "HCO3",
        ]
        missing = _missing_keys(s_list, required_phreeqc, detection_policy)
        if missing:
            missing_reports.append(
                "PHREEQC requires pH and major ions (Ca, Mg, Na, K, Cl, SO4, NO3, F, HCO3) "
                f"for all samples (missing: {missing})"
            )

    if config.isotope_enabled and config.lmwl_defined:
        missing = _missing_keys(
            s_list,
            [config.isotope_d18o_key, config.isotope_d2h_key],
            detection_policy,
        )
        if missing:
            missing_reports.append(
                "Isotope penalties require both "
                f"{config.isotope_d18o_key} and {config.isotope_d2h_key} "
                f"for all samples (missing: {missing})"
            )

    if config.nitrate_source_enabled:
        missing = _missing_keys(s_list, ["NO3"], detection_policy)
        if missing:
            missing_reports.append(
                f"Nitrate source requires NO3 for all samples (missing: {missing})"
            )

    if missing_reports:
        raise ValueError("; ".join(missing_reports))


def auto_disable_missing_modules(samples: object, config: Config) -> Config:
    """Disable feature flags when required inputs are missing across samples."""
    s_list = sample_list(samples)
    detection_policy = config.detection_limit_policy
    updates: Dict[str, object] = {}

    if config.phreeqc_enabled and not _any_numeric(s_list, "pH", detection_policy):
        updates["phreeqc_enabled"] = False

    if config.isotope_enabled and not _any_pair_numeric(
        s_list,
        config.isotope_d18o_key,
        config.isotope_d2h_key,
        detection_policy,
    ):
        updates["isotope_enabled"] = False

    if config.nitrate_source_enabled and not _any_numeric(
        s_list, "NO3", detection_policy
    ):
        updates["nitrate_source_enabled"] = False

    if updates:
        return replace(config, **updates)
    return config
