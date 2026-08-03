"""Truth-blind handling of repeated or model-derived head channels.

Groundwater benchmark rows often contain both ``head_meas`` and
``hydraulic_head``.  Those fields are not automatically two independent
pieces of evidence: they may measure the same latent head, or one may be a
model-derived value with structured bias.  This module makes that choice
explicit before an inference run and records an auditable consolidation
summary.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping


DEFAULT_HEAD_EVIDENCE_MODEL: Mapping[str, object] = {
    "primary_channel": "head_meas",
    "secondary_channel": "hydraulic_head",
    "primary_sigma_m": 0.10,
    "secondary_sigma_m": 0.10,
    "measurement_error_correlation": 0.0,
    "combination": "gls",
}


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _positive_model_float(model: Mapping[str, object], key: str, default: float) -> float:
    value = _finite(model.get(key))
    if value is None or value <= 0.0:
        return float(default)
    return float(value)


def _model_with_aliases(model: Mapping[str, object] | None) -> dict[str, object]:
    merged = dict(DEFAULT_HEAD_EVIDENCE_MODEL)
    if model:
        merged.update(model)
    primary = str(merged["primary_channel"])
    secondary = str(merged["secondary_channel"])
    merged["primary_channel"] = primary
    merged["secondary_channel"] = secondary
    merged["primary_sigma_m"] = _positive_model_float(
        merged,
        "primary_sigma_m",
        _positive_model_float(merged, f"{primary}_sigma_m", 0.10),
    )
    merged["secondary_sigma_m"] = _positive_model_float(
        merged,
        "secondary_sigma_m",
        _positive_model_float(merged, f"{secondary}_sigma_m", 0.10),
    )
    correlation = _finite(
        merged.get(
            "measurement_error_correlation",
            merged.get("correlation", 0.0),
        )
    )
    if correlation is None:
        correlation = 0.0
    merged["measurement_error_correlation"] = max(-0.99, min(0.99, correlation))
    combination = str(merged.get("combination", "gls"))
    if combination not in {"gls", "primary_only_with_discrepancy"}:
        raise ValueError(
            "Head evidence combination must be 'gls' or "
            "'primary_only_with_discrepancy'."
        )
    merged["combination"] = combination
    return merged


def _gls_pair(
    primary: float,
    secondary: float,
    sigma_primary: float,
    sigma_secondary: float,
    correlation: float,
) -> tuple[float, float]:
    """Return the latent-head GLS estimate and its standard deviation."""

    covariance = correlation * sigma_primary * sigma_secondary
    determinant = sigma_primary**2 * sigma_secondary**2 - covariance**2
    if determinant <= 0.0:
        raise ValueError("Head covariance matrix must be positive definite.")
    inv00 = sigma_secondary**2 / determinant
    inv11 = sigma_primary**2 / determinant
    inv01 = -covariance / determinant
    denominator = inv00 + inv11 + 2.0 * inv01
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise ValueError("Head covariance matrix gives no finite GLS precision.")
    estimate = (inv00 * primary + inv11 * secondary + inv01 * (primary + secondary)) / denominator
    sigma = math.sqrt(1.0 / denominator)
    return float(estimate), float(sigma)


def consolidate_head_evidence(
    observations: Iterable[Mapping[str, object]],
    *,
    model: Mapping[str, object] | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Consolidate correlated head channels without using generator truth.

    ``gls`` treats the two channels as repeated measurements of one latent
    head with the declared covariance.  ``primary_only_with_discrepancy``
    retains the primary channel and records the secondary-minus-primary
    residual without feeding that biased channel as independent evidence.
    The secondary raw value is removed from the returned inference rows in
    both modes, so downstream modules cannot count it a second time.
    """

    resolved = _model_with_aliases(model)
    primary_key = str(resolved["primary_channel"])
    secondary_key = str(resolved["secondary_channel"])
    sigma_primary = float(resolved["primary_sigma_m"])
    sigma_secondary = float(resolved["secondary_sigma_m"])
    correlation = float(resolved["measurement_error_correlation"])
    combination = str(resolved["combination"])

    output: list[dict[str, object]] = []
    combined_count = 0
    fallback_count = 0
    missing_count = 0
    discrepancy_values: list[float] = []
    effective_sigmas: list[float] = []

    for source_row in observations:
        row = dict(source_row)
        primary = _finite(row.get(primary_key))
        secondary = _finite(row.get(secondary_key))
        source_count = int(primary is not None) + int(secondary is not None)
        estimate: float | None = None
        estimate_sigma: float | None = None
        discrepancy: float | None = None

        if primary is not None and secondary is not None:
            if combination == "gls":
                estimate, estimate_sigma = _gls_pair(
                    primary,
                    secondary,
                    sigma_primary,
                    sigma_secondary,
                    correlation,
                )
                combined_count += 1
            else:
                estimate = primary
                estimate_sigma = sigma_primary
                discrepancy = float(secondary - primary)
                discrepancy_values.append(discrepancy)
        elif primary is not None:
            estimate = primary
            estimate_sigma = sigma_primary
            fallback_count += 1
        elif secondary is not None:
            estimate = secondary
            estimate_sigma = sigma_secondary
            fallback_count += 1
        else:
            missing_count += 1

        if estimate is not None:
            row[primary_key] = estimate
            row["head_evidence_sigma_m"] = estimate_sigma
            if estimate_sigma is not None:
                effective_sigmas.append(float(estimate_sigma))
        row.pop(secondary_key, None)
        row["head_evidence_model"] = combination
        row["head_evidence_source_count"] = source_count
        row["head_evidence_primary_channel"] = primary_key
        row["head_evidence_secondary_channel"] = secondary_key
        row["head_evidence_error_correlation"] = correlation
        if discrepancy is not None:
            row["head_evidence_secondary_minus_primary_m"] = discrepancy
        output.append(row)

    audit = {
        "model": resolved,
        "n_rows": len(output),
        "n_two_channel_rows": combined_count,
        "n_single_channel_rows": fallback_count,
        "n_missing_head_rows": missing_count,
        "n_discrepancy_rows": len(discrepancy_values),
        "effective_sigma_m": (
            sum(effective_sigmas) / len(effective_sigmas)
            if effective_sigmas
            else None
        ),
        "max_abs_secondary_minus_primary_m": max(
            (abs(value) for value in discrepancy_values),
            default=0.0,
        ),
        "secondary_channel_consumed_as_independent_evidence": False,
    }
    return output, audit


__all__ = [
    "DEFAULT_HEAD_EVIDENCE_MODEL",
    "consolidate_head_evidence",
]
