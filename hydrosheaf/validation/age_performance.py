"""Held-out, abstention-aware performance accounting for the age specialist.

Prediction is deliberately kept in :mod:`age_competent_baseline`.  This
module starts at the truth-release boundary and therefore may consume held-
out ages, case labels, mechanism labels, and missingness labels.  It never
feeds those values back into a predictor or calibrator.

The report is intentionally conservative: a metric cannot pass its gate when
the denominator is missing, when development and held-out cases overlap, or
when the requested calibration record is absent.  Abstention is counted in
the unconditional denominator and is never allowed to masquerade as low
selective risk.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from collections.abc import Iterable, Mapping


ABSTAIN = "abstain"
SELECT = "select"


def _finite(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _fingerprint(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class AgePerformanceThresholds:
    """Predeclared thresholds for a bounded age-specialist claim."""

    min_conditional_coverage: float = 0.90
    min_coverage_including_abstention: float = 0.80
    max_mae_years: float = 15.0
    max_mean_interval_width_years: float = 80.0
    max_selective_risk_years: float = 12.0
    selective_risk_acceptance: float = 0.50
    max_false_commitment_rate: float = 0.10
    min_selected_rate: float = 0.50
    min_rows: int = 4
    min_rows_per_stratum: int = 2

    def __post_init__(self) -> None:
        bounded = (
            ("min_conditional_coverage", self.min_conditional_coverage, 0.0, 1.0, False),
            ("min_coverage_including_abstention", self.min_coverage_including_abstention, 0.0, 1.0, False),
            ("selective_risk_acceptance", self.selective_risk_acceptance, 0.0, 1.0, False),
            ("min_selected_rate", self.min_selected_rate, 0.0, 1.0, False),
            ("max_mae_years", self.max_mae_years, 0.0, float("inf"), True),
            ("max_mean_interval_width_years", self.max_mean_interval_width_years, 0.0, float("inf"), True),
            ("max_selective_risk_years", self.max_selective_risk_years, 0.0, float("inf"), True),
            ("max_false_commitment_rate", self.max_false_commitment_rate, 0.0, 1.0, True),
        )
        for name, raw, lower, upper, inclusive_lower in bounded:
            value = _finite(raw)
            if value is None or (value < lower if inclusive_lower else value <= lower) or value > upper:
                raise ValueError(f"{name} is outside its valid range.")
            object.__setattr__(self, name, float(value))
        if int(self.min_rows) < 1 or int(self.min_rows_per_stratum) < 1:
            raise ValueError("minimum row counts must be positive.")
        object.__setattr__(self, "min_rows", int(self.min_rows))
        object.__setattr__(self, "min_rows_per_stratum", int(self.min_rows_per_stratum))

    def to_dict(self) -> dict[str, object]:
        return {
            "min_conditional_coverage": self.min_conditional_coverage,
            "min_coverage_including_abstention": self.min_coverage_including_abstention,
            "max_mae_years": self.max_mae_years,
            "max_mean_interval_width_years": self.max_mean_interval_width_years,
            "max_selective_risk_years": self.max_selective_risk_years,
            "selective_risk_acceptance": self.selective_risk_acceptance,
            "max_false_commitment_rate": self.max_false_commitment_rate,
            "min_selected_rate": self.min_selected_rate,
            "min_rows": self.min_rows,
            "min_rows_per_stratum": self.min_rows_per_stratum,
        }


def _truth_age(value: object) -> float | None:
    if isinstance(value, Mapping):
        for key in ("age_years", "truth_age_years", "value"):
            if key in value:
                return _finite(value[key])
        return None
    return _finite(value)


def _label(metadata: Mapping[str, object], *keys: str) -> str:
    for key in keys:
        value = metadata.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return "unspecified"


def _row_metrics(rows: list[dict[str, object]]) -> dict[str, object]:
    truth_rows = rows
    valid = [row for row in rows if bool(row["valid"])]
    selected = [row for row in valid if bool(row["selected"])]
    covered = [row for row in selected if bool(row["covered"])]
    errors = [float(row["error"]) for row in selected]
    widths = [float(row["width"]) for row in selected]
    all_valid_widths = [float(row["width"]) for row in valid]
    return {
        "n_truth": len(truth_rows),
        "n_valid_outputs": len(valid),
        "n_selected": len(selected),
        "n_abstain": sum(1 for row in rows if row.get("abstain")),
        "n_missing": sum(1 for row in rows if row.get("missing")),
        "n_invalid": sum(1 for row in rows if row.get("invalid")),
        "selection_rate": len(selected) / len(truth_rows) if truth_rows else float("nan"),
        "mae_years": sum(errors) / len(errors) if errors else float("nan"),
        "conditional_interval_coverage": len(covered) / len(selected) if selected else float("nan"),
        "coverage_including_abstention": len(covered) / len(truth_rows) if truth_rows else float("nan"),
        "mean_interval_width_years": sum(widths) / len(widths) if widths else float("nan"),
        "mean_interval_width_all_valid_outputs_years": (
            sum(all_valid_widths) / len(all_valid_widths) if all_valid_widths else float("nan")
        ),
        "false_commitment_rate": 1.0 - len(covered) / len(selected) if selected else float("nan"),
    }


def _selective_curve(rows: list[dict[str, object]], thresholds: AgePerformanceThresholds) -> list[dict[str, object]]:
    valid = [row for row in rows if bool(row["valid"])]
    ordered = sorted(valid, key=lambda row: (float(row["uncertainty"]), str(row["target_id"])))
    if not ordered:
        return []
    points: list[dict[str, object]] = []
    for requested in (0.25, thresholds.selective_risk_acceptance, 0.75, 1.0):
        if any(math.isclose(requested, point["requested_acceptance"]) for point in points):
            continue
        count = max(1, int(math.ceil(requested * len(ordered))))
        accepted = ordered[:count]
        errors = [float(row["error"]) for row in accepted]
        covered = sum(bool(row["covered"]) for row in accepted)
        points.append(
            {
                "requested_acceptance": float(requested),
                "accepted": count,
                "total_valid_outputs": len(ordered),
                "total_truth": len(rows),
                "acceptance_rate": count / len(rows) if rows else float("nan"),
                "mean_absolute_error_years": sum(errors) / len(errors),
                "interval_coverage": covered / count,
                "false_commitment_rate": 1.0 - covered / count,
            }
        )
    return points


def _gate(name: str, value: object, requirement: str, passed: bool | None) -> dict[str, object]:
    return {
        "name": name,
        "value": value,
        "requirement": requirement,
        "status": "PASS" if passed is True else "FAIL" if passed is False else "ABSTAIN",
    }


def _stratify(rows: list[dict[str, object]], key: str, thresholds: AgePerformanceThresholds) -> dict[str, object]:
    groups: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(str(row[key]), []).append(row)
    report: dict[str, object] = {}
    for group, group_rows in sorted(groups.items()):
        metrics = _row_metrics(group_rows)
        report[group] = {
            "metrics": metrics,
            "status": "PASSABLE_DENOMINATOR" if len(group_rows) >= thresholds.min_rows_per_stratum else "INSUFFICIENT_DENOMINATOR",
        }
    return report


def evaluate_age_performance(
    truth: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]],
    *,
    metadata: Mapping[str, Mapping[str, object]] | None = None,
    thresholds: AgePerformanceThresholds | None = None,
    use_calibrated: bool = True,
    development_case_ids: Iterable[object] = (),
    held_out_case_ids: Iterable[object] = (),
) -> dict[str, object]:
    """Evaluate age predictions after truth release.

    ``metadata`` is scoring-only and may contain ``case_id``, ``family`` or
    ``generator_family``, ``mechanism`` and ``missingness``.  The function
    requires explicit development/held-out case sets for a PASS claim; this
    prevents a pooled score from being mistaken for held-out validation.
    """

    cfg = thresholds or AgePerformanceThresholds()
    meta = metadata or {}
    rows: list[dict[str, object]] = []
    calibration_statuses: list[str] = []
    for raw_target, raw_truth in truth.items():
        target_id = str(raw_target)
        target_meta = meta.get(target_id, {})
        if not isinstance(target_meta, Mapping):
            target_meta = {}
        age = _truth_age(raw_truth)
        family = _label(target_meta, "family", "generator_family", "case_family")
        mechanism = _label(target_meta, "mechanism", "age_mechanism", "mechanism_family")
        missingness = _label(target_meta, "missingness", "missingness_pattern")
        base = {
            "target_id": target_id,
            "truth": age,
            "family": family,
            "mechanism": mechanism,
            "missingness": missingness,
            "case_id": _label(target_meta, "case_id"),
            "valid": False,
            "selected": False,
            "covered": False,
            "abstain": False,
            "missing": False,
            "invalid": False,
        }
        prediction = predictions.get(target_id)
        if age is None:
            base["invalid"] = True
            rows.append(base)
            continue
        if not isinstance(prediction, Mapping):
            base["missing"] = True
            rows.append(base)
            continue
        decision = str(prediction.get("decision", SELECT))
        if use_calibrated:
            estimate_key, lower_key, upper_key = (
                "calibrated_mean_age_years",
                "calibrated_age_low",
                "calibrated_age_high",
            )
        else:
            estimate_key, lower_key, upper_key = ("mean_age_years", "age_95_low", "age_95_high")
        estimate = _finite(prediction.get(estimate_key))
        lower = _finite(prediction.get(lower_key))
        upper = _finite(prediction.get(upper_key))
        # Calibrators intentionally leave abstained rows uncalibrated.  They
        # still count as abstentions in the unconditional denominator rather
        # than becoming invalid/missing rows.
        if decision == ABSTAIN and (estimate is None or lower is None or upper is None):
            estimate = _finite(prediction.get("mean_age_years"))
            lower = _finite(prediction.get("age_95_low"))
            upper = _finite(prediction.get("age_95_high"))
        if estimate is None or lower is None or upper is None or lower > upper or not lower <= estimate <= upper:
            if decision == ABSTAIN:
                base["abstain"] = True
                rows.append(base)
                continue
            base["invalid"] = True
            rows.append(base)
            continue
        if use_calibrated and decision != ABSTAIN:
            status = str(prediction.get("calibration_status", ""))
            calibration_statuses.append(status or "missing")
        uncertainty = _finite(prediction.get("uncertainty_years"))
        if uncertainty is None:
            uncertainty = max(upper - lower, 0.0) / 2.0
        base.update(
            {
                "valid": True,
                "selected": decision != ABSTAIN,
                "abstain": decision == ABSTAIN,
                "estimate": estimate,
                "lower": lower,
                "upper": upper,
                "width": upper - lower,
                "uncertainty": uncertainty,
                "error": abs(estimate - age),
                "covered": lower <= age <= upper,
            }
        )
        rows.append(base)

    overall = _row_metrics(rows)
    curve = _selective_curve(rows, cfg)
    selected_rate = _finite(overall["selection_rate"])
    selective_point = next(
        (point for point in curve if math.isclose(point["requested_acceptance"], cfg.selective_risk_acceptance)),
        None,
    )
    overall["selective_risk_at_acceptance_years"] = (
        None if selective_point is None else selective_point["mean_absolute_error_years"]
    )
    overall["selective_false_commitment_at_acceptance"] = (
        None if selective_point is None else selective_point["false_commitment_rate"]
    )
    development = {str(value) for value in development_case_ids}
    held_out = {str(value) for value in held_out_case_ids}
    observed_cases = {str(row["case_id"]) for row in rows if str(row["case_id"]) != "unspecified"}
    if not development or not held_out:
        held_out_gate = _gate(
            "held_out_case_disjoint",
            {"development_cases": sorted(development), "held_out_cases": sorted(held_out)},
            "explicit non-empty disjoint development and held-out case sets",
            None,
        )
    else:
        overlap = sorted(development & held_out)
        missing_labels = sorted(observed_cases - held_out)
        held_out_gate = _gate(
            "held_out_case_disjoint",
            {"overlap": overlap, "observed_cases": sorted(observed_cases), "outside_held_out": missing_labels},
            "overlap is empty and every labelled scoring case is held out",
            not overlap and not missing_labels,
        )

    if use_calibrated:
        calibration_gate = _gate(
            "development_fitted_calibration",
            sorted(set(calibration_statuses)),
            "all comparable selected outputs carry development_fitted calibration metadata",
            bool(calibration_statuses) and set(calibration_statuses) == {"development_fitted"},
        )
    else:
        calibration_gate = _gate(
            "development_fitted_calibration",
            "not_requested",
            "calibration gate not requested for raw diagnostic scoring",
            None,
        )

    gates = [
        held_out_gate,
        calibration_gate,
        _gate("mae", overall["mae_years"], f"<= {cfg.max_mae_years} years with n >= {cfg.min_rows}",
              None if int(overall["n_selected"]) < cfg.min_rows else _finite(overall["mae_years"]) is not None and float(overall["mae_years"]) <= cfg.max_mae_years),
        _gate("conditional_interval_coverage", overall["conditional_interval_coverage"], f">= {cfg.min_conditional_coverage} with n >= {cfg.min_rows}",
              None if int(overall["n_selected"]) < cfg.min_rows else _finite(overall["conditional_interval_coverage"]) is not None and float(overall["conditional_interval_coverage"]) >= cfg.min_conditional_coverage),
        _gate("coverage_including_abstention", overall["coverage_including_abstention"], f">= {cfg.min_coverage_including_abstention}",
              None if int(overall["n_truth"]) < cfg.min_rows else _finite(overall["coverage_including_abstention"]) is not None and float(overall["coverage_including_abstention"]) >= cfg.min_coverage_including_abstention),
        _gate("mean_interval_width", overall["mean_interval_width_years"], f"<= {cfg.max_mean_interval_width_years} years with n >= {cfg.min_rows}",
              None if int(overall["n_selected"]) < cfg.min_rows else _finite(overall["mean_interval_width_years"]) is not None and float(overall["mean_interval_width_years"]) <= cfg.max_mean_interval_width_years),
        _gate("selective_risk", None if selective_point is None else selective_point["mean_absolute_error_years"], f"<= {cfg.max_selective_risk_years} years at acceptance {cfg.selective_risk_acceptance}",
              None if selective_point is None else float(selective_point["mean_absolute_error_years"]) <= cfg.max_selective_risk_years),
        _gate("false_commitment", overall["false_commitment_rate"], f"<= {cfg.max_false_commitment_rate} on selected outputs",
              None if int(overall["n_selected"]) < cfg.min_rows else _finite(overall["false_commitment_rate"]) is not None and float(overall["false_commitment_rate"]) <= cfg.max_false_commitment_rate),
        _gate("nontrivial_selection", selected_rate, f">= {cfg.min_selected_rate}",
              None if selected_rate is None else selected_rate >= cfg.min_selected_rate),
    ]

    stratified = {
        "family": _stratify(rows, "family", cfg),
        "mechanism": _stratify(rows, "mechanism", cfg),
        "missingness": _stratify(rows, "missingness", cfg),
    }
    stratification_gate = _gate(
        "stratified_denominators",
        {
            dimension: {group: record["metrics"]["n_truth"] for group, record in groups.items()}
            for dimension, groups in stratified.items()
        },
        f"each declared non-unspecified stratum has n >= {cfg.min_rows_per_stratum}",
        all(
            int(record["metrics"]["n_truth"]) >= cfg.min_rows_per_stratum
            for groups in stratified.values()
            for group, record in groups.items()
            if group != "unspecified"
        ) and any(group != "unspecified" for groups in stratified.values() for group in groups),
    )
    gates.append(stratification_gate)
    declared_groups = [
        (dimension, group, record["metrics"])
        for dimension, groups in stratified.items()
        for group, record in groups.items()
        if group != "unspecified"
    ]
    stratified_failures = [
        f"{dimension}:{group}"
        for dimension, group, metrics in declared_groups
        if int(metrics["n_truth"]) >= cfg.min_rows_per_stratum
        and (
            int(metrics["n_selected"]) < cfg.min_rows_per_stratum
            or _finite(metrics["mae_years"]) is None
            or float(metrics["mae_years"]) > cfg.max_mae_years
            or _finite(metrics["conditional_interval_coverage"]) is None
            or float(metrics["conditional_interval_coverage"]) < cfg.min_conditional_coverage
            or _finite(metrics["coverage_including_abstention"]) is None
            or float(metrics["coverage_including_abstention"]) < cfg.min_coverage_including_abstention
            or _finite(metrics["false_commitment_rate"]) is None
            or float(metrics["false_commitment_rate"]) > cfg.max_false_commitment_rate
        )
    ]
    gates.append(
        _gate(
            "stratified_performance",
            {"failures": stratified_failures},
            "every adequately sized declared family/mechanism/missingness stratum meets age-risk and coverage thresholds",
            None if any(
                int(metrics["n_truth"]) < cfg.min_rows_per_stratum
                or int(metrics["n_selected"]) < cfg.min_rows_per_stratum
                for _, _, metrics in declared_groups
            ) else not stratified_failures,
        )
    )
    statuses = [str(gate["status"]) for gate in gates]
    claim_status = "FAIL" if "FAIL" in statuses else "ABSTAIN" if "ABSTAIN" in statuses else "PASS"
    return {
        "status": "scored" if rows else "no_comparable_outputs",
        "claim_status": claim_status,
        "metrics": overall,
        "selective_risk_curve": curve,
        "stratified": stratified,
        "gates": gates,
        "audit": {
            "truth_release_boundary": True,
            "calibration_requested": use_calibrated,
            "calibration_statuses_observed": sorted(set(calibration_statuses)),
            "thresholds": cfg.to_dict(),
            "prediction_keys_fingerprint": _fingerprint(sorted(predictions)),
        },
    }


__all__ = ["AgePerformanceThresholds", "evaluate_age_performance"]
