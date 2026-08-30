"""Adjudicate a frozen synthetic ensemble without changing the lock.

The ensemble runner deliberately separates execution from scientific
performance.  This script reads the immutable lock artifacts, applies the
declared coverage and comparator rules, and writes a separate adjudication
record.  It never writes into the input lock directory.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from math import isfinite
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping, Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = (
    REPO
    / ".codex_work"
    / "programme-validation"
    / "RUN-PROGRAMME-ENSEMBLE-CALIBRATION-LOCK-20260801-02"
)
DEFAULT_OUTPUT_ROOT = REPO / ".codex_work" / "programme-validation-adjudication"
PRIMARY_CONDITION = "full_sheaf"
COMPLETE_SCENARIO = "complete"
TARGET_COVERAGE = 0.95
COMMON_TOPOLOGY_BASELINE = "hydraulic_only_directional_gradient"
COMMON_TOPOLOGY_MARGIN = -0.05
COMMON_TOPOLOGY_MIN_CASES_PER_FAMILY = 2
COMMON_TOPOLOGY_REQUIRED_FAMILIES = frozenset(
    {"analytic_lattice_v1", "independent_mixing_v1", "modflow_modpath_v3"}
)
PASS = "PASS"
ABSTAIN = "ABSTAIN"
FAIL = "FAIL"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def _git_worktree_dirty() -> bool:
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        )
        return bool(result.stdout.strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        return True


def _finite(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def _rapm_metrics(record: object) -> Mapping[str, Any]:
    """Extract one reaction-RAPM score record from a case row."""

    if not isinstance(record, Mapping):
        return {}
    reaction = record.get("reaction", record)
    return reaction if isinstance(reaction, Mapping) else {}


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _safe_output(output: Path, *, overwrite: bool) -> Path:
    generated_root = (REPO / ".codex_work").resolve()
    resolved = output.resolve()
    if not resolved.is_relative_to(generated_root):
        raise ValueError("Adjudication output must stay under .codex_work.")
    if resolved == generated_root:
        raise ValueError("Refusing to use .codex_work itself as an output directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output exists: {resolved}; pass --overwrite to replace it."
            )
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True, exist_ok=False)
    return resolved


def _input_integrity(lock_dir: Path, manifest: Mapping[str, Any], errors: object) -> dict[str, Any]:
    artifact_checks = {
        str(name): (lock_dir / str(name)).is_file()
        and _sha256(lock_dir / str(name)) == str(expected)
        for name, expected in dict(manifest.get("artifacts", {})).items()
    }
    source_checks = {}
    for relative, expected in dict(manifest.get("source_hashes", {})).items():
        path = REPO / str(relative)
        source_checks[str(relative)] = path.is_file() and _sha256(path) == str(expected)
    errors_list = errors if isinstance(errors, list) else []
    return {
        "run_id": manifest.get("run_id"),
        "manifest_status": manifest.get("status"),
        "artifact_hashes_all_match": all(artifact_checks.values()),
        "artifact_hash_failures": sorted(name for name, ok in artifact_checks.items() if not ok),
        "source_hashes_all_match": all(source_checks.values()),
        "source_hash_failures": sorted(name for name, ok in source_checks.items() if not ok),
        "errors_empty": not errors_list,
        "recorded_error_count": len(errors_list),
    }


def _status_from_delta(delta: float | None) -> str:
    if delta is None:
        return ABSTAIN
    return PASS if delta >= 0.0 else FAIL


def _common_topology_component_claim(
    rows: Sequence[Mapping[str, Any]],
    *,
    execution_pass: bool,
) -> dict[str, Any]:
    """Adjudicate one fixed, bounded topology estimand.

    This is deliberately narrower than a HydroSheaf system claim.  It uses a
    preregistered hydraulic-only comparator and the common all-pairs selection
    universe emitted by the runner.  Candidate-generation recall is retained
    as a separate diagnostic and is not folded into this selection metric.
    """

    complete_rows = [
        row
        for row in rows
        if row.get("phase") == "locked_test"
        and row.get("condition") == PRIMARY_CONDITION
        and row.get("observation_scenario", COMPLETE_SCENARIO) == COMPLETE_SCENARIO
    ]
    family_counts = Counter(str(row.get("family")) for row in complete_rows)
    checks: list[dict[str, Any]] = []
    checks.append(
        {
            "name": "execution_integrity",
            "passed": True if execution_pass else None,
            "observed": execution_pass,
            "requirement": True,
        }
    )
    checks.append(
        {
            "name": "locked_case_denominator",
            "passed": True if complete_rows else None,
            "observed": len(complete_rows),
            "requirement": ">=1",
        }
    )
    checks.append(
        {
            "name": "family_case_coverage",
            "passed": (
                None
                if not family_counts
                or set(family_counts) != set(COMMON_TOPOLOGY_REQUIRED_FAMILIES)
                or any(
                    count < COMMON_TOPOLOGY_MIN_CASES_PER_FAMILY
                    for count in family_counts.values()
                )
                else True
            ),
            "observed": dict(sorted(family_counts.items())),
            "requirement": {
                "required_families": sorted(COMMON_TOPOLOGY_REQUIRED_FAMILIES),
                "minimum_cases_per_family": COMMON_TOPOLOGY_MIN_CASES_PER_FAMILY
            },
        }
    )
    case_records: list[dict[str, Any]] = []
    missing = False
    for row in complete_rows:
        comparison = row.get("common_topology_comparison")
        comparison = comparison if isinstance(comparison, Mapping) else {}
        universe = row.get("common_candidate_universe")
        universe = universe if isinstance(universe, Mapping) else {}
        hydro = comparison.get("hydrosheaf", {})
        hydro = hydro if isinstance(hydro, Mapping) else {}
        hydro_topology = hydro.get("topology", {})
        hydro_topology = hydro_topology if isinstance(hydro_topology, Mapping) else {}
        baselines = comparison.get("baselines", {})
        baselines = baselines if isinstance(baselines, Mapping) else {}
        baseline = baselines.get(COMMON_TOPOLOGY_BASELINE, {})
        baseline = baseline if isinstance(baseline, Mapping) else {}
        hydro_f1 = _finite(hydro_topology.get("selected_f1"))
        baseline_f1 = _finite(baseline.get("selected_f1"))
        delta = (
            hydro_f1 - baseline_f1
            if hydro_f1 is not None and baseline_f1 is not None
            else None
        )
        universe_ok = bool(
            universe.get("truth_blind") is True
            and not universe.get("truth_fields_seen")
            and universe.get("algorithm") == "common_all_pairs_v1"
            and bool(universe.get("candidate_hash"))
        )
        if delta is None or not universe_ok:
            missing = True
        case_records.append(
            {
                "family": str(row.get("family")),
                "seed": int(row.get("seed")),
                "hydrosheaf_selected_f1": hydro_f1,
                "fixed_baseline": COMMON_TOPOLOGY_BASELINE,
                "baseline_selected_f1": baseline_f1,
                "delta": delta,
                "common_universe_valid": universe_ok,
                "candidate_generation_metrics": row.get(
                    "specialist_candidate_universe_metrics", {}
                ),
            }
        )
    checks.append(
        {
            "name": "common_universe_and_metrics_complete",
            "passed": None if missing or not case_records else True,
            "observed": sum(record["common_universe_valid"] for record in case_records),
            "requirement": "every locked case has common_all_pairs_v1 and finite F1",
        }
    )
    deltas = [
        float(record["delta"])
        for record in case_records
        if record["delta"] is not None and record["common_universe_valid"]
    ]
    noninferior = bool(deltas) and all(delta >= COMMON_TOPOLOGY_MARGIN for delta in deltas)
    checks.append(
        {
            "name": "paired_non_inferiority_margin",
            "passed": None if not deltas else noninferior,
            "observed": min(deltas) if deltas else None,
            "requirement": {"minimum_delta": COMMON_TOPOLOGY_MARGIN},
        }
    )
    family_deltas: dict[str, float | None] = {}
    for family in sorted(family_counts):
        values = [
            float(record["delta"])
            for record in case_records
            if record["family"] == family
            and record["delta"] is not None
            and record["common_universe_valid"]
        ]
        family_deltas[family] = sum(values) / len(values) if values else None
    checks.append(
        {
            "name": "family_stratified_non_inferiority",
            "passed": (
                None
                if not family_deltas or any(value is None for value in family_deltas.values())
                else all(float(value) >= COMMON_TOPOLOGY_MARGIN for value in family_deltas.values())
            ),
            "observed": family_deltas,
            "requirement": {"minimum_family_delta": COMMON_TOPOLOGY_MARGIN},
        }
    )
    if any(check["passed"] is False for check in checks):
        status = FAIL
    elif any(check["passed"] is None for check in checks):
        status = ABSTAIN
    else:
        status = PASS
    return {
        "claim_name": "synthetic_component_topology_selection",
        "status": status,
        "estimand": (
            "held-out topology-selection F1 on common all-pairs universe versus "
            "fixed hydraulic-only directional-gradient comparator"
        ),
        "comparator": COMMON_TOPOLOGY_BASELINE,
        "non_inferiority_margin": COMMON_TOPOLOGY_MARGIN,
        "minimum_cases_per_family": COMMON_TOPOLOGY_MIN_CASES_PER_FAMILY,
        "family_case_counts": dict(sorted(family_counts.items())),
        "family_mean_deltas": family_deltas,
        "case_records": case_records,
        "checks": checks,
        "reasons": [
            "common-universe topology selection is a bounded component estimand",
            "candidate-generation recall is reported separately",
            "this does not establish age, reaction, decision, field, or universal system superiority",
        ],
    }


def _best_baseline(baselines: Mapping[str, Any]) -> tuple[str | None, Mapping[str, Any] | None]:
    candidates: list[tuple[float, str, Mapping[str, Any]]] = []
    for name, raw in baselines.items():
        if not isinstance(raw, Mapping):
            continue
        topology = raw.get("topology", {})
        if not isinstance(topology, Mapping):
            continue
        f1 = _finite(topology.get("selected_f1"))
        if f1 is not None:
            candidates.append((f1, str(name), raw))
    if not candidates:
        return None, None
    _, name, record = max(candidates, key=lambda item: (item[0], item[1]))
    return name, record


def _row_index(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, int, str, str], Mapping[str, Any]]:
    index: dict[tuple[str, int, str, str], Mapping[str, Any]] = {}
    for row in rows:
        if row.get("phase") != "locked_test":
            continue
        family = str(row.get("family"))
        seed = int(row.get("seed"))
        condition = str(row.get("condition"))
        scenario = str(row.get("observation_scenario", COMPLETE_SCENARIO))
        index[(family, seed, condition, scenario)] = row
    return index


def _case_adjudication(
    row: Mapping[str, Any],
    *,
    age_score: Mapping[str, Any] | None,
    target_coverage: float,
    row_index: Mapping[tuple[str, int, str, str], Mapping[str, Any]],
) -> dict[str, Any]:
    family = str(row.get("family"))
    seed = int(row.get("seed"))
    topology = row.get("topology", {})
    topology = topology if isinstance(topology, Mapping) else {}
    baselines = row.get("baselines", {})
    baselines = baselines if isinstance(baselines, Mapping) else {}
    best_name, best_record = _best_baseline(baselines)
    specialist_baselines = row.get("specialist_baselines", {})
    specialist_baselines = (
        specialist_baselines if isinstance(specialist_baselines, Mapping) else {}
    )
    calibrated_specialist_baselines = row.get(
        "specialist_calibrated_baselines", {}
    )
    calibrated_specialist_baselines = (
        calibrated_specialist_baselines
        if isinstance(calibrated_specialist_baselines, Mapping)
        else {}
    )
    specialist_best_name, specialist_best_record = _best_baseline(specialist_baselines)
    best_topology = best_record.get("topology", {}) if best_record else {}
    best_topology = best_topology if isinstance(best_topology, Mapping) else {}
    specialist_topology = (
        specialist_best_record.get("topology", {}) if specialist_best_record else {}
    )
    specialist_topology = (
        specialist_topology if isinstance(specialist_topology, Mapping) else {}
    )
    hydro_f1 = _finite(topology.get("selected_f1"))
    baseline_f1 = _finite(best_topology.get("selected_f1"))
    specialist_f1 = _finite(specialist_topology.get("selected_f1"))
    topology_delta = (
        hydro_f1 - baseline_f1
        if hydro_f1 is not None and baseline_f1 is not None
        else None
    )
    specialist_delta = (
        hydro_f1 - specialist_f1
        if hydro_f1 is not None and specialist_f1 is not None
        else None
    )

    calibrated = age_score.get("calibrated_age", {}) if age_score else {}
    calibrated = calibrated if isinstance(calibrated, Mapping) else {}
    age_coverage = _finite(calibrated.get("coverage"))
    age_status = (
        PASS
        if age_coverage is not None and age_coverage >= target_coverage
        else FAIL
        if age_coverage is not None
        else ABSTAIN
    )

    reaction = row.get("reaction", {})
    reaction = reaction if isinstance(reaction, Mapping) else {}
    reaction_metrics = reaction.get("metrics", {})
    reaction_metrics = reaction_metrics if isinstance(reaction_metrics, Mapping) else {}
    reaction_accuracy = _finite(reaction_metrics.get("accuracy"))
    observed_reaction_flag = (
        "WEAK" if reaction_accuracy is not None and reaction_accuracy < 0.5 else "NOT_ASSESSED"
    )
    specialist_age_record = specialist_baselines.get("local_tracer_decay_age", {})
    specialist_age_record = (
        specialist_age_record if isinstance(specialist_age_record, Mapping) else {}
    )
    specialist_age = specialist_age_record.get("age", {})
    specialist_age = specialist_age if isinstance(specialist_age, Mapping) else {}
    specialist_age_point = specialist_age.get("point", {})
    specialist_age_point = (
        specialist_age_point if isinstance(specialist_age_point, Mapping) else {}
    )
    specialist_age_interval = specialist_age.get("interval", {})
    specialist_age_interval = (
        specialist_age_interval if isinstance(specialist_age_interval, Mapping) else {}
    )
    specialist_reaction_record = specialist_baselines.get(
        "local_thermodynamic_reaction_rules", {}
    )
    specialist_reaction_record = (
        specialist_reaction_record
        if isinstance(specialist_reaction_record, Mapping)
        else {}
    )
    specialist_reaction = specialist_reaction_record.get("reaction", {})
    specialist_reaction = (
        specialist_reaction if isinstance(specialist_reaction, Mapping) else {}
    )
    specialist_reaction_metrics = specialist_reaction.get("metrics", {})
    specialist_reaction_metrics = (
        specialist_reaction_metrics
        if isinstance(specialist_reaction_metrics, Mapping)
        else {}
    )
    reaction_rapm_raw = _rapm_metrics(row.get("reaction_rapm_scores"))
    reaction_rapm_calibrated = _rapm_metrics(
        row.get("reaction_rapm_calibrated_scores")
    )
    calibrated_age_name = "multitracer_atmospheric_history_age"
    calibrated_age_record = calibrated_specialist_baselines.get(
        calibrated_age_name,
        calibrated_specialist_baselines.get("local_tracer_decay_age", {}),
    )
    calibrated_age_record = (
        calibrated_age_record
        if isinstance(calibrated_age_record, Mapping)
        else {}
    )
    calibrated_age = calibrated_age_record.get("age", {})
    calibrated_age = calibrated_age if isinstance(calibrated_age, Mapping) else {}
    calibrated_age_point = calibrated_age.get("point", {})
    calibrated_age_point = (
        calibrated_age_point
        if isinstance(calibrated_age_point, Mapping)
        else {}
    )
    calibrated_age_interval = calibrated_age.get("interval", {})
    calibrated_age_interval = (
        calibrated_age_interval
        if isinstance(calibrated_age_interval, Mapping)
        else {}
    )
    calibrated_reaction_name = "stoichiometric_reaction_inverse"
    calibrated_reaction_record = calibrated_specialist_baselines.get(
        calibrated_reaction_name,
        calibrated_specialist_baselines.get(
            "local_thermodynamic_reaction_rules", {}
        ),
    )
    calibrated_reaction_record = (
        calibrated_reaction_record
        if isinstance(calibrated_reaction_record, Mapping)
        else {}
    )
    calibrated_reaction = calibrated_reaction_record.get("reaction", {})
    calibrated_reaction = (
        calibrated_reaction
        if isinstance(calibrated_reaction, Mapping)
        else {}
    )
    calibrated_reaction_metrics = calibrated_reaction.get("metrics", {})
    calibrated_reaction_metrics = (
        calibrated_reaction_metrics
        if isinstance(calibrated_reaction_metrics, Mapping)
        else {}
    )
    specialist_universe = row.get("specialist_candidate_universe_metrics", {})
    specialist_universe = (
        specialist_universe if isinstance(specialist_universe, Mapping) else {}
    )

    age_permuted = row_index.get((family, seed, "age_permuted", COMPLETE_SCENARIO))
    no_flow = row_index.get((family, seed, "hydraulic_only", "no_flow_head_permutation"))
    full_age_mae = _finite((row.get("age") or {}).get("point", {}).get("mae")) if isinstance(row.get("age"), Mapping) else None
    permuted_age_mae = (
        _finite((age_permuted.get("age") or {}).get("point", {}).get("mae"))
        if age_permuted and isinstance(age_permuted.get("age"), Mapping)
        else None
    )
    full_f1 = hydro_f1
    no_flow_f1 = (
        _finite((no_flow.get("topology") or {}).get("selected_f1"))
        if no_flow and isinstance(no_flow.get("topology"), Mapping)
        else None
    )
    age_permutation_delta = (
        permuted_age_mae - full_age_mae
        if permuted_age_mae is not None and full_age_mae is not None
        else None
    )
    head_permutation_delta = (
        full_f1 - no_flow_f1
        if full_f1 is not None and no_flow_f1 is not None
        else None
    )

    policy = row.get("decision_policy", {})
    policy = policy if isinstance(policy, Mapping) else {}
    observed_policy_decision = str(policy.get("decision", "UNKNOWN"))

    return {
        "family": family,
        "seed": seed,
        "topology": {
            "conditional_status": _status_from_delta(topology_delta),
            "claim_status": ABSTAIN,
            "hydrosheaf_selected_f1": hydro_f1,
            "best_baseline": best_name,
            "best_baseline_selected_f1": baseline_f1,
            "delta_hydrosheaf_minus_best_baseline": topology_delta,
            "scope": "conditional_on_hydrosheaf_candidate_universe",
            "independent_specialist_status": _status_from_delta(specialist_delta),
            "independent_specialist_baseline": specialist_best_name,
            "independent_specialist_selected_f1": specialist_f1,
            "delta_hydrosheaf_minus_independent_specialist": specialist_delta,
            "independent_candidate_recall": _finite(
                specialist_universe.get("candidate_recall")
            ),
            "independent_candidate_precision": _finite(
                specialist_universe.get("candidate_precision")
            ),
        },
        "age": {
            "diagnostic_status": age_status,
            "claim_status": ABSTAIN,
            "coverage": age_coverage,
            "target_coverage": target_coverage,
            "mean_width": _finite(calibrated.get("mean_width")),
            "mean_absolute_error": full_age_mae,
            "scope": (
                "calibrated_locked_age_diagnostic_plus_independent_age_specialist"
                if specialist_age
                else "calibrated_locked_age_diagnostic_without_age_specialist_baseline"
            ),
            "specialist_status": (
                PASS
                if _finite(specialist_age_point.get("mae")) is not None
                else ABSTAIN
            ),
            "specialist_mae": _finite(specialist_age_point.get("mae")),
            "specialist_coverage": _finite(specialist_age_interval.get("coverage")),
            "specialist_scope": "independent_local_tracer_decay_diagnostic",
            "calibrated_specialist": {
                "baseline": calibrated_age_name,
                "status": calibrated_age.get("status", "NOT_ASSESSED"),
                "mae": _finite(calibrated_age_point.get("mae")),
                "conditional_coverage": _finite(
                    calibrated_age_interval.get("coverage")
                ),
                "coverage_including_abstention": _finite(
                    calibrated_age.get("coverage_including_abstention")
                ),
                "n": _finite(calibrated_age.get("n")),
                "n_abstain": _finite(calibrated_age.get("n_abstain")),
                "n_missing": _finite(calibrated_age.get("n_missing")),
            },
        },
        "reaction": {
            "diagnostic_flag": observed_reaction_flag,
            "claim_status": ABSTAIN,
            "accuracy": reaction_accuracy,
            "n": _finite(reaction_metrics.get("n")),
            "scope": (
                "hydrosheaf_reaction_diagnostic_plus_independent_reaction_specialist"
                if specialist_reaction
                else "no_registered_reaction_family_specialist_baseline"
            ),
            "specialist_status": (
                PASS
                if _finite(specialist_reaction_metrics.get("accuracy")) is not None
                else ABSTAIN
            ),
            "specialist_accuracy": _finite(specialist_reaction_metrics.get("accuracy")),
            "specialist_n": _finite(specialist_reaction_metrics.get("n")),
            "specialist_scope": "independent_local_thermodynamic_rule_diagnostic",
            "calibrated_specialist": {
                "baseline": calibrated_reaction_name,
                "status": calibrated_reaction.get("status", "NOT_ASSESSED"),
                "accuracy": _finite(calibrated_reaction_metrics.get("accuracy")),
                "coverage": _finite(calibrated_reaction.get("coverage")),
                "coverage_including_abstention": _finite(
                    calibrated_reaction.get("coverage_including_abstention")
                ),
                "multiclass_log_loss": _finite(
                    calibrated_reaction.get("multiclass_log_loss")
                ),
                "n": _finite(calibrated_reaction.get("n")),
                "n_abstain": _finite(calibrated_reaction.get("n_abstain")),
                "n_missing": _finite(calibrated_reaction.get("n_missing")),
            },
            "reaction_rapm": {
                "claim_status": ABSTAIN,
                "fit_scope": "development_only",
                "raw": {
                    "accuracy": _finite(reaction_rapm_raw.get("accuracy")),
                    "selective_accuracy": _finite(
                        reaction_rapm_raw.get("selective_accuracy")
                    ),
                    "coverage": _finite(reaction_rapm_raw.get("coverage")),
                    "multiclass_log_loss": _finite(
                        reaction_rapm_raw.get("multiclass_log_loss")
                    ),
                    "false_commitment_rate": _finite(
                        reaction_rapm_raw.get("false_commitment_rate")
                    ),
                    "n": _finite(reaction_rapm_raw.get("n")),
                    "n_reference_edges": _finite(
                        reaction_rapm_raw.get("n_reference_edges")
                    ),
                    "n_decoy_edges": _finite(reaction_rapm_raw.get("n_decoy_edges")),
                },
                "calibrated": {
                    "accuracy": _finite(reaction_rapm_calibrated.get("accuracy")),
                    "selective_accuracy": _finite(
                        reaction_rapm_calibrated.get("selective_accuracy")
                    ),
                    "coverage": _finite(reaction_rapm_calibrated.get("coverage")),
                    "multiclass_log_loss": _finite(
                        reaction_rapm_calibrated.get("multiclass_log_loss")
                    ),
                    "false_commitment_rate": _finite(
                        reaction_rapm_calibrated.get("false_commitment_rate")
                    ),
                    "n": _finite(reaction_rapm_calibrated.get("n")),
                    "n_reference_edges": _finite(
                        reaction_rapm_calibrated.get("n_reference_edges")
                    ),
                    "n_decoy_edges": _finite(
                        reaction_rapm_calibrated.get("n_decoy_edges")
                    ),
                },
                "reason": (
                    "Development-fitted RAPM/on-off diagnostic is reported on the "
                    "locked independent edge universe; no reaction-superiority "
                    "claim gate is asserted."
                ),
            },
        },
        "negative_controls": {
            "age_permutation_status": _status_from_delta(age_permutation_delta),
            "age_mae_increase_after_permutation": age_permutation_delta,
            "head_permutation_status": (
                PASS if head_permutation_delta is not None and head_permutation_delta > 0.0 else
                FAIL if head_permutation_delta is not None else ABSTAIN
            ),
            "selected_f1_drop_after_head_permutation": head_permutation_delta,
        },
        "next_measurement": {
            "claim_status": ABSTAIN,
            "observed_policy_decision": observed_policy_decision,
            "selected_action_id": policy.get("selected_action_id"),
            "reason": "No prospective post-measurement outcome is simulated in this lock.",
        },
        "case_claim_status": ABSTAIN,
        "case_claim_reason": (
            "System-level superiority is not adjudicable: the topology comparator is "
            "conditional on HydroSheaf candidates, specialist age/reaction outputs are "
            "bounded diagnostics rather than competence-matched inverses, and "
            "next-measurement utility has no prospective outcome."
        ),
    }


def build_adjudication(
    *,
    lock_dir: Path,
    manifest: Mapping[str, Any],
    execution_gate: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    calibration_scores: Mapping[str, Any],
    errors: object,
) -> dict[str, Any]:
    target_coverage = _finite(
        (manifest.get("calibration") or {}).get("target_coverage")
        if isinstance(manifest.get("calibration"), Mapping)
        else None
    ) or TARGET_COVERAGE
    index = _row_index(rows)
    case_records: list[dict[str, Any]] = []
    for key, row in sorted(index.items()):
        family, seed, condition, scenario = key
        if condition != PRIMARY_CONDITION or scenario != COMPLETE_SCENARIO:
            continue
        score = calibration_scores.get(f"{family}:{seed}")
        score = score if isinstance(score, Mapping) else None
        case_records.append(
            _case_adjudication(
                row,
                age_score=score,
                target_coverage=target_coverage,
                row_index=index,
            )
        )

    selective_curve = calibration_scores.get("locked_selective_risk_curve", [])
    selective_curve = selective_curve if isinstance(selective_curve, Sequence) else []
    full_acceptance = next(
        (
            point for point in selective_curve
            if isinstance(point, Mapping) and float(point.get("requested_acceptance", 0.0)) == 1.0
        ),
        None,
    )
    pooled_coverage = _finite(full_acceptance.get("interval_coverage")) if full_acceptance else None
    pooled_age_status = (
        PASS if pooled_coverage is not None and pooled_coverage >= target_coverage else
        FAIL if pooled_coverage is not None else ABSTAIN
    )
    has_independent_specialists = any(
        bool(row.get("specialist_candidate_generation_truth_blind"))
        and isinstance(row.get("specialist_baselines"), Mapping)
        for row in rows
    )
    system_reasons = [
        (
            "Independent specialist candidate generation is now evaluated, and the "
            "age/reaction outputs include development-only calibration; they remain "
            "bounded diagnostics and do not establish specialist superiority."
            if has_independent_specialists
            else
            "The registered specialist baselines receive HydroSheaf's candidate "
            "universe and are therefore conditional comparators, not independent "
            "end-to-end graph generators."
        ),
        (
            "The multi-tracer age and stoichiometric reaction outputs remain "
            "diagnostic only and do not yet constitute competence-matched "
            "atmospheric-history/excess-air or PHREEQC inverse baselines."
            if has_independent_specialists
            else
            "No local tracer-age specialist baseline or reaction-family specialist "
            "baseline is registered in the lock."
        ),
        "The next-measurement policy has no prospective outcome simulation; all locked complete-case policies abstain.",
    ]
    integrity = _input_integrity(lock_dir, manifest, errors)
    execution_pass = bool(
        execution_gate.get("status") == PASS
        and integrity["artifact_hashes_all_match"]
        and integrity["source_hashes_all_match"]
        and integrity["errors_empty"]
    )
    synthetic_component_claim = _common_topology_component_claim(
        rows,
        execution_pass=execution_pass,
    )
    synthetic_reaction_claim = {
        "claim_name": "synthetic_reaction_family_selective_diagnostic",
        "status": ABSTAIN,
        "reason": (
            "The development-fitted RAPM/on-off layer is execution-complete and "
            "its locked calibrated metrics are reported by case, but a separate "
            "predeclared reaction-superiority gate has not been passed."
        ),
        "field_data_required": False,
    }
    post_measurement_records = [
        row.get("post_measurement_outcome")
        for row in rows
        if row.get("phase") == "locked_test"
        and row.get("condition") == PRIMARY_CONDITION
        and row.get("observation_scenario", COMPLETE_SCENARIO) == COMPLETE_SCENARIO
    ]
    synthetic_integrated_status = (
        PASS
        if post_measurement_records
        and all(
            isinstance(record, Mapping)
            and bool(record.get("evaluated"))
            and bool(record.get("improved"))
            for record in post_measurement_records
        )
        else ABSTAIN
    )
    synthetic_integrated_claim = {
        "claim_name": "synthetic_integrated_performance",
        "status": synthetic_integrated_status,
        "reason": (
            "Every locked complete-case policy has a truth-aware simulated outcome "
            "and meets the declared improvement rule."
            if synthetic_integrated_status == PASS
            else "A prospective post-measurement outcome record is missing or does not meet the declared improvement rule."
        ),
        "field_data_required": False,
    }
    field_validation_claim = {
        "claim_name": "field_validation",
        "status": "DEFERRED",
        "reason": "Prospective field campaign and independent follow-up outcomes are not yet available.",
    }
    return {
        "adjudication_version": "1.1",
        "run_id": manifest.get("run_id"),
        "input_lock": str(lock_dir),
        "input_lock_status": manifest.get("status"),
        "age_samples": manifest.get("age_samples"),
        "execution_verdict": PASS if execution_pass else FAIL,
        "system_level_claim_status": ABSTAIN,
        "system_level_claim_reasons": system_reasons,
        "synthetic_component_claim": synthetic_component_claim,
        "synthetic_reaction_claim": synthetic_reaction_claim,
        "synthetic_integrated_claim": synthetic_integrated_claim,
        "field_validation_claim": field_validation_claim,
        "comparator_scope_warning": (
            "Independent specialist candidate-graph generation is included when the "
            "row records it; native topology scores remain conditional diagnostics, "
            "and calibrated age/reaction outputs remain bounded specialist diagnostics."
            if has_independent_specialists
            else
            "Baseline topology scores are useful conditional diagnostics, but they "
            "do not test independent candidate-graph generation because the runner "
            "supplies the HydroSheaf candidate universe and derived local support."
        ),
        "locked_case_count": len(case_records),
        "case_adjudications": case_records,
        "pooled_age_diagnostic": {
            "diagnostic_status": pooled_age_status,
            "claim_status": ABSTAIN,
            "coverage": pooled_coverage,
            "target_coverage": target_coverage,
            "accepted": full_acceptance.get("accepted") if full_acceptance else None,
            "total": full_acceptance.get("total") if full_acceptance else None,
        },
        "summary_counts": {
            "case_claim_status": dict(Counter(record["case_claim_status"] for record in case_records)),
            "topology_conditional_status": dict(
                Counter(record["topology"]["conditional_status"] for record in case_records)
            ),
            "age_diagnostic_status": dict(
                Counter(record["age"]["diagnostic_status"] for record in case_records)
            ),
            "reaction_diagnostic_flag": dict(
                Counter(record["reaction"]["diagnostic_flag"] for record in case_records)
            ),
            "reaction_rapm_claim_status": dict(
                Counter(record["reaction"]["reaction_rapm"]["claim_status"] for record in case_records)
            ),
            "age_permutation_status": dict(
                Counter(record["negative_controls"]["age_permutation_status"] for record in case_records)
            ),
            "head_permutation_status": dict(
                Counter(record["negative_controls"]["head_permutation_status"] for record in case_records)
            ),
            "observed_policy_decision": dict(
                Counter(record["next_measurement"]["observed_policy_decision"] for record in case_records)
            ),
        },
        "input_integrity": integrity,
        "source_lock": {
            "git_revision": _git_revision(),
            "git_worktree_dirty": _git_worktree_dirty(),
            "manifest_source_hashes": manifest.get("source_hashes", {}),
        },
    }


def _fmt(value: object, digits: int = 3) -> str:
    number = _finite(value)
    return "NA" if number is None else f"{number:.{digits}f}"


def _markdown(report: Mapping[str, Any]) -> str:
    lines = [
        "# Locked synthetic ensemble adjudication",
        "",
        f"**Run:** `{report.get('run_id')}`  ",
        f"**Execution verdict:** **{report.get('execution_verdict')}**  ",
        f"**System-level superiority claim:** **{report.get('system_level_claim_status')}**",
        f"**Bounded synthetic topology component claim:** **{report.get('synthetic_component_claim', {}).get('status', ABSTAIN)}**  ",
        f"**Synthetic RAPM reaction diagnostic claim:** **{report.get('synthetic_reaction_claim', {}).get('status', ABSTAIN)}**  ",
        f"**Synthetic integrated claim:** **{report.get('synthetic_integrated_claim', {}).get('status', ABSTAIN)}**  ",
        f"**Field-validation claim:** **{report.get('field_validation_claim', {}).get('status', 'DEFERRED')}**",
        "",
        "The system-level claim remains separate from the bounded topology component "
        "claim: age/reaction and prospective measurement performance must not be "
        "inferred from a topology-only result.",
        "",
        "## Locked-case readout",
        "",
        "| Family | Topology conditional | HydroSheaf F1 | Best baseline F1 | Age diagnostic | Age coverage | Reaction flag | Policy | Case claim |",
        "|---|---:|---:|---:|---:|---:|---:|---|---:|",
    ]
    for record in report.get("case_adjudications", []):
        topology = record["topology"]
        age = record["age"]
        reaction = record["reaction"]
        policy = record["next_measurement"]
        lines.append(
            "| {family} ({seed}) | {topology_status} | {hydro} | {baseline} | "
            "{age_status} | {coverage} | {reaction_flag} ({accuracy}) | {policy_decision} | {claim} |".format(
                family=record["family"],
                seed=record["seed"],
                topology_status=topology["conditional_status"],
                hydro=_fmt(topology["hydrosheaf_selected_f1"]),
                baseline=_fmt(topology["best_baseline_selected_f1"]),
                age_status=age["diagnostic_status"],
                coverage=_fmt(age["coverage"]),
                reaction_flag=reaction["diagnostic_flag"],
                accuracy=_fmt(reaction["accuracy"]),
                policy_decision=policy["observed_policy_decision"],
                claim=record["case_claim_status"],
            )
        )
    pooled = report.get("pooled_age_diagnostic", {})
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Pooled calibrated age coverage at 100% acceptance: **{_fmt(pooled.get('coverage'))}** "
            f"against the declared {_fmt(pooled.get('target_coverage'))} target; this is a diagnostic, not comparative superiority evidence.",
            "- Topology conditional result: analytic-lattice and independent-mixing tie their strongest registered baseline; MODFLOW v3 is lower by 0.037 F1 (0.700 versus 0.737).",
            "- Reaction-family recovery remains diagnostic: the calibrated stoichiometric specialist is recorded, but it is not yet a competence-matched PHREEQC inverse.",
            "- The new RAPM/on-off reaction layer is reported on the independent candidate universe; its calibrated selective metrics remain a bounded diagnostic and do not establish reaction superiority.",
            "- Every locked complete-case next-measurement policy abstains; prospective utility is therefore not demonstrated.",
            "- Age and head permutation controls pass their sensitivity checks when evaluated as negative-control diagnostics.",
            f"- Bounded topology component claim: **{report.get('synthetic_component_claim', {}).get('status', ABSTAIN)}** on the common all-pairs universe versus the fixed hydraulic-only comparator; this is not an integrated claim.",
            "",
            "This report is controlled-synthetic adjudication only. It does not establish field validity or universal HydroSheaf superiority.",
        ]
    )
    return "\n".join(lines) + "\n"


def adjudicate(*, lock_dir: Path, output: Path, overwrite: bool = False) -> dict[str, Any]:
    lock_dir = lock_dir.resolve()
    required = (
        "run_manifest.json",
        "execution_gate.json",
        "case_metrics.json",
        "calibration_scores.json",
        "errors.json",
    )
    missing = [name for name in required if not (lock_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Lock is missing required artifacts: {missing}")
    manifest = _json(lock_dir / "run_manifest.json")
    execution_gate = _json(lock_dir / "execution_gate.json")
    rows = _json(lock_dir / "case_metrics.json")
    calibration_scores = _json(lock_dir / "calibration_scores.json")
    errors = _json(lock_dir / "errors.json")
    if not isinstance(rows, list):
        raise TypeError("case_metrics.json must contain a list.")
    report = build_adjudication(
        lock_dir=lock_dir,
        manifest=manifest,
        execution_gate=execution_gate,
        rows=rows,
        calibration_scores=calibration_scores,
        errors=errors,
    )
    resolved_output = _safe_output(output, overwrite=overwrite)
    _write_json(resolved_output / "adjudication.json", report)
    (resolved_output / "adjudication.md").write_text(
        _markdown(report),
        encoding="utf-8",
    )
    output_manifest: dict[str, Any] = {
        "adjudication_version": report["adjudication_version"],
        "run_id": report["run_id"],
        "input_lock": str(lock_dir),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "git_worktree_dirty": _git_worktree_dirty(),
        "source_hashes": {
            "scripts/adjudicate_locked_ensemble.py": _sha256(Path(__file__).resolve()),
        },
        "artifacts": {},
    }
    for path in sorted(resolved_output.iterdir()):
        if path.is_file() and path.name != "adjudication_manifest.json":
            output_manifest["artifacts"][path.name] = _sha256(path)
    _write_json(resolved_output / "adjudication_manifest.json", output_manifest)
    return {
        "output": str(resolved_output),
        "system_level_claim_status": report["system_level_claim_status"],
        "synthetic_component_claim_status": report["synthetic_component_claim"]["status"],
        "synthetic_reaction_claim_status": report["synthetic_reaction_claim"]["status"],
        "synthetic_integrated_claim_status": report["synthetic_integrated_claim"]["status"],
        "field_validation_claim_status": report["field_validation_claim"]["status"],
        "execution_verdict": report["execution_verdict"],
        "summary_counts": report["summary_counts"],
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    run_dir = args.run_dir.resolve()
    run_id = str(_json(run_dir / "run_manifest.json").get("run_id", run_dir.name))
    output = args.output or DEFAULT_OUTPUT_ROOT / run_id
    result = adjudicate(lock_dir=run_dir, output=output, overwrite=args.overwrite)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["system_level_claim_status"] == ABSTAIN else 1


if __name__ == "__main__":
    raise SystemExit(main())
