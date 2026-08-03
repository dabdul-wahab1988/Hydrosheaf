"""Run development-fitted calibration on an independent synthetic ensemble.

The historical M7 v2 MODFLOW/MODPATH extension remains available for audit,
but the active MODFLOW/MODPATH family is the independent multilayer v3
generator.  The
analytic-lattice generator is a separate geometry and forward-model family.
Calibration truth is used only for the declared development split; locked
family/seed cases receive a frozen calibrator before their truth is scored.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Mapping, Sequence

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
M7_SCRIPTS = REPO / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
for import_path in (REPO, SCRIPTS, M7_SCRIPTS):
    if str(import_path) not in sys.path:
        sys.path.insert(0, str(import_path))

from hydrosheaf.validation import (  # noqa: E402
    AgeCalibrationObservation,
    ProgrammeDecision,
    ProgrammeRun,
    ProgrammeStage,
    DecisionKind,
    IdentifiabilityStatus,
    StageStatus,
    SELECT,
    apply_age_interval_calibrator,
    fit_age_interval_calibrator,
    score_age_posteriors,
    score_topology,
    score_calibrated_age_intervals,
    score_selective_risk,
    finite_numeric_mapping,
    probability_metrics,
    programme_stages_from_status,
    audit_generator_case,
    OBSERVATION_STRESS_SCENARIOS,
    apply_observation_stress,
    make_no_flow_control,
    make_permutation_control,
    default_baseline_registry,
    default_specialist_baseline_registry,
    score_age_baseline_outputs,
    score_reaction_baseline_outputs,
    normalise_reaction_family,
    IndependentCandidateEdge,
    generate_independent_candidate_universe,
    SpecialistAgeCalibrationObservation,
    SpecialistReactionCalibrationObservation,
    fit_specialist_age_calibrator,
    fit_specialist_reaction_calibrator,
    score_calibrated_specialist_age,
    score_calibrated_specialist_reaction,
    assess_age_gate,
    assess_integrated_gate,
    assess_kinetics_gate,
    assess_reaction_gate,
    aggregate_specialist_gates,
    run_kinetics_specialist_benchmark,
    AgePerformanceThresholds,
    evaluate_age_performance,
    ProspectiveMeasurementCase,
    ProspectivePolicy,
    evaluate_prospective_policies,
    select_random_measurement,
    select_declared_utility_measurement,
    select_specialist_measurement,
    ReactionRAPM,
    ReactionRAPMCalibrator,
    ReactionRAPMTrainingExample,
    cross_fitted_reaction_rapm_calibration_examples,
    fit_reaction_rapm,
    score_reaction_rapm_outputs,
    training_examples_from_observations,
    normalise_rapm_reaction_family,
    DiscreteModelObservation,
    fit_discrete_model_weights,
    score_locked_model_average,
    CandidateMeasurementAction,
    ScenarioBelief,
    select_next_measurement,
    consolidate_head_evidence,
    PreMeasurementPosteriorSummary,
    SyntheticObservationModel,
    simulate_post_measurement,
    DiscrepancyCalibrationObservation,
    apply_discrepancy_calibrator_to_record,
    fit_discrepancy_calibrator,
    score_locked_discrepancy,
)
from independent_lattice_generator import (  # noqa: E402
    IndependentLatticeAquifer,
    generate_independent_lattice,
)
from independent_mixing_generator import (  # noqa: E402
    IndependentMixingAquifer,
    generate_independent_mixing,
)
from independent_modflow_generator_v3 import (  # noqa: E402
    IndependentAquiferV3,
    generate_independent_aquifer_v3,
)
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402
from run_programme_synthetic_gate import (  # noqa: E402
    CONDITIONS,
    DEFAULT_BIN_DIR,
    DEFAULT_SIMULATOR_WORKSPACE,
    RUN_CONDITIONS,
    _build_discrepancy_reports,
    _condition_config,
    _git_revision,
    _git_worktree_dirty,
    _run_condition,
    _prepare_rows,
    _safe_output,
    _sha256,
    _write_json,
)

Case = IndependentAquiferV3 | IndependentLatticeAquifer | IndependentMixingAquifer

RUN_ID = "RUN-PROGRAMME-ENSEMBLE-SPECIALIST-CALIBRATION-20260801-06"
DEVELOPMENT_CASES = {
    "modflow_modpath_v3": (9101, 9102, 9103, 9104),
    "analytic_lattice_v1": (9301, 9302, 9303, 9304),
    "independent_mixing_v1": (9501, 9502, 9503, 9504),
}
LOCKED_CASES = {
    "modflow_modpath_v3": (9201, 9202),
    "analytic_lattice_v1": (9401, 9402),
    "independent_mixing_v1": (9601, 9602),
}
PROSPECTIVE_MINIMUM_LOCKED_CASES = 6
REQUIRED_GENERATOR_FAMILIES = {
    "modflow_modpath_v3",
    "analytic_lattice_v1",
    "independent_mixing_v1",
}
PRIMARY_GENERATOR_FAMILIES = {
    "analytic_lattice_v1",
    "independent_mixing_v1",
}
ALL_GENERATOR_FAMILIES = frozenset(REQUIRED_GENERATOR_FAMILIES)
BASELINE_REGISTRY = default_baseline_registry()
SPECIALIST_BASELINE_REGISTRY = default_specialist_baseline_registry()
PROGRAMME_SOURCE_FILES = (
    "hydrosheaf/validation/__init__.py",
    "hydrosheaf/validation/calibration.py",
    "hydrosheaf/validation/generator_critic.py",
    "hydrosheaf/validation/head_evidence.py",
    "hydrosheaf/validation/observation_stress.py",
    "hydrosheaf/validation/programme_contract.py",
    "hydrosheaf/validation/programme_gate.py",
    "hydrosheaf/validation/discrepancy.py",
    "hydrosheaf/validation/metrics.py",
    "hydrosheaf/validation/controls.py",
    "hydrosheaf/validation/baselines.py",
    "hydrosheaf/validation/specialist_baselines.py",
    "hydrosheaf/validation/specialist_candidate_generation.py",
    "hydrosheaf/validation/specialist_calibration.py",
    "hydrosheaf/validation/age_competent_baseline.py",
    "hydrosheaf/validation/age_performance.py",
    "hydrosheaf/validation/reaction_competent_baseline.py",
    "hydrosheaf/validation/reaction_rapm.py",
    "hydrosheaf/validation/synthetic_claims.py",
    "hydrosheaf/validation/model_averaging.py",
    "hydrosheaf/validation/decision_utility.py",
    "hydrosheaf/validation/kinetics_specialist_benchmark.py",
    "hydrosheaf/validation/performance_contract.py",
    "hydrosheaf/inference/network_fit.py",
    "hydrosheaf/config.py",
    "scripts/independent_lattice_generator.py",
    "scripts/independent_mixing_generator.py",
    "M7/m7_nonuniqueness_benchmark/scripts/independent_modflow_generator_v3.py",
    "scripts/run_programme_synthetic_gate.py",
    "scripts/run_programme_generator_ensemble.py",
    "scripts/adjudicate_locked_ensemble.py",
)
DEFAULT_OUTPUT = REPO / ".codex_work" / "programme-validation" / RUN_ID
M7_V3_SOURCE = REPO / "M7" / "m7_nonuniqueness_benchmark" / "scripts" / "independent_modflow_generator_v3.py"
GENERATOR_SOURCE_PATHS = {
    "modflow_modpath_v3": M7_V3_SOURCE,
    "analytic_lattice_v1": REPO / "scripts" / "independent_lattice_generator.py",
    "independent_mixing_v1": REPO / "scripts" / "independent_mixing_generator.py",
}


def _head_model_for_family(
    family: str,
    case: Case,
) -> dict[str, object]:
    """Return the preregistered head observation model for one family."""

    source_model = case.provenance.get("head_channel_covariance_model", {})
    model = dict(source_model) if isinstance(source_model, Mapping) else {}
    model.update(
        {
            "primary_channel": "head_meas",
            "secondary_channel": "hydraulic_head",
        }
    )
    if family == "analytic_lattice_v1":
        # The lattice secondary channel contains a structured local bias.  It
        # is retained as a discrepancy diagnostic, not as a second independent
        # measurement in the topology likelihood.
        model.update(
            {
                "primary_sigma_m": 0.06,
                "secondary_sigma_m": 0.14,
                "measurement_error_correlation": 0.0,
                "combination": "primary_only_with_discrepancy",
            }
        )
    elif family == "independent_mixing_v1":
        model.update(
            {
                "primary_sigma_m": 0.08,
                "secondary_sigma_m": 0.12,
                "measurement_error_correlation": 0.0,
                "combination": "primary_only_with_discrepancy",
            }
        )
    else:
        # v3 exposes two noisy measurements of the same multilayer MODFLOW
        # field.  They are combined with an explicit zero conditional-error
        # correlation, rather than counted as two unrelated priors.
        model.update(
            {
                "primary_sigma_m": 0.10,
                "secondary_sigma_m": 0.10,
                "measurement_error_correlation": 0.0,
                "combination": "gls",
            }
        )
    return model


def _prepare_observation_scenarios(
    family: str,
    case: Case,
) -> tuple[
    dict[str, tuple[dict[str, object], ...]],
    dict[str, object],
    dict[str, object],
]:
    """Build blind stress variants and the head-evidence audit record."""

    head_model = _head_model_for_family(family, case)
    consolidated, head_audit = consolidate_head_evidence(
        case.observations,
        model=head_model,
    )
    scenarios: dict[str, tuple[dict[str, object], ...]] = {}
    stress_records: dict[str, object] = {}
    for scenario_index, scenario in enumerate(OBSERVATION_STRESS_SCENARIOS):
        result = apply_observation_stress(
            consolidated,
            scenario,
            seed=int(case.seed) + 301_000 + scenario_index * 7_919,
        )
        scenarios[scenario] = result.observations
        stress_records[scenario] = result.to_dict()
    controls = (
        make_no_flow_control(
            consolidated,
            seed=int(case.seed) + 401_001,
            head_fields=("head_meas", "hydraulic_head"),
        ),
        make_permutation_control(
            consolidated,
            ("tritium_TU",),
            seed=int(case.seed) + 401_002,
            scenario="tracer_permutation",
        ),
        make_permutation_control(
            consolidated,
            ("SO4", "NO3", "Fe"),
            seed=int(case.seed) + 401_003,
            scenario="chemistry_permutation",
        ),
    )
    for control in controls:
        scenarios[control.scenario] = control.observations
        stress_records[control.scenario] = control.to_dict()
    return scenarios, head_audit, stress_records


def _edge_pair(edge: object) -> tuple[str, str]:
    if hasattr(edge, "u") and hasattr(edge, "v"):
        return str(getattr(edge, "u")), str(getattr(edge, "v"))
    if isinstance(edge, Mapping):
        return str(edge["u"]), str(edge["v"])
    values = tuple(edge)  # type: ignore[arg-type]
    if len(values) < 2:
        raise ValueError(f"Edge must contain at least two values: {edge!r}")
    return str(values[0]), str(values[1])


def _edge_attrs(edge: object) -> Mapping[str, object]:
    """Return declared edge metadata without assuming a concrete edge class."""

    if hasattr(edge, "attrs"):
        attrs = getattr(edge, "attrs")
        return attrs if isinstance(attrs, Mapping) else {}
    if isinstance(edge, Mapping):
        attrs = edge.get("attrs", {})
        return attrs if isinstance(attrs, Mapping) else {}
    return {}


def _finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _baseline_observations(
    rows: Sequence[Mapping[str, object]],
    candidates: Sequence[object],
) -> dict[str, object]:
    """Build only declared observed channels for specialist baselines."""

    by_id = {
        str(row.get("site_id")): row
        for row in rows
        if row.get("site_id") is not None
    }
    node_heads = {
        node_id: _finite_float(row.get("head_meas"))
        for node_id, row in by_id.items()
        if _finite_float(row.get("head_meas")) is not None
    }
    tracer_nodes: dict[str, dict[str, float]] = {}
    chemistry_nodes: dict[str, dict[str, float]] = {}
    tracer_fields = (
        "tritium_TU",
        "tritium_sigma_TU",
        "argon39_pmc",
        "argon39_sigma_pmc",
        "c14_pmc",
        "cfc11_pptv",
        "cfc11_sigma_pptv",
        "cfc12_pptv",
        "cfc12_sigma_pptv",
        "cfc113_pptv",
        "cfc113_sigma_pptv",
        "sf6_pptv",
        "sf6_sigma_pptv",
        "he4_ccpg",
        "h3_he3_TU",
        "sample_date",
    )
    chemistry_fields = (
        "Ca",
        "Cl",
        "F",
        "Fe",
        "HCO3",
        "K",
        "Mg",
        "Na",
        "NO3",
        "PO4",
        "SO4",
        "SiO2",
        "pH",
    )
    for node_id, row in by_id.items():
        tracer_nodes[node_id] = {
            field: value
            for field in tracer_fields
            if (value := _finite_float(row.get(field))) is not None
        }
        chemistry_nodes[node_id] = {
            field: value
            for field in chemistry_fields
            if (value := _finite_float(row.get(field))) is not None
        }
    tracer_history_nodes = {
        node_id: {
            **dict(features),
            "sample_year": (
                _finite_float(by_id[node_id].get("sample_date")) or 2025.0
            ),
        }
        for node_id, features in tracer_nodes.items()
    }
    hydraulic_edges: dict[str, dict[str, float]] = {}
    topology_edges: dict[str, dict[str, float]] = {}
    reaction_edges: dict[str, dict[str, object]] = {}
    for edge in candidates:
        source, target = _edge_pair(edge)
        source_head = _finite_float(by_id.get(source, {}).get("head_meas"))
        target_head = _finite_float(by_id.get(target, {}).get("head_meas"))
        features: dict[str, float] = {}
        if source_head is not None and target_head is not None:
            features["head_drop"] = source_head - target_head
        attrs = _edge_attrs(edge)
        distance = _finite_float(attrs.get("distance_km"))
        if distance is not None and distance > 0.0:
            features["distance"] = distance
        key = f"{source}->{target}"
        hydraulic_edges[key] = features
        support = _finite_float(
            attrs.get(
                "independent_support",
                attrs.get("p_uv", attrs.get("edge_confidence", attrs.get("flow_probability"))),
            )
        )
        topology_edges[key] = {"support": support if support is not None else 0.5}
        reaction_edges[key] = {
            "upstream": dict(chemistry_nodes.get(source, {})),
            "downstream": dict(chemistry_nodes.get(target, {})),
        }
    return {
        "hydraulic": {"node_heads": node_heads, "edges": hydraulic_edges},
        "topology": {"edges": topology_edges},
        "tracer_age": {"nodes": tracer_nodes},
        "tracer_age_history": {"nodes": tracer_history_nodes},
        "reaction_chemistry": {"edges": reaction_edges},
    }


def _infer_baseline_outputs(
    rows: Sequence[Mapping[str, object]],
    candidates: Sequence[object],
) -> dict[str, dict[str, object]]:
    """Run every registered baseline before any generator truth is consulted."""

    observations = _baseline_observations(rows, candidates)
    outputs: dict[str, dict[str, object]] = {}
    for spec in BASELINE_REGISTRY.specs():
        predictions = spec.score(candidates, observations)
        outputs[spec.name] = {
            "output_kind": "topology",
            "metadata": spec.to_audit_record(),
            "predictions": [prediction.to_audit_record() for prediction in predictions],
        }
    for spec in SPECIALIST_BASELINE_REGISTRY.specs():
        predictions = spec.predict(observations)
        outputs[spec.name] = {
            "output_kind": spec.output_kind,
            "metadata": spec.to_audit_record(),
            "predictions": [
                {"target": target, **dict(record)}
                for target, record in sorted(predictions.items())
            ],
        }
    return outputs


def _score_baseline_outputs(
    reference_edges: Sequence[Sequence[object]],
    candidates: Sequence[object],
    outputs: Mapping[str, Mapping[str, object]],
    *,
    true_ages_years: Mapping[str, object] | None = None,
    true_processes: Mapping[str, object] | None = None,
) -> dict[str, dict[str, object]]:
    """Score frozen baseline predictions only after inference is complete."""

    reference = {tuple(str(value) for value in edge[:2]) for edge in reference_edges}
    scored: dict[str, dict[str, object]] = {}
    candidate_pairs = [_edge_pair(edge) for edge in candidates]
    labels = [1.0 if pair in reference else 0.0 for pair in candidate_pairs]
    for name, record in outputs.items():
        prediction_rows = record.get("predictions", [])
        if not isinstance(prediction_rows, Sequence):
            continue
        output_kind = str(record.get("output_kind", "topology"))
        metadata = record.get("metadata", {})
        if output_kind == "age_interval":
            age_predictions = {
                str(prediction.get("target")): prediction
                for prediction in prediction_rows
                if isinstance(prediction, Mapping) and prediction.get("target") is not None
            }
            scored[name] = {
                "metadata": metadata,
                "age": score_age_baseline_outputs(
                    true_ages_years or {},
                    age_predictions,
                ),
            }
            continue
        if output_kind == "reaction_family":
            reaction_predictions = {
                str(prediction.get("target")): prediction
                for prediction in prediction_rows
                if isinstance(prediction, Mapping) and prediction.get("target") is not None
            }
            scored[name] = {
                "metadata": metadata,
                "reaction": score_reaction_baseline_outputs(
                    true_processes or {},
                    reaction_predictions,
                    candidate_edge_ids=tuple(
                        f"{pair[0]}->{pair[1]}" for pair in candidate_pairs
                    ),
                ),
            }
            continue
        selected = []
        probabilities = []
        for prediction in prediction_rows:
            if not isinstance(prediction, Mapping):
                continue
            edge = prediction.get("edge", [])
            if isinstance(edge, Sequence) and len(edge) >= 2:
                pair = (str(edge[0]), str(edge[1]))
                probabilities.append(_finite_float(prediction.get("probability")))
                if str(prediction.get("decision")) == "select":
                    selected.append(pair)
        usable_probabilities = [
            value if value is not None else float("nan") for value in probabilities
        ]
        scored[name] = {
            "metadata": metadata,
            "topology": score_topology(reference_edges, candidates, selected),
            "probability_calibration": dict(
                probability_metrics(labels, usable_probabilities, n_bins=10)
            ),
        }
    return scored


def _candidate_universe_metrics(
    reference_edges: Sequence[Sequence[object]],
    candidates: Sequence[object],
) -> dict[str, object]:
    """Report what the independent candidate generator could have recovered."""

    reference = {
        (str(edge[0]), str(edge[1]))
        for edge in reference_edges
        if len(edge) >= 2
    }
    candidate_pairs = {_edge_pair(edge) for edge in candidates}
    recovered = reference & candidate_pairs
    omitted = sorted(reference - candidate_pairs)
    return {
        "reference_edge_count": len(reference),
        "candidate_edge_count": len(candidate_pairs),
        "recovered_reference_edge_count": len(recovered),
        "omitted_reference_edge_count": len(omitted),
        "candidate_recall": (
            len(recovered) / len(reference) if reference else 0.0
        ),
        "candidate_precision": (
            len(recovered) / len(candidate_pairs) if candidate_pairs else 0.0
        ),
        "omitted_reference_edges": [list(edge) for edge in omitted],
    }


def _reaction_rapm_observations_for_case(
    case: Case,
) -> tuple[dict[str, object], dict[str, object]]:
    """Build the independent chemistry-edge universe used by reaction RAPM."""

    blind_rows = _prepare_rows(
        case,
        permute_age=False,
        observations=case.observations,
    )
    independent_universe = generate_independent_candidate_universe(
        blind_rows,
        max_neighbors=4,
        head_tie_tolerance_m=0.10,
    )
    return (
        _baseline_observations(blind_rows, independent_universe.edges),
        independent_universe.to_audit_record(),
    )


def _candidate_edge_ids_from_universe(
    universe: Mapping[str, object],
) -> tuple[str, ...]:
    """Return the frozen edge IDs represented by an independent universe."""

    raw_edges = universe.get("edges", ())
    if not isinstance(raw_edges, Sequence):
        return ()
    edge_ids: list[str] = []
    for record in raw_edges:
        if not isinstance(record, Mapping):
            continue
        edge = record.get("edge")
        if not isinstance(edge, Sequence) or isinstance(edge, (str, bytes)):
            continue
        if len(edge) < 2:
            continue
        edge_ids.append(f"{edge[0]}->{edge[1]}")
    return tuple(dict.fromkeys(edge_ids))


def _reaction_rapm_training_examples(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[ReactionRAPMTrainingExample, ...]:
    """Create fit-only labelled examples from development cases.

    Candidate generation remains independent and truth-blind.  Truth is used
    here only after the blind candidate/observation record has been created,
    and only for the declared development fit.
    """

    examples: list[ReactionRAPMTrainingExample] = []
    for family, case, _rows in cases:
        observations, _universe = _reaction_rapm_observations_for_case(case)
        examples.extend(
            training_examples_from_observations(
                observations,
                case.true_processes,
                case_id=f"{family}:{int(case.seed)}",
            )
        )
    return tuple(examples)


def _apply_reaction_rapm_model(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    model: ReactionRAPM,
) -> tuple[dict[str, dict[str, object]], bool]:
    """Apply a frozen RAPM model and score it only after prediction."""

    model_audit = model.to_audit_record()
    case_scores: dict[str, dict[str, object]] = {}
    records_complete = True
    for family, case, rows in cases:
        case_id = f"{family}:{int(case.seed)}"
        observations, universe = _reaction_rapm_observations_for_case(case)
        predictions = model.predict(observations)
        candidate_edge_ids = _candidate_edge_ids_from_universe(universe)
        raw_score = score_reaction_rapm_outputs(
            case.true_processes,
            predictions,
            candidate_edge_ids=candidate_edge_ids,
        )
        inference = {
            "output_kind": "reaction_family",
            "metadata": model_audit,
            "candidate_universe": universe,
            "candidate_edge_ids": list(candidate_edge_ids),
            "predictions": [
                {"target": target, **dict(record)}
                for target, record in sorted(predictions.items())
            ],
            "truth_blind": True,
            "truth_released_for_scoring": False,
        }
        for row in rows:
            if row.get("observation_scenario", "complete") != "complete":
                continue
            row["reaction_rapm_inference"] = inference
            row["reaction_rapm_truth_blind"] = True
            row["reaction_rapm_scores"] = {"reaction": raw_score}
        complete_rows = [
            row
            for row in rows
            if row.get("condition") == "full_sheaf"
            and row.get("observation_scenario", "complete") == "complete"
        ]
        if not complete_rows:
            records_complete = False
        case_scores[case_id] = {
            "family": family,
            "seed": int(case.seed),
            "candidate_universe": universe,
            "candidate_edge_ids": list(candidate_edge_ids),
            "raw_score": raw_score,
            "complete_case_recorded": bool(complete_rows),
        }
    return case_scores, records_complete


def _reaction_rapm_calibration_observations(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[SpecialistReactionCalibrationObservation, ...]:
    """Extract cross-fitted development logits for frozen temperature fitting.

    Each development case is predicted by a model fitted on the other
    development cases.  This prevents the temperature calibrator from seeing
    in-sample logits before it is frozen for the locked cases.
    """

    observations: list[SpecialistReactionCalibrationObservation] = []
    for held_out_index, (family, case, _rows) in enumerate(cases):
        training_cases = tuple(
            item for index, item in enumerate(cases) if index != held_out_index
        )
        if len(training_cases) < 2:
            continue
        try:
            cross_fitted_model = fit_reaction_rapm(
                _reaction_rapm_training_examples(training_cases),
                phase="development",
            )
        except (TypeError, ValueError, np.linalg.LinAlgError):
            continue
        blind_observations, _universe = _reaction_rapm_observations_for_case(case)
        predictions = cross_fitted_model.predict(blind_observations)
        case_id = f"{family}:{int(case.seed)}"
        for edge_id, prediction in predictions.items():
            logits = prediction.get("logits", {})
            if not edge_id or not isinstance(logits, Mapping) or not logits:
                continue
            truth_family = normalise_rapm_reaction_family(
                case.true_processes.get(edge_id, "none")
            )
            observations.append(
                SpecialistReactionCalibrationObservation(
                    case_id=case_id,
                    edge_id=edge_id,
                    truth_family=truth_family,
                    logits={str(key): float(value) for key, value in logits.items()},
                )
            )
    return tuple(observations)


def _apply_reaction_rapm_calibrator(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    calibrator: object,
) -> tuple[dict[str, dict[str, object]], bool]:
    """Apply a development-fitted reaction calibrator to locked cases."""

    scores: dict[str, dict[str, object]] = {}
    records_complete = True
    for family, case, rows in cases:
        case_id = f"{family}:{int(case.seed)}"
        row = _complete_case_row(rows)
        if not isinstance(row, Mapping):
            records_complete = False
            continue
        inference = row.get("reaction_rapm_inference", {})
        if not isinstance(inference, Mapping):
            records_complete = False
            continue
        prediction_rows = inference.get("predictions", ())
        if not isinstance(prediction_rows, Sequence):
            records_complete = False
            continue
        predictions = {
            str(prediction.get("target")): {
                key: value for key, value in dict(prediction).items() if key != "target"
            }
            for prediction in prediction_rows
            if isinstance(prediction, Mapping) and prediction.get("target") is not None
        }
        calibrated_predictions = calibrator.apply(predictions)
        truth = {
            str(edge_id): normalise_rapm_reaction_family(process)
            for edge_id, process in case.true_processes.items()
        }
        raw_candidate_edge_ids = inference.get("candidate_edge_ids", ())
        candidate_edge_ids = (
            tuple(str(edge_id) for edge_id in raw_candidate_edge_ids)
            if isinstance(raw_candidate_edge_ids, Sequence)
            else ()
        )
        calibrated_score = score_reaction_rapm_outputs(
            truth,
            calibrated_predictions,
            candidate_edge_ids=candidate_edge_ids,
        )
        metadata = dict(inference.get("metadata", {}))
        metadata["calibration"] = calibrator.to_dict()
        uncertainty = dict(metadata.get("uncertainty", {}))
        uncertainty.update(
            {
                "calibrated": True,
                "calibration_scope": "development_only",
                "calibration_status": "frozen_for_locked_test",
            }
        )
        metadata["uncertainty"] = uncertainty
        calibrated_inference = {
            "output_kind": "reaction_family",
            "metadata": metadata,
            "candidate_universe": inference.get("candidate_universe"),
            "candidate_edge_ids": list(candidate_edge_ids),
            "predictions": [
                {"target": target, **dict(prediction)}
                for target, prediction in sorted(calibrated_predictions.items())
            ],
            "truth_blind": True,
            "truth_released_for_scoring": False,
        }
        for target_row in rows:
            if target_row.get("observation_scenario", "complete") != "complete":
                continue
            target_row["reaction_rapm_calibrated_inference"] = calibrated_inference
            target_row["reaction_rapm_calibrated_scores"] = {
                "reaction": calibrated_score
            }
        scores[case_id] = {
            "family": family,
            "seed": int(case.seed),
            "reaction": calibrated_score,
        }
    return scores, records_complete


def _common_all_pairs_universe(
    rows: Sequence[Mapping[str, object]],
) -> tuple[tuple[IndependentCandidateEdge, ...], dict[str, object]]:
    """Build a truth-blind common evaluation universe for topology selection.

    Candidate generation recall is scored separately.  This all-pairs universe
    prevents the topology comparison from giving HydroSheaf and a specialist
    different opportunities to propose an edge while retaining each method's
    own candidate generator as a separate diagnostic.
    """

    nodes = tuple(
        sorted(
            {
                str(row.get("site_id"))
                for row in rows
                if row.get("site_id") is not None
            }
        )
    )
    edges = tuple(
        IndependentCandidateEdge(u, v, {"candidate_generator": "common_all_pairs_v1"})
        for u in nodes
        for v in nodes
        if u != v
    )
    payload = {
        "algorithm": "common_all_pairs_v1",
        "nodes": list(nodes),
        "edges": [list(edge.edge) for edge in edges],
        "input_channels": ["site_id"],
        "truth_blind": True,
    }
    candidate_hash = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    audit = {
        **payload,
        "version": "1.0",
        "candidate_hash": candidate_hash,
        "edge_count": len(edges),
        "candidate_universe_scope": "all_observed_node_pairs_for_selection_only",
        "truth_fields_seen": [],
    }
    return edges, audit


def _complete_case_row(
    rows: Sequence[Mapping[str, object]],
) -> Mapping[str, object] | None:
    return next(
        (
            row
            for row in rows
            if row.get("condition") == "full_sheaf"
            and row.get("observation_scenario", "complete") == "complete"
        ),
        None,
    )


def _age_selection_score(prediction: Mapping[str, object]) -> float | None:
    """Return the first finite declared age-uncertainty score."""

    for field_name in (
        "tracer_age_spread_years",
        "interval_width_years",
        "uncertainty_years",
    ):
        score = _finite_float(prediction.get(field_name))
        if score is not None:
            return score
    return None


def _specialist_calibration_observations(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[
    dict[str, list[SpecialistAgeCalibrationObservation]],
    dict[str, list[SpecialistReactionCalibrationObservation]],
]:
    """Extract development-only specialist predictions after blind inference."""

    age_observations: dict[str, list[SpecialistAgeCalibrationObservation]] = {}
    reaction_observations: dict[str, list[SpecialistReactionCalibrationObservation]] = {}
    for family, case, rows in cases:
        row = _complete_case_row(rows)
        if not isinstance(row, Mapping):
            continue
        raw_inference = row.get("specialist_baseline_inference", {})
        if not isinstance(raw_inference, Mapping):
            continue
        case_id = f"{family}:{int(case.seed)}"
        for spec in SPECIALIST_BASELINE_REGISTRY.specs():
            record = raw_inference.get(spec.name, {})
            if not isinstance(record, Mapping):
                continue
            prediction_rows = record.get("predictions", ())
            if not isinstance(prediction_rows, Sequence):
                continue
            if spec.output_kind == "age_interval":
                bucket = age_observations.setdefault(spec.name, [])
                for prediction in prediction_rows:
                    if not isinstance(prediction, Mapping):
                        continue
                    target_id = str(prediction.get("target", ""))
                    truth = case.true_ages_years.get(target_id)
                    if truth is None:
                        continue
                    required = ("mean_age_years", "age_95_low", "age_95_high")
                    if not all(field_name in prediction for field_name in required):
                        continue
                    selection_score = _age_selection_score(prediction)
                    bucket.append(
                        SpecialistAgeCalibrationObservation(
                            case_id=case_id,
                            target_id=target_id,
                            truth=float(truth),
                            estimate=float(prediction[required[0]]),
                            lower=float(prediction[required[1]]),
                            upper=float(prediction[required[2]]),
                            selection_score=selection_score,
                            selected=str(prediction.get("decision", SELECT)) != "abstain",
                        )
                    )
            elif spec.output_kind == "reaction_family":
                bucket = reaction_observations.setdefault(spec.name, [])
                for prediction in prediction_rows:
                    if not isinstance(prediction, Mapping):
                        continue
                    edge_id = str(prediction.get("target", ""))
                    logits = prediction.get("logits", prediction.get("raw_scores"))
                    if not isinstance(logits, Mapping):
                        probabilities = prediction.get("probabilities", {})
                        if isinstance(probabilities, Mapping):
                            logits = {
                                str(key): float(np.log(max(value, 1.0e-12)))
                                for key, value in (
                                    (key, _finite_float(value))
                                    for key, value in probabilities.items()
                                )
                                if value is not None
                            }
                    if not edge_id or not isinstance(logits, Mapping):
                        continue
                    truth_label = normalise_reaction_family(
                        str(case.true_processes.get(edge_id, "none"))
                    )
                    bucket.append(
                        SpecialistReactionCalibrationObservation(
                            case_id=case_id,
                            edge_id=edge_id,
                            truth_family=truth_label,
                            logits={str(key): float(value) for key, value in logits.items()},
                        )
                    )
    return age_observations, reaction_observations


def _fit_specialist_calibrators(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[
    dict[str, object],
    dict[str, object],
    dict[str, str],
]:
    age_observations, reaction_observations = _specialist_calibration_observations(cases)
    age_calibrators: dict[str, object] = {}
    reaction_calibrators: dict[str, object] = {}
    errors: dict[str, str] = {}
    for spec in SPECIALIST_BASELINE_REGISTRY.specs():
        if spec.output_kind == "age_interval":
            observations = age_observations.get(spec.name, [])
            try:
                age_calibrators[spec.name] = fit_specialist_age_calibrator(
                    observations,
                    target_coverage=0.95,
                    phase="development",
                    minimum_selection_rate=0.90,
                )
            except Exception as exc:
                errors[spec.name] = f"{type(exc).__name__}: {exc}"
        elif spec.output_kind == "reaction_family":
            observations = reaction_observations.get(spec.name, [])
            try:
                reaction_calibrators[spec.name] = fit_specialist_reaction_calibrator(
                    observations,
                    decision_threshold=0.60,
                    phase="development",
                )
            except Exception as exc:
                errors[spec.name] = f"{type(exc).__name__}: {exc}"
    return age_calibrators, reaction_calibrators, errors


def _apply_specialist_calibrators(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    age_calibrators: Mapping[str, object],
    reaction_calibrators: Mapping[str, object],
) -> tuple[dict[str, dict[str, object]], bool]:
    """Apply frozen calibrators to locked rows and score them after application."""

    locked_scores: dict[str, dict[str, object]] = {}
    records_complete = True
    for family, case, rows in cases:
        case_id = f"{family}:{int(case.seed)}"
        for row in rows:
            raw_inference = row.get("specialist_baseline_inference", {})
            if not isinstance(raw_inference, Mapping):
                records_complete = False
                continue
            calibrated_inference: dict[str, dict[str, object]] = {}
            calibrated_scores: dict[str, dict[str, object]] = {}
            for name, calibrator in {**age_calibrators, **reaction_calibrators}.items():
                raw_record = raw_inference.get(name)
                if not isinstance(raw_record, Mapping):
                    records_complete = False
                    continue
                prediction_rows = raw_record.get("predictions", ())
                if not isinstance(prediction_rows, Sequence):
                    records_complete = False
                    continue
                predictions = {
                    str(prediction.get("target")): dict(prediction)
                    for prediction in prediction_rows
                    if isinstance(prediction, Mapping) and prediction.get("target") is not None
                }
                calibrated_predictions = calibrator.apply(predictions)
                metadata = raw_record.get("metadata", {})
                metadata = dict(metadata) if isinstance(metadata, Mapping) else {}
                metadata["calibration"] = calibrator.to_dict()
                uncertainty = metadata.get("uncertainty", {})
                uncertainty = dict(uncertainty) if isinstance(uncertainty, Mapping) else {}
                uncertainty.update(
                    {
                        "calibrated": True,
                        "calibration_scope": "development_only",
                        "calibration_status": "frozen_for_locked_test",
                    }
                )
                metadata["uncertainty"] = uncertainty
                calibrated_inference[name] = {
                    "output_kind": raw_record.get("output_kind"),
                    "metadata": metadata,
                    "predictions": [
                        {"target": target, **dict(prediction)}
                        for target, prediction in sorted(calibrated_predictions.items())
                    ],
                }
                if name in age_calibrators:
                    calibrated_scores[name] = {
                        "metadata": metadata,
                        "age": score_calibrated_specialist_age(
                            case.true_ages_years,
                            calibrated_predictions,
                        ),
                    }
                else:
                    truth_by_edge = {
                        target: normalise_reaction_family(
                            str(case.true_processes.get(target, "none"))
                        )
                        for target in predictions
                    }
                    calibrated_scores[name] = {
                        "metadata": metadata,
                        "reaction": score_calibrated_specialist_reaction(
                            truth_by_edge,
                            calibrated_predictions,
                        ),
                    }
            row["specialist_calibrated_baseline_inference"] = calibrated_inference
            row["specialist_calibrated_baselines"] = calibrated_scores
            row["specialist_calibration_applied"] = True
            row["specialist_calibration_case_id"] = case_id
        locked_scores[case_id] = {
            "family": family,
            "seed": int(case.seed),
            "complete_case": next(
                (
                    row.get("specialist_calibrated_baselines", {})
                    for row in rows
                    if row.get("condition") == "full_sheaf"
                    and row.get("observation_scenario", "complete") == "complete"
                ),
                {},
            ),
        }
    return locked_scores, records_complete


def _specialist_record_complete(row: Mapping[str, object]) -> bool:
    """Fail closed unless the independent comparator record is structurally complete."""

    universe = row.get("specialist_candidate_universe")
    inference = row.get("specialist_baseline_inference")
    scored = row.get("specialist_baselines")
    if not (
        bool(row.get("specialist_candidate_generation_truth_blind"))
        and isinstance(universe, Mapping)
        and isinstance(inference, Mapping)
        and isinstance(scored, Mapping)
    ):
        return False
    if not (
        bool(universe.get("truth_blind"))
        and not universe.get("truth_fields_seen")
        and str(universe.get("algorithm", ""))
        and str(universe.get("version", ""))
        and str(universe.get("input_hash", ""))
        and str(universe.get("candidate_hash", ""))
        and isinstance(universe.get("edges"), Sequence)
    ):
        return False
    expected = set(BASELINE_REGISTRY.names()) | set(SPECIALIST_BASELINE_REGISTRY.names())
    return expected.issubset(set(inference)) and expected.issubset(set(scored))


def _hydrosheaf_probability_calibration(
    reference_edges: Sequence[Sequence[object]],
    candidates: Sequence[object],
) -> dict[str, object]:
    """Score the pipeline's emitted candidate probabilities, without tuning."""

    reference = {tuple(str(value) for value in edge[:2]) for edge in reference_edges}
    labels: list[float] = []
    probabilities: list[object] = []
    for edge in candidates:
        labels.append(1.0 if _edge_pair(edge) in reference else 0.0)
        attrs = _edge_attrs(edge)
        probabilities.append(
            attrs.get("p_uv", attrs.get("flow_probability", attrs.get("edge_confidence")))
        )
    return dict(probability_metrics(labels, probabilities, n_bins=10))


def _hydrosheaf_probability_predictions(
    candidates: Sequence[object],
) -> list[dict[str, object]]:
    """Serialise native candidate probabilities before truth is available."""

    predictions: list[dict[str, object]] = []
    for edge in candidates:
        attrs = _edge_attrs(edge)
        probability = _finite_float(
            attrs.get("p_uv", attrs.get("flow_probability", attrs.get("edge_confidence")))
        )
        probability = 0.5 if probability is None else min(1.0, max(0.0, probability))
        source, target = _edge_pair(edge)
        predictions.append(
            {
                "edge": [source, target],
                "probability": probability,
                "evidence_channel": "hydrosheaf_native_candidate_probability",
            }
        )
    return predictions


def _topology_model_observations(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    *,
    phase: str,
) -> tuple[DiscreteModelObservation, ...]:
    """Build a case-blocked, conditional candidate-universe model matrix."""

    observations: list[DiscreteModelObservation] = []
    baseline_names = tuple(spec.name for spec in BASELINE_REGISTRY.specs())
    for family, case, rows in cases:
        row = next(
            (
                item
                for item in rows
                if item.get("condition") == "full_sheaf"
                and item.get("observation_scenario", "complete") == "complete"
            ),
            None,
        )
        if not isinstance(row, Mapping):
            continue
        native_rows = row.get("hydrosheaf_topology_predictions", ())
        baseline_inference = row.get("baseline_inference", {})
        if not isinstance(native_rows, Sequence) or not isinstance(baseline_inference, Mapping):
            continue
        native_by_edge = {
            tuple(str(value) for value in item.get("edge", ())[:2]): item
            for item in native_rows
            if isinstance(item, Mapping)
            and isinstance(item.get("edge"), Sequence)
            and len(item.get("edge", ())) >= 2
        }
        baseline_by_name: dict[str, dict[tuple[str, str], Mapping[str, object]]] = {}
        for name in baseline_names:
            record = baseline_inference.get(name, {})
            prediction_rows = record.get("predictions", []) if isinstance(record, Mapping) else []
            baseline_by_name[name] = {
                tuple(str(value) for value in item.get("edge", ())[:2]): item
                for item in prediction_rows
                if isinstance(item, Mapping)
                and isinstance(item.get("edge"), Sequence)
                and len(item.get("edge", ())) >= 2
            }
        common_edges = set(native_by_edge)
        common_edges &= set.intersection(
            *(set(by_edge) for by_edge in baseline_by_name.values())
        ) if baseline_by_name else common_edges
        reference = {
            tuple(str(value) for value in edge[:2])
            for edge in case.true_edges
        }
        for edge in sorted(common_edges):
            probabilities: dict[str, dict[str, float]] = {}
            native_probability = _finite_float(native_by_edge[edge].get("probability"))
            native_probability = 0.5 if native_probability is None else native_probability
            probabilities["hydrosheaf_native"] = {
                "present": native_probability,
                "absent": 1.0 - native_probability,
            }
            for name in baseline_names:
                probability = _finite_float(baseline_by_name[name][edge].get("probability"))
                probability = 0.5 if probability is None else min(1.0, max(0.0, probability))
                probabilities[name] = {
                    "present": probability,
                    "absent": 1.0 - probability,
                }
            observations.append(
                DiscreteModelObservation(
                    case_id=f"{phase}:{family}:{case.seed}",
                    target_id=f"{edge[0]}->{edge[1]}",
                    truth="present" if edge in reference else "absent",
                    predictions=probabilities,
                    phase=phase,
                )
            )
    return tuple(observations)


def _measurement_scenarios(
    prior: Mapping[str, float],
    *,
    informativeness: float = 1.0,
) -> dict[str, ScenarioBelief]:
    """Return a pre-measurement scenario envelope for one action.

    ``informativeness`` is a truth-blind observation-space discriminability
    score.  It controls how strongly the declared measurement likelihoods
    separate support from no-support; it is not a realised label or hidden
    reference-edge lookup.
    """

    signal = min(1.0, max(0.0, float(informativeness)))
    nominal_support_present = 0.50 + 0.35 * signal
    nominal_support_absent = 0.50 - 0.25 * signal
    discrepant_support_present = 0.50 + 0.15 * signal
    discrepant_support_absent = 0.50 - 0.05 * signal

    return {
        "nominal": ScenarioBelief(
            prior=prior,
            outcome_likelihoods={
                "support": {
                    "present": nominal_support_present,
                    "absent": nominal_support_absent,
                },
                "no_support": {
                    "present": 1.0 - nominal_support_present,
                    "absent": 1.0 - nominal_support_absent,
                },
            },
        ),
        "discrepant": ScenarioBelief(
            prior=prior,
            outcome_likelihoods={
                "support": {
                    "present": discrepant_support_present,
                    "absent": discrepant_support_absent,
                },
                "no_support": {
                    "present": 1.0 - discrepant_support_present,
                    "absent": 1.0 - discrepant_support_absent,
                },
            },
        ),
    }


PROSPECTIVE_DECLARED_UTILITY_MODEL: dict[str, object] = {
    "version": "bayes_edge_classification_accuracy_v1",
    "measurement_separation": {
        "head": 0.55,
        "tracer_age": 0.45,
        "chemistry": 0.35,
    },
    "cost_penalty": 0.05,
    "score_definition": (
        "expected posterior Bayes classification accuracy minus cost penalty; "
        "not realised measurement evidence"
    ),
    "frozen_before_locked_scoring": True,
}
PROSPECTIVE_DECLARED_UTILITY_MODEL_HASH = hashlib.sha256(
    json.dumps(
        PROSPECTIVE_DECLARED_UTILITY_MODEL,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
).hexdigest()


def _declared_utility_contract_record() -> dict[str, object]:
    """Return the frozen, hashable declaration used by the Hydro policy."""

    return {
        **PROSPECTIVE_DECLARED_UTILITY_MODEL,
        "measurement_separation": dict(
            PROSPECTIVE_DECLARED_UTILITY_MODEL["measurement_separation"]
        ),
        "hash": PROSPECTIVE_DECLARED_UTILITY_MODEL_HASH,
    }


def _declared_expected_classification_utility(
    prior_probability: float,
    measurement: str,
    cost: float,
) -> float:
    """Return the pre-measurement Bayes classification utility declaration.

    The declaration is deliberately a different functional object from the
    post-selection simulator's evidence benefit.  It uses only the prior,
    a modality separation constant, and cost; it never consumes a realised
    signal, selected truth state, or generator-specific observation quality.
    """

    separation_by_measurement = PROSPECTIVE_DECLARED_UTILITY_MODEL[
        "measurement_separation"
    ]
    if not isinstance(separation_by_measurement, Mapping):
        raise TypeError("declared measurement separation must be a mapping")
    separation = float(separation_by_measurement[measurement])
    probability = min(0.95, max(0.05, float(prior_probability)))
    present_support = 0.50 + 0.50 * separation
    absent_support = 0.50 - 0.50 * separation
    support_mass = (
        probability * present_support + (1.0 - probability) * absent_support
    )
    no_support_mass = 1.0 - support_mass
    correct_mass = max(
        probability * present_support,
        (1.0 - probability) * absent_support,
    ) + max(
        probability * (1.0 - present_support),
        (1.0 - probability) * (1.0 - absent_support),
    )
    if not math.isfinite(correct_mass) or support_mass <= 0.0 or no_support_mass <= 0.0:
        raise ValueError("declared classification utility is not finite")
    cost_penalty = float(PROSPECTIVE_DECLARED_UTILITY_MODEL["cost_penalty"])
    return float(correct_mass - cost_penalty * float(cost))


def _measurement_actions(
    candidates: Sequence[object],
    *,
    probability_by_edge: Mapping[tuple[str, str], float] | None = None,
    informativeness_by_action: Mapping[str, float] | None = None,
) -> list[CandidateMeasurementAction]:
    """Build the declared truth-blind action universe used by all policies.

    ``probability_by_edge`` is an optional HydroSheaf prior for a common
    all-pairs universe.  It is an observation-only prior, never a released
    reference label; absent native edges retain the neutral 0.5 prior.
    """

    actions: list[CandidateMeasurementAction] = []
    edge_probabilities = probability_by_edge or {}
    action_informativeness = informativeness_by_action or {}
    # The policy knows only the prior edge probability, measurement modality,
    # and cost.  The held-out scorer evaluates the choice against separately
    # generated signals and a different evidence-benefit functional form.
    for edge in sorted(candidates, key=lambda item: _edge_pair(item)):
        source, target = _edge_pair(edge)
        attrs = getattr(edge, "attrs", {})
        raw_probability = _finite_float(
            attrs.get("p_uv", attrs.get("flow_probability", attrs.get("edge_confidence")))
        )
        probability = edge_probabilities.get((source, target), raw_probability)
        probability = min(0.95, max(0.05, probability if probability is not None else 0.5))
        prior = {"present": probability, "absent": 1.0 - probability}
        for measurement, cost in (("head", 1.0), ("tracer_age", 2.0), ("chemistry", 3.0)):
            action_id = f"{measurement}:{source}->{target}"
            informativeness = min(
                1.0,
                max(0.0, float(action_informativeness.get(action_id, 0.0))),
            )
            declared_expected_utility = _declared_expected_classification_utility(
                probability,
                measurement,
                cost,
            )
            actions.append(
                CandidateMeasurementAction(
                    action_id=action_id,
                    cost=cost,
                    scenarios=_measurement_scenarios(
                        prior,
                        informativeness=informativeness,
                    ),
                    metadata={
                        "measurement": measurement,
                        "target": f"{source}->{target}",
                        "prior_probability": probability,
                        "observation_informativeness": informativeness,
                        "declared_expected_utility": declared_expected_utility,
                        "declared_utility_model": PROSPECTIVE_DECLARED_UTILITY_MODEL[
                            "version"
                        ],
                        "declared_utility_model_hash": (
                            PROSPECTIVE_DECLARED_UTILITY_MODEL_HASH
                        ),
                        "model": "predeclared_binary_edge_support_v1",
                    },
                )
            )
    return actions


def _prospective_action_records(
    actions: Sequence[CandidateMeasurementAction],
) -> list[dict[str, object]]:
    """Serialise action declarations without serialising hidden truth."""

    return [
        {
            "action_id": action.action_id,
            "cost": float(action.cost),
            "feasible": bool(action.feasible),
            "prior_probability": float(
                dict(action.metadata).get("prior_probability", 0.5)
            ),
            "metadata": dict(action.metadata),
        }
        for action in sorted(actions, key=lambda item: item.action_id)
    ]


def _rebuild_prospective_actions(
    records: Sequence[Mapping[str, object]],
) -> tuple[CandidateMeasurementAction, ...]:
    """Rebuild the exact declared action objects from a locked audit record."""

    actions: list[CandidateMeasurementAction] = []
    for record in records:
        action_id = str(record.get("action_id", "")).strip()
        if not action_id or ":" not in action_id:
            continue
        raw_prior = _finite_float(record.get("prior_probability"))
        if raw_prior is None:
            metadata = record.get("metadata", {})
            raw_prior = _finite_float(
                metadata.get("prior_probability")
                if isinstance(metadata, Mapping)
                else None
            )
        prior_probability = min(0.95, max(0.05, raw_prior if raw_prior is not None else 0.5))
        metadata = record.get("metadata", {})
        metadata = dict(metadata) if isinstance(metadata, Mapping) else {}
        informativeness = _finite_float(
            metadata.get("observation_informativeness")
        )
        actions.append(
            CandidateMeasurementAction(
                action_id=action_id,
                cost=float(record.get("cost", 1.0)),
                feasible=bool(record.get("feasible", True)),
                scenarios=_measurement_scenarios(
                    {"present": prior_probability, "absent": 1.0 - prior_probability},
                    informativeness=(
                        0.0 if informativeness is None else informativeness
                    ),
                ),
                metadata=metadata,
            )
        )
    return tuple(sorted(actions, key=lambda item: item.action_id))


def _pair_observation_quality(
    rows: Sequence[Mapping[str, object]],
    source: str,
    target: str,
    measurement: str,
) -> float:
    """Compute a generator-realised, truth-blind discriminability score.

    This is deliberately an observation-space utility, not a hidden-truth
    likelihood.  The two node observations are already available before a
    new action is selected; truth enters only when the scorer assigns the
    action-specific present/absent state after selection.
    """

    by_id = {
        str(row.get("site_id")): row
        for row in rows
        if row.get("site_id") is not None
    }
    upstream = by_id.get(source, {})
    downstream = by_id.get(target, {})
    z_values: list[float] = []
    if measurement == "head":
        before = _finite_float(upstream.get("head_meas"))
        after = _finite_float(downstream.get("head_meas"))
        if before is not None and after is not None:
            sigma = max(
                0.05,
                _finite_float(upstream.get("head_sigma_m")) or 0.06,
                _finite_float(downstream.get("head_sigma_m")) or 0.06,
            )
            z_values.append(abs(before - after) / sigma)
    elif measurement == "tracer_age":
        tracer_scales = {
            "tritium_TU": 0.15,
            "argon39_pmc": 2.0,
            "c14_pmc": 0.8,
            "cfc11_pptv": 0.25,
            "cfc12_pptv": 0.25,
            "cfc113_pptv": 0.15,
            "sf6_pptv": 0.05,
        }
        for field, scale in tracer_scales.items():
            before = _finite_float(upstream.get(field))
            after = _finite_float(downstream.get(field))
            if before is not None and after is not None:
                z_values.append(abs(before - after) / max(scale, 1.0e-9))
    elif measurement == "chemistry":
        chemistry_fields = (
            "Ca", "Cl", "F", "Fe", "HCO3", "K", "Mg", "Na", "NO3", "PO4", "SO4", "SiO2", "pH"
        )
        for field in chemistry_fields:
            before = _finite_float(upstream.get(field))
            after = _finite_float(downstream.get(field))
            if before is None or after is None:
                continue
            if field == "pH":
                scale = 0.10
                z_values.append(abs(before - after) / scale)
            elif before >= 0.0 and after >= 0.0:
                z_values.append(abs(math.log1p(after) - math.log1p(before)) / 0.10)
    if not z_values:
        return 0.0
    median_z = float(np.median(np.asarray(z_values, dtype=float)))
    return float(np.clip(1.0 - math.exp(-median_z / 2.0), 0.0, 1.0))


def _prospective_realised_outcomes(
    rows: Sequence[Mapping[str, object]],
    actions: Sequence[CandidateMeasurementAction],
    reference_edges: Sequence[Sequence[object]],
    *,
    seed: int,
) -> tuple[dict[str, dict[str, float]], dict[str, str], dict[str, object]]:
    """Build a post-selection table from a state-conditioned measurement model.

    This simulator is deliberately independent of the observation-space
    discriminability used by the selector.  It generates counterfactual
    signals under each hidden edge state from a fixed measurement model,
    action-specific noise, censoring, and detection limits.  The hidden state
    and realised outcome are written only after all selectors have returned.
    This remains a controlled synthetic measurement simulator, not a field
    utility model or a second MODFLOW forward solve.
    """

    reference = {
        (str(edge[0]), str(edge[1]))
        for edge in reference_edges
        if len(edge) >= 2
    }
    benefits: dict[str, dict[str, float]] = {}
    states: dict[str, str] = {}
    outcome_records: dict[str, dict[str, object]] = {}
    state_mean_by_measurement = {
        "head": {"present": 0.82, "absent": 0.18},
        "tracer_age": {"present": 0.78, "absent": 0.22},
        "chemistry": {"present": 0.74, "absent": 0.26},
    }
    noise_sd_by_measurement = {
        "head": 0.16,
        "tracer_age": 0.20,
        "chemistry": 0.22,
    }
    detection_limit_by_measurement = {
        "head": 0.15,
        "tracer_age": 0.15,
        "chemistry": 0.15,
    }
    for action in actions:
        measurement, _, target = action.action_id.partition(":")
        if "->" not in target:
            continue
        source, destination = target.split("->", 1)
        action_seed = int.from_bytes(
            hashlib.sha256(
                f"prospective-independent-v2|{int(seed)}|{action.action_id}".encode(
                    "utf-8"
                )
            ).digest()[:8],
            "little",
        )
        rng = np.random.default_rng(action_seed)
        noise_sd = noise_sd_by_measurement.get(measurement, 0.20)
        detection_limit = detection_limit_by_measurement.get(measurement, 0.15)
        target_state = (
            "present" if (source, destination) in reference else "absent"
        )
        states[action.action_id] = target_state
        measurement_means = state_mean_by_measurement.get(
            measurement,
            {"present": 0.75, "absent": 0.25},
        )
        latent_noise = float(rng.normal(0.0, noise_sd * 0.50))
        measurement_noise = float(rng.normal(0.0, noise_sd))
        signal_by_state: dict[str, float] = {}
        detected_by_state: dict[str, bool] = {}
        evidence_by_state: dict[str, float] = {}
        for state in ("present", "absent"):
            latent_signal = float(
                np.clip(measurement_means[state] + latent_noise, 0.0, 1.0)
            )
            measured_signal = float(
                np.clip(latent_signal + measurement_noise, 0.0, 1.0)
            )
            detected = measured_signal >= detection_limit
            evidence = (
                float(np.clip(abs(measured_signal - 0.5) * 2.0, 0.0, 1.0))
                if detected
                else 0.0
            )
            signal_by_state[state] = measured_signal
            detected_by_state[state] = detected
            evidence_by_state[state] = evidence
        benefits[action.action_id] = {
            "present": float(1.0 + evidence_by_state["present"]),
            "absent": float(0.60 * evidence_by_state["absent"]),
        }
        realised_signal = signal_by_state[target_state]
        outcome_records[action.action_id] = {
            "measurement": measurement,
            "state_conditioned": True,
            "state_means": dict(measurement_means),
            "latent_noise": latent_noise,
            "measurement_noise": measurement_noise,
            "signal_by_state": signal_by_state,
            "detected_by_state": detected_by_state,
            "evidence_by_state": evidence_by_state,
            "realised_state": target_state,
            "measured_signal": realised_signal,
            "noise_sd": noise_sd,
            "detection_limit": detection_limit,
            "detected": detected_by_state[target_state],
            "rng_seed": action_seed,
        }
    audit = {
        "rule": "independent_state_conditioned_measurement_v2",
        "measurement_outcomes": outcome_records,
        "state_conditioned": True,
        "reuses_selection_quality": False,
        "noise_model": noise_sd_by_measurement,
        "detection_limits": detection_limit_by_measurement,
        "seed": int(seed),
        "truth_release": "action_specific_edge_state_only_after_selector_return",
        "not_a_field_utility_model": True,
    }
    return benefits, states, audit


def _specialist_action_scores(
    actions: Sequence[CandidateMeasurementAction],
    baseline_inference: Mapping[str, Mapping[str, object]],
) -> dict[str, float]:
    """Map the independent topology comparator onto the common action set."""

    preferred = baseline_inference.get("edge_local_topology_support", {})
    prediction_rows = preferred.get("predictions", ()) if isinstance(preferred, Mapping) else ()
    by_edge: dict[tuple[str, str], float] = {}
    if isinstance(prediction_rows, Sequence):
        for prediction in prediction_rows:
            if not isinstance(prediction, Mapping):
                continue
            edge = prediction.get("edge", ())
            if not isinstance(edge, Sequence) or len(edge) < 2:
                continue
            probability = _finite_float(prediction.get("probability"))
            if probability is not None:
                by_edge[(str(edge[0]), str(edge[1]))] = min(1.0, max(0.0, probability))
    return {
        action.action_id: by_edge.get(
            tuple(action.action_id.split(":", 1)[1].split("->", 1)),
            0.5,
        ) / max(float(action.cost), 1.0e-12)
        for action in actions
    }


def _decision_policy(candidates: Sequence[object]) -> dict[str, object]:
    """Select a robust one-step action from declared, truth-blind models."""

    actions = _measurement_actions(candidates)
    # The prospective benchmark evaluates the selected action against random
    # and specialist policies.  A zero threshold is the declared screening
    # rule here: actions with non-positive robust information gain still
    # abstain, while a positive but small gain is not silently discarded before
    # the independent outcome benchmark can score it.
    report = select_next_measurement(actions, min_utility_per_cost=0.0)
    return report.to_dict()


def _discrepancy_calibration_bundle(
    development_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    locked_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    discrepancy_reports: Sequence[Mapping[str, object]],
) -> tuple[
    object | None,
    dict[str, object],
    list[dict[str, object]],
    bool,
    str | None,
]:
    """Fit discrepancy dilation on development cases and score locked cases."""

    cases = {
        (family, int(case.seed)): case
        for family, case, _rows in (*development_cases, *locked_cases)
    }
    development_observations: list[DiscrepancyCalibrationObservation] = []
    locked_observations: list[DiscrepancyCalibrationObservation] = []
    development_reports: dict[str, Mapping[str, object]] = {}
    locked_reports: dict[str, Mapping[str, object]] = {}
    for report in discrepancy_reports:
        phase = str(report.get("phase", ""))
        family = str(report.get("family", ""))
        try:
            seed = int(report.get("seed"))
        except (TypeError, ValueError):
            continue
        case = cases.get((family, seed))
        target = str(report.get("target", ""))
        if case is None or not target.startswith("age:"):
            continue
        node = target.split(":", 1)[1]
        truth = case.true_ages_years.get(node)
        if truth is None:
            continue
        key = f"{family}:{seed}|{target}"
        try:
            observation = DiscrepancyCalibrationObservation(
                case_id=f"{family}:{seed}",
                target_id=key,
                truth=float(truth),
                estimate=float(report["estimate"]),
                lower=float(report["lower"]),
                upper=float(report["upper"]),
                phase=phase,
            )
        except (KeyError, TypeError, ValueError):
            continue
        if phase == "development":
            development_observations.append(observation)
            development_reports[key] = report
        elif phase == "locked_test":
            locked_observations.append(observation)
            locked_reports[key] = report

    if not development_observations or not locked_observations:
        return (
            None,
            {"status": "ABSTAIN_MISSING_DISCREPANCY_CALIBRATION_DATA"},
            [dict(report) for report in discrepancy_reports],
            False,
            "development and locked discrepancy observations are required",
        )
    try:
        calibrator = fit_discrepancy_calibrator(
            development_observations,
            target_coverage=0.95,
            phase="development",
        )
        calibrated_locked = {
            key: apply_discrepancy_calibrator_to_record(calibrator, report)
            for key, report in locked_reports.items()
        }
        raw_locked = {key: dict(report) for key, report in locked_reports.items()}
        scores = score_locked_discrepancy(
            locked_observations,
            raw_locked,
            calibrated_locked,
            target_coverage=calibrator.target_coverage,
        )
        calibrated_records: list[dict[str, object]] = []
        for report in discrepancy_reports:
            family = str(report.get("family", ""))
            try:
                seed = int(report.get("seed"))
            except (TypeError, ValueError):
                calibrated_records.append(dict(report))
                continue
            key = f"{family}:{seed}|{report.get('target', '')}"
            if key in development_reports or key in locked_reports:
                calibrated_records.append(
                    apply_discrepancy_calibrator_to_record(calibrator, report)
                )
            else:
                calibrated_records.append(dict(report))
        return calibrator, scores, calibrated_records, True, None
    except Exception as exc:
        return (
            None,
            {"status": "ABSTAIN_DISCREPANCY_CALIBRATION_ERROR", "error": f"{type(exc).__name__}: {exc}"},
            [dict(report) for report in discrepancy_reports],
            False,
            f"{type(exc).__name__}: {exc}",
        )


def _age_performance_payload(
    development_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    locked_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[dict[str, object], dict[str, dict[str, object]], tuple[str, ...], tuple[str, ...]]:
    """Collect unique locked age-specialist predictions for the age gate."""

    baseline_name = "independent_competence_matched_age_specialist"
    truth: dict[str, object] = {}
    predictions: dict[str, dict[str, object]] = {}
    metadata: dict[str, dict[str, object]] = {}
    development_ids = tuple(
        f"{family}:{int(case.seed)}" for family, case, _rows in development_cases
    )
    held_out_ids = tuple(
        f"{family}:{int(case.seed)}" for family, case, _rows in locked_cases
    )
    for family, case, rows in locked_cases:
        case_id = f"{family}:{int(case.seed)}"
        row = next(
            (
                item
                for item in rows
                if item.get("condition") == "full_sheaf"
                and item.get("observation_scenario", "complete") == "complete"
            ),
            None,
        )
        if not isinstance(row, Mapping):
            continue
        inference = row.get("specialist_calibrated_baseline_inference", {})
        if not isinstance(inference, Mapping):
            continue
        record = inference.get(baseline_name)
        if not isinstance(record, Mapping):
            continue
        raw_predictions = record.get("predictions", ())
        if not isinstance(raw_predictions, Sequence):
            continue
        for raw_prediction in raw_predictions:
            if not isinstance(raw_prediction, Mapping):
                continue
            target = raw_prediction.get("target")
            if target is None or str(target) not in case.true_ages_years:
                continue
            unique_target = f"{case_id}|{target}"
            prediction = dict(raw_prediction)
            predictions[unique_target] = prediction
            truth[unique_target] = case.true_ages_years[str(target)]
            metadata[unique_target] = {
                "case_id": case_id,
                "family": family,
                "generator_family": family,
                "mechanism": str(case.provenance.get("age_mechanism", family)),
                "missingness": str(row.get("observation_scenario", "complete")),
            }
    payload = {
        "baseline": baseline_name,
        "truth": truth,
        "predictions": predictions,
        "metadata": metadata,
        "truth_blind_inference": True,
        "calibration_fit_scope": "development_only",
    }
    return payload, metadata, development_ids, held_out_ids


def _age_point_interval_summary(
    truth: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]],
) -> dict[str, float | int]:
    """Return a small truth-release summary for the fixed age comparator."""

    errors: list[float] = []
    widths: list[float] = []
    covered = 0
    n_abstain = 0
    n_missing = 0
    n_invalid = 0
    for target, raw_truth in truth.items():
        prediction = predictions.get(str(target))
        if not isinstance(prediction, Mapping):
            n_missing += 1
            continue
        if str(prediction.get("decision", SELECT)) == "abstain":
            n_abstain += 1
            continue
        estimate = _finite_float(
            prediction.get("calibrated_mean_age_years", prediction.get("mean_age_years"))
        )
        lower = _finite_float(
            prediction.get("calibrated_age_low", prediction.get("age_95_low"))
        )
        upper = _finite_float(
            prediction.get("calibrated_age_high", prediction.get("age_95_high"))
        )
        reference = _finite_float(raw_truth)
        if None in (estimate, lower, upper, reference) or lower > upper:
            n_invalid += 1
            continue
        errors.append(abs(float(estimate) - float(reference)))
        widths.append(float(upper) - float(lower))
        covered += int(float(lower) <= float(reference) <= float(upper))
    count = len(errors)
    return {
        "n_truth": len(truth),
        "n_selected": count,
        "n_abstain": n_abstain,
        "n_missing": n_missing,
        "n_invalid": n_invalid,
        "selection_rate": count / len(truth) if truth else float("nan"),
        "mae_years": sum(errors) / count if count else float("nan"),
        "mean_interval_width_years": sum(widths) / len(widths) if widths else float("nan"),
        "coverage": covered / count if count else float("nan"),
        "coverage_including_abstention": covered / len(truth) if truth else float("nan"),
    }


def _age_performance_bundle(
    development_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    locked_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> tuple[dict[str, object], dict[str, object]]:
    """Evaluate the calibrated age comparator and its fixed-history control."""

    payload, metadata, development_ids, held_out_ids = _age_performance_payload(
        development_cases,
        locked_cases,
    )
    report = evaluate_age_performance(
        payload["truth"],  # type: ignore[arg-type]
        payload["predictions"],  # type: ignore[arg-type]
        metadata=metadata,
        thresholds=AgePerformanceThresholds(),
        use_calibrated=True,
        development_case_ids=development_ids,
        held_out_case_ids=held_out_ids,
    )
    baseline_truth: dict[str, object] = {}
    baseline_predictions: dict[str, dict[str, object]] = {}
    baseline_name = "multitracer_atmospheric_history_age"
    baseline_calibration_observations: list[SpecialistAgeCalibrationObservation] = []
    for family, case, rows in development_cases:
        row = _complete_case_row(rows)
        if not isinstance(row, Mapping):
            continue
        raw_inference = row.get("specialist_baseline_inference", {})
        if not isinstance(raw_inference, Mapping):
            continue
        record = raw_inference.get(baseline_name, {})
        if not isinstance(record, Mapping):
            continue
        raw_predictions = record.get("predictions", ())
        if not isinstance(raw_predictions, Sequence):
            continue
        case_id = f"{family}:{int(case.seed)}"
        for raw_prediction in raw_predictions:
            if not isinstance(raw_prediction, Mapping):
                continue
            target = raw_prediction.get("target")
            if target is None or str(target) not in case.true_ages_years:
                continue
            required = ("mean_age_years", "age_95_low", "age_95_high")
            if not all(field in raw_prediction for field in required):
                continue
            selection_score = _age_selection_score(raw_prediction)
            baseline_calibration_observations.append(
                SpecialistAgeCalibrationObservation(
                    case_id=case_id,
                    target_id=str(target),
                    truth=float(case.true_ages_years[str(target)]),
                    estimate=float(raw_prediction[required[0]]),
                    lower=float(raw_prediction[required[1]]),
                    upper=float(raw_prediction[required[2]]),
                    selection_score=selection_score,
                    selected=str(raw_prediction.get("decision", SELECT)) != "abstain",
                )
            )
    baseline_calibrator = fit_specialist_age_calibrator(
        baseline_calibration_observations,
        target_coverage=0.95,
        phase="development",
        minimum_selection_rate=0.90,
    )
    for family, case, rows in locked_cases:
        case_id = f"{family}:{int(case.seed)}"
        row = _complete_case_row(rows)
        if not isinstance(row, Mapping):
            continue
        raw_inference = row.get("specialist_baseline_inference", {})
        if not isinstance(raw_inference, Mapping):
            continue
        record = raw_inference.get(baseline_name, {})
        if not isinstance(record, Mapping):
            continue
        raw_predictions = record.get("predictions", ())
        if not isinstance(raw_predictions, Sequence):
            continue
        for raw_prediction in raw_predictions:
            if not isinstance(raw_prediction, Mapping):
                continue
            target = raw_prediction.get("target")
            if target is None or str(target) not in case.true_ages_years:
                continue
            key = f"{case_id}|{target}"
            baseline_truth[key] = case.true_ages_years[str(target)]
            calibrated_baseline = baseline_calibrator.apply(
                {str(target): dict(raw_prediction)}
            ).get(str(target))
            if isinstance(calibrated_baseline, Mapping):
                baseline_predictions[key] = dict(calibrated_baseline)
    baseline = _age_point_interval_summary(baseline_truth, baseline_predictions)
    baseline["calibrated"] = True
    baseline["calibration_fit_count"] = len(baseline_calibration_observations)
    baseline["calibration_case_count"] = len(baseline_calibrator.case_ids)
    baseline["calibrator"] = baseline_calibrator.to_dict()
    age_metrics = dict(report.get("metrics", {}))
    specialist_mae = _finite_float(age_metrics.get("mae_years"))
    baseline_mae = _finite_float(baseline.get("mae_years"))
    specialist_width = _finite_float(age_metrics.get("mean_interval_width_years"))
    baseline_width = _finite_float(baseline.get("mean_interval_width_years"))
    age_metrics.update(
        {
            "acceptance_rate": age_metrics.get("selection_rate"),
            "baseline_acceptance_rate": baseline.get("selection_rate"),
            "baseline_coverage_including_abstention": baseline.get(
                "coverage_including_abstention"
            ),
            "specialist_mae": specialist_mae,
            "baseline_mae": baseline_mae,
            "mean_width": specialist_width,
            "baseline_mean_width": baseline_width,
            "relative_width": (
                specialist_width / baseline_width
                if specialist_width is not None
                and baseline_width is not None
                and baseline_width > 0.0
                else None
            ),
            "selective_risk_years": age_metrics.get("selective_risk_at_acceptance_years"),
            "selective_risk": age_metrics.get("selective_risk_at_acceptance_years"),
            "held_out_generators": len(
                {
                    str(item.get("family"))
                    for item in metadata.values()
                    if isinstance(item, Mapping)
                }
            ) >= len(REQUIRED_GENERATOR_FAMILIES),
            "held_out_mechanisms": len(
                {
                    str(item.get("mechanism"))
                    for item in metadata.values()
                    if isinstance(item, Mapping)
                }
            ) >= 2,
            "competence_matched_baseline": bool(baseline_predictions),
            "calibrated": report.get("audit", {}).get("calibration_statuses_observed") == [
                "development_fitted"
            ]
            if isinstance(report.get("audit"), Mapping)
            else False,
            "family_stratified": bool(report.get("stratified"))
            and all(
                str(record.get("status")) == "PASSABLE_DENOMINATOR"
                for record in report.get("stratified", {}).get("family", {}).values()
                if isinstance(record, Mapping)
            ),
            "baseline_noninferior": (
                specialist_mae is not None
                and baseline_mae is not None
                and specialist_mae <= baseline_mae
            ),
            "acceptance_noninferior": (
                _finite_float(age_metrics.get("selection_rate")) is not None
                and _finite_float(baseline.get("selection_rate")) is not None
                and float(age_metrics["selection_rate"]) >= float(baseline["selection_rate"])
            ),
        }
    )
    return report, {
        "metrics": age_metrics,
        "baseline": baseline,
        "payload": payload,
    }


def _reaction_rapm_performance_bundle(
    locked_cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    *,
    calibrator: object | None,
) -> dict[str, object]:
    """Aggregate calibrated RAPM records across held-out generator families."""

    scores: list[Mapping[str, object]] = []
    baseline_records: dict[str, list[tuple[float, int]]] = {}
    baseline_case_counts: dict[str, int] = {}
    baseline_incomplete_counts: dict[str, int] = {}
    family_counts: dict[str, int] = {}
    classwise_complete = calibrator is not None
    candidate_baseline_names = {
        "local_thermodynamic_reaction_rules",
        "stoichiometric_reaction_inverse",
    }
    for family, _case, rows in locked_cases:
        row = _complete_case_row(rows)
        if not isinstance(row, Mapping):
            continue
        record = row.get("reaction_rapm_calibrated_scores", {})
        if not isinstance(record, Mapping):
            continue
        score = record.get("reaction")
        if isinstance(score, Mapping) and score.get("status") == "scored":
            scores.append(score)
            family_counts[family] = family_counts.get(family, 0) + 1
            reliability = score.get("classwise_reliability", {})
            classwise_complete = classwise_complete and isinstance(reliability, Mapping)
        baselines = row.get("specialist_baselines", {})
        if isinstance(baselines, Mapping):
            for name in candidate_baseline_names:
                candidate = baselines.get(name)
                candidate_score = candidate.get("reaction") if isinstance(candidate, Mapping) else None
                if isinstance(candidate_score, Mapping):
                    nested_metrics = candidate_score.get("metrics", {})
                    log_loss = _finite_float(
                        candidate_score.get("multiclass_log_loss")
                        if candidate_score.get("multiclass_log_loss") is not None
                        else nested_metrics.get("multiclass_log_loss")
                        if isinstance(nested_metrics, Mapping)
                        else None
                    )
                    baseline_case_counts[name] = baseline_case_counts.get(name, 0) + 1
                    missing = int(candidate_score.get("n_missing_outputs", 0) or 0)
                    complete = bool(candidate_score.get("outputs_complete", missing == 0))
                    if not complete:
                        baseline_incomplete_counts[name] = (
                            baseline_incomplete_counts.get(name, 0) + 1
                        )
                else:
                    log_loss = None
                if log_loss is not None:
                    n_baseline = int(
                        candidate_score.get("n", 0)
                        if isinstance(candidate_score, Mapping)
                        else 0
                    )
                    if n_baseline > 0 and bool(
                        candidate_score.get(
                            "outputs_complete",
                            int(candidate_score.get("n_missing_outputs", 0) or 0) == 0,
                        )
                    ):
                        baseline_records.setdefault(name, []).append(
                            (log_loss, n_baseline)
                        )
    n = sum(int(score.get("n", 0)) for score in scores)
    accepted = sum(int(score.get("n_accepted", 0)) for score in scores)
    accepted_correct = sum(
        int(round(float(score.get("selective_accuracy", 0.0)) * int(score.get("n_accepted", 0))))
        for score in scores
        if _finite_float(score.get("selective_accuracy")) is not None
    )
    decoys = sum(int(score.get("n_decoy_edges", 0)) for score in scores)
    false_commitments = sum(
        float(score.get("false_commitment_rate", 0.0)) * int(score.get("n_decoy_edges", 0))
        for score in scores
        if _finite_float(score.get("false_commitment_rate")) is not None
    )
    log_loss_values = [
        (_finite_float(score.get("multiclass_log_loss")), int(score.get("n", 0)))
        for score in scores
    ]
    finite_log_loss = [(value, count) for value, count in log_loss_values if value is not None and count]
    model_log_loss = (
        sum(float(value) * count for value, count in finite_log_loss)
        / sum(count for _value, count in finite_log_loss)
        if finite_log_loss
        else None
    )
    baseline_means: dict[str, float] = {}
    for name, records in baseline_records.items():
        if not records or baseline_incomplete_counts.get(name, 0) > 0:
            continue
        total = sum(count for _value, count in records)
        if total:
            baseline_means[name] = sum(value * count for value, count in records) / total
    strongest_baseline_name = (
        min(baseline_means, key=lambda name: (baseline_means[name], name))
        if baseline_means
        else None
    )
    baseline_log_loss = (
        baseline_means[strongest_baseline_name]
        if strongest_baseline_name is not None
        else None
    )
    selection_target_met = bool(
        calibrator is not None
        and bool(
            getattr(calibrator, "provenance", {}).get("threshold_tuning", {}).get(
                "target_met", False
            )
        )
    )
    max_classwise_ece = None
    if scores:
        eces = []
        for score in scores:
            reliability = score.get("classwise_reliability", {})
            if not isinstance(reliability, Mapping):
                continue
            for record in reliability.values():
                if isinstance(record, Mapping):
                    value = _finite_float(record.get("ece"))
                    if value is not None:
                        eces.append(value)
        if eces:
            max_classwise_ece = max(eces)
    return {
        "status": "scored" if scores else "no_comparable_outputs",
        "n": n,
        "coverage": accepted / n if n else None,
        "coverage_including_abstention": accepted / n if n else None,
        "n_missing_outputs": sum(
            int(score.get("n_missing_outputs", 0) or 0) for score in scores
        ),
        "outputs_complete": bool(scores)
        and all(int(score.get("n_missing_outputs", 0) or 0) == 0 for score in scores),
        "selective_risk": 1.0 - accepted_correct / accepted if accepted else None,
        "false_commitment_rate": false_commitments / decoys if decoys else None,
        "multiclass_log_loss": model_log_loss,
        "baseline_log_loss": baseline_log_loss,
        "baseline_log_losses": dict(sorted(baseline_means.items())),
        "strongest_baseline": strongest_baseline_name,
        "baseline_records_complete": bool(baseline_means),
        "held_out_generators": len(family_counts) >= len(REQUIRED_GENERATOR_FAMILIES),
        "competence_matched_baseline": bool(baseline_means),
        "calibrated": calibrator is not None,
        "classwise_calibrated": classwise_complete and max_classwise_ece is not None,
        "max_classwise_ece": max_classwise_ece,
        "selection_rule_target_met": selection_target_met,
        "mechanism_stratified": len(family_counts) >= 2,
        "family_counts": dict(sorted(family_counts.items())),
    }


def _score_decision_policy(
    policy: Mapping[str, object],
    reference_edges: Sequence[Sequence[object]],
) -> dict[str, object]:
    """Release truth only for a labelled synthetic target diagnostic.

    This is not a prospective information-gain result: no post-measurement
    observation is simulated.  It reports only whether the selected synthetic
    action targeted a reference edge, so the decision stage remains clearly
    bounded until prospective outcomes are available.
    """

    selected = policy.get("selected_action_id")
    reference = {
        (str(edge[0]), str(edge[1]))
        for edge in reference_edges
        if len(edge) >= 2
    }
    if not selected:
        return {
            "status": "abstained",
            "selected_target_is_reference_edge": None,
            "evaluation_scope": "synthetic_oracle_target_diagnostic_only",
        }
    action_text = str(selected)
    target_text = action_text.split(":", 1)[1] if ":" in action_text else ""
    if "->" not in target_text:
        return {
            "status": "invalid_selected_target",
            "selected_target_is_reference_edge": None,
            "evaluation_scope": "synthetic_oracle_target_diagnostic_only",
        }
    source, target = target_text.split("->", 1)
    return {
        "status": "scored_after_truth_release",
        "selected_action_id": action_text,
        "selected_target_is_reference_edge": (source, target) in reference,
        "evaluation_scope": "synthetic_oracle_target_diagnostic_only",
        "not_prospective_field_utility": True,
    }


_POST_MEASUREMENT_LIKELIHOODS: dict[str, dict[str, dict[str, float]]] = {
    "head": {
        "present": {"support": 0.85, "no_support": 0.15},
        "absent": {"support": 0.25, "no_support": 0.75},
    },
    "tracer_age": {
        "present": {"support": 0.80, "no_support": 0.20},
        "absent": {"support": 0.20, "no_support": 0.80},
    },
    "chemistry": {
        "present": {"support": 0.78, "no_support": 0.22},
        "absent": {"support": 0.22, "no_support": 0.78},
    },
}
_POST_MEASUREMENT_DISCREPANT_LIKELIHOODS: dict[str, dict[str, dict[str, float]]] = {
    "head": {
        "present": {"support": 0.65, "no_support": 0.35},
        "absent": {"support": 0.45, "no_support": 0.55},
    },
    "tracer_age": {
        "present": {"support": 0.62, "no_support": 0.38},
        "absent": {"support": 0.42, "no_support": 0.58},
    },
    "chemistry": {
        "present": {"support": 0.60, "no_support": 0.40},
        "absent": {"support": 0.40, "no_support": 0.60},
    },
}


def _post_measurement_outcome_evaluation(
    policy: Mapping[str, object],
    candidates: Sequence[object],
    reference_edges: Sequence[Sequence[object]],
) -> dict[str, object]:
    """Score realised synthetic outcome utility after a blind policy action."""

    selected = policy.get("selected_action_id")
    if not selected:
        return {
            "status": "abstained",
            "evaluated": False,
            "improved": False,
            "reason": "policy_abstained_no_post_measurement_score",
            "truth_released_after_policy": False,
            "field_validation": "deferred",
        }
    action_id = str(selected)
    measurement, _, _ = action_id.partition(":")
    target_text = action_id.split(":", 1)[1] if ":" in action_id else ""
    if measurement not in _POST_MEASUREMENT_LIKELIHOODS or "->" not in target_text:
        return {
            "status": "invalid_selected_action",
            "evaluated": False,
            "improved": False,
            "reason": "selected_action_has_no_declared_observation_model",
            "truth_released_after_policy": False,
            "field_validation": "deferred",
        }
    candidate_probability: float | None = None
    for candidate in candidates:
        if _edge_pair(candidate) != tuple(target_text.split("->", 1)):
            continue
        attrs = _edge_attrs(candidate)
        candidate_probability = _finite_float(
            attrs.get("p_uv", attrs.get("flow_probability", attrs.get("edge_confidence")))
        )
        break
    if candidate_probability is None:
        candidate_probability = 0.5
    candidate_probability = min(0.95, max(0.05, candidate_probability))
    reference = {
        (str(edge[0]), str(edge[1]))
        for edge in reference_edges
        if len(edge) >= 2
    }
    true_hypothesis = "present" if tuple(target_text.split("->", 1)) in reference else "absent"
    posterior = PreMeasurementPosteriorSummary(
        posterior={"present": candidate_probability, "absent": 1.0 - candidate_probability},
        decision="MEASURE",
        action=action_id,
    )
    scenario_records: dict[str, object] = {}
    for scenario_name, likelihoods in (
        ("nominal", _POST_MEASUREMENT_LIKELIHOODS),
        ("discrepant", _POST_MEASUREMENT_DISCREPANT_LIKELIHOODS),
    ):
        model = SyntheticObservationModel(
            likelihoods_by_action={action_id: likelihoods[measurement]},
            target_metric="posterior_entropy",
        )
        evaluation = simulate_post_measurement(
            posterior,
            model,
            true_hypothesis,
            minimum_improvement=0.0,
        )
        scenario_records[scenario_name] = evaluation.to_dict()
    robust_improved = all(
        bool(record.get("improved"))
        for record in scenario_records.values()
        if isinstance(record, Mapping)
    )
    return {
        "status": "evaluated",
        "evaluated": True,
        "improved": robust_improved,
        "reason": (
            "target_metric_improved_under_all_declared_scenarios"
            if robust_improved
            else "target_metric_failed_at_least_one_declared_scenario"
        ),
        "action_id": action_id,
        "target": target_text,
        "measurement": measurement,
        "true_hypothesis": true_hypothesis,
        "pre_measurement_probability": candidate_probability,
        "scenarios": scenario_records,
        "truth_released_after_policy": True,
        "policy_truth_blind": True,
        "field_validation": "deferred",
    }


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object is not JSON serialisable: {type(value).__name__}")


def _generate_case(
    family: str,
    seed: int,
    *,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
) -> Case:
    if family == "modflow_modpath_v3":
        return generate_independent_aquifer_v3(
            seed,
            simulator_workspace / f"modflow_modpath_v3_{seed}",
            mf6_executable,
            mp7_executable,
        )
    if family == "analytic_lattice_v1":
        return generate_independent_lattice(
            seed,
            simulator_workspace / f"analytic_lattice_v1_{seed}",
        )
    if family == "independent_mixing_v1":
        return generate_independent_mixing(
            seed,
            simulator_workspace / f"independent_mixing_v1_{seed}",
        )
    raise ValueError(f"Unknown generator family: {family}")


def _run_case(
    family: str,
    phase: str,
    case: Case,
    *,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    age_samples: int,
    observation_scenarios: Mapping[str, Sequence[Mapping[str, object]]],
    head_evidence_audit: Mapping[str, object],
) -> tuple[
    list[dict[str, object]],
    list[ProgrammeStage],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    rows: list[dict[str, object]] = []
    stages: list[ProgrammeStage] = []
    errors: list[dict[str, object]] = []
    condition_rows: dict[str, dict[str, object]] = {}
    plans = [(condition, "complete") for condition in RUN_CONDITIONS]
    plans.extend(
        ("full_sheaf", scenario)
        for scenario in OBSERVATION_STRESS_SCENARIOS
        if scenario != "complete"
    )
    plans.extend(
        (
            "hydraulic_only" if scenario == "no_flow_head_permutation" else "full_sheaf",
            scenario,
        )
        for scenario in (
            "no_flow_head_permutation",
            "tracer_permutation",
            "chemistry_permutation",
        )
    )
    model = head_evidence_audit.get("model", {})
    head_sigma_m = None
    raw_sigma = head_evidence_audit.get("effective_sigma_m")
    try:
        head_sigma_m = float(raw_sigma) if raw_sigma is not None else None
    except (TypeError, ValueError):
        head_sigma_m = None
    if head_sigma_m is None and isinstance(model, Mapping):
        raw_sigma = model.get("primary_sigma_m")
        try:
            head_sigma_m = float(raw_sigma) if raw_sigma is not None else None
        except (TypeError, ValueError):
            head_sigma_m = None
    for condition, observation_scenario in plans:
        try:
            scenario_rows = observation_scenarios[observation_scenario]
            # Run comparator inference from the same blind rows and candidate
            # universe before the HydroSheaf result is scored.  The generator
            # truth enters only in the scoring calls below.
            blind_rows = _prepare_rows(
                case,
                permute_age=condition == "age_permuted",
                observations=scenario_rows,
            )
            blind_config = _condition_config(condition, head_sigma_m=head_sigma_m)
            blind_candidates = infer_edges(
                blind_rows,
                method="probabilistic",
                config=blind_config,
            )
            baseline_inference = _infer_baseline_outputs(
                blind_rows,
                blind_candidates,
            )
            independent_universe = generate_independent_candidate_universe(
                blind_rows,
                max_neighbors=4,
                head_tie_tolerance_m=0.10,
            )
            specialist_candidates = independent_universe.edges
            specialist_baseline_inference = _infer_baseline_outputs(
                blind_rows,
                specialist_candidates,
            )
            common_candidates, common_universe = _common_all_pairs_universe(blind_rows)
            common_baseline_inference = _infer_baseline_outputs(
                blind_rows,
                common_candidates,
            )
            hydrosheaf_prior_by_edge = {
                _edge_pair(edge): (
                    _finite_float(
                        _edge_attrs(edge).get(
                            "p_uv",
                            _edge_attrs(edge).get(
                                "flow_probability",
                                _edge_attrs(edge).get("edge_confidence"),
                            ),
                        )
                    )
                    or 0.5
                )
                for edge in blind_candidates
            }
            prospective_informativeness_by_action = {
                f"{measurement}:{source}->{target}": _pair_observation_quality(
                    blind_rows,
                    source,
                    target,
                    measurement,
                )
                for source, target in sorted(
                    {
                        _edge_pair(edge)
                        for edge in common_candidates
                    }
                )
                for measurement in ("head", "tracer_age", "chemistry")
            }
            prospective_actions = tuple(
                _measurement_actions(
                    common_candidates,
                    probability_by_edge=hydrosheaf_prior_by_edge,
                    informativeness_by_action=prospective_informativeness_by_action,
                )
            )
            # These selectors are run before the synthetic reference edge set
            # is read.  The resulting reports are the truth-blind decisions;
            # the realised outcome table is constructed only afterwards.
            prospective_hydrosheaf_report = select_declared_utility_measurement(
                prospective_actions,
                {
                    action.action_id: float(
                        dict(action.metadata).get(
                            "declared_expected_utility", -1.0e308
                        )
                    )
                    for action in prospective_actions
                },
                min_utility_per_cost=-1.0e308,
            )
            prospective_specialist_scores = _specialist_action_scores(
                prospective_actions,
                common_baseline_inference,
            )
            prospective_actions = tuple(
                CandidateMeasurementAction(
                    action_id=action.action_id,
                    cost=action.cost,
                    feasible=action.feasible,
                    scenarios=action.scenarios,
                    metadata={
                        **dict(action.metadata),
                        "specialist_score": prospective_specialist_scores[action.action_id],
                    },
                )
                for action in prospective_actions
            )
            prospective_specialist_report = select_specialist_measurement(
                prospective_actions,
                prospective_specialist_scores,
            )
            prospective_random_report = select_random_measurement(
                prospective_actions,
                seed=int(case.seed),
            )
            prospective_benefits, prospective_states, prospective_outcome_audit = (
                _prospective_realised_outcomes(
                    blind_rows,
                    prospective_actions,
                    case.true_edges,
                    seed=int(case.seed),
                )
            )
            prospective_case_record = {
                "case_id": f"{family}:{int(case.seed)}",
                "actions": _prospective_action_records(prospective_actions),
                "specialist_scores": prospective_specialist_scores,
                "pre_selection_reports": {
                    "hydrosheaf": prospective_hydrosheaf_report.to_dict(),
                    "specialist": prospective_specialist_report.to_dict(),
                    "random": prospective_random_report.to_dict(),
                },
                "benefit_by_action_and_state": prospective_benefits,
                "true_state_by_action": prospective_states,
                "truth_blind_selection": True,
                "candidate_universe_scope": "common_all_pairs_v1",
                "declared_utility_contract": _declared_utility_contract_record(),
                "outcome_audit": prospective_outcome_audit,
            }
            decision_policy = _decision_policy(blind_candidates)
            result = _run_condition(
                case,
                condition,
                simulator_workspace=simulator_workspace,
                mf6_executable=mf6_executable,
                mp7_executable=mp7_executable,
                age_samples=age_samples,
                observations=scenario_rows,
                observation_scenario=observation_scenario,
                head_sigma_m=head_sigma_m,
                head_evidence_audit=head_evidence_audit,
            )
            hydrosheaf_selected_edges = [
                tuple(str(value) for value in edge[:2])
                for edge in result.get("hydrosheaf_selected_edges", ())
                if isinstance(edge, Sequence) and len(edge) >= 2
            ]
            common_baseline_scores = _score_baseline_outputs(
                case.true_edges,
                common_candidates,
                common_baseline_inference,
                true_ages_years=case.true_ages_years,
                true_processes=case.true_processes,
            )
            common_topology_comparison = {
                "candidate_universe": common_universe,
                "hydrosheaf": {
                    "topology": score_topology(
                        case.true_edges,
                        common_candidates,
                        hydrosheaf_selected_edges,
                    )
                },
                "baselines": {
                    name: record["topology"]
                    for name, record in common_baseline_scores.items()
                    if isinstance(record, Mapping)
                    and isinstance(record.get("topology"), Mapping)
                },
                "comparison_scope": (
                    "selection_on_common_all_pairs_universe; candidate_generation "
                    "scored separately"
                ),
            }
            post_measurement_outcome = _post_measurement_outcome_evaluation(
                decision_policy,
                blind_candidates,
                case.true_edges,
            )
            result = {
                **result,
                "family": family,
                "phase": phase,
                "baseline_inference": baseline_inference,
                "hydrosheaf_topology_predictions": _hydrosheaf_probability_predictions(
                    blind_candidates
                ),
                "baselines": _score_baseline_outputs(
                    case.true_edges,
                    blind_candidates,
                    baseline_inference,
                    true_ages_years=case.true_ages_years,
                    true_processes=case.true_processes,
                ),
                "specialist_candidate_generation_truth_blind": True,
                "specialist_candidate_universe": independent_universe.to_audit_record(),
                "specialist_baseline_inference": specialist_baseline_inference,
                "specialist_baselines": _score_baseline_outputs(
                    case.true_edges,
                    specialist_candidates,
                    specialist_baseline_inference,
                    true_ages_years=case.true_ages_years,
                    true_processes=case.true_processes,
                ),
                "common_candidate_universe": common_universe,
                "common_baseline_inference": common_baseline_inference,
                "common_baselines": common_baseline_scores,
                "common_topology_comparison": common_topology_comparison,
                "specialist_candidate_universe_metrics": _candidate_universe_metrics(
                    case.true_edges,
                    specialist_candidates,
                ),
                "hydrosheaf_topology_probability": _hydrosheaf_probability_calibration(
                    case.true_edges,
                    blind_candidates,
                ),
                "decision_policy": decision_policy,
                "decision_policy_truth_blind": True,
                "decision_policy_evaluation": _score_decision_policy(
                    decision_policy,
                    case.true_edges,
                ),
                "post_measurement_outcome": post_measurement_outcome,
                "prospective_decision_case": prospective_case_record,
            }
            rows.append(result)
            condition_rows[f"{observation_scenario}:{condition}"] = result
            if observation_scenario == "complete":
                condition_rows[condition] = result
            stages.extend(
                programme_stages_from_status(
                    result["stage_status"],
                    name_prefix=(
                        f"{phase}:{family}:{case.seed}:"
                        f"{observation_scenario}:{condition}:"
                    ),
                )
            )
            stages.append(
                ProgrammeStage(
                    name=(
                        f"{phase}:{family}:{case.seed}:"
                        f"{observation_scenario}:{condition}:chemistry_reaction"
                    ),
                    status=StageStatus.COMPLETED,
                    detail=result["reaction_stage"]["detail"],
                )
            )
            stages.append(
                ProgrammeStage(
                    name=(
                        f"{phase}:{family}:{case.seed}:"
                        f"{observation_scenario}:{condition}:specialist_comparators"
                    ),
                    status=StageStatus.COMPLETED,
                    detail=(
                        "Independent geometry/head candidate universe generated "
                        f"with {len(specialist_candidates)} directed edges; "
                        "specialist age and reaction outputs scored separately."
                    ),
                )
            )
        except Exception as exc:  # preserve failed stages in the evidence record
            detail = f"{type(exc).__name__}: {exc}"
            errors.append(
                {
                    "phase": phase,
                    "family": family,
                    "seed": int(case.seed),
                    "condition": condition,
                    "observation_scenario": observation_scenario,
                    "error": detail,
                }
            )
            stages.append(
                ProgrammeStage(
                    name=(
                        f"{phase}:{family}:{case.seed}:"
                        f"{observation_scenario}:{condition}:runner"
                    ),
                    status=StageStatus.FAILED,
                    detail=detail,
                )
            )
    reports = _build_discrepancy_reports(int(case.seed), condition_rows)
    for report in reports:
        report["family"] = family
        report["phase"] = phase
    stages.append(
        ProgrammeStage(
            name=f"{phase}:{family}:{case.seed}:model_discrepancy",
            status=StageStatus.COMPLETED if reports else StageStatus.FAILED,
            detail=f"reported {len(reports)} age scenario comparisons",
        )
    )
    return rows, stages, errors, reports


def _development_observations(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
) -> list[AgeCalibrationObservation]:
    observations: list[AgeCalibrationObservation] = []
    for family, case, rows in cases:
        full_sheaf = next(
            (
                row
                for row in rows
                if row.get("condition") == "full_sheaf"
                and row.get("observation_scenario", "complete") == "complete"
            ),
            None,
        )
        if not full_sheaf:
            continue
        predictions = full_sheaf.get("age_predictions", {})
        if not isinstance(predictions, Mapping):
            continue
        for node, truth in case.true_ages_years.items():
            prediction = predictions.get(node)
            if not isinstance(prediction, Mapping):
                continue
            observations.append(
                AgeCalibrationObservation(
                    f"{family}:{node}",
                    truth=float(truth),
                    estimate=float(prediction["mean_age_years"]),
                    lower=float(prediction["age_95_low"]),
                    upper=float(prediction["age_95_high"]),
                )
            )
    return observations


def _locked_calibration_outputs(
    cases: Sequence[tuple[str, Case, list[dict[str, object]]]],
    calibrator: object,
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    scores: dict[str, object] = {}
    calibrated_predictions: dict[str, dict[str, object]] = {}
    all_selective_rows: list[dict[str, object]] = []
    for family, case, rows in cases:
        full_sheaf = next(
            (
                row
                for row in rows
                if row.get("condition") == "full_sheaf"
                and row.get("observation_scenario", "complete") == "complete"
            ),
            None,
        )
        if not full_sheaf:
            continue
        raw_predictions = full_sheaf.get("age_predictions", {})
        if not isinstance(raw_predictions, Mapping):
            continue
        calibrated = apply_age_interval_calibrator(calibrator, raw_predictions)
        calibrated_predictions[f"{family}:{case.seed}"] = calibrated
        truth = {str(node): float(value) for node, value in case.true_ages_years.items()}
        raw_score = score_age_posteriors(truth, raw_predictions)
        calibrated_score = score_calibrated_age_intervals(truth, calibrated)
        scores[f"{family}:{case.seed}"] = {
            "family": family,
            "seed": int(case.seed),
            "raw_age": raw_score,
            "calibrated_age": calibrated_score,
        }
        for node, true_value in truth.items():
            prediction = calibrated.get(node)
            if not prediction:
                continue
            low = float(prediction["calibrated_age_low"])
            high = float(prediction["calibrated_age_high"])
            all_selective_rows.append(
                {
                    "target_id": f"{family}:{case.seed}:{node}",
                    "truth": true_value,
                    "estimate": float(prediction["mean_age_years"]),
                    "lower": low,
                    "upper": high,
                    "uncertainty": high - low,
                }
            )
    curve = score_selective_risk(all_selective_rows)
    scores["locked_selective_risk_curve"] = [point.to_dict() for point in curve]
    scores["locked_total_age_predictions"] = len(all_selective_rows)
    return scores, calibrated_predictions


def run_ensemble_benchmark(
    *,
    run_id: str = RUN_ID,
    output: Path,
    simulator_workspace: Path,
    mf6_executable: Path,
    mp7_executable: Path,
    age_samples: int = 400,
    overwrite: bool = False,
    families: Sequence[str] | None = None,
) -> dict[str, object]:
    if int(age_samples) < 100:
        raise ValueError("age_samples must be at least 100 for a smoke gate.")
    output = _safe_output(Path(output), overwrite=overwrite)
    simulator_workspace = Path(simulator_workspace)
    mf6_executable = Path(mf6_executable).resolve()
    mp7_executable = Path(mp7_executable).resolve()
    selected_families = tuple(
        sorted(ALL_GENERATOR_FAMILIES if families is None else {str(item) for item in families})
    )
    unknown_families = set(selected_families) - ALL_GENERATOR_FAMILIES
    if unknown_families:
        raise ValueError(f"Unknown generator families: {sorted(unknown_families)}")
    if not selected_families:
        raise ValueError("At least one generator family must be selected.")
    if "modflow_modpath_v3" in selected_families and (
        not mf6_executable.exists() or not mp7_executable.exists()
    ):
        raise FileNotFoundError(
            "Both official mf6.exe and mp7.exe are required when the "
            "modflow_modpath_v3 family is selected."
        )
    simulator_workspace.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    all_stages: list[ProgrammeStage] = []
    errors: list[dict[str, object]] = []
    discrepancy_reports: list[dict[str, object]] = []
    generator_records: dict[str, Mapping[str, object]] = {}
    critic_records: dict[str, dict[str, object]] = {}
    observation_stress_records: dict[str, object] = {}
    development_cases: list[tuple[str, Case, list[dict[str, object]]]] = []
    locked_cases: list[tuple[str, Case, list[dict[str, object]]]] = []

    selected_development_cases = {
        family: DEVELOPMENT_CASES[family]
        for family in selected_families
    }
    selected_locked_cases = {
        family: LOCKED_CASES[family]
        for family in selected_families
    }
    for phase, case_registry, destination in (
        ("development", selected_development_cases, development_cases),
        ("locked_test", selected_locked_cases, locked_cases),
    ):
        for family, seeds in case_registry.items():
            for seed in seeds:
                try:
                    case = _generate_case(
                        family,
                        seed,
                        simulator_workspace=simulator_workspace,
                        mf6_executable=mf6_executable,
                        mp7_executable=mp7_executable,
                    )
                    provenance = dict(case.provenance)
                    if provenance.get("imports_hydrosheaf") is not False:
                        raise RuntimeError(
                            f"Independent generator import audit failed for {family}."
                        )
                    if family == "modflow_modpath_v3":
                        provenance["generator_source_sha256"] = _sha256(M7_V3_SOURCE)
                    observation_scenarios, head_audit, stress_records = (
                        _prepare_observation_scenarios(family, case)
                    )
                    provenance["head_channel_covariance_model"] = dict(
                        head_audit["model"]
                    )
                    provenance["head_evidence_audit"] = dict(head_audit)
                    provenance["observation_stress"] = dict(stress_records)
                    case_key = f"{phase}:{family}:{seed}"
                    generator_records[case_key] = provenance
                    observation_stress_records[case_key] = stress_records
                    repeat_case = _generate_case(
                        family,
                        seed,
                        simulator_workspace=simulator_workspace,
                        mf6_executable=mf6_executable,
                        mp7_executable=mp7_executable,
                    )
                    alternate_case = _generate_case(
                        family,
                        int(seed) + 1,
                        simulator_workspace=simulator_workspace,
                        mf6_executable=mf6_executable,
                        mp7_executable=mp7_executable,
                    )
                    critic = audit_generator_case(
                        family,
                        case,
                        source_path=GENERATOR_SOURCE_PATHS[family],
                        repeat_case=repeat_case,
                        alternate_case=alternate_case,
                        observation_scenarios=observation_scenarios,
                        provenance_override=provenance,
                        covariance_consumed=True,
                    )
                    critic_key = f"{phase}:{family}:{seed}"
                    critic_records[critic_key] = critic.to_dict()
                    all_stages.append(
                        ProgrammeStage(
                            name=f"{critic_key}:generator_critic",
                            status=StageStatus.COMPLETED,
                            detail=(
                                f"verdict={critic.verdict}; "
                                f"major={len(critic.majors)}; blocker={len(critic.blockers)}"
                            ),
                        )
                    )
                    rows, stages, case_errors, case_reports = _run_case(
                        family,
                        phase,
                        case,
                        simulator_workspace=simulator_workspace,
                        mf6_executable=mf6_executable,
                        mp7_executable=mp7_executable,
                        age_samples=age_samples,
                        observation_scenarios=observation_scenarios,
                        head_evidence_audit=head_audit,
                    )
                    destination.append((family, case, rows))
                    all_rows.extend(rows)
                    all_stages.extend(stages)
                    errors.extend(case_errors)
                    discrepancy_reports.extend(case_reports)
                except Exception as exc:
                    detail = f"{type(exc).__name__}: {exc}"
                    errors.append(
                        {
                            "phase": phase,
                            "family": family,
                            "seed": int(seed),
                            "condition": "generation",
                            "error": detail,
                        }
                    )
                    all_stages.append(
                        ProgrammeStage(
                            name=f"{phase}:{family}:{seed}:generation",
                            status=StageStatus.FAILED,
                            detail=detail,
                        )
                    )

    (
        discrepancy_calibrator,
        discrepancy_calibration_scores,
        discrepancy_reports_calibrated,
        discrepancy_calibration_records_complete,
        discrepancy_calibration_error,
    ) = _discrepancy_calibration_bundle(
        development_cases,
        locked_cases,
        discrepancy_reports,
    )
    development_observations = _development_observations(development_cases)
    calibrator = fit_age_interval_calibrator(
        development_observations,
        target_coverage=0.95,
        phase="development",
    )
    locked_scores, calibrated_predictions = _locked_calibration_outputs(
        locked_cases,
        calibrator,
    )
    model_averaging_error: str | None = None
    model_weight_fit = None
    model_averaging_scores: dict[str, object] = {"status": "not_available", "n": 0}
    development_model_observations = _topology_model_observations(
        development_cases,
        phase="development",
    )
    locked_model_observations = _topology_model_observations(
        locked_cases,
        phase="locked_test",
    )
    try:
        model_weight_fit = fit_discrete_model_weights(development_model_observations)
        model_averaging_scores = score_locked_model_average(
            locked_model_observations,
            model_weight_fit,
        )
    except Exception as exc:  # preserve an explicit failed validation stage
        model_averaging_error = f"{type(exc).__name__}: {exc}"
        errors.append(
            {
                "phase": "programme",
                "condition": "model_averaging",
                "error": model_averaging_error,
            }
        )
    model_averaging_records_complete = bool(
        model_weight_fit is not None
        and model_weight_fit.fit_scope == "development_only"
        and bool(development_model_observations)
        and bool(locked_model_observations)
        and model_averaging_scores.get("status") == "scored"
    )
    all_stages.append(
        ProgrammeStage(
            name="programme:model_averaging",
            status=(
                StageStatus.COMPLETED
                if model_averaging_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"case-blocked log-score fit on {len(development_model_observations)} "
                f"development targets; scored {len(locked_model_observations)} locked targets"
                if model_averaging_records_complete
                else model_averaging_error or "No complete model-averaging matrix was produced."
            ),
        )
    )
    reaction_rapm_model: ReactionRAPM | None = None
    reaction_rapm_fit_error: str | None = None
    reaction_rapm_development_scores: dict[str, dict[str, object]] = {}
    reaction_rapm_locked_scores: dict[str, dict[str, object]] = {}
    reaction_rapm_calibrator: object | None = None
    reaction_rapm_calibration_error: str | None = None
    reaction_rapm_calibration_observation_count = 0
    reaction_rapm_records_complete = False
    reaction_rapm_calibration_records_complete = False
    try:
        reaction_rapm_examples = _reaction_rapm_training_examples(development_cases)
        reaction_rapm_model = fit_reaction_rapm(
            reaction_rapm_examples,
            phase="development",
        )
        reaction_rapm_development_scores, development_rapm_records = (
            _apply_reaction_rapm_model(development_cases, reaction_rapm_model)
        )
        reaction_rapm_locked_scores, locked_rapm_records = _apply_reaction_rapm_model(
            locked_cases,
            reaction_rapm_model,
        )
        reaction_rapm_records_complete = bool(
            development_rapm_records and locked_rapm_records
        )
    except Exception as exc:
        reaction_rapm_fit_error = f"{type(exc).__name__}: {exc}"
        errors.append(
            {
                "phase": "programme",
                "condition": "reaction_rapm_fit",
                "error": reaction_rapm_fit_error,
            }
        )
    all_stages.append(
        ProgrammeStage(
            name="programme:reaction_rapm",
            status=(
                StageStatus.COMPLETED
                if reaction_rapm_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"fitted on {reaction_rapm_model.fit_record_count} development chemistry edges "
                f"across {len(reaction_rapm_model.fit_case_ids)} cases; "
                f"locked cases scored={len(reaction_rapm_locked_scores)}"
                if reaction_rapm_model is not None and reaction_rapm_records_complete
                else reaction_rapm_fit_error or "No complete reaction RAPM records were produced."
            ),
        )
    )
    if reaction_rapm_model is not None and reaction_rapm_records_complete:
        try:
            reaction_rapm_calibration_examples = (
                cross_fitted_reaction_rapm_calibration_examples(
                    reaction_rapm_examples,
                    config=reaction_rapm_model.config,
                )
            )
            reaction_rapm_calibration_observation_count = len(
                reaction_rapm_calibration_examples
            )
            reaction_rapm_calibrator = ReactionRAPMCalibrator.fit(
                reaction_rapm_calibration_examples,
                classes=reaction_rapm_model.config.classes,
                decision_threshold=0.60,
                decision_margin=0.10,
                tune_selection_rule=True,
                target_coverage=0.25,
                max_selective_risk=0.40,
                phase="development",
            )
            calibrated_locked_scores, calibration_records = (
                _apply_reaction_rapm_calibrator(
                    locked_cases,
                    reaction_rapm_calibrator,
                )
            )
            reaction_rapm_locked_scores = {
                **reaction_rapm_locked_scores,
                **{
                    case_id: {
                        **reaction_rapm_locked_scores.get(case_id, {}),
                        "calibrated": value,
                    }
                    for case_id, value in calibrated_locked_scores.items()
                },
            }
            reaction_rapm_calibration_records_complete = bool(
                calibration_records and calibrated_locked_scores
            )
        except Exception as exc:
            reaction_rapm_calibration_error = f"{type(exc).__name__}: {exc}"
            errors.append(
                {
                    "phase": "programme",
                    "condition": "reaction_rapm_calibration",
                    "error": reaction_rapm_calibration_error,
                }
            )
    else:
        reaction_rapm_calibration_error = "Reaction RAPM model was not available."
    all_stages.append(
        ProgrammeStage(
            name="programme:reaction_rapm_calibration",
            status=(
                StageStatus.COMPLETED
                if reaction_rapm_calibration_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"fitted frozen temperature scaling on "
                f"{reaction_rapm_calibration_observation_count} development predictions"
                if reaction_rapm_calibration_records_complete
                else reaction_rapm_calibration_error
                or "No complete reaction RAPM calibration records were produced."
            ),
        )
    )
    specialist_age_calibrators: dict[str, object] = {}
    specialist_reaction_calibrators: dict[str, object] = {}
    specialist_calibration_errors: dict[str, str] = {}
    locked_specialist_calibration_scores: dict[str, dict[str, object]] = {}
    specialist_calibration_records_complete = False
    try:
        (
            specialist_age_calibrators,
            specialist_reaction_calibrators,
            specialist_calibration_errors,
        ) = _fit_specialist_calibrators(development_cases)
        if specialist_calibration_errors:
            raise RuntimeError(
                json.dumps(specialist_calibration_errors, sort_keys=True)
            )
        (
            locked_specialist_calibration_scores,
            specialist_calibration_records_complete,
        ) = _apply_specialist_calibrators(
            locked_cases,
            specialist_age_calibrators,
            specialist_reaction_calibrators,
        )
        if not specialist_calibration_records_complete:
            raise RuntimeError("Locked specialist calibration records are incomplete.")
    except Exception as exc:
        if not specialist_calibration_errors:
            specialist_calibration_errors["runner"] = f"{type(exc).__name__}: {exc}"
        errors.append(
            {
                "phase": "programme",
                "condition": "specialist_calibration",
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
    all_stages.append(
        ProgrammeStage(
            name="programme:specialist_calibration",
            status=(
                StageStatus.COMPLETED
                if specialist_calibration_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"fitted {len(specialist_age_calibrators)} age and "
                f"{len(specialist_reaction_calibrators)} reaction calibrators "
                f"on {len(development_cases)} development cases"
                if specialist_calibration_records_complete
                else json.dumps(specialist_calibration_errors, sort_keys=True)
            ),
        )
    )
    kinetics_benchmark_report = None
    kinetics_benchmark_error: str | None = None
    kinetics_benchmark_records_complete = False
    try:
        kinetics_benchmark_report = run_kinetics_specialist_benchmark()
        kinetics_benchmark_records_complete = bool(
            kinetics_benchmark_report
            and kinetics_benchmark_report.specialist_score.n_cases > 0
            and kinetics_benchmark_report.specialist_score.status == "SCORED"
        )
    except Exception as exc:
        kinetics_benchmark_error = f"{type(exc).__name__}: {exc}"
        errors.append(
            {
                "phase": "programme",
                "condition": "m8_kinetics_benchmark",
                "error": kinetics_benchmark_error,
            }
        )
    all_stages.append(
        ProgrammeStage(
            name="programme:m8_kinetics_benchmark",
            status=(
                StageStatus.COMPLETED
                if kinetics_benchmark_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                "deterministic multi-regime M8 kinetics benchmark scored"
                if kinetics_benchmark_records_complete
                else kinetics_benchmark_error
                or "No complete M8 kinetics benchmark records were produced."
            ),
        )
    )
    age_performance_report: dict[str, object] = {
        "status": "not_available",
        "claim_status": "ABSTAIN",
    }
    age_performance_bundle: dict[str, object] = {"metrics": {}}
    age_performance_error: str | None = None
    try:
        age_performance_report, age_performance_bundle = _age_performance_bundle(
            development_cases,
            locked_cases,
        )
    except Exception as exc:
        age_performance_error = f"{type(exc).__name__}: {exc}"
        errors.append(
            {
                "phase": "programme",
                "condition": "age_performance",
                "error": age_performance_error,
            }
        )
    age_performance_records_complete = bool(
        age_performance_report.get("status") == "scored"
        and isinstance(age_performance_bundle.get("metrics"), Mapping)
    )
    age_gate = assess_age_gate(
        age_performance_bundle.get("metrics", {}),
        max_selective_risk_years=12.0,
    )
    all_stages.append(
        ProgrammeStage(
            name="programme:age_performance",
            status=(
                StageStatus.COMPLETED
                if age_performance_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"held-out calibrated age report scored; claim_gate={age_gate.status}"
                if age_performance_records_complete
                else age_performance_error or "No complete age performance record was produced."
            ),
        )
    )
    reaction_performance_metrics = _reaction_rapm_performance_bundle(
        locked_cases,
        calibrator=reaction_rapm_calibrator,
    )
    reaction_gate = assess_reaction_gate(reaction_performance_metrics)
    reaction_performance_records_complete = bool(
        reaction_performance_metrics.get("status") == "scored"
    )
    all_stages.append(
        ProgrammeStage(
            name="programme:reaction_performance",
            status=(
                StageStatus.COMPLETED
                if reaction_performance_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"held-out RAPM reaction records aggregated; claim_gate={reaction_gate.status}"
                if reaction_performance_records_complete
                else "No complete reaction RAPM performance record was produced."
            ),
        )
    )
    if kinetics_benchmark_report is not None:
        kinetics_score = kinetics_benchmark_report.specialist_score
        kinetic_regimes = {
            str(observation.regime)
            for observation in kinetics_benchmark_report.dataset.observations
            if observation.case_id in kinetics_benchmark_report.dataset.locked_case_ids
        }
        locked_predictions = tuple(
            prediction
            for prediction in kinetics_benchmark_report.specialist_predictions
            if prediction.case_id in kinetics_benchmark_report.dataset.locked_case_ids
        )
        locked_observations = tuple(
            observation
            for observation in kinetics_benchmark_report.dataset.observations
            if observation.case_id in kinetics_benchmark_report.dataset.locked_case_ids
        )
        area_informed_count = sum(
            observation.surface_area_measurement is not None
            for observation in locked_observations
        )
        identified_area_informed_count = sum(
            prediction.k_a_identifiable for prediction in locked_predictions
        )
        identifiability_given_area = (
            identified_area_informed_count / area_informed_count
            if area_informed_count
            else None
        )
        parameter_errors = [
            value
            for value in (
                kinetics_score.k_log_rmse_identified,
                kinetics_score.surface_area_log_rmse_identified,
            )
            if value is not None and math.isfinite(float(value))
        ]
        kinetics_performance_metrics: dict[str, object] = {
            "interval_coverage": kinetics_score.effective_rate_interval_coverage,
            # Separate k/A recovery is assessed conditional on the declared
            # independent surface-area measurement.  The unconditional rate is
            # retained separately because missing area is a structural
            # confounding condition, not a failed prediction.
            "identifiability_rate": identifiability_given_area,
            "identifiability_rate_overall": kinetics_score.identified_case_rate,
            "surface_area_informed_case_count": area_informed_count,
            "parameter_abstention_rate": kinetics_score.parameter_abstention_rate,
            "selective_risk": kinetics_score.false_commitment_rate,
            "selective_error": kinetics_score.selective_effective_rate_log_rmse,
            "predictive_rmse": kinetics_score.predictive_rmse,
            "parameter_error": max(parameter_errors) if parameter_errors else kinetics_score.effective_rate_log_rmse,
            "held_out_kinetic_regimes": len(kinetic_regimes) >= 3,
            "competence_matched_baseline": (
                kinetics_benchmark_report.comparator_status
                == "COMPETENCE_MATCHED_DIAGNOSTIC_ONLY"
                and kinetics_benchmark_report.comparator_score.status == "SCORED"
            ),
            "calibrated": (
                kinetics_benchmark_report.interval_calibrator.fit_scope
                == "development_only"
                and kinetics_benchmark_report.interval_calibrator.calibration_status
                == "development_truth_only_locked_apply"
            ),
            "ka_confounded_reported": bool(
                kinetics_score.parameter_abstention_rate > 0.0
                or kinetics_score.surface_area_log_rmse_identified is None
            ),
            "transport_stratified": len(kinetic_regimes) >= 3,
        }
    else:
        kinetics_performance_metrics = {}
    kinetics_gate = assess_kinetics_gate(kinetics_performance_metrics)
    kinetics_performance_records_complete = bool(
        kinetics_benchmark_records_complete
        and kinetics_benchmark_report is not None
    )
    all_stages.append(
        ProgrammeStage(
            name="programme:m8_kinetics_performance_gate",
            status=(
                StageStatus.COMPLETED
                if kinetics_performance_records_complete
                else StageStatus.FAILED
            ),
            detail=(
                f"M8 held-out metrics recorded; claim_gate={kinetics_gate.status}"
                if kinetics_performance_records_complete
                else "No complete M8 performance record was produced."
            ),
        )
    )
    all_scores_finite = True
    for value in locked_scores.values():
        if isinstance(value, Mapping):
            for nested in value.values():
                if isinstance(nested, Mapping):
                    for metric in nested.values():
                        if isinstance(metric, Mapping):
                            numeric = {
                                key: item
                                for key, item in metric.items()
                                if isinstance(item, (int, float)) and not isinstance(item, bool)
                            }
                            all_scores_finite = all_scores_finite and finite_numeric_mapping(numeric)

    required_case_count = sum(len(values) for values in selected_development_cases.values()) + sum(
        len(values) for values in selected_locked_cases.values()
    )
    case_count_complete = len(generator_records) == required_case_count
    all_rows_truth_blind = bool(all_rows) and all(bool(row.get("truth_blind")) for row in all_rows)
    all_rows_complete = bool(all_rows) and all(bool(row.get("stages_complete")) for row in all_rows)
    baseline_records_complete = bool(all_rows) and all(
        isinstance(row.get("baseline_inference"), Mapping)
        and isinstance(row.get("baselines"), Mapping)
        for row in all_rows
    )
    specialist_records_complete = bool(all_rows) and all(
        _specialist_record_complete(row)
        and isinstance(row.get("specialist_candidate_universe_metrics"), Mapping)
        for row in all_rows
    )
    locked_rows = [row for row in all_rows if row.get("phase") == "locked_test"]
    specialist_calibration_records_complete = bool(
        specialist_calibration_records_complete
        and locked_rows
        and all(
            bool(row.get("specialist_calibration_applied"))
            and isinstance(row.get("specialist_calibrated_baseline_inference"), Mapping)
            and isinstance(row.get("specialist_calibrated_baselines"), Mapping)
            for row in locked_rows
        )
    )
    decision_policy_records_complete = bool(all_rows) and all(
        bool(row.get("decision_policy_truth_blind"))
        and isinstance(row.get("decision_policy"), Mapping)
        for row in all_rows
    )
    post_measurement_records_complete = bool(all_rows) and all(
        isinstance(row.get("post_measurement_outcome"), Mapping)
        and bool(row.get("post_measurement_outcome", {}).get("policy_truth_blind", True))
        for row in all_rows
    )
    locked_complete_outcomes = [
        row.get("post_measurement_outcome")
        for row in all_rows
        if row.get("phase") == "locked_test"
        and row.get("condition") == "full_sheaf"
        and row.get("observation_scenario", "complete") == "complete"
    ]
    post_measurement_claim_gate = bool(
        locked_complete_outcomes
        and all(
            isinstance(record, Mapping)
            and bool(record.get("evaluated"))
            and bool(record.get("improved"))
            for record in locked_complete_outcomes
        )
    )
    prospective_case_objects: list[ProspectiveMeasurementCase] = []
    prospective_case_record_count = 0
    prospective_case_build_errors: list[str] = []
    prospective_outcomes_independent = True
    for family, case, rows in locked_cases:
        row = _complete_case_row(rows)
        record = row.get("prospective_decision_case") if isinstance(row, Mapping) else None
        if not isinstance(record, Mapping):
            prospective_case_build_errors.append(f"{family}:{case.seed}:missing_record")
            prospective_outcomes_independent = False
            continue
        raw_actions = record.get("actions", ())
        raw_benefits = record.get("benefit_by_action_and_state", {})
        raw_states = record.get("true_state_by_action", {})
        if not (
            isinstance(raw_actions, Sequence)
            and isinstance(raw_benefits, Mapping)
            and isinstance(raw_states, Mapping)
            and bool(record.get("truth_blind_selection"))
        ):
            prospective_case_build_errors.append(f"{family}:{case.seed}:malformed_record")
            prospective_outcomes_independent = False
            continue
        outcome_audit = record.get("outcome_audit", {})
        if not (
            isinstance(outcome_audit, Mapping)
            and outcome_audit.get("rule")
            == "independent_state_conditioned_measurement_v2"
            and bool(outcome_audit.get("state_conditioned"))
            and outcome_audit.get("reuses_selection_quality") is False
        ):
            prospective_outcomes_independent = False
        action_records = [item for item in raw_actions if isinstance(item, Mapping)]
        actions = _rebuild_prospective_actions(action_records)
        try:
            prospective_case_objects.append(
                ProspectiveMeasurementCase(
                    case_id=f"{family}:{int(case.seed)}",
                    actions=actions,
                    benefit_by_action_and_state={
                        str(action_id): dict(by_state)
                        for action_id, by_state in raw_benefits.items()
                        if isinstance(by_state, Mapping)
                    },
                    true_state_by_action={
                        str(action_id): str(state)
                        for action_id, state in raw_states.items()
                    },
                )
            )
            prospective_case_record_count += 1
        except Exception as exc:
            prospective_case_build_errors.append(
                f"{family}:{case.seed}:{type(exc).__name__}:{exc}"
            )
            prospective_outcomes_independent = False
    prospective_decision_records_complete = bool(
        prospective_case_objects
        and len(prospective_case_objects) == len(locked_cases)
        and not prospective_case_build_errors
    )
    prospective_decision_benchmark: dict[str, object] = {
        "status": "ABSTAIN",
        "claim_status": "ABSTAIN",
        "records_complete": prospective_decision_records_complete,
        "reason": "prospective benchmark not yet evaluated",
        "case_count": prospective_case_record_count,
    }
    if prospective_decision_records_complete:
        try:
            benchmark = evaluate_prospective_policies(
                prospective_case_objects,
                [
                    ProspectivePolicy(
                        "hydrosheaf",
                        lambda actions: select_declared_utility_measurement(
                            actions,
                            {
                                action.action_id: float(
                                    dict(action.metadata).get(
                                        "declared_expected_utility", -1.0e308
                                    )
                                )
                                for action in actions
                            },
                            min_utility_per_cost=-1.0e308,
                        ),
                    ),
                    ProspectivePolicy(
                        "random",
                        lambda actions: select_random_measurement(actions, seed=0),
                        scoring_mode="uniform_action_expectation",
                    ),
                    ProspectivePolicy(
                        "specialist",
                        lambda actions: select_specialist_measurement(
                            actions,
                            {
                                action.action_id: float(
                                    dict(action.metadata).get("specialist_score", 0.5)
                                )
                                for action in actions
                            },
                        ),
                    ),
                ],
                cost_penalty=float(PROSPECTIVE_DECLARED_UTILITY_MODEL["cost_penalty"]),
                required_policy_ids=("random", "specialist"),
                calibration_sufficient=(
                    prospective_case_record_count >= PROSPECTIVE_MINIMUM_LOCKED_CASES
                ),
            )
            prospective_decision_benchmark = {
                **benchmark.to_dict(),
                "records_complete": True,
                "evidence_contract": {
                    "minimum_locked_case_count": PROSPECTIVE_MINIMUM_LOCKED_CASES,
                    "observed_locked_case_count": prospective_case_record_count,
                    "truth_used_only_after_selection": True,
                },
                "truth_blind_selection": True,
                "truth_released_only_for_post_selection_scoring": True,
                "claim_scope": "controlled_synthetic_generator_realised_decision_utility_only",
                "candidate_universe_scope": "common_all_pairs_v1",
                "declared_utility_contract": _declared_utility_contract_record(),
                "outcome_simulator": "independent_state_conditioned_measurement_v2",
                "outcomes_independent_of_selection_quality": True,
                "random_policy_evaluation": (
                    "uniform_action_expectation_over_all_feasible_actions"
                ),
            }
        except Exception as exc:
            detail = f"{type(exc).__name__}: {exc}"
            prospective_decision_benchmark = {
                "status": "ABSTAIN",
                "claim_status": "ABSTAIN",
                "records_complete": False,
                "reason": detail,
                "case_count": prospective_case_record_count,
            }
            prospective_case_build_errors.append(detail)
    elif prospective_case_build_errors:
        prospective_decision_benchmark["reason"] = "; ".join(prospective_case_build_errors)
    model_averaging_converged = bool(
        model_weight_fit is not None and bool(model_weight_fit.converged)
    )
    prospective_policies = prospective_decision_benchmark.get("policies", {})
    prospective_policies = (
        prospective_policies if isinstance(prospective_policies, Mapping) else {}
    )
    hydro_policy_score = prospective_policies.get("hydrosheaf", {})
    random_policy_score = prospective_policies.get("random", {})
    specialist_policy_score = prospective_policies.get("specialist", {})
    prospective_pairwise = prospective_decision_benchmark.get("pairwise", {})
    prospective_pairwise = (
        prospective_pairwise if isinstance(prospective_pairwise, Mapping) else {}
    )
    hydro_random_pair = prospective_pairwise.get("hydrosheaf_vs_random", {})
    hydro_specialist_pair = prospective_pairwise.get("hydrosheaf_vs_specialist", {})
    hydro_random_pair = (
        hydro_random_pair if isinstance(hydro_random_pair, Mapping) else {}
    )
    hydro_specialist_pair = (
        hydro_specialist_pair
        if isinstance(hydro_specialist_pair, Mapping)
        else {}
    )
    source_snapshot_recorded = bool(
        PROGRAMME_SOURCE_FILES
        and all((REPO / relative).exists() for relative in PROGRAMME_SOURCE_FILES)
        and all(
            isinstance(_sha256(REPO / relative), str)
            and len(_sha256(REPO / relative)) == 64
            for relative in PROGRAMME_SOURCE_FILES
        )
    )
    prospective_outcomes_complete = bool(
        prospective_decision_benchmark.get("status") == "SCORED"
        and all(
            isinstance(score, Mapping) and score.get("status") == "SCORED"
            for score in (
                hydro_policy_score,
                random_policy_score,
                specialist_policy_score,
            )
        )
    )
    integrated_gate_metrics: dict[str, object] = {
        "model_averaging_converged": model_averaging_converged,
        "discrepancy_calibrated": (
            True
            if discrepancy_calibrator is not None
            and discrepancy_calibration_records_complete
            and discrepancy_calibrator.fit_scope == "development_only"
            and discrepancy_calibration_scores.get("status") == "scored"
            else None
        ),
        "prospective_outcomes_complete": prospective_outcomes_complete,
        "prospective_evidence_sufficient": bool(
            prospective_decision_benchmark.get("calibration_sufficient")
        ),
        "prospective_case_count": prospective_decision_benchmark.get("case_count"),
        "prospective_outcomes_independent": prospective_outcomes_independent,
        "prospective_random_baseline_valid": bool(
            isinstance(random_policy_score, Mapping)
            and random_policy_score.get("evaluation_mode")
            == "uniform_action_expectation"
            and prospective_decision_benchmark.get("random_policy_evaluation")
            == "uniform_action_expectation_over_all_feasible_actions"
        ),
        "paired_uncertainty_available": bool(
            hydro_random_pair.get("paired_uncertainty_available")
            and hydro_specialist_pair.get("paired_uncertainty_available")
        ),
        "paired_random_delta_ci_low": hydro_random_pair.get(
            "paired_delta_ci_low"
        ),
        "paired_random_delta_ci_high": hydro_random_pair.get(
            "paired_delta_ci_high"
        ),
        "paired_specialist_delta_ci_low": hydro_specialist_pair.get(
            "paired_delta_ci_low"
        ),
        "paired_specialist_delta_ci_high": hydro_specialist_pair.get(
            "paired_delta_ci_high"
        ),
        "source_snapshot_recorded": source_snapshot_recorded,
        "raw_benchmark_verified": bool(
            prospective_decision_records_complete
            and prospective_decision_benchmark.get("status") == "SCORED"
            and prospective_decision_benchmark.get("candidate_universe_scope")
            == "common_all_pairs_v1"
        ),
        "hydrosheaf_mean_utility_per_cost": (
            hydro_policy_score.get("mean_cost_adjusted_utility")
            if isinstance(hydro_policy_score, Mapping)
            else None
        ),
        "random_mean_utility_per_cost": (
            random_policy_score.get("mean_cost_adjusted_utility")
            if isinstance(random_policy_score, Mapping)
            else None
        ),
        "strongest_specialist_mean_utility_per_cost": (
            specialist_policy_score.get("mean_cost_adjusted_utility")
            if isinstance(specialist_policy_score, Mapping)
            else None
        ),
        "calibration_degradation": (
            discrepancy_calibration_scores.get("calibration_degradation")
            if discrepancy_calibration_scores.get("status") == "scored"
            else None
        ),
        "observed_false_commitment_rate": (
            discrepancy_calibration_scores.get("calibrated", {}).get(
                "false_commitment_rate"
            )
            if isinstance(discrepancy_calibration_scores.get("calibrated"), Mapping)
            else None
        ),
        "no_material_calibration_degradation": (
            float(discrepancy_calibration_scores["calibration_degradation"]) <= 0.05
            if discrepancy_calibration_scores.get("status") == "scored"
            and _finite_float(discrepancy_calibration_scores.get("calibration_degradation")) is not None
            else None
        ),
        "false_commitment_controlled": (
            float(discrepancy_calibration_scores["calibrated"]["false_commitment_rate"]) <= 0.10
            if isinstance(discrepancy_calibration_scores.get("calibrated"), Mapping)
            and _finite_float(
                discrepancy_calibration_scores.get("calibrated", {}).get(
                    "false_commitment_rate"
                )
            ) is not None
            else None
        ),
    }
    integrated_gate = assess_integrated_gate(
        integrated_gate_metrics,
        minimum_prospective_case_count=PROSPECTIVE_MINIMUM_LOCKED_CASES,
    )
    specialist_claim_summary = aggregate_specialist_gates(
        age=age_gate,
        reaction=reaction_gate,
        kinetics=kinetics_gate,
        integrated=integrated_gate,
    )
    present_families = {key.split(":", 2)[1] for key in generator_records}
    selected_generator_families_present = set(selected_families).issubset(
        present_families
    )
    required_generator_families_present = REQUIRED_GENERATOR_FAMILIES.issubset(
        present_families
    )
    selected_modflow_case_count = sum(
        len(registry.get("modflow_modpath_v3", ()))
        for registry in (selected_development_cases, selected_locked_cases)
    )
    modflow_provenance_records = [
        record
        for key, record in generator_records.items()
        if ":modflow_modpath_v3:" in key
    ]
    external_solver_execution_gate = (
        selected_modflow_case_count == 0
        or (
            len(modflow_provenance_records) == selected_modflow_case_count
            and all(
                all(
                    record.get(field)
                    for field in (
                        "mf6_executable",
                        "mf6_sha256",
                        "mf6_version",
                        "mp7_executable",
                        "mp7_sha256",
                        "mp7_version",
                    )
                )
                for record in modflow_provenance_records
            )
        )
    )
    discrepancy_reports_finite = bool(discrepancy_reports) and all(
        finite_numeric_mapping(
            {
                "estimate": report["estimate"],
                "lower": report["lower"],
                "upper": report["upper"],
                "scenario_range": report["scenario_range"],
                "relative_disagreement": report["relative_disagreement"],
            }
        )
        for report in discrepancy_reports
    )
    critic_complete = len(critic_records) == required_case_count
    critic_gate = bool(critic_records) and all(
        bool(record.get("critic_gate")) for record in critic_records.values()
    )
    critic_blocker_count = sum(
        int(record.get("summary", {}).get("BLOCKER", 0))
        for record in critic_records.values()
    )
    critic_major_count = sum(
        int(record.get("summary", {}).get("MAJOR", 0))
        for record in critic_records.values()
    )
    expected_stress_record_count = required_case_count
    observation_stress_complete = len(observation_stress_records) == expected_stress_record_count
    required_stress_scenarios = set(OBSERVATION_STRESS_SCENARIOS) - {"complete"}
    observation_stress_has_missingness = bool(observation_stress_records) and all(
        all(
            int(dict(scenario).get("missing_count", 0)) > 0
            or int(dict(scenario).get("censored_count", 0)) > 0
            for scenario_name, scenario in dict(records).items()
            if scenario_name in required_stress_scenarios
        )
        for records in observation_stress_records.values()
        if isinstance(records, Mapping)
    )
    primary_critic_records = {
        key: record
        for key, record in critic_records.items()
        if any(f":{family_name}:" in key for family_name in PRIMARY_GENERATOR_FAMILIES)
    }
    selected_primary_families = set(selected_families) & PRIMARY_GENERATOR_FAMILIES
    primary_critic_gate = (
        not selected_primary_families
        or bool(primary_critic_records)
    ) and all(
        bool(record.get("critic_gate")) for record in primary_critic_records.values()
    )
    execution_gate = bool(
        not errors
        and case_count_complete
        and selected_generator_families_present
        and external_solver_execution_gate
        and all_rows_truth_blind
        and all_rows_complete
        and baseline_records_complete
        and specialist_records_complete
        and specialist_calibration_records_complete
        and reaction_rapm_records_complete
        and reaction_rapm_calibration_records_complete
        and kinetics_benchmark_records_complete
        and age_performance_records_complete
        and reaction_performance_records_complete
        and kinetics_performance_records_complete
        and decision_policy_records_complete
        and post_measurement_records_complete
        and prospective_decision_records_complete
        and critic_complete
        and critic_gate
        and observation_stress_complete
        and observation_stress_has_missingness
        and calibrator.fit_count == len(development_observations)
        and bool(locked_scores.get("locked_total_age_predictions"))
        and all_scores_finite
        and discrepancy_reports_finite
        and discrepancy_calibration_records_complete
        and model_averaging_records_complete
    )
    integrated_performance_gate = bool(
        execution_gate
        and required_generator_families_present
        and integrated_gate.status == "PASS"
    )

    programme_decisions: list[ProgrammeDecision] = []
    for row in all_rows:
        if row.get("condition") != "full_sheaf" or row.get("observation_scenario") != "complete":
            continue
        policy = row.get("decision_policy")
        if not isinstance(policy, Mapping):
            continue
        selected_id = policy.get("selected_action_id")
        audit_rows = policy.get("audit_records", [])
        selected_record = next(
            (
                item
                for item in audit_rows
                if isinstance(item, Mapping) and item.get("action_id") == selected_id
            ),
            None,
        ) if isinstance(audit_rows, Sequence) else None
        if policy.get("decision") == "MEASURE" and isinstance(selected_record, Mapping):
            metadata = selected_record.get("metadata", {})
            metadata = metadata if isinstance(metadata, Mapping) else {}
            measurement = str(metadata.get("measurement", "measurement"))
            target = str(metadata.get("target", selected_id))
            programme_decisions.append(
                ProgrammeDecision(
                    decision_kind=DecisionKind.ACTION,
                    identifiability=IdentifiabilityStatus.PARTIALLY_IDENTIFIED,
                    reason=(
                        "Selected by the truth-blind robust information-gain-per-cost "
                        "policy under declared nominal/discrepant scenarios."
                    ),
                    measurement=measurement,
                    target=target,
                    cost=float(selected_record["cost"]),
                    expected_utility=float(selected_record["robust_utility_per_cost"]),
                    scenario_count=len(selected_record.get("scenario_utilities", ())),
                    provenance={
                        "family": row.get("family"),
                        "seed": row.get("seed"),
                        "decision_policy_truth_blind": True,
                        "prospective_outcome_validation": "deferred",
                    },
                )
            )
        else:
            programme_decisions.append(
                ProgrammeDecision(
                    decision_kind=DecisionKind.ABSTAIN,
                    identifiability=IdentifiabilityStatus.UNKNOWN,
                    reason=(
                        "No feasible action cleared the predeclared robust utility "
                        "threshold."
                    ),
                    scenario_count=0,
                    provenance={
                        "family": row.get("family"),
                        "seed": row.get("seed"),
                        "decision_policy_truth_blind": True,
                        "prospective_outcome_validation": "deferred",
                    },
                )
            )
    if not programme_decisions:
        programme_decisions.append(
            ProgrammeDecision(
                decision_kind=DecisionKind.ABSTAIN,
                identifiability=IdentifiabilityStatus.UNKNOWN,
                reason="No complete full-sheaf policy record was available.",
                scenario_count=0,
                provenance={"prospective_outcome_validation": "deferred"},
            )
        )

    all_stages.append(
        ProgrammeStage(
            name="programme:calibration",
            status=StageStatus.COMPLETED if execution_gate else StageStatus.FAILED,
            detail=(
                f"fitted on {calibrator.fit_count} development predictions and "
                "applied before locked-test scoring"
            ),
        )
    )
    all_stages.append(
        ProgrammeStage(
            name="programme:next_measurement",
            status=StageStatus.COMPLETED if execution_gate else StageStatus.FAILED,
            detail=(
                "Truth-blind one-step policy selected or abstained with a "
                "predeclared cost threshold; declared synthetic outcomes are "
                "scored after the policy decision."
            ),
        )
    )
    programme_run = ProgrammeRun(
        run_id=run_id,
        generator="independent_generator_ensemble",
        generator_independent=True,
        stages=tuple(all_stages),
        decisions=tuple(programme_decisions),
        claim_boundary=(
            "development-fitted calibration on held-out controlled-synthetic "
            "generator families with declared post-measurement outcome scoring; "
            "not field validation or universal superiority"
        ),
        truth_released_for_scoring=False,
    )

    _write_json(output / "case_metrics.json", all_rows)
    _write_json(
        output / "specialist_comparators.json",
        [
            {
                "family": row.get("family"),
                "phase": row.get("phase"),
                "seed": row.get("seed"),
                "condition": row.get("condition"),
                "observation_scenario": row.get("observation_scenario"),
                "candidate_generation_truth_blind": row.get(
                    "specialist_candidate_generation_truth_blind"
                ),
                "candidate_universe": row.get("specialist_candidate_universe"),
                "candidate_universe_metrics": row.get(
                    "specialist_candidate_universe_metrics"
                ),
                "baseline_inference": row.get("specialist_baseline_inference"),
                "baseline_scores": row.get("specialist_baselines"),
                "reaction_rapm_inference": row.get("reaction_rapm_inference"),
                "reaction_rapm_scores": row.get("reaction_rapm_scores"),
                "reaction_rapm_calibrated_inference": row.get(
                    "reaction_rapm_calibrated_inference"
                ),
                "reaction_rapm_calibrated_scores": row.get(
                    "reaction_rapm_calibrated_scores"
                ),
                "reaction_rapm_truth_blind": row.get("reaction_rapm_truth_blind"),
            }
            for row in all_rows
            if "specialist_candidate_universe" in row
        ],
    )
    _write_json(
        output / "specialist_calibration.json",
        {
            "fit_scope": "development_only",
            "age_calibrators": {
                name: calibrator.to_dict()
                for name, calibrator in specialist_age_calibrators.items()
            },
            "reaction_calibrators": {
                name: calibrator.to_dict()
                for name, calibrator in specialist_reaction_calibrators.items()
            },
            "locked_scores": locked_specialist_calibration_scores,
            "errors": specialist_calibration_errors,
            "records_complete": specialist_calibration_records_complete,
            "reaction_rapm_calibrator": (
                reaction_rapm_calibrator.to_dict()
                if reaction_rapm_calibrator is not None
                else None
            ),
            "reaction_rapm_calibration_observation_count": (
                reaction_rapm_calibration_observation_count
            ),
            "reaction_rapm_calibration_error": reaction_rapm_calibration_error,
            "reaction_rapm_records_complete": reaction_rapm_records_complete,
            "reaction_rapm_calibration_records_complete": (
                reaction_rapm_calibration_records_complete
            ),
            "reaction_rapm_calibration_prediction_protocol": (
                "development_case_cross_fitted"
            ),
            "claim_status": "calibrated_locked_diagnostic_not_superiority",
        },
    )
    _write_json(output / "discrepancy_reports.json", discrepancy_reports)
    _write_json(output / "discrepancy_reports_calibrated.json", discrepancy_reports_calibrated)
    _write_json(
        output / "discrepancy_calibration.json",
        {
            "calibrator": (
                discrepancy_calibrator.to_dict()
                if discrepancy_calibrator is not None
                else None
            ),
            "scores": discrepancy_calibration_scores,
            "records_complete": discrepancy_calibration_records_complete,
            "error": discrepancy_calibration_error,
            "claim_status": "development_fitted_locked_diagnostic_not_superiority",
        },
    )
    _write_json(
        output / "age_performance.json",
        {
            "report": age_performance_report,
            "bundle": age_performance_bundle,
            "gate": age_gate.to_dict(),
        },
    )
    _write_json(
        output / "reaction_performance.json",
        {
            "metrics": reaction_performance_metrics,
            "gate": reaction_gate.to_dict(),
        },
    )
    _write_json(
        output / "kinetics_benchmark.json",
        {
            "report": (
                kinetics_benchmark_report.as_dict()
                if kinetics_benchmark_report is not None
                else None
            ),
            "metrics": kinetics_performance_metrics,
            "gate": kinetics_gate.to_dict(),
        },
    )
    _write_json(
        output / "prospective_decision_benchmark.json",
        {
            "benchmark": prospective_decision_benchmark,
            "case_build_errors": prospective_case_build_errors,
            "records_complete": prospective_decision_records_complete,
            "claim_boundary": (
                "generator-realised observation-separation utility with action-specific "
                "truth released only after truth-blind selector return; not field utility"
            ),
        },
    )
    _write_json(output / "performance_gates.json", specialist_claim_summary)
    _write_json(output / "generator_provenance.json", generator_records)
    _write_json(output / "generator_critic.json", critic_records)
    _write_json(output / "observation_stress.json", observation_stress_records)
    _write_json(output / "stage_status.json", [stage.to_dict() for stage in all_stages])
    _write_json(output / "programme_record.json", programme_run.to_dict())
    _write_json(output / "calibrator.json", calibrator.to_dict())
    _write_json(output / "calibration_scores.json", locked_scores)
    _write_json(output / "calibrated_predictions.json", calibrated_predictions)
    _write_json(
        output / "model_averaging.json",
        {
            "fit": model_weight_fit.to_dict() if model_weight_fit is not None else None,
            "development_observation_count": len(development_model_observations),
            "locked_observation_count": len(locked_model_observations),
            "candidate_universe_scope": "hydrosheaf_native_candidates_only",
            "score": model_averaging_scores,
            "error": model_averaging_error,
            "claim_status": "conditional_diagnostic_not_superiority",
        },
    )
    _write_json(
        output / "reaction_rapm.json",
        {
            "model": (
                reaction_rapm_model.to_dict()
                if reaction_rapm_model is not None
                else None
            ),
            "development_scores": reaction_rapm_development_scores,
            "locked_scores": reaction_rapm_locked_scores,
            "calibrator": (
                reaction_rapm_calibrator.to_dict()
                if reaction_rapm_calibrator is not None
                else None
            ),
            "calibration_observation_count": reaction_rapm_calibration_observation_count,
            "calibration_prediction_protocol": "development_case_cross_fitted",
            "fit_error": reaction_rapm_fit_error,
            "calibration_error": reaction_rapm_calibration_error,
            "records_complete": reaction_rapm_records_complete,
            "calibration_records_complete": reaction_rapm_calibration_records_complete,
            "claim_status": "development_fitted_locked_diagnostic_not_superiority",
        },
    )
    _write_json(
        output / "baseline_registry.json",
        list(BASELINE_REGISTRY.audit_table())
        + list(SPECIALIST_BASELINE_REGISTRY.audit_table()),
    )
    _write_json(
        output / "decision_policy.json",
        [
            {
                "family": row.get("family"),
                "phase": row.get("phase"),
                "seed": row.get("seed"),
                "condition": row.get("condition"),
                "observation_scenario": row.get("observation_scenario"),
                "decision_policy": row.get("decision_policy"),
                "truth_blind": row.get("decision_policy_truth_blind"),
            }
            for row in all_rows
            if "decision_policy" in row
        ],
    )
    _write_json(
        output / "decision_policy_evaluation.json",
        [
            {
                "family": row.get("family"),
                "phase": row.get("phase"),
                "seed": row.get("seed"),
                "condition": row.get("condition"),
                "observation_scenario": row.get("observation_scenario"),
                "evaluation": row.get("decision_policy_evaluation"),
            }
            for row in all_rows
            if "decision_policy_evaluation" in row
        ],
    )
    _write_json(
        output / "post_measurement_outcomes.json",
        [
            {
                "family": row.get("family"),
                "phase": row.get("phase"),
                "seed": row.get("seed"),
                "condition": row.get("condition"),
                "observation_scenario": row.get("observation_scenario"),
                "evaluation": row.get("post_measurement_outcome"),
            }
            for row in all_rows
            if "post_measurement_outcome" in row
        ],
    )
    _write_json(output / "errors.json", errors)
    _write_json(
        output / "execution_gate.json",
        {
            "age_samples": int(age_samples),
            "status": "PASS" if execution_gate else "FAIL",
            "gate_scope": f"selected synthetic generator families: {list(selected_families)}",
            "checks": {
                "selected_generator_families_present": selected_generator_families_present,
                "external_solver_execution_gate": external_solver_execution_gate,
                "required_generator_families_present": required_generator_families_present,
                "case_count_complete": case_count_complete,
                "truth_blind_inference_rows": all_rows_truth_blind,
                "all_required_stages_completed": all_rows_complete,
                "baseline_records_complete": baseline_records_complete,
                "specialist_records_complete": specialist_records_complete,
                "specialist_calibration_records_complete": specialist_calibration_records_complete,
                "reaction_rapm_records_complete": reaction_rapm_records_complete,
                "reaction_rapm_calibration_records_complete": (
                    reaction_rapm_calibration_records_complete
                ),
                "age_performance_records_complete": age_performance_records_complete,
                "reaction_performance_records_complete": reaction_performance_records_complete,
                "kinetics_performance_records_complete": kinetics_performance_records_complete,
                "decision_policy_records_complete": decision_policy_records_complete,
                "post_measurement_records_complete": post_measurement_records_complete,
                "prospective_decision_records_complete": prospective_decision_records_complete,
                "calibrator_fitted_on_development": calibrator.fit_scope == "development_only",
                "locked_scores_finite": all_scores_finite,
                "discrepancy_reports_finite": discrepancy_reports_finite,
                "discrepancy_calibration_records_complete": discrepancy_calibration_records_complete,
                "discrepancy_calibrator_fit_scope_development_only": bool(
                    discrepancy_calibrator is not None
                    and discrepancy_calibrator.fit_scope == "development_only"
                ),
                "model_averaging_records_complete": model_averaging_records_complete,
                "model_averaging_fit_scope_development_only": bool(
                    model_weight_fit is not None
                    and model_weight_fit.fit_scope == "development_only"
                ),
                "generator_critic_complete": critic_complete,
                "generator_critic_gate": critic_gate,
                "primary_generator_critic_gate": primary_critic_gate,
                "observation_stress_complete": observation_stress_complete,
                "observation_stress_has_missingness_or_censoring": observation_stress_has_missingness,
                "no_recorded_errors": not errors,
                "model_averaging_converged": model_averaging_converged,
                "post_measurement_claim_gate": post_measurement_claim_gate,
                "age_claim_gate": age_gate.status,
                "reaction_claim_gate": reaction_gate.status,
                "kinetics_claim_gate": kinetics_gate.status,
                "integrated_claim_gate": integrated_gate.status,
            },
            "programme_workflow_status": (
                "INCOMPLETE_CRITIC_REVIEW" if not critic_gate else "SCOPED_EXECUTION"
            ),
            "critic_status": "PASS" if critic_gate else "REVIEW_REQUIRED",
            "critic_major_count": critic_major_count,
            "critic_blocker_count": critic_blocker_count,
            "primary_critic_status": "PASS" if primary_critic_gate else "REVIEW_REQUIRED",
            "missing_programme_stages": (
                [] if post_measurement_records_complete else ["prospective_measurement_outcome"]
            ),
            "performance_gate": (
                "PASS" if integrated_performance_gate
                else "DEFERRED_UNTIL_CONVERGED_MODEL_AVERAGING_AND_REALIZED_SYNTHETIC_OUTCOME_IMPROVEMENT"
            ),
            "field_validation": "DEFERRED_UNTIL_FIELD_DATA_READY",
            "synthetic_component_claim": "adjudicated_separately_on_locked_common-universe_metrics",
            "synthetic_integrated_claim": (
                "PASS" if integrated_performance_gate else "ABSTAIN"
            ),
            "specialist_claim_status": specialist_claim_summary["claim_status"],
        },
    )

    manifest: dict[str, object] = {
        "run_id": run_id,
        "status": (
            "PASS" if integrated_performance_gate
            else "PASS_SCOPED_EXECUTION" if execution_gate
            else "FAIL"
        ),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "git_worktree_dirty": _git_worktree_dirty(),
        "source_hashes": {
            relative: _sha256(REPO / relative) for relative in PROGRAMME_SOURCE_FILES
        },
        "python": sys.version,
        "platform": platform.platform(),
        "age_samples": int(age_samples),
        "development_cases": {
            family: list(seeds) for family, seeds in selected_development_cases.items()
        },
        "locked_cases": {
            family: list(seeds) for family, seeds in selected_locked_cases.items()
        },
        "selected_families": list(selected_families),
        "generator_families": sorted(present_families),
        "conditions": list(CONDITIONS),
        "discrepancy_scenarios": ["age_history_tropical"],
        "calibration": calibrator.to_dict(),
        "discrepancy_report_count": len(discrepancy_reports),
        "discrepancy_calibration": {
            "records_complete": discrepancy_calibration_records_complete,
            "fit_scope": (
                discrepancy_calibrator.fit_scope
                if discrepancy_calibrator is not None
                else None
            ),
            "fit_data_hash": (
                discrepancy_calibrator.fit_data_hash
                if discrepancy_calibrator is not None
                else None
            ),
            "locked_status": discrepancy_calibration_scores.get("status"),
        },
        "model_averaging": {
            "development_observation_count": len(development_model_observations),
            "locked_observation_count": len(locked_model_observations),
            "records_complete": model_averaging_records_complete,
            "candidate_universe_scope": "hydrosheaf_native_candidates_only",
            "converged": model_averaging_converged,
        },
        "specialist_candidate_generation": {
            "algorithm": "independent_geometry_head_knn_v1",
            "candidate_universe_scope": "blind_observation_geometry_and_heads",
            "records_complete": specialist_records_complete,
            "age_and_reaction_outputs": "scored_on_independent_universe_or_all_observed_nodes_as_declared",
        },
        "specialist_calibration": {
            "fit_scope": "development_only",
            "age_calibrator_count": len(specialist_age_calibrators),
            "reaction_calibrator_count": len(specialist_reaction_calibrators),
            "records_complete": specialist_calibration_records_complete,
            "claim_status": "calibrated_locked_diagnostic_not_superiority",
        },
        "reaction_rapm": {
            "fit_scope": "development_only",
            "calibration_prediction_protocol": "development_case_cross_fitted",
            "records_complete": reaction_rapm_records_complete,
            "calibration_records_complete": reaction_rapm_calibration_records_complete,
            "fit_case_count": (
                len(reaction_rapm_model.fit_case_ids)
                if reaction_rapm_model is not None
                else 0
            ),
            "fit_record_count": (
                reaction_rapm_model.fit_record_count
                if reaction_rapm_model is not None
                else 0
            ),
            "cv_log_loss": (
                reaction_rapm_model.cv_log_loss
                if reaction_rapm_model is not None
                else None
            ),
            "ridge_lambda": (
                reaction_rapm_model.ridge_lambda
                if reaction_rapm_model is not None
                else None
            ),
            "on_off_weight": (
                reaction_rapm_model.on_off_weight
                if reaction_rapm_model is not None
                else None
            ),
            "claim_status": "locked_diagnostic_not_superiority",
        },
        "generator_critic_gate": critic_gate,
        "primary_generator_critic_gate": primary_critic_gate,
        "generator_critic_major_count": critic_major_count,
        "generator_critic_blocker_count": critic_blocker_count,
        "external_solver_execution_gate": external_solver_execution_gate,
        "integrated_performance_gate": integrated_performance_gate,
        "specialist_claim_status": specialist_claim_summary["claim_status"],
        "specialist_claim_gates": specialist_claim_summary["gates"],
        "age_performance_records_complete": age_performance_records_complete,
        "reaction_performance_records_complete": reaction_performance_records_complete,
        "kinetics_performance_records_complete": kinetics_performance_records_complete,
        "prospective_decision_records_complete": prospective_decision_records_complete,
        "post_measurement_records_complete": post_measurement_records_complete,
        "post_measurement_claim_gate": post_measurement_claim_gate,
        "model_averaging_converged": model_averaging_converged,
        "observation_stress_scenarios": list(OBSERVATION_STRESS_SCENARIOS),
        "observation_stress_complete": observation_stress_complete,
        "observation_stress_has_missingness_or_censoring": observation_stress_has_missingness,
        "claim_boundary": (
            "held-out controlled-synthetic calibration and truth-blind one-step "
            "policy diagnostics with realized-outcome scoring when actions are "
            "selected; no field-validation or universal-superiority claim"
        ),
        "errors": errors,
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)  # type: ignore[index]
    _write_json(output / "run_manifest.json", manifest)
    return {
        "manifest": manifest,
        "execution_gate": execution_gate,
        "integrated_performance_gate": integrated_performance_gate,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-id",
        default=RUN_ID,
        help="immutable run identifier written to the programme and manifest records",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--simulator-workspace", type=Path, default=DEFAULT_SIMULATOR_WORKSPACE)
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument("--age-samples", type=int, default=400)
    parser.add_argument("--quick", action="store_true", help="use one case per family")
    parser.add_argument(
        "--families",
        nargs="+",
        choices=sorted(ALL_GENERATOR_FAMILIES),
        help=(
            "generator families to run; omit for all families. Analytic-only "
            "runs are scoped execution checks, not the full integrated gate."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    global DEVELOPMENT_CASES, LOCKED_CASES
    if args.quick:
        DEVELOPMENT_CASES = {
            "modflow_modpath_v3": (9101,),
            "analytic_lattice_v1": (9301,),
            "independent_mixing_v1": (9501,),
        }
        LOCKED_CASES = {
            "modflow_modpath_v3": (9201,),
            "analytic_lattice_v1": (9401,),
            "independent_mixing_v1": (9601,),
        }
    output = args.output
    if args.quick and output == DEFAULT_OUTPUT:
        output = DEFAULT_OUTPUT.parent / f"{RUN_ID}-SMOKE"
    result = run_ensemble_benchmark(
        run_id=args.run_id,
        output=output,
        simulator_workspace=args.simulator_workspace,
        mf6_executable=args.bin_dir / "mf6.exe",
        mp7_executable=args.bin_dir / "mp7.exe",
        age_samples=args.age_samples,
        overwrite=args.overwrite,
        families=args.families,
    )
    print(json.dumps(result, indent=2, sort_keys=True, default=_json_default))
    return 0 if bool(result["execution_gate"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
