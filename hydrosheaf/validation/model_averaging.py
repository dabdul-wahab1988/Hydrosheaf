"""Development-only probabilistic model averaging for validation outputs.

This module combines discrete predictive distributions with a proper log score;
it does not average point estimates or interval endpoints.  Weights are fitted
with equal contribution from each case (rather than each node), then frozen
before locked-test predictions are evaluated.  Material component disagreement
is retained in the audit record and suppresses a reportable aggregate decision.

The fitted weights are forecast-combination weights, not posterior model
probabilities.  They are useful for the controlled-synthetic validation
programme but do not establish a Bayesian model probability or field validity.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from typing import Any


_EPSILON = 1.0e-12


def _finite(value: object, *, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return result


def _json_hash(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class DiscreteModelObservation:
    """One truth-labelled target used only by development fitting or scoring."""

    case_id: str
    target_id: str
    truth: str
    predictions: Mapping[str, Mapping[str, float]]
    phase: str = "development"

    def __post_init__(self) -> None:
        case_id = str(self.case_id).strip()
        target_id = str(self.target_id).strip()
        truth = str(self.truth).strip()
        phase = str(self.phase).strip()
        if not case_id or not target_id or not truth:
            raise ValueError("Model observations require case_id, target_id, and truth.")
        if not phase:
            raise ValueError("Model observation phase must be non-empty.")
        normalised = _normalise_prediction_matrix(self.predictions)
        for model, values in normalised.items():
            if truth not in values:
                raise ValueError(
                    f"Truth outcome {truth!r} is absent from model {model!r}."
                )
        object.__setattr__(self, "case_id", case_id)
        object.__setattr__(self, "target_id", target_id)
        object.__setattr__(self, "truth", truth)
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "predictions", normalised)

    @property
    def model_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self.predictions))

    def to_dict(self, *, include_truth: bool = True) -> dict[str, Any]:
        record: dict[str, Any] = {
            "case_id": self.case_id,
            "target_id": self.target_id,
            "phase": self.phase,
            "predictions": {
                model: dict(self.predictions[model]) for model in self.model_ids
            },
        }
        if include_truth:
            record["truth"] = self.truth
        return record


@dataclass(frozen=True)
class ModelWeightFit:
    """Immutable development-only fit record for forecast-combination weights."""

    model_ids: tuple[str, ...]
    weights: Mapping[str, float]
    prior_weights: Mapping[str, float]
    case_ids: tuple[str, ...]
    observation_count: int
    score_rule: str
    fit_scope: str
    weight_floor: float
    iterations: int
    converged: bool
    fit_data_hash: str
    fit_log_score: float
    max_iterations: int = 0
    convergence_tolerance: float = 0.0
    gradient_tolerance: float = 0.0
    final_gradient_norm: float = 0.0
    kkt_residual: float = 0.0
    last_weight_change: float = 0.0
    simplex_residual: float = 0.0
    objective_initial: float = 0.0
    objective_delta: float = 0.0
    objective_trace: tuple[float, ...] = ()
    convergence_status: str = "ABSTAIN"
    convergence_reason: str = "not evaluated"

    def __post_init__(self) -> None:
        model_ids = tuple(str(model) for model in self.model_ids)
        if not model_ids or len(model_ids) != len(set(model_ids)):
            raise ValueError("A model fit requires unique model IDs.")
        if self.fit_scope != "development_only":
            raise ValueError("Model weights may only be fitted on development data.")
        floor = _finite(self.weight_floor, name="weight_floor")
        if floor < 0.0 or floor * len(model_ids) >= 1.0:
            raise ValueError("weight_floor must be non-negative and leave simplex mass.")
        weights = _normalise_weights(self.weights, model_ids, name="weights")
        priors = _normalise_weights(self.prior_weights, model_ids, name="prior_weights")
        if any(value < floor - 1.0e-10 for value in weights.values()):
            raise ValueError("Fitted weights must honour weight_floor.")
        if int(self.observation_count) <= 0 or not self.case_ids:
            raise ValueError("A model fit requires observations and cases.")
        if not self.score_rule:
            raise ValueError("A model fit requires a score rule.")
        if len(str(self.fit_data_hash)) != 64:
            raise ValueError("fit_data_hash must be a SHA-256 hex digest.")
        score = _finite(self.fit_log_score, name="fit_log_score")
        max_iterations = int(self.max_iterations)
        if max_iterations < 0:
            raise ValueError("max_iterations must be non-negative.")
        convergence_tolerance = _finite(
            self.convergence_tolerance, name="convergence_tolerance"
        )
        gradient_tolerance = _finite(
            self.gradient_tolerance, name="gradient_tolerance"
        )
        final_gradient_norm = _finite(
            self.final_gradient_norm, name="final_gradient_norm"
        )
        kkt_residual = _finite(self.kkt_residual, name="kkt_residual")
        last_weight_change = _finite(
            self.last_weight_change, name="last_weight_change"
        )
        simplex_residual = _finite(self.simplex_residual, name="simplex_residual")
        objective_initial = _finite(
            self.objective_initial, name="objective_initial"
        )
        objective_delta = _finite(self.objective_delta, name="objective_delta")
        trace = tuple(_finite(value, name="objective_trace") for value in self.objective_trace)
        status = str(self.convergence_status).strip().upper()
        if status not in {"CONVERGED", "ABSTAIN"}:
            raise ValueError("convergence_status must be CONVERGED or ABSTAIN.")
        if status == "CONVERGED" and not self.converged:
            raise ValueError("A CONVERGED status requires converged=True.")
        if status == "ABSTAIN" and self.converged:
            raise ValueError("An ABSTAIN status requires converged=False.")
        if any(value < 0.0 for value in (convergence_tolerance, gradient_tolerance,
                                         final_gradient_norm, kkt_residual,
                                         last_weight_change, simplex_residual)):
            raise ValueError("Convergence diagnostics cannot be negative.")
        object.__setattr__(self, "model_ids", model_ids)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "prior_weights", priors)
        object.__setattr__(self, "case_ids", tuple(sorted(str(case) for case in self.case_ids)))
        object.__setattr__(self, "weight_floor", floor)
        object.__setattr__(self, "observation_count", int(self.observation_count))
        object.__setattr__(self, "fit_log_score", score)
        object.__setattr__(self, "max_iterations", max_iterations)
        object.__setattr__(self, "convergence_tolerance", convergence_tolerance)
        object.__setattr__(self, "gradient_tolerance", gradient_tolerance)
        object.__setattr__(self, "final_gradient_norm", final_gradient_norm)
        object.__setattr__(self, "kkt_residual", kkt_residual)
        object.__setattr__(self, "last_weight_change", last_weight_change)
        object.__setattr__(self, "simplex_residual", simplex_residual)
        object.__setattr__(self, "objective_initial", objective_initial)
        object.__setattr__(self, "objective_delta", objective_delta)
        object.__setattr__(self, "objective_trace", trace)
        object.__setattr__(self, "convergence_status", status)
        object.__setattr__(self, "convergence_reason", str(self.convergence_reason))

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_ids": list(self.model_ids),
            "weights": dict(self.weights),
            "prior_weights": dict(self.prior_weights),
            "case_ids": list(self.case_ids),
            "observation_count": self.observation_count,
            "score_rule": self.score_rule,
            "fit_scope": self.fit_scope,
            "weight_floor": self.weight_floor,
            "iterations": self.iterations,
            "converged": self.converged,
            "fit_data_hash": self.fit_data_hash,
            "fit_log_score": self.fit_log_score,
            "max_iterations": self.max_iterations,
            "convergence_tolerance": self.convergence_tolerance,
            "gradient_tolerance": self.gradient_tolerance,
            "final_gradient_norm": self.final_gradient_norm,
            "kkt_residual": self.kkt_residual,
            "last_weight_change": self.last_weight_change,
            "simplex_residual": self.simplex_residual,
            "objective_initial": self.objective_initial,
            "objective_delta": self.objective_delta,
            "objective_trace": list(self.objective_trace),
            "convergence_status": self.convergence_status,
            "convergence_reason": self.convergence_reason,
        }


@dataclass(frozen=True)
class ModelAveragedDiscreteReport:
    """Audited model mixture with an explicit reportability decision."""

    target_id: str
    mixture_probabilities: Mapping[str, float]
    model_probabilities: Mapping[str, Mapping[str, float]]
    weights: Mapping[str, float]
    max_pairwise_total_variation: float
    disagreement_threshold: float
    reportable: bool
    decision: str
    identifiability: str
    compatible_outcomes: tuple[str, ...]
    reason: str
    fit_data_hash: str
    weighted_pairwise_total_variation: float = 0.0
    disagreement_status: str = "UNKNOWN"
    fit_converged: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_id": self.target_id,
            "mixture_probabilities": dict(self.mixture_probabilities),
            "model_probabilities": {
                model: dict(values) for model, values in self.model_probabilities.items()
            },
            "weights": dict(self.weights),
            "max_pairwise_total_variation": self.max_pairwise_total_variation,
            "disagreement_threshold": self.disagreement_threshold,
            "reportable": self.reportable,
            "decision": self.decision,
            "identifiability": self.identifiability,
            "compatible_outcomes": list(self.compatible_outcomes),
            "reason": self.reason,
            "fit_data_hash": self.fit_data_hash,
            "weighted_pairwise_total_variation": self.weighted_pairwise_total_variation,
            "disagreement_status": self.disagreement_status,
            "fit_converged": self.fit_converged,
        }


def fit_discrete_model_weights(
    observations: Iterable[DiscreteModelObservation],
    *,
    model_ids: Iterable[str] | None = None,
    weight_floor: float = 0.05,
    max_iterations: int = 2_000,
    tolerance: float = 1.0e-10,
    gradient_tolerance: float | None = None,
) -> ModelWeightFit:
    """Fit case-blocked log-score weights using development observations only."""

    records = tuple(observations)
    if not records:
        raise ValueError("At least one development observation is required.")
    if any(record.phase != "development" for record in records):
        raise ValueError("Locked or field observations cannot be used for fitting.")
    ids = tuple(sorted(str(model) for model in (model_ids or records[0].model_ids)))
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("model_ids must contain unique models.")
    for record in records:
        if record.model_ids != ids:
            raise ValueError(
                f"Model matrix is incomplete for {record.case_id}/{record.target_id}."
            )
    floor = _finite(weight_floor, name="weight_floor")
    if floor < 0.0 or floor * len(ids) >= 1.0:
        raise ValueError("weight_floor must be non-negative and leave simplex mass.")
    if max_iterations < 1 or tolerance <= 0.0:
        raise ValueError("max_iterations and tolerance must be positive.")
    gradient_tol = (
        max(1.0e-8, math.sqrt(tolerance))
        if gradient_tolerance is None
        else _finite(gradient_tolerance, name="gradient_tolerance")
    )
    if gradient_tol <= 0.0:
        raise ValueError("gradient_tolerance must be positive.")

    by_case: dict[str, list[DiscreteModelObservation]] = defaultdict(list)
    seen_targets: set[tuple[str, str]] = set()
    for record in records:
        key = (record.case_id, record.target_id)
        if key in seen_targets:
            raise ValueError(f"Duplicate case/target observation: {key!r}")
        seen_targets.add(key)
        by_case[record.case_id].append(record)
    case_ids = tuple(sorted(by_case))
    prior = {model: 1.0 / len(ids) for model in ids}
    remaining = 1.0 - floor * len(ids)
    simplex = {model: 1.0 / len(ids) for model in ids}

    def weights_from_simplex(values: Mapping[str, float]) -> dict[str, float]:
        return {model: floor + remaining * values[model] for model in ids}

    def objective(values: Mapping[str, float]) -> float:
        weights = weights_from_simplex(values)
        case_scores = []
        for case_id in case_ids:
            rows = by_case[case_id]
            row_scores = []
            for row in rows:
                mixture = sum(
                    weights[model] * row.predictions[model][row.truth]
                    for model in ids
                )
                row_scores.append(math.log(max(_EPSILON, mixture)))
            case_scores.append(sum(row_scores) / len(row_scores))
        return sum(case_scores) / len(case_scores)

    def gradient(values: Mapping[str, float]) -> dict[str, float]:
        weights = weights_from_simplex(values)
        case_gradients: list[dict[str, float]] = []
        for case_id in case_ids:
            rows = by_case[case_id]
            local = {model: 0.0 for model in ids}
            for row in rows:
                mixture = sum(
                    weights[model] * row.predictions[model][row.truth]
                    for model in ids
                )
                denominator = max(_EPSILON, mixture)
                for model in ids:
                    local[model] += (
                        remaining * row.predictions[model][row.truth] / denominator
                    )
            case_gradients.append(
                {model: value / len(rows) for model, value in local.items()}
            )
        return {
            model: sum(local[model] for local in case_gradients) / len(case_gradients)
            for model in ids
        }

    current_score = objective(simplex)
    objective_initial = current_score
    objective_trace = [current_score]
    step = 0.5
    converged = False
    iterations = 0
    last_weight_change = 0.0
    convergence_reason = "max_iterations_reached"

    def kkt_residual(values: Mapping[str, float], gradients: Mapping[str, float]) -> float:
        free = [model for model in ids if values[model] > 1.0e-12]
        reference = max(gradients[model] for model in free or ids)
        residuals = []
        for model in ids:
            if values[model] > 1.0e-12:
                residuals.append(abs(gradients[model] - reference))
            else:
                residuals.append(max(0.0, gradients[model] - reference))
        return max(residuals, default=0.0)

    for iteration in range(1, max_iterations + 1):
        iterations = iteration
        grad = gradient(simplex)
        accepted = False
        trial_step = step
        next_simplex = simplex
        next_score = current_score
        for _ in range(60):
            exponent = {
                model: math.log(max(_EPSILON, simplex[model]))
                + trial_step * grad[model]
                for model in ids
            }
            maximum = max(exponent.values())
            unnormalised = {
                model: math.exp(value - maximum) for model, value in exponent.items()
            }
            total = sum(unnormalised.values())
            candidate = {model: value / total for model, value in unnormalised.items()}
            candidate_score = objective(candidate)
            if candidate_score + max(1.0e-15, tolerance * 1.0e-3) >= current_score:
                next_simplex = candidate
                next_score = candidate_score
                accepted = True
                break
            trial_step *= 0.5
        if not accepted:
            convergence_reason = "line_search_stalled"
            break
        change = sum(abs(next_simplex[model] - simplex[model]) for model in ids)
        last_weight_change = change
        simplex = next_simplex
        current_score = next_score
        objective_trace.append(current_score)
        step = min(8.0, max(1.0e-6, trial_step * 1.25))
        residual = kkt_residual(simplex, gradient(simplex))
        if change <= tolerance and residual <= gradient_tol:
            converged = True
            convergence_reason = "weight_change_and_kkt_residual_within_tolerance"
            break

    # The multiplicative mirror-descent path above is retained as a portable
    # first pass, but it can stagnate near a simplex face when model scores are
    # almost tied: the objective change becomes tiny while the KKT residual is
    # still material.  Resolve that state with the repository's declared
    # SciPy dependency using the exact same case-blocked objective and
    # analytic gradient.  A fit is still marked ABSTAIN unless the optimizer
    # reports success and the post-fit KKT residual satisfies the frozen gate.
    if not converged and max_iterations > 1:
        try:
            from scipy.optimize import minimize, root

            model_index = {model: index for index, model in enumerate(ids)}
            start = [float(simplex[model]) for model in ids]
            callback_trace: list[float] = []

            def vector_objective(values: object) -> float:
                vector = [float(value) for value in values]  # type: ignore[union-attr]
                return -objective(
                    {model: vector[index] for model, index in model_index.items()}
                )

            def vector_gradient(values: object) -> list[float]:
                vector = [float(value) for value in values]  # type: ignore[union-attr]
                local = gradient(
                    {model: vector[index] for model, index in model_index.items()}
                )
                return [-local[model] for model in ids]

            def callback(values: object) -> None:
                vector = [float(value) for value in values]  # type: ignore[union-attr]
                callback_trace.append(-vector_objective(vector))

            result = minimize(
                vector_objective,
                start,
                jac=vector_gradient,
                method="SLSQP",
                bounds=[(0.0, 1.0) for _ in ids],
                constraints={
                    "type": "eq",
                    "fun": lambda values: sum(float(value) for value in values) - 1.0,
                    "jac": lambda values: [1.0 for _ in values],
                },
                callback=callback,
                options={
                    "maxiter": int(max_iterations),
                    "ftol": max(float(tolerance), 1.0e-12),
                    "disp": False,
                },
            )
            vector = [float(value) for value in result.x]
            total = sum(vector)
            if total > 0.0 and all(math.isfinite(value) and value >= 0.0 for value in vector):
                simplex = {
                    model: max(0.0, vector[index] / total)
                    for model, index in model_index.items()
                }
                current_score = objective(simplex)
                if callback_trace:
                    objective_trace.extend(callback_trace)
                last_weight_change = sum(
                    abs(simplex[model] - start[index])
                    for model, index in model_index.items()
                )
                iterations = int(result.nit)
                final_candidate_kkt = kkt_residual(simplex, gradient(simplex))
                converged = bool(result.success) and final_candidate_kkt <= gradient_tol
                convergence_reason = (
                    "scipy_slsqp_success_with_kkt_residual_within_tolerance"
                    if converged
                    else f"scipy_slsqp_not_converged:{str(result.message)}"
                )

            # SLSQP can report success on this deliberately flat objective
            # before the absolute KKT residual reaches the strict audit
            # tolerance.  Solve the interior simplex stationarity equations
            # directly in log-ratio coordinates; this preserves the simplex
            # exactly and does not relax the convergence gate.
            if not converged and all(simplex[model] > 0.0 for model in ids):
                def logits_to_simplex(logits: object) -> dict[str, float]:
                    raw = [float(value) for value in logits]  # type: ignore[union-attr]
                    extended = [*raw, 0.0]
                    maximum = max(extended)
                    exponentials = [math.exp(value - maximum) for value in extended]
                    total = sum(exponentials)
                    values = [value / total for value in exponentials]
                    return {
                        model: values[index]
                        for model, index in model_index.items()
                    }

                reference_model = ids[-1]

                def stationarity(logits: object) -> list[float]:
                    candidate = logits_to_simplex(logits)
                    candidate_gradient = gradient(candidate)
                    return [
                        candidate_gradient[model] - candidate_gradient[reference_model]
                        for model in ids[:-1]
                    ]

                start_logits = [
                    math.log(max(simplex[model], _EPSILON) / max(simplex[reference_model], _EPSILON))
                    for model in ids[:-1]
                ]
                root_result = root(
                    stationarity,
                    start_logits,
                    method="hybr",
                    options={"xtol": min(1.0e-10, math.sqrt(tolerance))},
                )
                if bool(root_result.success):
                    candidate_simplex = logits_to_simplex(root_result.x)
                    candidate_weights = weights_from_simplex(candidate_simplex)
                    candidate_kkt = kkt_residual(
                        candidate_simplex,
                        gradient(candidate_simplex),
                    )
                    if (
                        candidate_kkt <= gradient_tol
                        and all(value >= floor - 1.0e-10 for value in candidate_weights.values())
                    ):
                        simplex = candidate_simplex
                        fitted_candidate_score = objective(simplex)
                        objective_trace.append(fitted_candidate_score)
                        current_score = fitted_candidate_score
                        last_weight_change = sum(
                            abs(simplex[model] - start[index])
                            for model, index in model_index.items()
                        )
                        iterations = max(iterations, int(getattr(root_result, "nfev", 0)))
                        converged = True
                        convergence_reason = "interior_stationarity_root_within_tolerance"
        except Exception as exc:  # preserve the explicit fail-closed status
            convergence_reason = f"fallback_optimizer_error:{type(exc).__name__}"

    fitted = weights_from_simplex(simplex)
    final_gradient = gradient(simplex)
    final_gradient_norm = math.sqrt(
        sum(value * value for value in final_gradient.values())
    )
    final_kkt_residual = kkt_residual(simplex, final_gradient)
    if not converged and final_kkt_residual <= gradient_tol and last_weight_change <= tolerance:
        converged = True
        convergence_reason = "final_kkt_residual_within_tolerance"
    if converged and convergence_reason == "max_iterations_reached":
        convergence_reason = "convergence_diagnostics_within_tolerance"
    fit_payload = [record.to_dict(include_truth=True) for record in sorted(records, key=lambda item: (item.case_id, item.target_id))]
    return ModelWeightFit(
        model_ids=ids,
        weights=fitted,
        prior_weights=prior,
        case_ids=case_ids,
        observation_count=len(records),
        score_rule="case_blocked_mean_log_score",
        fit_scope="development_only",
        weight_floor=floor,
        iterations=iterations,
        converged=converged,
        fit_data_hash=_json_hash(fit_payload),
        fit_log_score=current_score,
        max_iterations=max_iterations,
        convergence_tolerance=tolerance,
        gradient_tolerance=gradient_tol,
        final_gradient_norm=final_gradient_norm,
        kkt_residual=final_kkt_residual,
        last_weight_change=last_weight_change,
        simplex_residual=abs(sum(fitted.values()) - 1.0),
        objective_initial=objective_initial,
        objective_delta=current_score - objective_initial,
        objective_trace=tuple(objective_trace),
        convergence_status="CONVERGED" if converged else "ABSTAIN",
        convergence_reason=convergence_reason,
    )


def apply_discrete_model_average(
    target_id: str,
    predictions: Mapping[str, Mapping[str, float]],
    fit: ModelWeightFit,
    *,
    disagreement_threshold: float = 0.25,
    report_threshold: float = 0.60,
) -> ModelAveragedDiscreteReport:
    """Apply frozen weights and suppress the aggregate on material disagreement."""

    target = str(target_id).strip()
    if not target:
        raise ValueError("target_id must be non-empty.")
    threshold = _finite(disagreement_threshold, name="disagreement_threshold")
    report_threshold = _finite(report_threshold, name="report_threshold")
    if threshold < 0.0 or not 0.0 < report_threshold <= 1.0:
        raise ValueError("Invalid disagreement or report threshold.")
    normalised_predictions = _normalise_prediction_matrix(predictions)
    if tuple(sorted(normalised_predictions)) != fit.model_ids:
        raise ValueError(
            f"Prediction model IDs {tuple(sorted(normalised_predictions))!r} do not match "
            f"fit {fit.model_ids!r}."
        )
    outcomes = tuple(sorted(next(iter(normalised_predictions.values())).keys()))
    mixture = {
        outcome: sum(
            fit.weights[model] * normalised_predictions[model][outcome]
            for model in fit.model_ids
        )
        for outcome in outcomes
    }
    max_tv = 0.0
    weighted_tv = 0.0
    for index, left in enumerate(fit.model_ids):
        for right in fit.model_ids[index + 1 :]:
            pairwise_tv = 0.5 * sum(
                abs(
                    normalised_predictions[left][outcome]
                    - normalised_predictions[right][outcome]
                )
                for outcome in outcomes
            )
            max_tv = max(max_tv, pairwise_tv)
            weighted_tv += fit.weights[left] * fit.weights[right] * pairwise_tv
    disagreement_status = (
        "MODEL_DISAGREEMENT" if max_tv >= threshold else "NO_MATERIAL_DISAGREEMENT"
    )
    reportable = fit.converged and max_tv < threshold
    top = max(mixture, key=mixture.get)
    compatible = tuple(
        outcome for outcome in outcomes if mixture[outcome] >= 1.0 / len(outcomes)
    )
    if not fit.converged:
        decision = "ABSTAIN"
        identifiability = "FIT_NOT_CONVERGED"
        reason = (
            "Frozen model-averaging weights failed the numerical convergence "
            f"contract ({fit.convergence_status}; {fit.convergence_reason}; "
            f"KKT residual={fit.kkt_residual:.6g})."
        )
        compatible = ()
    elif not reportable:
        decision = "ABSTAIN"
        identifiability = "MODEL_DISAGREEMENT"
        reason = (
            "Component predictive distributions disagree beyond the frozen "
            f"total-variation threshold ({max_tv:.6g} >= {threshold:.6g}); "
            "the mixture is retained for audit but is not reportable."
        )
        compatible = ()
    elif mixture[top] >= report_threshold:
        decision = "SET_REPORT"
        identifiability = "PARTIALLY_IDENTIFIED"
        reason = "Frozen development weights produced a reportable mixture."
    else:
        decision = "ABSTAIN"
        identifiability = "UNKNOWN"
        reason = "The weighted mixture did not clear the report threshold."
    return ModelAveragedDiscreteReport(
        target_id=target,
        mixture_probabilities=mixture,
        model_probabilities={model: dict(normalised_predictions[model]) for model in fit.model_ids},
        weights=dict(fit.weights),
        max_pairwise_total_variation=max_tv,
        disagreement_threshold=threshold,
        reportable=reportable,
        decision=decision,
        identifiability=identifiability,
        compatible_outcomes=compatible,
        reason=reason,
        fit_data_hash=fit.fit_data_hash,
        weighted_pairwise_total_variation=weighted_tv,
        disagreement_status=disagreement_status,
        fit_converged=fit.converged,
    )


def score_locked_model_average(
    observations: Iterable[DiscreteModelObservation],
    fit: ModelWeightFit,
    *,
    disagreement_threshold: float = 0.25,
    report_threshold: float = 0.60,
) -> dict[str, Any]:
    """Score a frozen model average on locked observations only."""

    records = tuple(observations)
    if not records:
        return {"status": "not_available", "n": 0, "held_out_scoring": True}
    if any(record.phase != "locked_test" for record in records):
        raise ValueError("Locked scoring accepts locked_test observations only.")
    log_scores: list[float] = []
    brier_scores: list[float] = []
    disagreement = 0
    reports = []
    per_model_log_scores = {model: [] for model in fit.model_ids}
    per_model_brier_scores = {model: [] for model in fit.model_ids}
    disagreement_status_counts: dict[str, int] = {}
    for record in records:
        report = apply_discrete_model_average(
            record.target_id,
            record.predictions,
            fit,
            disagreement_threshold=disagreement_threshold,
            report_threshold=report_threshold,
        )
        reports.append(report)
        probability = max(_EPSILON, report.mixture_probabilities[record.truth])
        log_scores.append(math.log(probability))
        brier_scores.append(
            sum(
                (report.mixture_probabilities[outcome] - (1.0 if outcome == record.truth else 0.0)) ** 2
                for outcome in report.mixture_probabilities
            )
        )
        for model in fit.model_ids:
            model_probability = max(_EPSILON, record.predictions[model][record.truth])
            per_model_log_scores[model].append(math.log(model_probability))
            per_model_brier_scores[model].append(
                sum(
                    (record.predictions[model][outcome]
                     - (1.0 if outcome == record.truth else 0.0)) ** 2
                    for outcome in record.predictions[model]
                )
            )
        disagreement_status_counts[report.disagreement_status] = (
            disagreement_status_counts.get(report.disagreement_status, 0) + 1
        )
        disagreement += int(not report.reportable)
    return {
        "status": "scored" if fit.converged else "ABSTAIN_FIT_NOT_CONVERGED",
        "held_out_scoring": True,
        "fit_converged": fit.converged,
        "fit_convergence_status": fit.convergence_status,
        "n": len(records),
        "mean_log_score": sum(log_scores) / len(log_scores),
        "mean_brier_score": sum(brier_scores) / len(brier_scores),
        "disagreement_rate": disagreement / len(records),
        "reportable_rate": 1.0 - disagreement / len(records),
        "per_model_mean_log_score": {
            model: sum(values) / len(values)
            for model, values in per_model_log_scores.items()
        },
        "per_model_mean_brier_score": {
            model: sum(values) / len(values)
            for model, values in per_model_brier_scores.items()
        },
        "disagreement_status_counts": disagreement_status_counts,
        "reports": [report.to_dict() for report in reports],
        "fit_data_hash": fit.fit_data_hash,
    }


def _normalise_weights(
    values: Mapping[str, float], model_ids: tuple[str, ...], *, name: str
) -> dict[str, float]:
    if set(values) != set(model_ids):
        raise ValueError(f"{name} must contain exactly the fitted model IDs.")
    result = {model: _finite(values[model], name=f"{name}[{model}]") for model in model_ids}
    if any(value < 0.0 for value in result.values()):
        raise ValueError(f"{name} cannot contain negative values.")
    total = sum(result.values())
    if total <= 0.0 or not math.isclose(total, 1.0, rel_tol=1.0e-8, abs_tol=1.0e-8):
        raise ValueError(f"{name} must sum to 1.")
    return result


def _normalise_prediction_matrix(
    predictions: Mapping[str, Mapping[str, float]] | object,
) -> dict[str, dict[str, float]]:
    if not isinstance(predictions, Mapping) or not predictions:
        raise ValueError("At least one model prediction is required.")
    normalised: dict[str, dict[str, float]] = {}
    for model_id, probabilities in predictions.items():
        model = str(model_id).strip()
        if not model:
            raise ValueError("Model IDs must be non-empty.")
        if not isinstance(probabilities, Mapping) or not probabilities:
            raise ValueError(f"Model {model!r} needs a non-empty distribution.")
        values: dict[str, float] = {}
        for outcome, probability in probabilities.items():
            outcome_id = str(outcome).strip()
            value = _finite(probability, name=f"{model}[{outcome_id}]")
            if not outcome_id or value < 0.0:
                raise ValueError(f"Invalid probability for {model!r}/{outcome!r}.")
            values[outcome_id] = value
        total = sum(values.values())
        if total <= 0.0 or not math.isclose(total, 1.0, rel_tol=1.0e-8, abs_tol=1.0e-8):
            raise ValueError(f"Probabilities for model {model!r} must sum to 1.")
        normalised[model] = values
    outcome_sets = {frozenset(values) for values in normalised.values()}
    if len(outcome_sets) != 1:
        raise ValueError("All models must use the same outcome set per target.")
    return normalised


__all__ = [
    "DiscreteModelObservation",
    "ModelAveragedDiscreteReport",
    "ModelWeightFit",
    "apply_discrete_model_average",
    "fit_discrete_model_weights",
    "score_locked_model_average",
]
