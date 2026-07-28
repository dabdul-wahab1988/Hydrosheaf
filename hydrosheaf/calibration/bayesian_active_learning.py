"""Bayesian value-of-information design for groundwater measurements.

The module ranks explicit measurement actions against a discrete ensemble of
scientific hypotheses.  An action is represented by one or more predictive
Gaussian scenarios, its resource cost, and its feasibility.  Ranking uses
expected reduction in Shannon entropy, calculated by deterministic
Gauss--Hermite quadrature.  Scenario robustness, cost adjustment, abstention,
posterior updating, and non-redundant batch selection are part of the public
contract.

The implementation is deliberately domain-neutral: a hypothesis may represent
a complete flow topology, a reaction mechanism, or a parameter regime.  The
caller owns the forward predictions and observation-error model.  Consequently
the routine is an experimental-design engine, not evidence that any supplied
model is scientifically correct.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import math
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp, ndtri
from scipy.stats import qmc


_PROBABILITY_FLOOR = 1.0e-300


@dataclass(frozen=True)
class PredictiveScenario:
    """Conditional Gaussian predictions for one model-discrepancy scenario.

    ``means[h]`` and ``standard_deviations[h]`` define the predictive
    distribution of the proposed observation under hypothesis ``h``.
    Scenario weights are used for model-averaged expected information and
    posterior updating; the worst scenario is retained separately for robust
    acquisition.
    """

    name: str
    means: Sequence[float]
    standard_deviations: Sequence[float] | float
    weight: float = 1.0


@dataclass(frozen=True)
class MeasurementOption:
    """A concrete candidate measurement that the workflow can execute."""

    option_id: str
    measurement_type: str
    target_id: str
    cost: float
    scenarios: Sequence[PredictiveScenario]
    feasible: bool = True
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class AcquisitionConfig:
    """Controls for robust, cost-aware acquisition and abstention."""

    quadrature_order: int = 21
    robustness_weight: float = 0.5
    decision_weight: float = 0.0
    cost_exponent: float = 1.0
    minimum_expected_information_gain: float = 1.0e-4
    batch_qmc_samples: int = 2048
    random_seed: int = 20260728


def _normalise_probabilities(
    values: Sequence[float],
    *,
    name: str,
    minimum_size: int = 2,
) -> np.ndarray:
    probabilities = np.asarray(values, dtype=float)
    if probabilities.ndim != 1 or probabilities.size < int(minimum_size):
        raise ValueError(
            f"{name} must contain at least {int(minimum_size)} probabilities."
        )
    if not np.all(np.isfinite(probabilities)) or np.any(probabilities < 0.0):
        raise ValueError(f"{name} must be finite and non-negative.")
    total = float(probabilities.sum())
    if total <= 0.0:
        raise ValueError(f"{name} must have positive total mass.")
    return probabilities / total


def _validate_config(config: AcquisitionConfig) -> None:
    if int(config.quadrature_order) < 5:
        raise ValueError("quadrature_order must be at least 5.")
    if not 0.0 <= float(config.robustness_weight) <= 1.0:
        raise ValueError("robustness_weight must lie in [0, 1].")
    if not 0.0 <= float(config.decision_weight) <= 1.0:
        raise ValueError("decision_weight must lie in [0, 1].")
    if float(config.cost_exponent) < 0.0:
        raise ValueError("cost_exponent must be non-negative.")
    if float(config.minimum_expected_information_gain) < 0.0:
        raise ValueError("minimum_expected_information_gain must be non-negative.")
    if int(config.batch_qmc_samples) < 64:
        raise ValueError("batch_qmc_samples must be at least 64.")


def _scenario_arrays(
    scenario: PredictiveScenario,
    n_hypotheses: int,
) -> tuple[np.ndarray, np.ndarray]:
    means = np.asarray(scenario.means, dtype=float)
    standard_deviations = np.asarray(scenario.standard_deviations, dtype=float)
    if standard_deviations.ndim == 0:
        standard_deviations = np.full(n_hypotheses, float(standard_deviations))
    if means.shape != (n_hypotheses,) or standard_deviations.shape != (
        n_hypotheses,
    ):
        raise ValueError(
            f"Scenario {scenario.name!r} must provide one mean and standard "
            "deviation per hypothesis."
        )
    if not np.all(np.isfinite(means)):
        raise ValueError(f"Scenario {scenario.name!r} has non-finite means.")
    if not np.all(np.isfinite(standard_deviations)) or np.any(
        standard_deviations <= 0.0
    ):
        raise ValueError(
            f"Scenario {scenario.name!r} standard deviations must be finite "
            "and strictly positive."
        )
    if not math.isfinite(float(scenario.weight)) or float(scenario.weight) <= 0.0:
        raise ValueError(f"Scenario {scenario.name!r} weight must be positive.")
    return means, standard_deviations


def _validate_option(option: MeasurementOption, n_hypotheses: int) -> None:
    if not str(option.option_id).strip():
        raise ValueError("option_id must be non-empty.")
    if not str(option.measurement_type).strip():
        raise ValueError(f"Option {option.option_id!r} has no measurement_type.")
    if not str(option.target_id).strip():
        raise ValueError(f"Option {option.option_id!r} has no target_id.")
    if not math.isfinite(float(option.cost)) or float(option.cost) <= 0.0:
        raise ValueError(f"Option {option.option_id!r} cost must be positive.")
    if not option.scenarios:
        raise ValueError(f"Option {option.option_id!r} requires a predictive scenario.")
    names = [scenario.name for scenario in option.scenarios]
    if len(set(names)) != len(names):
        raise ValueError(f"Option {option.option_id!r} scenario names must be unique.")
    for scenario in option.scenarios:
        _scenario_arrays(scenario, n_hypotheses)


def shannon_entropy(probabilities: Sequence[float]) -> float:
    """Return Shannon entropy in nats for a discrete probability vector."""

    p = _normalise_probabilities(probabilities, name="probabilities")
    positive = p > 0.0
    return float(-np.sum(p[positive] * np.log(p[positive])))


def expected_information_gain(
    prior_probabilities: Sequence[float],
    means: Sequence[float],
    standard_deviations: Sequence[float] | float,
    *,
    quadrature_order: int = 21,
) -> float:
    """Expected KL gain for a scalar Gaussian observation, in nats.

    The expectation is evaluated by Gauss--Hermite quadrature under every
    hypothesis and then weighted by the current posterior probability.  This
    is deterministic and avoids a nested Monte Carlo bias for the scalar
    acquisition used by :func:`rank_measurement_options`.
    """

    prior = _normalise_probabilities(prior_probabilities, name="prior_probabilities")
    scenario = PredictiveScenario(
        name="evaluation",
        means=means,
        standard_deviations=standard_deviations,
    )
    mu, sigma = _scenario_arrays(scenario, len(prior))
    # If several scientific hypotheses make exactly the same predictive
    # statement, the observation cannot distinguish within that equivalence
    # class. Collapsing identical (mean, sigma) pairs and summing their prior
    # mass is therefore exact for mutual information and avoids repeatedly
    # integrating hundreds of topology particles that encode only two local
    # edge states for a particular action.
    predictive_classes, inverse = np.unique(
        np.column_stack((mu, sigma)),
        axis=0,
        return_inverse=True,
    )
    if len(predictive_classes) < len(prior):
        prior = np.bincount(inverse, weights=prior).astype(float)
        prior /= prior.sum()
        mu = predictive_classes[:, 0]
        sigma = predictive_classes[:, 1]
    if len(prior) == 1:
        return 0.0
    if int(quadrature_order) < 5:
        raise ValueError("quadrature_order must be at least 5.")

    nodes, weights = hermgauss(int(quadrature_order))
    quadrature_weights = weights / math.sqrt(math.pi)
    log_prior = np.log(np.clip(prior, _PROBABILITY_FLOOR, None))
    prior_entropy = shannon_entropy(prior)
    expected_posterior_entropy = 0.0

    for true_index, true_probability in enumerate(prior):
        observations = mu[true_index] + math.sqrt(2.0) * sigma[true_index] * nodes
        residuals = (observations[:, None] - mu[None, :]) / sigma[None, :]
        log_likelihood = (
            -0.5 * residuals**2
            - np.log(sigma[None, :])
            - 0.5 * math.log(2.0 * math.pi)
        )
        log_posterior = log_likelihood + log_prior[None, :]
        log_posterior -= logsumexp(log_posterior, axis=1, keepdims=True)
        posterior = np.exp(log_posterior)
        entropy = -np.sum(
            posterior * np.where(posterior > 0.0, log_posterior, 0.0),
            axis=1,
        )
        expected_posterior_entropy += float(
            true_probability * np.dot(quadrature_weights, entropy)
        )

    # Numerical quadrature may produce a negative value at round-off scale.
    return float(max(0.0, prior_entropy - expected_posterior_entropy))


def expected_brier_risk_reduction(
    prior_probabilities: Sequence[float],
    means: Sequence[float],
    standard_deviations: Sequence[float] | float,
    decision_values: Sequence[Sequence[float]],
    *,
    quadrature_order: int = 21,
) -> float:
    """Expected reduction in mean squared Bayes risk for decision targets.

    ``decision_values[h, j]`` is the value of decision target ``j`` under
    hypothesis ``h``. For binary topology indicators this is the expected
    reduction in posterior Brier risk. The calculation collapses hypotheses
    with identical predictive likelihoods while retaining their conditional
    target means, so the reduction is exact for the declared Gaussian model.
    """

    prior = _normalise_probabilities(prior_probabilities, name="prior_probabilities")
    values = np.asarray(decision_values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if values.ndim != 2 or values.shape[0] != len(prior) or values.shape[1] < 1:
        raise ValueError(
            "decision_values must have one finite row per hypothesis and at "
            "least one decision target."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("decision_values must be finite.")
    scenario = PredictiveScenario(
        name="decision_evaluation",
        means=means,
        standard_deviations=standard_deviations,
    )
    mu, sigma = _scenario_arrays(scenario, len(prior))
    if int(quadrature_order) < 5:
        raise ValueError("quadrature_order must be at least 5.")

    predictive_classes, inverse = np.unique(
        np.column_stack((mu, sigma)), axis=0, return_inverse=True
    )
    raw_class_prior = np.bincount(inverse, weights=prior).astype(float)
    active_classes = np.flatnonzero(raw_class_prior > 1.0e-15)
    n_classes = len(active_classes)
    if n_classes <= 1:
        return 0.0
    class_prior = raw_class_prior[active_classes]
    class_prior /= class_prior.sum()
    class_decisions = np.zeros((n_classes, values.shape[1]), dtype=float)
    for compact_index, class_index in enumerate(active_classes):
        mask = inverse == class_index
        conditional = prior[mask] / raw_class_prior[class_index]
        class_decisions[compact_index] = conditional @ values[mask]
    prior_decision = class_prior @ class_decisions

    class_means = predictive_classes[active_classes, 0]
    class_sds = predictive_classes[active_classes, 1]
    nodes, weights = hermgauss(int(quadrature_order))
    quadrature_weights = weights / math.sqrt(math.pi)
    log_class_prior = np.log(np.clip(class_prior, _PROBABILITY_FLOOR, None))
    expected_reduction = 0.0
    for true_index, true_probability in enumerate(class_prior):
        observations = (
            class_means[true_index]
            + math.sqrt(2.0) * class_sds[true_index] * nodes
        )
        residuals = (
            observations[:, None] - class_means[None, :]
        ) / class_sds[None, :]
        log_likelihood = (
            -0.5 * residuals**2
            - np.log(class_sds[None, :])
            - 0.5 * math.log(2.0 * math.pi)
        )
        log_posterior = log_likelihood + log_class_prior[None, :]
        log_posterior -= logsumexp(log_posterior, axis=1, keepdims=True)
        posterior_decision = np.exp(log_posterior) @ class_decisions
        squared_shift = np.mean((posterior_decision - prior_decision) ** 2, axis=1)
        expected_reduction += float(
            true_probability * np.dot(quadrature_weights, squared_shift)
        )
    return float(max(0.0, expected_reduction))


def _prior_decision_risk(prior: np.ndarray, decision_values: np.ndarray) -> float:
    mean = prior @ decision_values
    return float(np.mean(prior @ (decision_values**2) - mean**2))


def _scenario_eigs(
    prior: np.ndarray,
    option: MeasurementOption,
    config: AcquisitionConfig,
    decision_values: np.ndarray | None,
    prior_decision_risk: float,
) -> tuple[
    list[dict[str, float | str]],
    float,
    float,
    float,
    float,
    float,
    float,
]:
    values: list[dict[str, float | str]] = []
    scenario_weights = []
    eigs = []
    risk_reductions = []
    utilities = []
    prior_entropy = shannon_entropy(prior)
    for scenario in option.scenarios:
        means, standard_deviations = _scenario_arrays(scenario, len(prior))
        eig = expected_information_gain(
            prior,
            means,
            standard_deviations,
            quadrature_order=config.quadrature_order,
        )
        risk_reduction = (
            expected_brier_risk_reduction(
                prior,
                means,
                standard_deviations,
                decision_values,
                quadrature_order=config.quadrature_order,
            )
            if decision_values is not None
            else 0.0
        )
        information_fraction = eig / max(prior_entropy, 1.0e-15)
        decision_fraction = risk_reduction / max(prior_decision_risk, 1.0e-15)
        utility = (
            (1.0 - config.decision_weight) * information_fraction
            + config.decision_weight * decision_fraction
        )
        values.append(
            {
                "scenario": scenario.name,
                "expected_information_gain": eig,
                "expected_brier_risk_reduction": risk_reduction,
                "normalised_information_utility": information_fraction,
                "normalised_decision_utility": decision_fraction,
                "combined_utility": utility,
            }
        )
        scenario_weights.append(float(scenario.weight))
        eigs.append(eig)
        risk_reductions.append(risk_reduction)
        utilities.append(utility)
    weights = _normalise_probabilities(
        scenario_weights,
        name="scenario_weights",
        minimum_size=1,
    )
    mean_eig = float(np.dot(weights, np.asarray(eigs, dtype=float)))
    worst_eig = float(min(eigs))
    robust_eig = float(
        (1.0 - config.robustness_weight) * mean_eig
        + config.robustness_weight * worst_eig
    )
    mean_risk_reduction = float(np.dot(weights, np.asarray(risk_reductions)))
    worst_risk_reduction = float(min(risk_reductions))
    mean_utility = float(np.dot(weights, np.asarray(utilities)))
    worst_utility = float(min(utilities))
    robust_utility = float(
        (1.0 - config.robustness_weight) * mean_utility
        + config.robustness_weight * worst_utility
    )
    return (
        values,
        mean_eig,
        worst_eig,
        robust_eig,
        mean_risk_reduction,
        worst_risk_reduction,
        robust_utility,
    )


def rank_measurement_options(
    hypothesis_ids: Sequence[str],
    prior_probabilities: Sequence[float],
    options: Sequence[MeasurementOption],
    *,
    decision_values: Sequence[Sequence[float]] | None = None,
    config: AcquisitionConfig | None = None,
) -> dict[str, Any]:
    """Rank feasible actions by scenario-robust EIG per unit cost.

    The returned status is ``ACTIONABLE`` only when the best action exceeds the
    configured minimum robust expected information gain.  Ties are resolved by
    ``option_id`` so recommendations are stable across processes and platforms.
    """

    active_config = config or AcquisitionConfig()
    _validate_config(active_config)
    ids = [str(value) for value in hypothesis_ids]
    if len(ids) != len(set(ids)) or len(ids) < 2:
        raise ValueError("hypothesis_ids must contain at least two unique values.")
    prior = _normalise_probabilities(
        prior_probabilities, name="prior_probabilities"
    )
    if len(prior) != len(ids):
        raise ValueError("hypothesis_ids and prior_probabilities lengths differ.")
    decision_array: np.ndarray | None = None
    prior_decision_risk = 0.0
    if decision_values is not None:
        decision_array = np.asarray(decision_values, dtype=float)
        if decision_array.ndim == 1:
            decision_array = decision_array[:, None]
        if (
            decision_array.ndim != 2
            or decision_array.shape[0] != len(prior)
            or decision_array.shape[1] < 1
            or not np.all(np.isfinite(decision_array))
        ):
            raise ValueError(
                "decision_values must be finite with one row per hypothesis."
            )
        prior_decision_risk = _prior_decision_risk(prior, decision_array)
    elif active_config.decision_weight > 0.0:
        raise ValueError("decision_values are required when decision_weight is positive.")

    ranked: list[dict[str, Any]] = []
    excluded: list[dict[str, str]] = []
    seen_options: set[str] = set()
    for option in options:
        _validate_option(option, len(prior))
        if option.option_id in seen_options:
            raise ValueError(f"Duplicate option_id: {option.option_id}")
        seen_options.add(option.option_id)
        if not option.feasible:
            excluded.append({"option_id": option.option_id, "reason": "infeasible"})
            continue
        (
            scenario_values,
            mean_eig,
            worst_eig,
            robust_eig,
            mean_risk_reduction,
            worst_risk_reduction,
            robust_utility,
        ) = _scenario_eigs(
            prior,
            option,
            active_config,
            decision_array,
            prior_decision_risk,
        )
        cost_adjusted = (
            robust_utility / float(option.cost) ** active_config.cost_exponent
        )
        ranked.append(
            {
                "option_id": option.option_id,
                "measurement_type": option.measurement_type,
                "target_id": option.target_id,
                "cost": float(option.cost),
                "expected_information_gain": mean_eig,
                "worst_case_information_gain": worst_eig,
                "robust_information_gain": robust_eig,
                "expected_brier_risk_reduction": mean_risk_reduction,
                "worst_case_brier_risk_reduction": worst_risk_reduction,
                "robust_combined_utility": robust_utility,
                "cost_adjusted_score": float(cost_adjusted),
                "scenario_scores": scenario_values,
                "metadata": dict(option.metadata),
            }
        )

    ranked.sort(key=lambda row: (-row["cost_adjusted_score"], row["option_id"]))
    for index, row in enumerate(ranked, start=1):
        row["rank"] = index

    best_eig = ranked[0]["robust_information_gain"] if ranked else 0.0
    actionable = bool(
        ranked and best_eig >= active_config.minimum_expected_information_gain
    )
    return {
        "status": "ACTIONABLE" if actionable else "ABSTAIN",
        "selected_option_id": ranked[0]["option_id"] if actionable else None,
        "prior_entropy": shannon_entropy(prior),
        "prior_decision_risk": prior_decision_risk,
        "rankings": ranked,
        "excluded_options": excluded,
        "config": {
            "quadrature_order": active_config.quadrature_order,
            "robustness_weight": active_config.robustness_weight,
            "decision_weight": active_config.decision_weight,
            "cost_exponent": active_config.cost_exponent,
            "minimum_expected_information_gain": (
                active_config.minimum_expected_information_gain
            ),
        },
        "claim_guardrail": (
            "A recommendation is conditional on the supplied hypothesis ensemble, "
            "predictive scenarios, observation-error model, costs, and feasibility."
        ),
    }


def update_hypothesis_posterior(
    hypothesis_ids: Sequence[str],
    prior_probabilities: Sequence[float],
    option: MeasurementOption,
    observed_value: float,
) -> dict[str, Any]:
    """Update hypothesis weights using scenario-averaged Gaussian likelihoods."""

    ids = [str(value) for value in hypothesis_ids]
    prior = _normalise_probabilities(
        prior_probabilities, name="prior_probabilities"
    )
    if len(ids) != len(prior) or len(ids) != len(set(ids)):
        raise ValueError("Hypothesis IDs and probabilities must align and be unique.")
    _validate_option(option, len(prior))
    if not math.isfinite(float(observed_value)):
        raise ValueError("observed_value must be finite.")

    scenario_weights = _normalise_probabilities(
        [scenario.weight for scenario in option.scenarios],
        name="scenario_weights",
        minimum_size=1,
    )
    scenario_log_likelihoods = []
    for scenario in option.scenarios:
        means, standard_deviations = _scenario_arrays(scenario, len(prior))
        residual = (float(observed_value) - means) / standard_deviations
        scenario_log_likelihoods.append(
            -0.5 * residual**2
            - np.log(standard_deviations)
            - 0.5 * math.log(2.0 * math.pi)
        )
    stacked = np.vstack(scenario_log_likelihoods)
    log_likelihood = logsumexp(
        stacked + np.log(scenario_weights)[:, None], axis=0
    )
    log_prior = np.log(np.clip(prior, _PROBABILITY_FLOOR, None))
    log_joint = log_prior + log_likelihood
    log_evidence = float(logsumexp(log_joint))
    posterior = np.exp(log_joint - log_evidence)
    positive = posterior > 0.0
    realised_kl = float(
        np.sum(
            posterior[positive]
            * (np.log(posterior[positive]) - log_prior[positive])
        )
    )
    return {
        "hypothesis_ids": ids,
        "posterior_probabilities": posterior.tolist(),
        "prior_entropy": shannon_entropy(prior),
        "posterior_entropy": shannon_entropy(posterior),
        "realised_information_gain": realised_kl,
        "log_predictive_density": log_evidence,
        "observed_value": float(observed_value),
        "option_id": option.option_id,
    }


def _stable_seed(base_seed: int, labels: Sequence[str]) -> int:
    payload = "\x1f".join(map(str, labels)).encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return (int(base_seed) + int.from_bytes(digest[:4], "big")) % (2**32)


def _validate_batch_scenarios(options: Sequence[MeasurementOption]) -> list[str]:
    if not options:
        return []
    names = [scenario.name for scenario in options[0].scenarios]
    for option in options[1:]:
        if [scenario.name for scenario in option.scenarios] != names:
            raise ValueError(
                "All batch options must expose the same ordered scenario names."
            )
    return names


def _joint_eig_qmc(
    prior: np.ndarray,
    options: Sequence[MeasurementOption],
    scenario_index: int,
    *,
    samples: int,
    seed: int,
) -> float:
    if not options:
        return 0.0
    exponent = int(math.ceil(math.log2(max(64, int(samples)))))
    n_samples = 2**exponent
    sampler = qmc.Sobol(d=len(options) + 1, scramble=True, seed=int(seed))
    unit = np.clip(sampler.random_base2(exponent), 1.0e-12, 1.0 - 1.0e-12)
    cumulative = np.cumsum(prior)
    true_indices = np.searchsorted(cumulative, unit[:, 0], side="right")
    true_indices = np.minimum(true_indices, len(prior) - 1)
    normal_draws = ndtri(unit[:, 1:])

    log_posterior = np.broadcast_to(
        np.log(np.clip(prior, _PROBABILITY_FLOOR, None)),
        (n_samples, len(prior)),
    ).copy()
    for option_index, option in enumerate(options):
        scenario = option.scenarios[scenario_index]
        means, standard_deviations = _scenario_arrays(scenario, len(prior))
        observations = (
            means[true_indices]
            + standard_deviations[true_indices] * normal_draws[:, option_index]
        )
        residuals = (observations[:, None] - means[None, :]) / standard_deviations[
            None, :
        ]
        log_posterior += (
            -0.5 * residuals**2
            - np.log(standard_deviations[None, :])
            - 0.5 * math.log(2.0 * math.pi)
        )
    log_posterior -= logsumexp(log_posterior, axis=1, keepdims=True)
    posterior = np.exp(log_posterior)
    posterior_entropy = -np.sum(
        posterior * np.where(posterior > 0.0, log_posterior, 0.0), axis=1
    )
    return float(max(0.0, shannon_entropy(prior) - np.mean(posterior_entropy)))


def _robust_joint_eig(
    prior: np.ndarray,
    options: Sequence[MeasurementOption],
    config: AcquisitionConfig,
) -> tuple[float, float, float, list[dict[str, float | str]]]:
    names = _validate_batch_scenarios(options)
    if not names:
        return 0.0, 0.0, 0.0, []
    scenario_weights = _normalise_probabilities(
        [scenario.weight for scenario in options[0].scenarios],
        name="scenario_weights",
        minimum_size=1,
    )
    eigs = []
    rows: list[dict[str, float | str]] = []
    labels = sorted(option.option_id for option in options)
    for scenario_index, name in enumerate(names):
        eig = _joint_eig_qmc(
            prior,
            options,
            scenario_index,
            samples=config.batch_qmc_samples,
            seed=_stable_seed(config.random_seed + scenario_index, labels),
        )
        eigs.append(eig)
        rows.append({"scenario": name, "joint_expected_information_gain": eig})
    mean_eig = float(np.dot(scenario_weights, np.asarray(eigs)))
    worst_eig = float(min(eigs))
    robust_eig = float(
        (1.0 - config.robustness_weight) * mean_eig
        + config.robustness_weight * worst_eig
    )
    return mean_eig, worst_eig, robust_eig, rows


def select_measurement_batch(
    hypothesis_ids: Sequence[str],
    prior_probabilities: Sequence[float],
    options: Sequence[MeasurementOption],
    *,
    batch_size: int,
    budget: float | None = None,
    config: AcquisitionConfig | None = None,
) -> dict[str, Any]:
    """Greedily maximise robust joint EIG for a non-redundant batch.

    Joint information is estimated with a deterministic scrambled Sobol rule.
    Each addition is scored by its marginal joint information gain, which
    prevents duplicate or near-duplicate measurements from being treated as
    independent value.
    """

    active_config = config or AcquisitionConfig()
    _validate_config(active_config)
    ids = [str(value) for value in hypothesis_ids]
    prior = _normalise_probabilities(
        prior_probabilities, name="prior_probabilities"
    )
    if len(ids) != len(prior) or len(ids) != len(set(ids)):
        raise ValueError("Hypothesis IDs and probabilities must align and be unique.")
    if int(batch_size) < 1:
        raise ValueError("batch_size must be positive.")
    if budget is not None and (not math.isfinite(float(budget)) or budget <= 0.0):
        raise ValueError("budget must be positive when provided.")

    available = []
    seen = set()
    for option in options:
        _validate_option(option, len(prior))
        if option.option_id in seen:
            raise ValueError(f"Duplicate option_id: {option.option_id}")
        seen.add(option.option_id)
        if option.feasible:
            available.append(option)
    _validate_batch_scenarios(available)

    selected: list[MeasurementOption] = []
    selection_rows: list[dict[str, Any]] = []
    total_cost = 0.0
    current_joint = 0.0
    for _ in range(min(int(batch_size), len(available))):
        candidates = []
        for option in available:
            if option in selected:
                continue
            proposed_cost = total_cost + float(option.cost)
            if budget is not None and proposed_cost > float(budget) + 1.0e-12:
                continue
            mean_eig, worst_eig, robust_eig, scenarios = _robust_joint_eig(
                prior, [*selected, option], active_config
            )
            marginal = max(0.0, robust_eig - current_joint)
            score = marginal / float(option.cost) ** active_config.cost_exponent
            candidates.append(
                (
                    -score,
                    option.option_id,
                    option,
                    mean_eig,
                    worst_eig,
                    robust_eig,
                    marginal,
                    scenarios,
                )
            )
        if not candidates:
            break
        (
            _,
            _,
            best,
            mean_eig,
            worst_eig,
            robust_eig,
            marginal,
            scenarios,
        ) = min(candidates, key=lambda value: (value[0], value[1]))
        if marginal < active_config.minimum_expected_information_gain:
            break
        selected.append(best)
        total_cost += float(best.cost)
        current_joint = robust_eig
        selection_rows.append(
            {
                "batch_rank": len(selected),
                "option_id": best.option_id,
                "measurement_type": best.measurement_type,
                "target_id": best.target_id,
                "cost": float(best.cost),
                "cumulative_cost": total_cost,
                "marginal_robust_information_gain": marginal,
                "joint_expected_information_gain": mean_eig,
                "joint_worst_case_information_gain": worst_eig,
                "joint_robust_information_gain": robust_eig,
                "scenario_scores": scenarios,
            }
        )

    return {
        "status": "ACTIONABLE" if selected else "ABSTAIN",
        "selected_option_ids": [option.option_id for option in selected],
        "selections": selection_rows,
        "total_cost": total_cost,
        "joint_robust_information_gain": current_joint,
        "prior_entropy": shannon_entropy(prior),
        "claim_guardrail": (
            "Batch value is conditional on the supplied scenario family and "
            "conditional-independence observation model."
        ),
    }


__all__ = [
    "AcquisitionConfig",
    "MeasurementOption",
    "PredictiveScenario",
    "expected_information_gain",
    "expected_brier_risk_reduction",
    "rank_measurement_options",
    "select_measurement_batch",
    "shannon_entropy",
    "update_hypothesis_posterior",
]
