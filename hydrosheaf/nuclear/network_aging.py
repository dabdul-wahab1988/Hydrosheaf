"""Topology-aware Bayesian groundwater dating.

The model in this module deliberately separates two questions:

* did the numerical sampler converge; and
* do the supplied tracers identify age?

A converged posterior is not automatically an identified posterior.  This is
especially important for post-bomb tritium and for old water near a detection
floor.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
import warnings
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import networkx as nx
import numpy as np

from hydrosheaf.nuclear.input_history import InputHistory, get_input_history
from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.nuclear.nuclides import Nuclide, TRITIUM
from ..log import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class TracerObservationSet:
    """One tracer panel observed on some or all network nodes."""

    name: str
    nuclide: Nuclide
    observations: Mapping[str, float]
    sigmas: Mapping[str, float]
    sample_date: float
    input_history: Optional[InputHistory] = None
    model_type: str = "PFM"
    reference_concentration: float = 100.0


def _load_network_aging_dependencies() -> Tuple[Any, Any, Any]:
    """Load optional Bayesian dependencies lazily."""

    try:
        import nutpie
        import pymc as pm
        import pytensor.tensor as pt
    except Exception as exc:
        raise ImportError(
            "Network aging requires pymc, pytensor, and nutpie with writable "
            "runtime directories."
        ) from exc
    return pm, pt, nutpie


def _linear_interp(x, x_ref, y_ref, pt_module: Any):
    """Differentiable linear interpolation on a regular lookup grid."""

    x_min = x_ref[0]
    x_max = x_ref[-1]
    x_clipped = pt_module.clip(x, x_min, x_max)
    dx = x_ref[1] - x_ref[0]
    idx_float = (x_clipped - x_min) / dx
    idx_floor = pt_module.floor(idx_float).astype("int64")
    idx_ceil = idx_floor + 1
    n_points = x_ref.shape[0]
    idx_floor = pt_module.clip(idx_floor, 0, n_points - 2)
    idx_ceil = pt_module.clip(idx_ceil, 1, n_points - 1)
    y_floor = y_ref[idx_floor]
    y_ceil = y_ref[idx_ceil]
    w_ceil = idx_float - idx_floor
    return y_floor * (1.0 - w_ceil) + y_ceil * w_ceil


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^0-9A-Za-z_]+", "_", value).strip("_")
    return safe or "tracer"


def _forward_grid(
    tracer: TracerObservationSet,
    age_grid: np.ndarray,
) -> np.ndarray:
    """Build a tracer response grid without embedding it in PyTensor."""

    decay = np.log(2.0) / tracer.nuclide.half_life_years
    history = tracer.input_history
    if history is None and tracer.nuclide.symbol == "3H":
        history = get_input_history("global")
    if history is None:
        return tracer.reference_concentration * np.exp(-decay * age_grid)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Recharge date .* exceeds input history start",
        )
        return np.asarray(
            [
                convolve_input(
                    tracer.sample_date,
                    float(age),
                    history.years,
                    history.values,
                    decay,
                    model_type=tracer.model_type,
                )
                for age in age_grid
            ],
            dtype=float,
        )


def _parent_mixing_weights(
    graph: nx.DiGraph,
    node: str,
    parents: Sequence[str],
) -> np.ndarray:
    """Resolve fixed, observable parent-flow weights.

    User-provided ``flow_fraction`` values are normalised.  When they are not
    available, equal mixing is an explicit assumption rather than an additional
    weakly identified latent variable.
    """

    raw = np.asarray(
        [
            max(
                0.0,
                float(
                    graph.get_edge_data(parent, node, default={}).get(
                        "flow_fraction", 1.0
                    )
                ),
            )
            for parent in parents
        ],
        dtype=float,
    )
    if not np.isfinite(raw).all() or raw.sum() <= 0.0:
        raw = np.ones(len(parents), dtype=float)
    return raw / raw.sum()


def _count_likelihood_modes(
    probabilities: np.ndarray,
    age_grid: np.ndarray,
) -> int:
    """Count material, separated modes in a normalised age likelihood."""

    if probabilities.size < 3 or not np.isfinite(probabilities).all():
        return 0
    threshold = float(np.max(probabilities)) * 0.05
    candidates = [
        idx
        for idx in range(1, probabilities.size - 1)
        if probabilities[idx] >= threshold
        and probabilities[idx] >= probabilities[idx - 1]
        and probabilities[idx] >= probabilities[idx + 1]
    ]
    if probabilities[0] >= threshold and probabilities[0] >= probabilities[1]:
        candidates.append(0)
    if probabilities[-1] >= threshold and probabilities[-1] >= probabilities[-2]:
        candidates.append(probabilities.size - 1)
    selected: list[int] = []
    for idx in sorted(candidates, key=lambda item: probabilities[item], reverse=True):
        if all(abs(age_grid[idx] - age_grid[other]) >= 5.0 for other in selected):
            selected.append(idx)
    return len(selected)


def _local_likelihood_diagnostics(
    node: str,
    tracers: Sequence[TracerObservationSet],
    grids: Sequence[np.ndarray],
    age_grid: np.ndarray,
) -> Dict[str, object]:
    """Diagnose tracer-only ambiguity independently of the network prior."""

    log_likelihood = np.zeros_like(age_grid, dtype=float)
    n_observed = 0
    for tracer, concentration_grid in zip(tracers, grids):
        if node not in tracer.observations:
            continue
        value = max(0.0, float(tracer.observations[node]))
        sigma = max(
            1e-6,
            float(tracer.sigmas.get(node, max(0.5, 0.1 * value))),
        )
        log_likelihood -= 0.5 * ((value - concentration_grid) / sigma) ** 2
        n_observed += 1
    if n_observed == 0:
        return {
            "n_tracers_observed": 0,
            "likelihood_mode_count": 0,
            "likelihood_age_95_width_years": None,
            "tracer_identifiable": False,
        }

    relative = np.exp(log_likelihood - float(np.max(log_likelihood)))
    total = float(relative.sum())
    if total <= 0.0 or not math.isfinite(total):
        return {
            "n_tracers_observed": n_observed,
            "likelihood_mode_count": 0,
            "likelihood_age_95_width_years": None,
            "tracer_identifiable": False,
        }
    probabilities = relative / total
    cdf = np.cumsum(probabilities)
    low = float(age_grid[int(np.searchsorted(cdf, 0.025))])
    high = float(age_grid[min(len(age_grid) - 1, int(np.searchsorted(cdf, 0.975)))])
    modes = _count_likelihood_modes(probabilities, age_grid)
    width = high - low
    # A 100-year support window is intentionally permissive.  Results wider
    # than this are not useful as groundwater "dates" even if MCMC converges.
    identifiable = bool(modes <= 1 and width <= 100.0)
    return {
        "n_tracers_observed": n_observed,
        "likelihood_mode_count": int(modes),
        "likelihood_age_95_width_years": float(width),
        "tracer_identifiable": identifiable,
    }


def _student_t_logpdf(
    residual: np.ndarray,
    sigma: float,
    degrees_of_freedom: float,
) -> np.ndarray:
    nu = float(degrees_of_freedom)
    scale = max(1e-12, float(sigma))
    constant = (
        math.lgamma((nu + 1.0) / 2.0)
        - math.lgamma(nu / 2.0)
        - 0.5 * math.log(nu * math.pi)
        - math.log(scale)
    )
    return constant - 0.5 * (nu + 1.0) * np.log1p(
        (np.asarray(residual, dtype=float) / scale) ** 2 / nu
    )


def _infer_independent_grid(
    *,
    nodes: Sequence[str],
    tracers: Sequence[TracerObservationSet],
    concentration_grids: Sequence[np.ndarray],
    age_grid: np.ndarray,
    n_samples: int,
    n_chains: int,
    random_seed: int,
    root_age_prior_median_years: float,
    root_age_prior_log_sigma: float,
    input_scale: float,
    degrees_of_freedom: float,
) -> Dict[str, Dict[str, float]]:
    """Exact quadrature for a factorised, edge-free age posterior."""

    if any(tracer.nuclide.symbol == "14C" for tracer in tracers):
        raise ValueError(
            "The independent grid sampler does not marginalise the 14C A0 "
            "correction; use sampler='smc' for radiocarbon panels."
        )
    positive_age = np.maximum(age_grid, 1e-9)
    log_sigma = max(float(root_age_prior_log_sigma), 1e-6)
    log_prior = (
        -np.log(positive_age * log_sigma * math.sqrt(2.0 * math.pi))
        - 0.5
        * (
            (
                np.log(positive_age)
                - math.log(max(float(root_age_prior_median_years), 1e-6))
            )
            / log_sigma
        )
        ** 2
    )
    log_prior[age_grid <= 0.0] = -np.inf
    probabilities_by_node: Dict[str, np.ndarray] = {}
    local_diagnostics: Dict[str, Dict[str, object]] = {}
    for node in nodes:
        log_posterior = log_prior.copy()
        for tracer_index, (tracer, grid) in enumerate(
            zip(tracers, concentration_grids)
        ):
            if node not in tracer.observations:
                continue
            expected = grid * (input_scale if tracer_index == 0 else 1.0)
            value = float(tracer.observations[node])
            sigma = max(
                1e-6,
                float(tracer.sigmas.get(node, max(0.5, 0.1 * max(value, 0.0)))),
            )
            log_posterior += _student_t_logpdf(
                value - expected,
                sigma,
                degrees_of_freedom,
            )
        relative = np.exp(log_posterior - np.nanmax(log_posterior))
        probabilities = relative / relative.sum()
        probabilities_by_node[node] = probabilities
        local_diagnostics[node] = _local_likelihood_diagnostics(
            node, tracers, concentration_grids, age_grid
        )

    draws = np.empty((int(n_chains), int(n_samples), len(nodes)), dtype=float)
    for chain in range(int(n_chains)):
        rng = np.random.default_rng(int(random_seed) + 104729 * chain)
        for node_index, node in enumerate(nodes):
            draws[chain, :, node_index] = rng.choice(
                age_grid,
                size=int(n_samples),
                replace=True,
                p=probabilities_by_node[node],
            )

    import arviz as az

    idata = az.from_dict(
        {"posterior": {"ages": draws}},
        coords={"node": list(nodes)},
        dims={"ages": ["node"]},
    )
    r_hat_values = np.asarray(
        az.rhat(idata, var_names=["ages"], method="rank")["ages"],
        dtype=float,
    )
    ess_bulk_values = np.asarray(
        az.ess(idata, var_names=["ages"], method="bulk")["ages"],
        dtype=float,
    )
    ess_tail_values = np.asarray(
        az.ess(idata, var_names=["ages"], method="tail")["ages"],
        dtype=float,
    )
    mcse_values = np.asarray(
        az.mcse(idata, var_names=["ages"], method="mean")["ages"],
        dtype=float,
    )
    flat_ages = draws.reshape(-1, len(nodes))
    posterior_predictive: Dict[str, Dict[str, float]] = {}
    from scipy.stats import t as student_t

    multiplier = float(student_t.ppf(0.975, df=degrees_of_freedom))
    for tracer_index, (tracer, grid) in enumerate(zip(tracers, concentration_grids)):
        observed_nodes = [node for node in nodes if node in tracer.observations]
        if not observed_nodes:
            continue
        predicted = np.column_stack(
            [
                np.interp(
                    flat_ages[:, nodes.index(node)],
                    age_grid,
                    grid,
                )
                * (input_scale if tracer_index == 0 else 1.0)
                for node in observed_nodes
            ]
        )
        values = np.asarray(
            [tracer.observations[node] for node in observed_nodes], dtype=float
        )
        errors = np.asarray(
            [
                max(
                    1e-6,
                    tracer.sigmas.get(
                        node, max(0.5, 0.1 * max(tracer.observations[node], 0.0))
                    ),
                )
                for node in observed_nodes
            ],
            dtype=float,
        )
        posterior_predictive[tracer.name] = {
            "standardized_rmse": float(
                np.sqrt(np.mean(((values - predicted.mean(axis=0)) / errors) ** 2))
            ),
            "observed_95_interval_coverage": float(
                np.mean(
                    (
                        values
                        >= np.quantile(predicted, 0.025, axis=0) - multiplier * errors
                    )
                    & (
                        values
                        <= np.quantile(predicted, 0.975, axis=0) + multiplier * errors
                    )
                )
            ),
            "n_observations": int(len(values)),
        }

    diagnostic: Dict[str, object] = {
        "inferred_input_scale": float(input_scale),
        "input_scale_mode": "fixed",
        "inferred_velocity_m_y": None,
        "velocity_mode": "not_applicable_edge_free",
        "model_type": tracers[0].model_type,
        "n_tracer_panels": len(tracers),
        "n_tune": 0,
        "target_accept": None,
        "sampler": "exact_grid",
        "likelihood": "student_t",
        "likelihood_degrees_of_freedom": float(degrees_of_freedom),
        "age_parameterization": "independent_lognormal_grid_quadrature",
        "age_r_hat_max": float(np.nanmax(r_hat_values)),
        "age_ess_bulk_min": float(np.nanmin(ess_bulk_values)),
        "age_ess_tail_min": float(np.nanmin(ess_tail_values)),
        "age_mcse_mean_max": float(np.nanmax(mcse_values)),
        "diagnostic_method": "rank_normalized_split_rhat_bulk_tail_ess",
        "divergences": 0,
        "posterior_predictive": posterior_predictive,
        "n_nodes_tracer_identifiable": int(
            sum(
                bool(value["tracer_identifiable"])
                for value in local_diagnostics.values()
            )
        ),
        "all_nodes_tracer_identifiable": bool(
            all(
                bool(value["tracer_identifiable"])
                for value in local_diagnostics.values()
            )
        ),
        "edge_age_order": {},
    }
    diagnostic["converged"] = bool(
        float(diagnostic["age_r_hat_max"]) <= 1.01
        and float(diagnostic["age_ess_bulk_min"]) >= 400.0
        and float(diagnostic["age_ess_tail_min"]) >= 400.0
    )
    results: Dict[str, Dict[str, float]] = {
        "_diagnostics": diagnostic  # type: ignore[dict-item]
    }
    for node_index, node in enumerate(nodes):
        node_ages = flat_ages[:, node_index]
        low = float(np.percentile(node_ages, 2.5))
        high = float(np.percentile(node_ages, 97.5))
        local = local_diagnostics[node]
        results[node] = {
            "mean_age_years": float(np.mean(node_ages)),
            "median_age_years": float(np.median(node_ages)),
            "std_age_years": float(np.std(node_ages)),
            "p_modern": float(np.mean(node_ages < 60.0)),
            "age_95_low": low,
            "age_95_high": high,
            "posterior_95_width_years": high - low,
            "n_tracers_observed": float(local["n_tracers_observed"]),
            "likelihood_mode_count": float(local["likelihood_mode_count"]),
            "likelihood_age_95_width_years": (
                float(local["likelihood_age_95_width_years"])
                if local["likelihood_age_95_width_years"] is not None
                else float("nan")
            ),
            "tracer_identifiable": float(bool(local["tracer_identifiable"])),
        }
    return results


def infer_network_ages_bayesian(
    graph: nx.DiGraph,
    node_observations: Dict[str, float],
    node_sigmas: Dict[str, float],
    sample_date: float,
    nuclide: Nuclide = TRITIUM,
    input_hist: Optional[InputHistory] = None,
    n_samples: int = 1000,
    n_chains: int = 4,
    random_seed: int = 42,
    model_type: str = "PFM",
    velocity_prior_mu: float = 10.0,
    velocity_prior_sigma: float = 2.5,
    input_scale_prior_mu: float = 1.0,
    input_scale_prior_sigma: float = 0.25,
    c14_initial_prior_mu: float = 85.0,
    c14_initial_prior_sigma: float = 10.0,
    *,
    additional_tracers: Optional[Sequence[TracerObservationSet]] = None,
    infer_velocity: bool = False,
    infer_input_scale: bool = False,
    n_tune: Optional[int] = None,
    target_accept: float = 0.95,
    root_age_prior_median_years: float = 10.0,
    root_age_prior_log_sigma: float = 1.0,
    age_increment_prior_scale_years: float = 15.0,
    max_age_years: float = 500.0,
    sampler: str = "smc",
    likelihood_degrees_of_freedom: float = 4.0,
) -> Dict[str, Dict[str, float]]:
    """Infer node-age posteriors under a directed acyclic mixing network.

    Ages are generated in topological order.  Each downstream mean age equals
    the flow-weighted mean of routed parent ages plus a positive local
    residence-time increment.  This removes the discontinuous hinge potential
    and the redundant independent-age parameterisation used by the earlier
    implementation.

    Velocity and atmospheric input scaling are fixed by default because tracer
    concentrations alone generally do not identify them separately from age.
    Set ``infer_velocity`` or ``infer_input_scale`` only when the supplied priors
    represent independent information.
    """

    if not isinstance(graph, nx.DiGraph):
        raise TypeError("graph must be a networkx.DiGraph.")
    if graph.number_of_nodes() == 0:
        raise ValueError("Network aging requires at least one node.")
    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError(
            "Topology-generative Bayesian aging requires a directed acyclic "
            "graph; resolve cycles or sample topology before dating."
        )
    if n_samples < 2 or n_chains < 1:
        raise ValueError("n_samples must be at least 2 and n_chains at least 1.")
    sampler_normalized = str(sampler).strip().lower()
    if sampler_normalized not in {"grid", "smc", "nuts"}:
        raise ValueError("sampler must be 'grid', 'smc', or 'nuts'.")
    if float(likelihood_degrees_of_freedom) <= 2.0:
        raise ValueError(
            "likelihood_degrees_of_freedom must exceed 2 for finite variance."
        )
    if sampler_normalized == "nuts" and not 0.5 < float(target_accept) < 1.0:
        raise ValueError("target_accept must lie strictly between 0.5 and 1.0.")

    nodes = [str(node) for node in graph.nodes()]
    if len(set(nodes)) != graph.number_of_nodes():
        raise ValueError("String-normalised node identifiers must be unique.")
    graph_string = nx.relabel_nodes(graph, {node: str(node) for node in graph.nodes()})
    node_idx = {node: index for index, node in enumerate(nodes)}
    topological_nodes = list(nx.topological_sort(graph_string))
    roots = [node for node in topological_nodes if graph_string.in_degree(node) == 0]
    non_roots = [node for node in topological_nodes if node not in roots]

    primary = TracerObservationSet(
        name=nuclide.symbol,
        nuclide=nuclide,
        observations=dict(node_observations),
        sigmas=dict(node_sigmas),
        sample_date=float(sample_date),
        input_history=input_hist,
        model_type=model_type,
    )
    tracers = [primary, *(additional_tracers or ())]
    if not any(
        str(node) in {str(key) for key in tracer.observations}
        for tracer in tracers
        for node in nodes
    ):
        raise ValueError("At least one tracer observation must match a graph node.")
    # Normalise external mapping keys once.
    tracers = [
        TracerObservationSet(
            name=tracer.name,
            nuclide=tracer.nuclide,
            observations={
                str(key): float(value)
                for key, value in tracer.observations.items()
                if value is not None and np.isfinite(float(value))
            },
            sigmas={
                str(key): float(value)
                for key, value in tracer.sigmas.items()
                if value is not None and np.isfinite(float(value))
            },
            sample_date=float(tracer.sample_date),
            input_history=tracer.input_history,
            model_type=tracer.model_type,
            reference_concentration=float(tracer.reference_concentration),
        )
        for tracer in tracers
    ]

    max_age = max(10.0, float(max_age_years))
    age_grid = np.linspace(0.0, max_age, max(501, int(max_age) + 1))
    concentration_grids = [_forward_grid(tracer, age_grid) for tracer in tracers]
    if sampler_normalized == "grid":
        if graph_string.number_of_edges() != 0:
            raise ValueError(
                "sampler='grid' is exact only for edge-free independent-node "
                "models; use sampler='smc' for a network."
            )
        if infer_velocity or infer_input_scale:
            raise ValueError(
                "sampler='grid' requires fixed velocity and input scaling."
            )
        return _infer_independent_grid(
            nodes=nodes,
            tracers=tracers,
            concentration_grids=concentration_grids,
            age_grid=age_grid,
            n_samples=n_samples,
            n_chains=n_chains,
            random_seed=random_seed,
            root_age_prior_median_years=root_age_prior_median_years,
            root_age_prior_log_sigma=root_age_prior_log_sigma,
            input_scale=input_scale_prior_mu,
            degrees_of_freedom=likelihood_degrees_of_freedom,
        )

    pm, pt, nutpie = _load_network_aging_dependencies()
    velocity_mu = max(float(velocity_prior_mu), 1e-6)
    velocity_sigma = max(float(velocity_prior_sigma), 1e-6)
    tune = (
        0
        if sampler_normalized == "smc"
        else (int(n_tune) if n_tune is not None else max(500, int(n_samples)))
    )

    logger.info(
        "Starting topology-generative Bayesian dating: %d nodes, %d edges, "
        "%d tracer panels, draws=%d, tune=%d.",
        len(nodes),
        graph_string.number_of_edges(),
        len(tracers),
        n_samples,
        tune,
    )

    with pm.Model() as model:
        if infer_input_scale:
            input_scale = pm.LogNormal(
                "input_scale",
                mu=np.log(max(float(input_scale_prior_mu), 1e-6)),
                sigma=max(float(input_scale_prior_sigma), 1e-6),
            )
        else:
            input_scale = pt.as_tensor_variable(float(input_scale_prior_mu))

        if infer_velocity:
            velocity_sigma_log = float(np.log1p(velocity_sigma / velocity_mu))
            velocity = pm.LogNormal(
                "velocity",
                mu=np.log(velocity_mu),
                sigma=velocity_sigma_log,
            )
        else:
            velocity = pt.as_tensor_variable(velocity_mu)

        root_ages = pm.LogNormal(
            "root_ages",
            mu=np.log(max(float(root_age_prior_median_years), 1e-3)),
            sigma=max(float(root_age_prior_log_sigma), 1e-3),
            shape=len(roots),
        )
        if non_roots:
            age_increments = pm.HalfNormal(
                "age_increments",
                sigma=max(float(age_increment_prior_scale_years), 1e-3),
                shape=len(non_roots),
            )
        else:
            age_increments = None

        age_expressions: Dict[str, Any] = {}
        for root_index, node in enumerate(roots):
            age_expressions[node] = root_ages[root_index]
        for increment_index, node in enumerate(non_roots):
            parents = list(graph_string.predecessors(node))
            weights = _parent_mixing_weights(graph_string, node, parents)
            routed = []
            for parent, weight in zip(parents, weights):
                edge_data = graph_string.get_edge_data(parent, node, default={})
                length_m = max(0.0, float(edge_data.get("length_m", 1.0)))
                routed.append(
                    float(weight) * (age_expressions[parent] + length_m / velocity)
                )
            parent_mean = sum(routed)
            age_expressions[node] = parent_mean + age_increments[increment_index]
        ages = pm.Deterministic(
            "ages", pt.stack([age_expressions[node] for node in nodes])
        )

        has_c14 = any(tracer.nuclide.symbol == "14C" for tracer in tracers)
        if has_c14:
            c14_a0_mu = pm.TruncatedNormal(
                "c14_a0_mu",
                mu=float(c14_initial_prior_mu),
                sigma=max(float(c14_initial_prior_sigma), 1e-3),
                lower=1.0,
                upper=120.0,
            )
            c14_a0_local = pm.TruncatedNormal(
                "c14_a0_local",
                mu=c14_a0_mu,
                sigma=5.0,
                lower=1.0,
                upper=120.0,
                shape=len(nodes),
            )
        else:
            c14_a0_local = None

        grid_x_pt = pt.as_tensor_variable(age_grid)
        for tracer_index, (tracer, concentration_grid) in enumerate(
            zip(tracers, concentration_grids)
        ):
            observed_nodes = [node for node in nodes if node in tracer.observations]
            if not observed_nodes:
                continue
            indices = np.asarray([node_idx[node] for node in observed_nodes], dtype=int)
            values = np.asarray(
                [max(0.0, tracer.observations[node]) for node in observed_nodes],
                dtype=float,
            )
            errors = np.asarray(
                [
                    max(
                        1e-6,
                        tracer.sigmas.get(
                            node, max(0.5, 0.1 * tracer.observations[node])
                        ),
                    )
                    for node in observed_nodes
                ],
                dtype=float,
            )
            expected = _linear_interp(
                ages[indices],
                grid_x_pt,
                pt.as_tensor_variable(concentration_grid),
                pt,
            )
            if tracer.nuclide.symbol == "14C":
                expected = expected * (c14_a0_local[indices] / 100.0)
            if tracer_index == 0 and infer_input_scale:
                expected = expected * input_scale
            variable_name = f"{tracer_index}_{_safe_name(tracer.name)}"
            pm.Deterministic(f"expected_{variable_name}", expected)
            pm.StudentT(
                f"obs_likelihood_{variable_name}",
                nu=float(likelihood_degrees_of_freedom),
                mu=expected,
                sigma=errors,
                observed=values,
            )

        if sampler_normalized == "smc":
            trace = pm.sample_smc(
                draws=int(n_samples),
                chains=int(n_chains),
                cores=1,
                random_seed=[
                    int(random_seed) + 104729 * chain for chain in range(int(n_chains))
                ],
                compute_convergence_checks=False,
                progressbar=False,
                return_inferencedata=True,
            )
        else:
            compiled = nutpie.compile_pymc_model(
                model,
                default_initialization_strategy="prior",
            )
            trace = nutpie.sample(
                compiled,
                draws=int(n_samples),
                tune=tune,
                chains=int(n_chains),
                seed=int(random_seed),
                target_accept=float(target_accept),
                save_warmup=False,
                progress_bar=False,
            )

    results: Dict[str, Dict[str, float]] = {}
    flat_ages = np.asarray(trace.posterior["ages"]).reshape(-1, len(nodes))
    diag: Dict[str, object] = {
        "inferred_input_scale": (
            float(np.asarray(trace.posterior["input_scale"]).mean())
            if "input_scale" in trace.posterior
            else float(input_scale_prior_mu)
        ),
        "input_scale_mode": "inferred" if infer_input_scale else "fixed",
        "inferred_velocity_m_y": (
            float(np.asarray(trace.posterior["velocity"]).mean())
            if "velocity" in trace.posterior
            else velocity_mu
        ),
        "velocity_mode": "inferred" if infer_velocity else "fixed",
        "model_type": model_type,
        "n_tracer_panels": len(tracers),
        "n_tune": tune,
        "target_accept": float(target_accept),
        "sampler": sampler_normalized,
        "likelihood": "student_t",
        "likelihood_degrees_of_freedom": float(likelihood_degrees_of_freedom),
        "age_parameterization": "dag_parent_mixture_plus_positive_increment",
    }

    try:
        import arviz as az

        r_hat_values = np.asarray(
            az.rhat(trace, var_names=["ages"], method="rank")["ages"],
            dtype=float,
        )
        ess_bulk_values = np.asarray(
            az.ess(trace, var_names=["ages"], method="bulk")["ages"],
            dtype=float,
        )
        ess_tail_values = np.asarray(
            az.ess(trace, var_names=["ages"], method="tail")["ages"],
            dtype=float,
        )
        mcse_mean_values = np.asarray(
            az.mcse(trace, var_names=["ages"], method="mean")["ages"],
            dtype=float,
        )
        diag.update(
            {
                "age_r_hat_max": float(np.nanmax(r_hat_values)),
                "age_ess_bulk_min": float(np.nanmin(ess_bulk_values)),
                "age_ess_tail_min": float(np.nanmin(ess_tail_values)),
                "age_mcse_mean_max": float(np.nanmax(mcse_mean_values)),
                "diagnostic_method": ("rank_normalized_split_rhat_bulk_tail_ess"),
            }
        )
    except Exception as exc:
        logger.warning("Could not compute Bayesian-age diagnostics: %s", exc)
        diag.update(
            {
                "age_r_hat_max": None,
                "age_ess_bulk_min": None,
                "age_ess_tail_min": None,
                "age_mcse_mean_max": None,
                "diagnostic_method": "unavailable",
            }
        )

    sample_stats = getattr(trace, "sample_stats", None)
    divergences = None
    if sample_stats is not None:
        for key in ("diverging", "divergence"):
            if key in sample_stats:
                divergences = int(np.asarray(sample_stats[key]).sum())
                break
    if sampler_normalized == "smc" and divergences is None:
        divergences = 0
    diag["divergences"] = divergences
    diag["converged"] = bool(
        diag.get("age_r_hat_max") is not None
        and float(diag["age_r_hat_max"]) <= 1.01
        and float(diag.get("age_ess_bulk_min") or 0.0) >= 400.0
        and float(diag.get("age_ess_tail_min") or 0.0) >= 400.0
        and int(divergences or 0) == 0
    )
    if "c14_a0_mu" in trace.posterior:
        diag["inferred_c14_a0_global"] = float(
            np.asarray(trace.posterior["c14_a0_mu"]).mean()
        )

    posterior_predictive: Dict[str, Dict[str, float]] = {}
    for tracer_index, tracer in enumerate(tracers):
        variable_name = f"expected_{tracer_index}_{_safe_name(tracer.name)}"
        if variable_name not in trace.posterior:
            continue
        expected_draws = np.asarray(trace.posterior[variable_name]).reshape(
            -1,
            np.asarray(trace.posterior[variable_name]).shape[-1],
        )
        observed_nodes = [node for node in nodes if node in tracer.observations]
        values = np.asarray(
            [tracer.observations[node] for node in observed_nodes], dtype=float
        )
        errors = np.asarray(
            [
                max(
                    1e-6,
                    tracer.sigmas.get(node, max(0.5, 0.1 * tracer.observations[node])),
                )
                for node in observed_nodes
            ],
            dtype=float,
        )
        means = expected_draws.mean(axis=0)
        # 97.5% Student-t quantile for the default df=4 is 2.776.  Use SciPy
        # for non-default degrees of freedom when available.
        try:
            from scipy.stats import t as student_t

            predictive_multiplier = float(
                student_t.ppf(0.975, df=float(likelihood_degrees_of_freedom))
            )
        except Exception:
            predictive_multiplier = 2.776
        predictive_low = (
            np.quantile(expected_draws, 0.025, axis=0) - predictive_multiplier * errors
        )
        predictive_high = (
            np.quantile(expected_draws, 0.975, axis=0) + predictive_multiplier * errors
        )
        posterior_predictive[tracer.name] = {
            "standardized_rmse": float(
                np.sqrt(np.mean(((values - means) / errors) ** 2))
            ),
            "observed_95_interval_coverage": float(
                np.mean((values >= predictive_low) & (values <= predictive_high))
            ),
            "n_observations": int(len(values)),
        }
    diag["posterior_predictive"] = posterior_predictive

    local_diagnostics = {
        node: _local_likelihood_diagnostics(
            node, tracers, concentration_grids, age_grid
        )
        for node in nodes
    }
    diag["n_nodes_tracer_identifiable"] = int(
        sum(bool(value["tracer_identifiable"]) for value in local_diagnostics.values())
    )
    diag["all_nodes_tracer_identifiable"] = bool(
        all(bool(value["tracer_identifiable"]) for value in local_diagnostics.values())
    )

    edge_order: Dict[str, Dict[str, float]] = {}
    velocity_estimate = float(diag["inferred_velocity_m_y"])
    for u, v, edge_data in graph_string.edges(data=True):
        length_m = max(0.0, float(edge_data.get("length_m", 1.0)))
        margin = length_m / max(velocity_estimate, 1e-9)
        differences = flat_ages[:, node_idx[v]] - flat_ages[:, node_idx[u]] - margin
        edge_order[f"{u}->{v}"] = {
            "p_downstream_older": float(np.mean(differences > 0.0)),
            "mean_age_margin_years": float(np.mean(differences)),
        }
    diag["edge_age_order"] = edge_order
    logger.info("Bayesian network dating diagnostics: %s", diag)
    results["_diagnostics"] = diag  # type: ignore[assignment]

    for index, node in enumerate(nodes):
        node_ages = flat_ages[:, index]
        low = float(np.percentile(node_ages, 2.5))
        high = float(np.percentile(node_ages, 97.5))
        local = local_diagnostics[node]
        results[node] = {
            "mean_age_years": float(np.mean(node_ages)),
            "median_age_years": float(np.median(node_ages)),
            "std_age_years": float(np.std(node_ages)),
            "p_modern": float(np.mean(node_ages < 60.0)),
            "age_95_low": low,
            "age_95_high": high,
            "posterior_95_width_years": high - low,
            "n_tracers_observed": float(local["n_tracers_observed"]),
            "likelihood_mode_count": float(local["likelihood_mode_count"]),
            "likelihood_age_95_width_years": (
                float(local["likelihood_age_95_width_years"])
                if local["likelihood_age_95_width_years"] is not None
                else float("nan")
            ),
            "tracer_identifiable": float(bool(local["tracer_identifiable"])),
        }

    return results
