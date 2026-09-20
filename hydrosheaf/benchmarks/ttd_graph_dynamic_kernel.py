"""Independent synthetic benchmark generator for dynamic state-dependent graph TTDs.

This module generates controlled-synthetic truth-sealed benchmarks across:
- stationary, slowly drifting, abrupt shift, and seasonal phase kernels;
- missing edges, reversed topology, randomized graphs;
- sparse sampling, high noise, unmodelled local recharge, and wrong forcing;
- tracer dropout and tracer conflict scenarios.

Strict truth-blindness:
Truth arrays reside exclusively in DynamicKernelTruth with a SHA-256 cryptographic
commitment. The inference function receives only DynamicKernelObservations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Any, Mapping, Optional

import numpy as np

from ..models.ttd_losses import wasserstein_ttd_w1
from ..nuclear.dynamic_kernel_inversion import DynamicTTDRecovery
from ..nuclear.graph_tracer_forward import convolve_causal_dynamic_kernel, edge_id_str


DYNAMIC_SCENARIOS = (
    "stationary",
    "slowly_drifting",
    "abrupt_shift",
    "seasonal_phase_varying",
    "wrong_forcing",
    "unmodelled_local_recharge",
    "missing_edge",
    "reversed_graph",
    "random_graph",
    "sparse_sampling",
    "high_observation_noise",
    "tracer_dropout",
    "tracer_conflict",
)


@dataclass(frozen=True)
class DynamicKernelBenchmarkConfig:
    seed: int = 20260917
    scenario: str = "seasonal_phase_varying"
    n_steps: int = 208
    step_days: float = 7.0
    season_period_steps: int = 52
    max_lag_steps: int = 18
    observation_probability: float = 0.85
    observation_noise_std: float = 0.035
    training_fraction: float = 0.70
    young_water_cutoff_steps: int = 4

    def __post_init__(self) -> None:
        if self.scenario not in DYNAMIC_SCENARIOS:
            raise ValueError(f"Unknown scenario {self.scenario!r}, must be one of {DYNAMIC_SCENARIOS}")
        if self.step_days <= 0.0 or self.n_steps < 2 or self.max_lag_steps < 1:
            raise ValueError("step_days, n_steps, and max_lag_steps must be positive")
        if not 0.0 < self.training_fraction < 1.0:
            raise ValueError("training_fraction must lie strictly between zero and one")


@dataclass(frozen=True)
class DynamicKernelTruth:
    """Sealed ground-truth container. Never provided to inference routines."""

    case_id: str
    scenario: str
    time_steps: np.ndarray
    true_forcing: np.ndarray
    node_signals: Mapping[str, np.ndarray]
    local_inputs: Mapping[str, np.ndarray]
    true_edge_kernels: Mapping[str, np.ndarray]  # edge_id -> (T, L)
    local_recharge_fractions: Mapping[str, np.ndarray]
    edge_mixing_fractions: Mapping[str, np.ndarray]
    young_water_fractions: Mapping[str, np.ndarray]
    mean_transit_times: Mapping[str, np.ndarray]
    sealed_commitment: str


@dataclass(frozen=True)
class DynamicKernelObservations:
    """Truth-blind observation bundle supplied to inverse solvers."""

    case_id: str
    scenario: str
    time_steps: np.ndarray
    regional_forcing: np.ndarray
    node_observations: Mapping[str, np.ndarray]  # With NaNs for unobserved/holdout steps
    observation_masks: Mapping[str, np.ndarray]
    candidate_edges: tuple[tuple[str, str], ...]
    training_end_step: int
    sealed_commitment: str
    metadata: Mapping[str, Any]
    local_inputs: Mapping[str, np.ndarray] = field(default_factory=dict)
    observation_commitment: str = ""


def _array_digest(values: np.ndarray) -> str:
    """Stable digest for a numeric array used in the truth commitment."""
    arr = np.asarray(values, dtype=np.float64)
    h = hashlib.sha256()
    h.update(str(arr.shape).encode("ascii"))
    h.update(arr.tobytes(order="C"))
    return h.hexdigest()


def generate_dynamic_kernel_case(
    config: Optional[DynamicKernelBenchmarkConfig] = None,
) -> tuple[DynamicKernelTruth, DynamicKernelObservations]:
    """Generate a reproducible, cryptographically sealed dynamic benchmark case."""
    cfg = config or DynamicKernelBenchmarkConfig()
    rng = np.random.default_rng(cfg.seed)

    T = cfg.n_steps
    times = np.arange(T, dtype=float)
    L = cfg.max_lag_steps + 1
    lags = np.arange(L, dtype=float) * cfg.step_days

    # Topology: Recharge R -> A -> B -> C, and direct bypass A -> C
    true_edges = (("R", "A"), ("A", "B"), ("B", "C"), ("A", "C"))
    # 1. Regional forcing: sinusoidal seasonal signal + stochastic events
    w0 = 2.0 * np.pi / cfg.season_period_steps
    regional_forcing = np.sin(w0 * times) + 0.35 * np.cos(2 * w0 * times) + 0.15 * rng.standard_normal(T)

    # For wrong_forcing scenario, provide corrupted forcing to observations
    if cfg.scenario == "wrong_forcing":
        obs_forcing = np.roll(regional_forcing, 12) + 0.4 * rng.standard_normal(T)
    else:
        obs_forcing = np.array(regional_forcing, copy=True)

    # 2. Build time-varying edge kernels h_uv(tau, t)
    kernels: dict[str, np.ndarray] = {}

    for u, v in true_edges:
        eid = edge_id_str(u, v)
        k_mat = np.zeros((T, L), dtype=float)

        if cfg.scenario == "stationary":
            # Constant mean lag
            c_lag = (4.0 if eid != "A->C" else 8.0) * cfg.step_days
            weights = np.exp(-0.5 * ((lags - c_lag) / 1.5) ** 2)
            weights /= np.sum(weights)
            k_mat[:] = weights
        elif cfg.scenario == "slowly_drifting":
            # Linear drift in mean lag
            c_lags = np.linspace(3.0, 7.0, T) * cfg.step_days
            for t in range(T):
                w = np.exp(-0.5 * ((lags - c_lags[t]) / 1.5) ** 2)
                k_mat[t] = w / np.sum(w)
        elif cfg.scenario == "abrupt_shift":
            # Step change at midpoint
            mid = T // 2
            for t in range(T):
                c = (3.5 if t < mid else 7.5) * cfg.step_days
                w = np.exp(-0.5 * ((lags - c) / 1.5) ** 2)
                k_mat[t] = w / np.sum(w)
        else:
            # Seasonal phase variation: mean lag oscillates with season
            phase = w0 * times + (0.5 if eid == "A->B" else 0.0)
            c_lags = (4.0 + 2.0 * np.sin(phase)) * cfg.step_days
            for t in range(T):
                w = np.exp(-0.5 * ((lags - c_lags[t]) / 1.4) ** 2)
                k_mat[t] = w / np.sum(w)

        kernels[eid] = k_mat

    # 3. Mixing fractions
    rho_dict: dict[str, np.ndarray] = {}
    pi_dict: dict[str, np.ndarray] = {}

    # Node R: pure source
    rho_dict["R"] = np.ones(T)
    # Node A: 20% local recharge, 80% from R
    rho_dict["A"] = np.full(T, 0.20)
    pi_dict["R->A"] = np.full(T, 0.80)
    # Node B: 15% local recharge, 85% from A
    rho_dict["B"] = np.full(T, 0.15)
    pi_dict["A->B"] = np.full(T, 0.85)
    # Node C: 10% local recharge, 60% from B, 30% from bypass A
    rho_dict["C"] = np.full(T, 0.10)
    pi_dict["B->C"] = np.full(T, 0.60)
    pi_dict["A->C"] = np.full(T, 0.30)

    # Forward simulate latent signals.  Local recharge histories are explicit
    # observed inputs when the declared experiment makes them available.
    node_signals: dict[str, np.ndarray] = {}
    node_signals["R"] = regional_forcing

    # Node A
    convolved_ra = convolve_causal_dynamic_kernel(
        node_signals["R"], kernels["R->A"], lags, step_days=cfg.step_days
    )
    local_a = 0.5 * np.cos(w0 * times) + 0.1 * rng.standard_normal(T)
    node_signals["A"] = rho_dict["A"] * local_a + pi_dict["R->A"] * convolved_ra

    # Node B
    convolved_ab = convolve_causal_dynamic_kernel(
        node_signals["A"], kernels["A->B"], lags, step_days=cfg.step_days
    )
    local_b = 0.3 * np.sin(w0 * times + 1.0) + 0.1 * rng.standard_normal(T)
    node_signals["B"] = rho_dict["B"] * local_b + pi_dict["A->B"] * convolved_ab

    # Node C (multi-path mixing)
    convolved_bc = convolve_causal_dynamic_kernel(
        node_signals["B"], kernels["B->C"], lags, step_days=cfg.step_days
    )
    convolved_ac = convolve_causal_dynamic_kernel(
        node_signals["A"], kernels["A->C"], lags, step_days=cfg.step_days
    )
    local_c = 0.2 * np.cos(w0 * times - 0.5) + 0.1 * rng.standard_normal(T)
    node_signals["C"] = rho_dict["C"] * local_c + pi_dict["B->C"] * convolved_bc + pi_dict["A->C"] * convolved_ac

    # Calculate truth metrics
    young_mask = lags <= cfg.young_water_cutoff_steps * cfg.step_days
    fy_truth: dict[str, np.ndarray] = {}
    mean_age_truth: dict[str, np.ndarray] = {}
    for eid, k_mat in kernels.items():
        fy_truth[eid] = np.sum(k_mat[:, young_mask], axis=-1)
        mean_age_truth[eid] = k_mat @ lags

    # 4. Generate noisy, masked observations
    train_end = int(T * cfg.training_fraction)
    noise_std = cfg.observation_noise_std * (3.0 if cfg.scenario == "high_observation_noise" else 1.0)
    obs_prob = cfg.observation_probability * (0.35 if cfg.scenario == "sparse_sampling" else 1.0)

    node_obs: dict[str, np.ndarray] = {}
    node_masks: dict[str, np.ndarray] = {}

    for n in ("A", "B", "C"):
        clean_sig = node_signals[n]
        noise = rng.normal(0.0, noise_std, T)
        noisy = clean_sig + noise
        mask = rng.random(T) < obs_prob
        # Keep future holdout
        obs_array = np.array(noisy, copy=True)
        obs_array[~mask] = np.nan
        node_obs[n] = obs_array
        node_masks[n] = mask

    # Local input histories are part of the visible observation contract unless
    # the declared stressor is specifically that recharge is unmodelled.  They
    # are not hidden scenario labels: the estimator only receives this mapping.
    local_inputs_truth = {"A": local_a, "B": local_b, "C": local_c}
    local_inputs_observed: dict[str, np.ndarray] = {}
    if cfg.scenario != "unmodelled_local_recharge":
        for node, series in local_inputs_truth.items():
            local_inputs_observed[node] = np.asarray(series, dtype=float).copy()

    # Candidate edges based on scenario
    if cfg.scenario == "reversed_graph":
        cand_edges = (("A", "R"), ("B", "A"), ("C", "B"), ("C", "A"))
    elif cfg.scenario == "random_graph":
        cand_edges = (("R", "C"), ("C", "A"), ("B", "R"))
    elif cfg.scenario == "missing_edge":
        cand_edges = (("R", "A"), ("A", "B"), ("B", "C"))  # Missing bypass A->C
    else:
        cand_edges = true_edges

    # Create cryptographic commitment
    truth_summary = {
        "scenario": cfg.scenario,
        "seed": cfg.seed,
        "arrays": {
            "forcing": _array_digest(regional_forcing),
            "node_signals": {node: _array_digest(values) for node, values in sorted(node_signals.items())},
            "local_inputs": {node: _array_digest(values) for node, values in sorted(local_inputs_truth.items())},
            "edge_kernels": {eid: _array_digest(values) for eid, values in sorted(kernels.items())},
            "edge_mixing": {eid: _array_digest(values) for eid, values in sorted(pi_dict.items())},
            "local_fractions": {node: _array_digest(values) for node, values in sorted(rho_dict.items())},
        },
    }
    commit = hashlib.sha256(json.dumps(truth_summary, sort_keys=True).encode()).hexdigest()
    case_id = f"DYN-{cfg.scenario[:6].upper()}-{cfg.seed}"

    truth = DynamicKernelTruth(
        case_id=case_id,
        scenario=cfg.scenario,
        time_steps=times,
        true_forcing=regional_forcing,
        node_signals=node_signals,
        local_inputs=local_inputs_truth,
        true_edge_kernels=kernels,
        local_recharge_fractions=rho_dict,
        edge_mixing_fractions=pi_dict,
        young_water_fractions=fy_truth,
        mean_transit_times=mean_age_truth,
        sealed_commitment=commit,
    )

    observations = DynamicKernelObservations(
        case_id=case_id,
        scenario=cfg.scenario,
        time_steps=times,
        regional_forcing=obs_forcing,
        node_observations=node_obs,
        observation_masks=node_masks,
        candidate_edges=cand_edges,
        training_end_step=train_end,
        sealed_commitment=commit,
        metadata={
            "max_lag_steps": cfg.max_lag_steps,
            "lag_grid_days": lags.tolist(),
            "young_water_cutoff_days": cfg.young_water_cutoff_steps * cfg.step_days,
            "season_period_steps": cfg.season_period_steps,
            "step_days": cfg.step_days,
            "training_end_step": train_end,
        },
        local_inputs=local_inputs_observed,
        observation_commitment=hashlib.sha256(
            json.dumps(
                {
                    "case_id": case_id,
                    "forcing": _array_digest(obs_forcing),
                    "node_observations": {
                        node: _array_digest(values) for node, values in sorted(node_obs.items())
                    },
                    "observation_masks": {
                        node: _array_digest(values.astype(np.uint8))
                        for node, values in sorted(node_masks.items())
                    },
                    "local_inputs": {
                        node: _array_digest(values) for node, values in sorted(local_inputs_observed.items())
                    },
                    "candidate_edges": list(cand_edges),
                    "training_end_step": train_end,
                },
                sort_keys=True,
            ).encode()
        ).hexdigest(),
    )

    return truth, observations


def evaluate_dynamic_kernel_benchmark(
    truth: DynamicKernelTruth,
    observations: DynamicKernelObservations,
    recoveries: Mapping[str, DynamicTTDRecovery],
) -> dict[str, Any]:
    """Score recoveries against sealed truth under pre-registered benchmark metrics."""
    if truth.sealed_commitment != observations.sealed_commitment:
        raise ValueError("Cryptographic commitment mismatch between truth and observations!")

    per_edge_metrics: dict[str, Any] = {}
    justified_recoveries = 0
    correct_abstentions = 0
    unsupported_estimates = 0
    false_abstentions = 0
    w1_errors: list[float] = []
    fy_errors: list[float] = []
    coverage_checks: list[bool] = []

    # The scorer knows the generator scenario; the estimator does not.  Only
    # scenarios with an actual observation/topology corruption are treated as
    # negative controls.  ``tracer_conflict`` is not a single-tracer dynamic
    # corruption and must not be turned into an oracle abstention.
    topology_corrupted = observations.scenario in {
        "reversed_graph",
        "random_graph",
    }

    lags = np.asarray(
        observations.metadata.get(
            "lag_grid_days",
            np.arange(observations.metadata["max_lag_steps"] + 1, dtype=float)
            * float(observations.metadata.get("step_days", 1.0)),
        ),
        dtype=float,
    )

    for eid, rec in recoveries.items():
        forcing_corrupted_for_edge = observations.scenario == "wrong_forcing" and eid.startswith("R->")
        local_recharge_unavailable = observations.scenario == "unmodelled_local_recharge"
        edge_is_corrupted = topology_corrupted or forcing_corrupted_for_edge or local_recharge_unavailable
        if rec.status == "ABSTAIN":
            if edge_is_corrupted or eid not in truth.true_edge_kernels:
                correct_abstentions += 1
                diag = "CORRECT_ABSTENTION"
            else:
                false_abstentions += 1
                diag = "FALSE_ABSTENTION"
            per_edge_metrics[eid] = {"status": "ABSTAIN", "diagnostic": diag, "reasons": list(rec.reason_codes)}
        else:
            # Status is RECOVERED
            if edge_is_corrupted or eid not in truth.true_edge_kernels:
                unsupported_estimates += 1
                diag = "UNSUPPORTED_ESTIMATE"
            else:
                justified_recoveries += 1
                diag = "JUSTIFIED_RECOVERY"

            # Compute W1 distance vs true kernel
            if eid in truth.true_edge_kernels and rec.estimated_kernel is not None:
                true_k = truth.true_edge_kernels[eid]
                est_k = rec.estimated_kernel
                w1 = float(np.mean([wasserstein_ttd_w1(est_k[t], true_k[t], lags) for t in range(len(true_k))]))
                w1_errors.append(w1)

                # Fy error & coverage
                true_fy = float(np.mean(truth.young_water_fractions[eid]))
                if rec.young_water_fraction is not None:
                    fy_errors.append(abs(rec.young_water_fraction - true_fy))
                if rec.young_water_interval is not None:
                    cov = rec.young_water_interval[0] <= true_fy <= rec.young_water_interval[1]
                    coverage_checks.append(cov)

            per_edge_metrics[eid] = {
                "status": "RECOVERED",
                "diagnostic": diag,
                "young_water_fraction": rec.young_water_fraction,
                "mean_age": rec.mean_age,
                "forecast_r2": rec.forecast_r2,
                "forecast_rmse": rec.forecast_rmse,
            }

    summary = {
        "case_id": truth.case_id,
        "scenario": truth.scenario,
        "commitment_verified": True,
        "truth_commitment": truth.sealed_commitment,
        "observation_commitment": observations.observation_commitment,
        "justified_recoveries": justified_recoveries,
        "correct_abstentions": correct_abstentions,
        "unsupported_estimates": unsupported_estimates,
        "false_abstentions": false_abstentions,
        "mean_w1_distance": float(np.mean(w1_errors)) if w1_errors else None,
        "mean_fy_absolute_error": float(np.mean(fy_errors)) if fy_errors else None,
        "empirical_coverage": float(np.mean(coverage_checks)) if coverage_checks else None,
        "per_edge": per_edge_metrics,
    }
    return summary
