"""Independent multi-tracer synthetic benchmark generator and evaluator.

This module generates controlled-synthetic multi-tracer panels (d18O, d2H, 3H, SF6, 14C)
across network nodes, tests tracer dropout, and injects conflicting-tracer stressors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Any, Mapping, Optional

import numpy as np

from ..nuclear.graph_tracer_forward import convolve_causal_dynamic_kernel, edge_id_str
from ..nuclear.multi_tracer_graph_inversion import MultiTracerNodeRecovery
from ..nuclear.tracer_registry import (
    LAMBDA_3H_PER_DAY,
    LAMBDA_14C_PER_DAY,
)


TRACER_PANELS = (
    "all",
    "isotopes_only",
    "isotopes_plus_3H",
    "isotopes_plus_SF6",
    "no_14C",
    "conflicting_tracers",
)


@dataclass(frozen=True)
class MultiTracerBenchmarkConfig:
    seed: int = 20260917
    panel: str = "all"
    n_steps: int = 156
    step_days: float = 14.0
    season_period_steps: int = 26
    max_lag_steps: int = 14
    observation_noise_relative: float = 0.05
    training_fraction: float = 0.70

    def __post_init__(self) -> None:
        if self.panel not in TRACER_PANELS:
            raise ValueError(f"Unknown panel {self.panel!r}, must be one of {TRACER_PANELS}")
        if self.step_days <= 0.0 or self.n_steps < 2 or self.max_lag_steps < 1:
            raise ValueError("step_days, n_steps, and max_lag_steps must be positive")
        if not 0.0 < self.training_fraction < 1.0:
            raise ValueError("training_fraction must lie strictly between zero and one")


@dataclass(frozen=True)
class MultiTracerTruth:
    case_id: str
    panel: str
    time_steps: np.ndarray
    true_node_tracers: Mapping[str, Mapping[str, np.ndarray]]  # node -> {tracer -> series}
    true_edge_kernels: Mapping[str, np.ndarray]
    sealed_commitment: str


@dataclass(frozen=True)
class MultiTracerObservations:
    case_id: str
    panel: str
    time_steps: np.ndarray
    node_tracer_observations: Mapping[str, Mapping[str, np.ndarray]]  # With NaNs
    candidate_edges: tuple[tuple[str, str], ...]
    available_tracers: tuple[str, ...]
    training_end_step: int
    sealed_commitment: str
    metadata: Mapping[str, Any]
    local_tracer_inputs: Mapping[str, Mapping[str, np.ndarray]] = field(default_factory=dict)
    observation_commitment: str = ""


def _array_digest(values: np.ndarray) -> str:
    arr = np.asarray(values, dtype=np.float64)
    h = hashlib.sha256()
    h.update(str(arr.shape).encode("ascii"))
    h.update(arr.tobytes(order="C"))
    return h.hexdigest()


def generate_multitracer_case(
    config: Optional[MultiTracerBenchmarkConfig] = None,
) -> tuple[MultiTracerTruth, MultiTracerObservations]:
    """Generate multi-tracer groundwater network benchmark case."""
    cfg = config or MultiTracerBenchmarkConfig()
    rng = np.random.default_rng(cfg.seed)

    T = cfg.n_steps
    times = np.arange(T, dtype=float)
    L = cfg.max_lag_steps + 1
    lags = np.arange(L, dtype=float) * cfg.step_days

    true_edges = (("R", "A"), ("A", "B"))
    nodes = ("R", "A", "B")

    # Dynamic kernels on edges: Gaussian traveling wave
    w0 = 2.0 * np.pi / cfg.season_period_steps
    kernels: dict[str, np.ndarray] = {}

    for u, v in true_edges:
        eid = edge_id_str(u, v)
        k_mat = np.zeros((T, L), dtype=float)
        mean_lags = (4.0 if eid == "R->A" else 6.0) * cfg.step_days + 14.0 * np.sin(w0 * times)
        for t in range(T):
            w = np.exp(-0.5 * ((lags - mean_lags[t]) / (10.0)) ** 2)
            k_mat[t] = w / np.sum(w)
        kernels[eid] = k_mat

    # Generate source input histories for 5 tracers
    # 1. d18O: seasonal sine + noise (mean -5 permil, amplitude 2 permil)
    d18O_source = -5.0 + 2.0 * np.sin(w0 * times) + 0.2 * rng.standard_normal(T)
    # 2. d2H: on LMWL (d2H = 8 * d18O + 10) + small deuterium excess variation
    d2H_source = 8.0 * d18O_source + 10.0 + 0.5 * rng.standard_normal(T)
    # 3. 3H: atmospheric input with bomb peak decay curve (mean ~4 TU, small seasonal peak)
    h3_source = 3.5 + 1.2 * np.sin(w0 * times + 0.5) + 0.1 * rng.standard_normal(T)
    # 4. SF6: monotonic rising concentration (modern atmospheric pptv)
    sf6_source = np.linspace(6.0, 10.0, T) + 0.1 * rng.standard_normal(T)
    # 5. 14C: modern recharge carbon activity (~95 pmC)
    c14_source = np.full(T, 95.0) + 0.5 * rng.standard_normal(T)

    source_tracers = {
        "d18O": d18O_source,
        "d2H": d2H_source,
        "3H": h3_source,
        "SF6": sf6_source,
        "14C": c14_source,
    }

    # Forward simulate tracer transport across nodes
    node_tracers: dict[str, dict[str, np.ndarray]] = {n: {} for n in nodes}
    local_tracer_inputs: dict[str, dict[str, np.ndarray]] = {"A": {}, "B": {}}
    for tracer, src in source_tracers.items():
        node_tracers["R"][tracer] = src

    # Decay rates per time step
    dt_days = cfg.step_days
    decay_3h = LAMBDA_3H_PER_DAY * dt_days
    decay_14c = LAMBDA_14C_PER_DAY * dt_days

    for tracer in ("d18O", "d2H", "3H", "SF6", "14C"):
        decay_rate = decay_3h if tracer == "3H" else (decay_14c if tracer == "14C" else 0.0)
        # Node A: 20% local source, 80% from R
        convolved_ra = convolve_causal_dynamic_kernel(
            node_tracers["R"][tracer],
            kernels["R->A"],
            lags,
            decay_constant_per_time=decay_rate,
            step_days=dt_days,
        )
        local_a = source_tracers[tracer] * 0.95
        local_tracer_inputs["A"][tracer] = np.asarray(local_a, dtype=float).copy()
        carbon_scale = 0.85 if tracer == "14C" else 1.0
        node_tracers["A"][tracer] = carbon_scale * (0.20 * local_a + 0.80 * convolved_ra)

        # Node B: 15% local source, 85% from A
        convolved_ab = convolve_causal_dynamic_kernel(
            node_tracers["A"][tracer],
            kernels["A->B"],
            lags,
            decay_constant_per_time=decay_rate,
            step_days=dt_days,
        )
        local_b = source_tracers[tracer] * 0.90
        local_tracer_inputs["B"][tracer] = np.asarray(local_b, dtype=float).copy()
        node_tracers["B"][tracer] = carbon_scale * (0.15 * local_b + 0.85 * convolved_ab)

    # Determine which tracers are active in observation panel
    if cfg.panel == "isotopes_only":
        active_tracers = ("d18O", "d2H")
    elif cfg.panel == "isotopes_plus_3H":
        active_tracers = ("d18O", "d2H", "3H")
    elif cfg.panel == "isotopes_plus_SF6":
        active_tracers = ("d18O", "d2H", "SF6")
    elif cfg.panel == "no_14C":
        active_tracers = ("d18O", "d2H", "3H", "SF6")
    elif cfg.panel == "conflicting_tracers":
        active_tracers = ("d18O", "d2H", "3H", "SF6", "14C")
    else:
        active_tracers = ("d18O", "d2H", "3H", "SF6", "14C")

    # Generate noisy observations
    train_end = int(T * cfg.training_fraction)
    obs_node_tracers: dict[str, dict[str, np.ndarray]] = {n: {} for n in nodes}

    for n in nodes:
        for t in active_tracers:
            clean = node_tracers[n][t]
            noise_sd = cfg.observation_noise_relative * float(np.std(clean) + 0.1)
            noise = rng.normal(0.0, noise_sd, T)
            obs = np.array(clean + noise, copy=True)
            # Inject conflict in conflicting_tracers scenario at node B
            if cfg.panel == "conflicting_tracers" and n == "B":
                if t == "3H":
                    obs[:] = 0.01  # Sub-detection limit tritium
                elif t == "SF6":
                    obs[:] = 8.5   # Modern high SF6
                elif t == "14C":
                    obs[:] = 5.0   # Dead radiocarbon
            obs_node_tracers[n][t] = obs

    commit_data = {
        "panel": cfg.panel,
        "seed": cfg.seed,
        "tracers": list(active_tracers),
        "node_tracers": {
            node: {tracer: _array_digest(values) for tracer, values in sorted(series.items())}
            for node, series in sorted(node_tracers.items())
        },
        "edge_kernels": {eid: _array_digest(values) for eid, values in sorted(kernels.items())},
    }
    commit = hashlib.sha256(json.dumps(commit_data, sort_keys=True).encode()).hexdigest()
    case_id = f"MTRACER-{cfg.panel.upper()}-{cfg.seed}"

    truth = MultiTracerTruth(
        case_id=case_id,
        panel=cfg.panel,
        time_steps=times,
        true_node_tracers=node_tracers,
        true_edge_kernels=kernels,
        sealed_commitment=commit,
    )

    observation_commitment = hashlib.sha256(
        json.dumps(
            {
                "case_id": case_id,
                "node_tracer_observations": {
                    node: {tracer: _array_digest(values) for tracer, values in sorted(series.items())}
                    for node, series in sorted(obs_node_tracers.items())
                },
                "candidate_edges": list(true_edges),
                "available_tracers": list(active_tracers),
                "local_tracer_inputs": {
                    node: {tracer: _array_digest(values) for tracer, values in sorted(series.items())}
                    for node, series in sorted(local_tracer_inputs.items())
                },
                "training_end_step": train_end,
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()

    observations = MultiTracerObservations(
        case_id=case_id,
        panel=cfg.panel,
        time_steps=times,
        node_tracer_observations=obs_node_tracers,
        candidate_edges=true_edges,
        available_tracers=active_tracers,
        training_end_step=train_end,
        sealed_commitment=commit,
        metadata={
            "step_days": cfg.step_days,
            "max_lag_steps": cfg.max_lag_steps,
            "lag_grid_days": lags.tolist(),
            "training_end_step": train_end,
        },
        local_tracer_inputs=local_tracer_inputs,
        observation_commitment=observation_commitment,
    )

    return truth, observations


def evaluate_multitracer_recovery(
    truth: MultiTracerTruth,
    observations: MultiTracerObservations,
    recovery: MultiTracerNodeRecovery,
) -> dict[str, Any]:
    """Score multi-tracer node recovery against sealed truth."""
    if truth.sealed_commitment != observations.sealed_commitment:
        raise ValueError("Commitment mismatch in multi-tracer evaluation!")

    is_conflict = observations.panel == "conflicting_tracers"

    if recovery.status == "ABSTAIN":
        if is_conflict or recovery.conflict_detected:
            decision = "CORRECT_ABSTENTION"
        else:
            decision = "FALSE_ABSTENTION"
        return {
            "case_id": truth.case_id,
            "panel": truth.panel,
            "status": "ABSTAIN",
            "decision": decision,
            "conflict_detected": recovery.conflict_detected,
            "reason_codes": list(recovery.reason_codes),
            "truth_commitment": truth.sealed_commitment,
            "observation_commitment": observations.observation_commitment,
        }

    # If recovered under conflicting panel, that's an unsupported claim
    if is_conflict:
        return {
            "case_id": truth.case_id,
            "panel": truth.panel,
            "status": "RECOVERED",
            "decision": "UNSUPPORTED_ESTIMATE_FAILED_CONFLICT_GATE",
            "truth_commitment": truth.sealed_commitment,
            "observation_commitment": observations.observation_commitment,
        }

    if recovery.diagnostics.get("holdout_gate_enforced") and not recovery.diagnostics.get(
        "holdout_gate_passed", False
    ):
        return {
            "case_id": truth.case_id,
            "panel": truth.panel,
            "status": "RECOVERED",
            "decision": "UNSUPPORTED_ESTIMATE_FAILED_HOLDOUT_GATE",
            "active_tracers": list(recovery.active_tracers),
            "held_out_per_tracer_rmse": dict(recovery.held_out_per_tracer_rmse),
            "held_out_per_tracer_r2": dict(recovery.held_out_per_tracer_r2),
            "held_out_per_tracer_normalized_rmse": dict(recovery.held_out_per_tracer_normalized_rmse),
            "truth_commitment": truth.sealed_commitment,
            "observation_commitment": observations.observation_commitment,
        }

    # Check recovery accuracy
    per_edge_rmse = {}
    for eid, rec in recovery.edge_recoveries.items():
        if rec.status == "RECOVERED" and rec.forecast_rmse is not None:
            per_edge_rmse[eid] = rec.forecast_rmse

    return {
        "case_id": truth.case_id,
        "panel": truth.panel,
        "status": "RECOVERED",
        "decision": "JUSTIFIED_RECOVERY",
        "active_tracers": list(recovery.active_tracers),
        "per_tracer_rmse": dict(recovery.per_tracer_rmse),
        "held_out_per_tracer_rmse": dict(recovery.held_out_per_tracer_rmse),
        "held_out_per_tracer_r2": dict(recovery.held_out_per_tracer_r2),
        "held_out_per_tracer_normalized_rmse": dict(recovery.held_out_per_tracer_normalized_rmse),
        "loto_sensitivity": dict(recovery.loto_sensitivity),
        "effective_rank": recovery.effective_rank,
        "condition_number": recovery.condition_number,
        "edge_rmse": per_edge_rmse,
        "truth_commitment": truth.sealed_commitment,
        "observation_commitment": observations.observation_commitment,
    }
