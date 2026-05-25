"""Causal discovery layer for groundwater flow-direction support.

Guarded modes (in order of data richness):

1. ``insufficient_data``
   Default for single-snapshot data. Adds no penalty, records reason.

2. ``temporal_partial_corr``
   Uses repeated/time-series samples to compute partial correlation
   between upstream and downstream chemistry controlling for shared
   drivers (seasonality, recharge events).

3. ``lagged_tracer_support``
   Tests whether upstream tracer history at lag tau predicts downstream
   tracer concentration better than the reverse direction.

4. ``common_driver_screen``
   Flags edges where chemistry similarity could be explained by shared
   recharge, lithology, or a common external source rather than direct
   flow connectivity.

Design principle: this is NOT a default topology engine. Static
one-sample-per-well data is insufficient for true PC/FCI algorithms.
The module is guarded and provides diagnostic scores only when
sufficient temporal or multilevel data is available.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np

from ..config import Config
from ..log import get_logger

logger = get_logger("causal.discovery")


@dataclass
class CausalResult:
    """Per-edge causal assessment result."""

    support_score: float = 0.0
    confounded_score: float = 0.0
    p_value: Optional[float] = None
    method: str = "insufficient_data"
    n_observations: int = 0
    status: str = "insufficient_data"
    diagnostics: Dict[str, Any] = field(default_factory=dict)


def _has_temporal_data(
    sample: Mapping[str, object],
) -> bool:
    """Check if a sample has temporal data (time series or repeated measures)."""
    temporal_keys = {"date", "timestamp", "sample_date", "time", "datetime"}
    has_time = any(k in sample for k in temporal_keys)

    # Check for repeated measures at the same site
    if "replicates" in sample or "n_measurements" in sample:
        return True

    return has_time


def _prepare_temporal_vectors(
    upstream_samples: Sequence[Mapping[str, object]],
    downstream_samples: Sequence[Mapping[str, object]],
    ion_order: List[str],
    config: Config,
) -> Optional[Dict[str, np.ndarray]]:
    """Prepare aligned temporal vectors for causal analysis.

    Returns None if temporal data is insufficient.
    """
    if len(upstream_samples) < 2 or len(downstream_samples) < 2:
        return None

    min_obs = int(getattr(config, "causal_min_observations", 5))
    if len(upstream_samples) < min_obs or len(downstream_samples) < min_obs:
        return None

    # Try to align by timestamp
    from ..data.schema import vector_from_sample

    def _get_vectors(samples):
        vectors = []
        timestamps = []
        for s in samples:
            val, _ = vector_from_sample(
                s, ion_order,
                config.missing_policy, config.detection_limit_policy,
            )
            if val is not None:
                vectors.append(np.array(val, dtype=float))
                # Try to extract timestamp
                ts = s.get("date") or s.get("timestamp") or s.get("sample_date")
                timestamps.append(ts)
        return vectors, timestamps

    up_vecs, up_ts = _get_vectors(upstream_samples)
    down_vecs, down_ts = _get_vectors(downstream_samples)

    if len(up_vecs) < min_obs or len(down_vecs) < min_obs:
        return None

    # If timestamps are available and match, align them
    if up_ts and down_ts and all(t is not None for t in up_ts) and all(t is not None for t in down_ts):
        # Build time-indexed matrices
        try:
            up_idx = {str(t): v for t, v in zip(up_ts, up_vecs)}
            down_idx = {str(t): v for t, v in zip(down_ts, down_vecs)}
            common_times = sorted(set(up_idx.keys()) & set(down_idx.keys()))
            if len(common_times) >= min_obs:
                up_aligned = np.array([up_idx[t] for t in common_times])
                down_aligned = np.array([down_idx[t] for t in common_times])
                return {"up": up_aligned, "down": down_aligned, "n": len(common_times)}
        except Exception:
            pass

    # No alignment possible — use summary statistics
    return None


def _temporal_partial_correlation(
    up_vectors: np.ndarray,
    down_vectors: np.ndarray,
) -> Dict[str, Any]:
    """Compute partial correlation between upstream and downstream chemistry.

    Uses the mean species vector over time as a proxy when aligned
    timestamps are not available.
    """
    n_up, n_dim = up_vectors.shape
    n_down, _ = down_vectors.shape

    # Simple correlation between mean vectors
    up_mean = np.mean(up_vectors, axis=0)
    down_mean = np.mean(down_vectors, axis=0)

    # Correlation across species
    if np.std(up_mean) > 1e-12 and np.std(down_mean) > 1e-12:
        corr = np.corrcoef(up_mean, down_mean)[0, 1]
    else:
        corr = 0.0

    # Partial correlation controlling for "time" (index as proxy)
    if n_up >= 5 and n_down >= 5:
        time_idx = np.arange(min(n_up, n_down))
        try:
            from scipy.stats import pearsonr

            # Correlate per-species trends
            correlations = []
            p_values = []
            for d in range(n_dim):
                u_slice = up_vectors[: len(time_idx), d]
                d_slice = down_vectors[: len(time_idx), d]
                if np.std(u_slice) > 1e-12 and np.std(d_slice) > 1e-12:
                    r, p = pearsonr(u_slice, d_slice)
                    if np.isfinite(r):
                        correlations.append(r)
                        p_values.append(p)

            avg_corr = float(np.mean(correlations)) if correlations else 0.0
            avg_p = float(np.mean(p_values)) if p_values else 1.0
        except Exception:
            avg_corr = float(np.clip(corr, -1, 1))
            avg_p = 1.0
    else:
        avg_corr = float(np.clip(corr, -1, 1))
        avg_p = 1.0

    return {"partial_correlation": avg_corr, "p_value": avg_p, "n_pairs": min(n_up, n_down)}


def _lagged_tracer_support(
    upstream_history: np.ndarray,
    downstream_history: np.ndarray,
    tracer_idx: int = 5,  # Default: Cl
    max_lag: int = 3,
) -> Dict[str, Any]:
    """Test whether upstream tracer history predicts downstream tracer.

    For each lag tau, computes correlation between upstream[t - tau]
    and downstream[t]. Returns the maximum correlation and optimal lag.
    """
    n = min(upstream_history.shape[0], downstream_history.shape[0])
    if n <= max_lag + 1:
        return {"lag_support": 0.0, "optimal_lag": 0, "max_correlation": 0.0}

    if tracer_idx >= upstream_history.shape[1]:
        return {"lag_support": 0.0, "optimal_lag": 0, "max_correlation": 0.0}

    u_tracer = upstream_history[:, tracer_idx]
    d_tracer = downstream_history[:, tracer_idx]

    best_corr = 0.0
    best_lag = 0

    for lag in range(1, max_lag + 1):
        # Upstream at t-lag, downstream at t
        u_lagged = u_tracer[: n - lag]
        d_current = d_tracer[lag:n]

        if np.std(u_lagged) > 1e-12 and np.std(d_current) > 1e-12:
            corr = np.corrcoef(u_lagged, d_current)[0, 1]
            if np.isfinite(corr) and corr > best_corr:
                best_corr = corr
                best_lag = lag

    # Compare with reverse direction (downstream predicting upstream)
    reverse_corr = 0.0
    for lag in range(1, max_lag + 1):
        d_lagged = d_tracer[: n - lag]
        u_current = u_tracer[lag:n]
        if np.std(d_lagged) > 1e-12 and np.std(u_current) > 1e-12:
            corr = np.corrcoef(d_lagged, u_current)[0, 1]
            if np.isfinite(corr):
                reverse_corr = max(reverse_corr, corr)

    # Support = forward_corr - reverse_corr (positive favors u->v direction)
    support = max(0.0, best_corr - reverse_corr)

    return {
        "lag_support": float(support),
        "optimal_lag": best_lag,
        "max_correlation": float(best_corr),
        "reverse_correlation": float(reverse_corr),
    }


def _common_driver_screen(
    upstream: Mapping[str, object],
    downstream: Mapping[str, object],
    config: Config,
) -> Dict[str, Any]:
    """Screen for common external drivers that could explain chemistry similarity.

    Checks:
    1. Same aquifer unit (lithology)
    2. Similar recharge elevation (topography)
    3. Similar well depth (sampling depth)
    """
    flags = []

    # Check aquifer unit
    aquifer_key = getattr(config, "edge_aquifer_key", "aquifer_unit")
    u_aq = upstream.get(aquifer_key)
    v_aq = downstream.get(aquifer_key)
    same_aquifer = False
    if u_aq is not None and v_aq is not None and str(u_aq) == str(v_aq):
        same_aquifer = True
        flags.append("shared_aquifer")

    # Check screen depth similarity
    screen_key = getattr(config, "edge_screen_depth_key", "screen_depth")
    u_sd = _safe_float(upstream.get(screen_key))
    v_sd = _safe_float(downstream.get(screen_key))
    similar_depth = False
    if u_sd is not None and v_sd is not None:
        if abs(u_sd - v_sd) < 10.0:  # Within 10m
            similar_depth = True
            flags.append("similar_screen_depth")

    # Check elevation proximity
    elev_key = getattr(config, "edge_elevation_key", "elevation")
    u_elev = _safe_float(upstream.get(elev_key))
    v_elev = _safe_float(downstream.get(elev_key))
    similar_elev = False
    if u_elev is not None and v_elev is not None:
        if abs(u_elev - v_elev) < 20.0:  # Within 20m
            similar_elev = True
            flags.append("similar_elevation")

    n_flags = len(flags)
    confound_score = min(1.0, n_flags / 3.0) if n_flags > 1 else 0.0

    return {
        "common_driver_flags": flags,
        "confound_score": confound_score,
        "same_aquifer": same_aquifer,
        "similar_depth": similar_depth,
        "similar_elevation": similar_elev,
    }


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def assess_edge_causality(
    upstream_sample: Mapping[str, object],
    downstream_sample: Mapping[str, object],
    upstream_history: Optional[Sequence[Mapping[str, object]]] = None,
    downstream_history: Optional[Sequence[Mapping[str, object]]] = None,
    config: Optional[Config] = None,
) -> CausalResult:
    """Assess causal direction support for a candidate flow edge.

    Parameters
    ----------
    upstream_sample : mapping
        Single-sample data for the upstream well.
    downstream_sample : mapping
        Single-sample data for the downstream well.
    upstream_history : sequence of mappings, optional
        Time-series data for the upstream well.
    downstream_history : sequence of mappings, optional
        Time-series data for the downstream well.
    config : Config, optional
        Hydrosheaf configuration.

    Returns
    -------
    CausalResult
        Per-edge causal assessment.
    """
    cfg = config or Config()
    min_obs = int(getattr(cfg, "causal_min_observations", 5))

    # Always run common driver screen (works with static data)
    driver_info = _common_driver_screen(upstream_sample, downstream_sample, cfg)

    has_up_history = upstream_history is not None and len(upstream_history) >= min_obs
    has_down_history = downstream_history is not None and len(downstream_history) >= min_obs

    if not has_up_history or not has_down_history:
        # insufficient_data: record driver info for diagnostics but do NOT
        # apply a confounded_score penalty — static data cannot distinguish
        # causal flow from shared-driver explanations.
        return CausalResult(
            support_score=0.0,
            confounded_score=0.0,
            p_value=None,
            method="insufficient_data",
            n_observations=max(
                len(upstream_history) if upstream_history else 0,
                len(downstream_history) if downstream_history else 0,
            ),
            status="insufficient_data",
            diagnostics={
                "reason": "Insufficient temporal data for causal discovery",
                "min_required": min_obs,
                **driver_info,
            },
        )

    # Temporal data available — run partial correlation
    ion_order = list(cfg.ion_order)
    temporal = _prepare_temporal_vectors(
        list(upstream_history), list(downstream_history), ion_order, cfg
    )

    if temporal is not None and temporal["n"] >= min_obs:
        corr_info = _temporal_partial_correlation(temporal["up"], temporal["down"])

        # Lagged tracer support
        lag_info = _lagged_tracer_support(
            temporal["up"], temporal["down"],
            tracer_idx=ion_order.index("Cl") if "Cl" in ion_order else 5,
        )

        support = float(corr_info.get("partial_correlation", 0.0))
        p_val = float(corr_info.get("p_value", 1.0))

        return CausalResult(
            support_score=max(0.0, support),
            confounded_score=driver_info["confound_score"],
            p_value=p_val,
            method="temporal_partial_corr",
            n_observations=temporal["n"],
            status="temporal_partial_corr",
            diagnostics={
                "partial_correlation": support,
                "p_value": p_val,
                **lag_info,
                **driver_info,
            },
        )

    # Temporal data exists but can't be aligned — use correlation of means
    from ..data.schema import vector_from_sample

    u_vals, _ = vector_from_sample(
        upstream_sample, ion_order, cfg.missing_policy, cfg.detection_limit_policy,
    )
    v_vals, _ = vector_from_sample(
        downstream_sample, ion_order, cfg.missing_policy, cfg.detection_limit_policy,
    )

    if u_vals is not None and v_vals is not None:
        u_arr = np.array(u_vals, dtype=float)
        v_arr = np.array(v_vals, dtype=float)
        if np.std(u_arr) > 1e-12 and np.std(v_arr) > 1e-12:
            corr = np.corrcoef(u_arr, v_arr)[0, 1]
            support = max(0.0, float(corr))
        else:
            support = 0.0
    else:
        support = 0.0

    return CausalResult(
        support_score=support,
        confounded_score=driver_info["confound_score"],
        p_value=None,
        method="temporal_partial_corr",
        n_observations=max(
            len(upstream_history) if upstream_history else 0,
            len(downstream_history) if downstream_history else 0,
        ),
        status="temporal_partial_corr",
        diagnostics={
            "partial_correlation": support,
            "note": "Temporal data present but not time-aligned; using mean correlation",
            **driver_info,
        },
    )


def compute_causal_support(
    upstream_sample: Mapping[str, object],
    downstream_sample: Mapping[str, object],
    upstream_history: Optional[Sequence[Mapping[str, object]]] = None,
    downstream_history: Optional[Sequence[Mapping[str, object]]] = None,
    config: Optional[Config] = None,
) -> Dict[str, object]:
    """Convenience wrapper returning a dict of causal fields for edge attrs.

    Returns
    -------
    dict with keys matching the EdgeResult causal fields.
    """
    result = assess_edge_causality(
        upstream_sample, downstream_sample,
        upstream_history, downstream_history, config,
    )
    return {
        "causal_support_score": result.support_score,
        "causal_confounded_score": result.confounded_score,
        "causal_p_value": result.p_value,
        "causal_method": result.method,
        "causal_n_observations": result.n_observations,
        "causal_status": result.status,
    }
