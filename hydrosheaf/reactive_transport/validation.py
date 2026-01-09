"""
Forward validation of inverse model results using reactive transport.
"""

from typing import Dict, List, Optional

from . import KineticParameters, ReactiveTransportResult, ValidationSummary
from .kinetic_phreeqc import build_kinetic_block, run_phreeqc_kinetic
from .metrics import (
    check_thermodynamic_consistency,
    compute_consistency_metrics,
    compute_per_ion_metrics,
)


def validate_edge_forward(
    edge_result: "EdgeResult",  # type: ignore
    x_u: List[float],
    x_v_observed: List[float],
    config: "Config",  # type: ignore
    kinetic_params: Optional[Dict[str, KineticParameters]] = None,
    residence_time_days: Optional[float] = None,
) -> ReactiveTransportResult:
    """
    Run forward validation for a single edge.

    Parameters
    ----------
    edge_result : EdgeResult
        Inverse model result with extents and transport params
    x_u : List[float]
        Upstream observed composition
    x_v_observed : List[float]
        Downstream observed composition
    config : Config
        Hydrosheaf configuration
    kinetic_params : Optional[Dict]
        Rate law parameters
    residence_time_days : Optional[float]
        Residence time (if None, estimate from distance/velocity)

    Returns
    -------
    ReactiveTransportResult
        Forward simulation with consistency metrics

    Mathematical Implementation
    ---------------------------
    1. Apply transport to upstream:
       x_transport = T(γ) · x_u
       (using evaporation_affine or mixing_affine from transport.py)

    2. Build kinetic block from inverse extents

    3. Run PHREEQC kinetic with x_transport as initial

    4. Extract final composition x_v^fwd

    5. Compute metrics:
       RMSE = sqrt(mean((x_v^fwd - x_v^obs)²))
       NSE = 1 - sum((x_v^fwd - x_v^obs)²) / sum((x_v^obs - mean(x_v^obs))²)
       PBIAS = 100 * sum(x_v^fwd - x_v^obs) / sum(x_v^obs)

    6. Check thermodynamic consistency:
       For each mineral k:
         If ξ_k > 0 (dissolution): require SI_k(0) < τ (was undersaturated)
         If ξ_k < 0 (precip): require SI_k(0) > -τ (was supersaturated)
    """
    from ..models.transport import evaporation_affine, mixing_affine

    # Use residence time if provided, otherwise default
    if residence_time_days is None:
        residence_time_days = getattr(config, "rt_default_residence_time", 30.0)

    # Apply transport to get post-transport composition
    if edge_result.transport_model == "evap" and edge_result.gamma is not None:
        x_transport = evaporation_affine(x_u, edge_result.gamma)
    elif edge_result.transport_model == "mix" and edge_result.f is not None:
        endmember_id = edge_result.endmember_id or list(config.mixing_endmembers.keys())[0]
        endmember = config.mixing_endmembers[endmember_id]
        x_transport = mixing_affine(x_u, endmember, edge_result.f)
    else:
        x_transport = x_u  # No transport

    # Build kinetic block from reaction extents
    reaction_labels = edge_result.z_labels
    extents = edge_result.z_extents

    if not reaction_labels or not extents:
        # No reactions, just return transport result
        result = ReactiveTransportResult(
            edge_id=edge_result.edge_id,
            simulator="none",
            inverse_extents=[],
            inverse_residence_time_days=residence_time_days,
            forward_x_v=x_transport,
        )

        # Compute metrics
        metrics = compute_consistency_metrics(x_transport, x_v_observed, config.weights)
        result.rmse = metrics["rmse"]
        result.nse = metrics["nse"]
        result.pbias = metrics["pbias"]
        result.thermodynamic_consistent = True

        return result

    # Build PHREEQC kinetic input
    kinetics_block = build_kinetic_block(
        reaction_labels=reaction_labels,
        extents=extents,
        residence_time_days=residence_time_days,
        kinetic_params=kinetic_params,
        temperature_c=config.temp_default_c,
    )

    # Prepare initial solution dict
    initial_solution = {}
    for i, ion in enumerate(config.ion_order):
        if i < len(x_transport):
            initial_solution[ion] = x_transport[i]

    # Run PHREEQC kinetic simulation
    phreeqc_result = run_phreeqc_kinetic(
        initial_solution=initial_solution,
        kinetics_block=kinetics_block,
        residence_time_days=residence_time_days,
        config=config,
    )

    # Create result object
    result = ReactiveTransportResult(
        edge_id=edge_result.edge_id,
        simulator="phreeqc_kinetic",
        inverse_extents=extents,
        inverse_residence_time_days=residence_time_days,
        forward_x_v=[],
    )

    if phreeqc_result["success"]:
        # Extract final composition
        x_v_forward = phreeqc_result["final_composition"]
        result.forward_x_v = x_v_forward

        # Store SI trajectory
        result.forward_si_trajectory = phreeqc_result["si_series"]

        # Compute consistency metrics
        if x_v_forward:
            metrics = compute_consistency_metrics(x_v_forward, x_v_observed, config.weights)
            result.rmse = metrics["rmse"]
            result.nse = metrics["nse"]
            result.pbias = metrics["pbias"]

            # Per-ion metrics
            per_ion = compute_per_ion_metrics(x_v_forward, x_v_observed, config.ion_order)
            result.per_ion_error = [per_ion[ion]["error"] for ion in config.ion_order if ion in per_ion]
            result.per_ion_bias = [per_ion[ion]["rel_error"] for ion in config.ion_order if ion in per_ion]

        # Check thermodynamic consistency
        is_consistent, violations = check_thermodynamic_consistency(
            extents=extents,
            si_initial=edge_result.si_u,
            reaction_labels=reaction_labels,
            si_threshold=config.si_threshold_tau,
        )
        result.thermodynamic_consistent = is_consistent

    else:
        # Forward simulation failed
        result.rmse = float("inf")
        result.nse = float("-inf")
        result.pbias = 100.0
        result.thermodynamic_consistent = False

    return result


def validate_network_forward(
    edge_results: List["EdgeResult"],  # type: ignore
    samples: Dict[str, Dict[str, float]],
    config: "Config",  # type: ignore
    kinetic_params: Optional[Dict[str, KineticParameters]] = None,
) -> ValidationSummary:
    """
    Run forward validation for all edges in network.

    Parameters
    ----------
    edge_results : List[EdgeResult]
        All inverse model results
    samples : Dict[str, Dict]
        Sample data keyed by sample_id
    config : Config
        Configuration
    kinetic_params : Optional[Dict]
        Kinetic parameters

    Returns
    -------
    ValidationSummary
        Network-level validation statistics

    Implementation
    --------------
    1. For each edge_result:
       a. Extract x_u, x_v from samples
       b. Call validate_edge_forward
       c. Store result

    2. Aggregate metrics:
       mean_rmse = mean([r.rmse for r in results])
       mean_nse = mean([r.nse for r in results])

    3. Flag inconsistent edges (NSE < 0 or RMSE > threshold)

    4. Return summary
    """
    validation_results = {}
    rmse_values = []
    nse_values = []
    n_consistent = 0

    rmse_threshold = getattr(config, "rt_consistency_rmse_threshold", 1.0)
    nse_threshold = getattr(config, "rt_consistency_nse_threshold", 0.5)

    for edge_result in edge_results:
        # Skip if skipped in inverse model
        if edge_result.skipped_reason:
            continue

        # Get upstream and downstream samples
        u_sample = samples.get(edge_result.u)
        v_sample = samples.get(edge_result.v)

        if not u_sample or not v_sample:
            continue

        # Extract concentration vectors
        x_u = [u_sample.get(ion, 0.0) for ion in config.ion_order]
        x_v = [v_sample.get(ion, 0.0) for ion in config.ion_order]

        # Validate edge
        try:
            result = validate_edge_forward(
                edge_result=edge_result,
                x_u=x_u,
                x_v_observed=x_v,
                config=config,
                kinetic_params=kinetic_params,
            )

            validation_results[edge_result.edge_id] = result

            # Collect metrics
            if result.rmse < float("inf"):
                rmse_values.append(result.rmse)
            if result.nse > float("-inf"):
                nse_values.append(result.nse)

            # Check consistency
            if (
                result.rmse < rmse_threshold
                and result.nse > nse_threshold
                and result.thermodynamic_consistent
            ):
                n_consistent += 1

        except Exception as e:
            # Skip failed validations
            print(f"Warning: Validation failed for edge {edge_result.edge_id}: {e}")
            continue

    # Compute summary statistics
    mean_rmse = sum(rmse_values) / len(rmse_values) if rmse_values else float("inf")
    mean_nse = sum(nse_values) / len(nse_values) if nse_values else float("-inf")

    # Identify inconsistent edges
    inconsistent = [
        edge_id
        for edge_id, result in validation_results.items()
        if result.rmse >= rmse_threshold or result.nse <= nse_threshold or not result.thermodynamic_consistent
    ]

    summary = ValidationSummary(
        n_edges_validated=len(validation_results),
        n_edges_consistent=n_consistent,
        mean_rmse=mean_rmse,
        mean_nse=mean_nse,
        inconsistent_edges=inconsistent,
        edge_results=validation_results,
    )

    return summary
