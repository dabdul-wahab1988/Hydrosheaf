"""Data validation and requirement checks.

The public pipeline has two different kinds of optional configuration:

* an explicit ``True`` or ``False`` is a user override; and
* ``None`` means *automatic* for capability-gated modules.

Keeping those states distinct lets the default configuration turn optional
components on when their evidence is available, while avoiding a late solver
failure when the corresponding observations are absent.  The resolver below
is intentionally conservative: it only decides whether a module has the
minimum information needed to run; it never fabricates a value or imputes a
missing observation.
"""
from dataclasses import replace
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from .parsing import sample_list
from .schema import parse_numeric


def _any_numeric(
    samples: Sequence[Mapping[str, object]],
    key: str,
    detection_policy: str,
) -> bool:
    """Check if any sample has a valid numeric value for the given key."""
    for sample in samples:
        if parse_numeric(sample.get(key), detection_policy) is not None:
            return True
    return False


def _any_pair_numeric(
    samples: Sequence[Mapping[str, object]],
    key_a: str,
    key_b: str,
    detection_policy: str,
) -> bool:
    """Check if any sample has valid numeric values for both keys."""
    for sample in samples:
        if (
            parse_numeric(sample.get(key_a), detection_policy) is not None
            and parse_numeric(sample.get(key_b), detection_policy) is not None
        ):
            return True
    return False


def _missing_keys(
    samples: Sequence[Mapping[str, object]],
    keys: Sequence[str],
    detection_policy: str,
) -> List[str]:
    """Identify sample IDs missing required numeric keys."""
    missing_ids: List[str] = []
    for sample in samples:
        missing = False
        for key in keys:
            if parse_numeric(sample.get(key), detection_policy) is None:
                missing = True
                break
        if missing:
            sample_id = sample.get("site_id") or sample.get("sample_id") or "unknown"
            missing_ids.append(str(sample_id))
    return missing_ids


def _any_numeric_alias(
    samples: Sequence[Mapping[str, object]],
    keys: Sequence[str],
    detection_policy: str,
) -> bool:
    """Return whether any row contains a finite numeric value under aliases."""
    return any(
        parse_numeric(sample.get(key), detection_policy) is not None
        for sample in samples
        for key in keys
    )


def _any_all_numeric(
    samples: Sequence[Mapping[str, object]],
    keys: Sequence[str],
    detection_policy: str,
) -> bool:
    """Return whether one row contains finite values for every requested key."""
    return any(
        all(parse_numeric(sample.get(key), detection_policy) is not None for key in keys)
        for sample in samples
    )


def _sample_id_map(
    samples: Sequence[Mapping[str, object]],
) -> Dict[str, Mapping[str, object]]:
    """Build a stable sample lookup for capability checks."""
    result: Dict[str, Mapping[str, object]] = {}
    for sample in samples:
        node_id = sample.get("site_id") or sample.get("sample_id")
        if node_id is not None:
            result[str(node_id)] = sample
    return result


def _has_complete_core_edge_panel(
    samples: Sequence[Mapping[str, object]],
    config: Config,
    candidate_edges: Optional[Sequence[object]],
) -> bool:
    """Return whether at least one candidate edge has a complete core panel."""
    if not candidate_edges or not config.ion_order:
        return False
    sample_map = _sample_id_map(samples)
    for edge in candidate_edges:
        upstream = sample_map.get(str(getattr(edge, "u", "")))
        downstream = sample_map.get(str(getattr(edge, "v", "")))
        if upstream is None or downstream is None:
            continue
        if all(
            parse_numeric(upstream.get(ion), config.detection_limit_policy) is not None
            and parse_numeric(downstream.get(ion), config.detection_limit_policy)
            is not None
            for ion in config.ion_order
        ):
            return True
    return False


def _has_complete_edge_measurement(
    samples: Sequence[Mapping[str, object]],
    candidate_edges: Optional[Sequence[object]],
    keys: Sequence[str],
    detection_policy: str,
) -> bool:
    """Return whether one candidate edge has all requested endpoint fields."""
    if not candidate_edges or not keys:
        return False
    sample_map = _sample_id_map(samples)
    for edge in candidate_edges:
        upstream = sample_map.get(str(getattr(edge, "u", "")))
        downstream = sample_map.get(str(getattr(edge, "v", "")))
        if upstream is None or downstream is None:
            continue
        if all(
            parse_numeric(upstream.get(key), detection_policy) is not None
            and parse_numeric(downstream.get(key), detection_policy) is not None
            for key in keys
        ):
            return True
    return False


def _has_topology_signal(
    samples: Sequence[Mapping[str, object]],
    config: Config,
    candidate_edges: Optional[Sequence[object]],
) -> bool:
    """Check for evidence that can support a topology posterior.

    A topology posterior needs more than node identifiers. It can use an
    explicit edge-inclusion prior (normally written by probabilistic edge
    generation) or a complete chemical panel whose global sheaf residual
    supplies the graph cost. Head/elevation and geometry alone activate the
    hydraulic diagnostics, but do not give a manually supplied edge a
    non-neutral topology prior; in that case MCMC would only resample an
    arbitrary candidate list and is therefore auto-disabled.
    """
    if not candidate_edges:
        return False

    detection_policy = config.detection_limit_policy
    edge_signal_keys = (
        "edge_confidence",
        "prior_edge_probability",
        "p_uv",
    )
    for edge in candidate_edges:
        attrs = getattr(edge, "attrs", {}) or {}
        if any(
            parse_numeric(attrs.get(key), detection_policy) is not None
            for key in edge_signal_keys
        ):
            return True

    return False


def _has_hydraulic_signal(
    samples: Sequence[Mapping[str, object]],
    config: Config,
) -> bool:
    """Return whether at least two nodes have a head/elevation observation."""
    detection_policy = config.detection_limit_policy
    keys = tuple(
        dict.fromkeys(
            (
                getattr(config, "hydraulic_hodge_head_key", "hydraulic_head"),
                "hydraulic_head",
                "head",
                "head_meas",
                "water_level",
                getattr(config, "edge_elevation_key", "elevation"),
                "elevation",
            )
        )
    )
    count = sum(
        _any_numeric_alias([sample], keys, detection_policy) for sample in samples
    )
    return count >= 2


def _has_geophysics_signal(
    samples: Sequence[Mapping[str, object]],
    config: Config,
    candidate_edges: Optional[Sequence[object]],
) -> bool:
    """Return whether configured or observed geophysical data are available."""
    if getattr(config, "geophysics_bedrock_elevation_raster", None):
        return True
    if getattr(config, "geophysics_ert_profiles", None):
        return True

    sample_keys = (
        "prob_geophys",
        "geophysical_barrier",
        "barrier_detected",
        "apparent_resistivity",
        "apparent_resistivity_ohm_m",
        "resistivity",
        "conductivity",
        "formation_conductivity",
        "ert_conductivity",
        "phi_sNMR",
        "effective_porosity",
        "t2_star_ms",
        "k_geophys_m_day",
        "k_geophys_std_m_day",
    )
    policy = config.detection_limit_policy
    for sample in samples:
        for key in sample_keys:
            value = sample.get(key)
            if isinstance(value, bool) and value:
                return True
            if parse_numeric(value, policy) is not None:
                return True

    for edge in candidate_edges or ():
        attrs = getattr(edge, "attrs", {}) or {}
        for key in sample_keys:
            value = attrs.get(key)
            if isinstance(value, bool) and value:
                return True
            if parse_numeric(value, policy) is not None:
                return True
    return False


def _has_geophysics_subsignal(
    samples: Sequence[Mapping[str, object]],
    config: Config,
    candidate_edges: Optional[Sequence[object]],
    kind: str,
) -> bool:
    """Return whether a specific geophysical submodule has its inputs."""
    policy = config.detection_limit_policy
    if kind == "barrier":
        keys = ("geophysical_barrier", "barrier_detected", "dyke_crosses_to")
    elif kind == "conductivity":
        keys = (
            "conductivity",
            "formation_conductivity",
            "ert_conductivity",
            "apparent_resistivity",
            "apparent_resistivity_ohm_m",
        )
    elif kind == "snmr":
        keys = ("phi_sNMR", "effective_porosity", "t2_star_ms", "k_geophys_m_day")
    else:
        return False

    for sample in samples:
        for key in keys:
            value = sample.get(key)
            if isinstance(value, bool) and value:
                return True
            if parse_numeric(value, policy) is not None:
                return True
    for edge in candidate_edges or ():
        attrs = getattr(edge, "attrs", {}) or {}
        for key in keys:
            value = attrs.get(key)
            if isinstance(value, bool) and value:
                return True
            if parse_numeric(value, policy) is not None:
                return True
    return False


def _resolve_auto_flag(
    config: Config,
    updates: Dict[str, object],
    decisions: Dict[str, Dict[str, object]],
    flag: str,
    available: bool,
    reason: str,
) -> bool:
    """Resolve a tri-state optional flag and record an auditable decision."""
    raw = getattr(config, flag, None)
    if raw is None:
        enabled = bool(available)
        updates[flag] = enabled
        decisions[flag] = {
            "status": "auto_enabled" if enabled else "auto_disabled",
            "configured": "auto",
            "enabled": enabled,
            "reason": "required inputs available" if enabled else reason,
        }
        return enabled

    enabled = bool(raw)
    decisions[flag] = {
        "status": "enabled" if enabled else "disabled",
        "configured": enabled,
        "enabled": enabled,
        "reason": (
            "explicit override; required inputs were not checked"
            if enabled and not available
            else "explicitly disabled"
            if not enabled
            else "explicitly enabled"
        ),
    }
    return enabled


def _deferred_auto_decision(reason: str) -> Dict[str, object]:
    """Record an automatic decision that needs a candidate-edge universe."""
    return {
        "status": "deferred",
        "configured": "auto",
        "enabled": None,
        "reason": reason,
    }


def resolve_optional_modules(
    samples: object,
    config: Config,
    *,
    candidate_edges: Optional[Sequence[object]] = None,
) -> Tuple[Config, Dict[str, Dict[str, object]]]:
    """Resolve automatic optional-module defaults against available inputs.

    The automatic defaults cover the optional sheaf stack: isotope and
    chloride evidence terms, age evidence, Bayesian topology, cohomology,
    hydraulic Hodge diagnostics, and the geophysical master/submodule
    switches.  An explicit boolean remains an override, which is useful for
    controlled tests and deliberately strict analyses.

    Returns
    -------
    (Config, decisions)
        A copied configuration with automatic flags resolved to booleans and
        a machine-readable decision record suitable for pipeline diagnostics.
    """
    s_list = sample_list(samples)
    policy = config.detection_limit_policy
    updates: Dict[str, object] = {}
    decisions: Dict[str, Dict[str, object]] = {}

    age_keys = [
        "mean_age_years",
        "age_years",
        "age",
        "3H",
        "tritium",
        "H3",
        "H3_TU",
        "tritium_TU",
        "Tritium",
        str(getattr(config, "residence_time_tracer", "3H")),
    ]
    age_available = _any_numeric_alias(s_list, tuple(dict.fromkeys(age_keys)), policy)
    _resolve_auto_flag(
        config,
        updates,
        decisions,
        "sheaf_age_enabled",
        age_available,
        "no numeric age or configured residence-time tracer observation",
    )

    topology_signal = _has_topology_signal(s_list, config, candidate_edges)
    core_edge_panel = _has_complete_core_edge_panel(s_list, config, candidate_edges)
    if candidate_edges is None:
        # Do not turn an unresolved ``None`` into an explicit False before a
        # caller has built its candidate universe.  This is important for
        # infer_edges(..., method="probabilistic_sheaf"), where edge-level
        # geometry is the evidence that can activate topology.
        decisions["topology_posterior_enabled"] = _deferred_auto_decision(
            "candidate edges are required to evaluate topology evidence"
        )
        decisions["sheaf_cohomology_enabled"] = _deferred_auto_decision(
            "candidate edges are required to evaluate endpoint panel completeness"
        )
    else:
        topology_available = bool(candidate_edges) and (
            topology_signal or core_edge_panel
        )
        _resolve_auto_flag(
            config,
            updates,
            decisions,
            "topology_posterior_enabled",
            topology_available,
            "no candidate edges with head, geometry, edge-prior, or complete chemical evidence",
        )

        cohomology_available = bool(candidate_edges) and core_edge_panel
        _resolve_auto_flag(
            config,
            updates,
            decisions,
            "sheaf_cohomology_enabled",
            cohomology_available,
            "no candidate edge has a complete configured core chemical panel at both endpoints",
        )

        isotope_available = _has_complete_edge_measurement(
            s_list,
            candidate_edges,
            (config.isotope_d18o_key, config.isotope_d2h_key),
            policy,
        )
        _resolve_auto_flag(
            config,
            updates,
            decisions,
            "sheaf_isotope_enabled",
            isotope_available,
            "no candidate edge has both isotope endpoints",
        )
        cl_available = _has_complete_edge_measurement(
            s_list,
            candidate_edges,
            ("Cl",),
            policy,
        )
        _resolve_auto_flag(
            config,
            updates,
            decisions,
            "sheaf_cl_enabled",
            cl_available,
            "no candidate edge has chloride at both endpoints",
        )

    hydraulic_available = _has_hydraulic_signal(s_list, config)
    _resolve_auto_flag(
        config,
        updates,
        decisions,
        "hydraulic_hodge_enabled",
        hydraulic_available,
        "fewer than two nodes have hydraulic head or elevation evidence",
    )

    geophysics_available = _has_geophysics_signal(s_list, config, candidate_edges)
    geophysics_enabled = _resolve_auto_flag(
        config,
        updates,
        decisions,
        "geophysics_enabled",
        geophysics_available,
        "no geophysical profile, raster, node field, or edge field is available",
    )
    submodule_specs = (
        ("geophysics_barrier_enabled", "barrier", "no barrier/dyke geophysical field is available"),
        (
            "geophysics_conductivity_conditioning",
            "conductivity",
            "no conductivity/resistivity field is available",
        ),
        (
            "geophysics_ec_crossval_enabled",
            "conductivity",
            "no conductivity/resistivity field is available",
        ),
        (
            "geophysics_fluid_ec_enabled",
            "conductivity",
            "no endpoint conductivity field is available",
        ),
        (
            "geophysics_snmr_k_enabled",
            "snmr",
            "no sNMR porosity, T2*, or direct geophysical K field is available",
        ),
    )
    for flag, kind, reason in submodule_specs:
        available = geophysics_enabled and _has_geophysics_subsignal(
            s_list, config, candidate_edges, kind
        )
        _resolve_auto_flag(config, updates, decisions, flag, available, reason)

    sheaf_available = bool(candidate_edges) and (
        topology_signal
        or core_edge_panel
        or _any_numeric_alias(s_list, ("Cl",), policy)
        or _any_numeric_alias(s_list, tuple(dict.fromkeys(age_keys)), policy)
        or _any_pair_numeric(
            s_list,
            config.isotope_d18o_key,
            config.isotope_d2h_key,
            policy,
        )
    )
    decisions["sheaf_refinement"] = {
        "status": "auto_enabled" if sheaf_available else "auto_disabled",
        "configured": "auto",
        "enabled": sheaf_available,
        "reason": (
            "candidate edges and at least one sheaf evidence channel are available"
            if sheaf_available
            else "no candidate edges or sheaf evidence channel is available"
        ),
    }

    if not updates:
        resolved = config
    else:
        resolved = replace(config, **updates)
    return resolved, decisions


def validate_required_inputs(samples: object, config: Config) -> None:
    """Raise ValueError if required inputs for enabled modules are missing."""
    s_list = sample_list(samples)
    detection_policy = config.detection_limit_policy
    missing_reports: List[str] = []

    if config.phreeqc_enabled:
        required_phreeqc = [
            "pH",
            "Ca",
            "Mg",
            "Na",
            "K",
            "Cl",
            "SO4",
            "NO3",
            "F",
            "HCO3",
        ]
        missing = _missing_keys(s_list, required_phreeqc, detection_policy)
        if missing:
            missing_reports.append(
                "PHREEQC requires pH and major ions (Ca, Mg, Na, K, Cl, SO4, NO3, F, HCO3) "
                f"for all samples (missing: {missing})"
            )

    if config.isotope_enabled and config.lmwl_defined:
        missing = _missing_keys(
            s_list,
            [config.isotope_d18o_key, config.isotope_d2h_key],
            detection_policy,
        )
        if missing:
            missing_reports.append(
                "Isotope penalties require both "
                f"{config.isotope_d18o_key} and {config.isotope_d2h_key} "
                f"for all samples (missing: {missing})"
            )

    if config.nitrate_source_enabled:
        missing = _missing_keys(s_list, ["NO3"], detection_policy)
        if missing:
            missing_reports.append(
                f"Nitrate source requires NO3 for all samples (missing: {missing})"
            )

    if missing_reports:
        raise ValueError("; ".join(missing_reports))


def auto_disable_missing_modules(samples: object, config: Config) -> Config:
    """Disable feature flags when required inputs are missing across samples."""
    s_list = sample_list(samples)
    detection_policy = config.detection_limit_policy
    updates: Dict[str, object] = {}

    if config.phreeqc_enabled and not _any_numeric(s_list, "pH", detection_policy):
        updates["phreeqc_enabled"] = False

    if config.isotope_enabled and not _any_pair_numeric(
        s_list,
        config.isotope_d18o_key,
        config.isotope_d2h_key,
        detection_policy,
    ):
        updates["isotope_enabled"] = False

    if config.nitrate_source_enabled and not _any_numeric(
        s_list, "NO3", detection_policy
    ):
        updates["nitrate_source_enabled"] = False

    if updates:
        return replace(config, **updates)
    return config
