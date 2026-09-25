"""
3D network construction and edge inference.
"""

import csv
import json
import logging
import math
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from .distance import classify_edge_type, compute_3d_distance, compute_screen_overlap
from .layers import assign_layers_to_nodes, get_aquitard_probability
from .constraints import check_hydraulic_feasibility
from .types_3d import Edge3D, LayeredAquiferSystem, Network3D, Node3D
from ..graph.head_inference import (
    infer_heads_bayesian_linear,
    infer_heads_bayesian_mcmc,
)

logger = logging.getLogger(__name__)


def _coerce_geophysical_row(row: Mapping[str, Any]) -> Dict[str, Any]:
    """Normalize one JSON/CSV geophysical observation row."""
    normalized: Dict[str, Any] = {}
    for key, value in row.items():
        if value is None:
            normalized[str(key)] = None
            continue
        if isinstance(value, str):
            stripped = value.strip()
            if stripped.lower() in {"true", "yes", "y", "1"}:
                normalized[str(key)] = True
                continue
            if stripped.lower() in {"false", "no", "n", "0"}:
                normalized[str(key)] = False
                continue
            try:
                number = float(stripped)
            except ValueError:
                normalized[str(key)] = value
            else:
                normalized[str(key)] = number if math.isfinite(number) else None
        else:
            normalized[str(key)] = value
    return normalized


def load_geophysical_observations(path: str | Path) -> Dict[str, Dict[str, Any]]:
    """Load node-keyed ERT/TEM/sNMR observations from JSON or CSV.

    JSON may be either ``{"node_id": {...}}`` or
    ``{"observations": [{"node_id": "W1", ...}]}``. CSV requires a
    ``node_id`` column and may use ``well_id`` or ``sample_id`` as aliases.
    The adapter deliberately does not infer coordinates or units; those must
    be declared by the source columns and validated at the node boundary.
    """
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)

    rows: List[Mapping[str, Any]] = []
    if source.suffix.lower() == ".json":
        with source.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, Mapping) and isinstance(payload.get("observations"), list):
            rows = [row for row in payload["observations"] if isinstance(row, Mapping)]
        elif isinstance(payload, list):
            rows = [row for row in payload if isinstance(row, Mapping)]
        elif isinstance(payload, Mapping):
            rows = [
                dict(value, node_id=key)
                for key, value in payload.items()
                if isinstance(value, Mapping)
            ]
    elif source.suffix.lower() in {".csv", ".tsv"}:
        delimiter = "\t" if source.suffix.lower() == ".tsv" else ","
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter=delimiter))
    else:
        raise ValueError(
            f"Unsupported geophysical observation format {source.suffix!r}; use JSON or CSV."
        )

    observations: Dict[str, Dict[str, Any]] = {}
    for raw_row in rows:
        row = _coerce_geophysical_row(raw_row)
        node_id = (
            row.get("node_id")
            or row.get("well_id")
            or row.get("sample_id")
            or row.get("id")
        )
        if node_id is None:
            continue
        row["_geophysics_source"] = str(source)
        observations[str(node_id)] = row
    return observations


def attach_geophysical_observations(
    nodes: Sequence[Node3D],
    observations: Mapping[str, Mapping[str, Any]],
) -> List[Node3D]:
    """Return nodes with validated node-keyed geophysical attributes attached."""
    attached: List[Node3D] = []
    for node in nodes:
        attrs = dict(node.attrs or {})
        row = observations.get(str(node.node_id))
        if row is not None:
            attrs.update(dict(row))
        attached.append(replace(node, attrs=attrs))
    return attached


def attach_bedrock_elevation_raster(
    nodes: Sequence[Node3D],
    path: str | Path,
    config: Config,
) -> List[Node3D]:
    """Sample a GeoTIFF bedrock-elevation raster at node coordinates.

    Raster sampling is deliberately explicit about coordinate reference
    systems.  ``geophysics_node_crs`` declares the CRS of ``Node3D.x/y`` and
    the raster must carry its own CRS.  Without both declarations, a sampled
    elevation would be numerically plausible but geospatially untrustworthy.
    ``rasterio`` remains an optional dependency; enabling this adapter without
    it fails with an actionable error rather than silently dropping the layer.
    """
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    node_crs = getattr(config, "geophysics_node_crs", None)
    if not node_crs:
        raise ValueError(
            "geophysics_node_crs must be set when sampling a bedrock raster."
        )
    try:
        import rasterio
        from rasterio.warp import transform
    except ImportError as exc:  # pragma: no cover - depends on optional install
        raise RuntimeError(
            "GeoTIFF geophysics input requires the optional 'rasterio' package."
        ) from exc

    attached: List[Node3D] = []
    with rasterio.open(source) as dataset:
        if dataset.crs is None:
            raise ValueError(f"Bedrock raster has no CRS: {source}")
        xs = [float(node.x) for node in nodes]
        ys = [float(node.y) for node in nodes]
        transformed_x, transformed_y = transform(
            node_crs,
            dataset.crs,
            xs,
            ys,
        )
        sampled = dataset.sample(zip(transformed_x, transformed_y))
        for node, values in zip(nodes, sampled):
            attrs = dict(node.attrs or {})
            value = values[0] if len(values) else None
            if value is not None and dataset.nodata is not None:
                if float(value) == float(dataset.nodata):
                    value = None
            if value is not None and math.isfinite(float(value)):
                attrs["bedrock_elevation_m"] = float(value)
                attrs["_geophysics_bedrock_source"] = str(source)
            attached.append(replace(node, attrs=attrs))
    return attached


def attach_configured_geophysics(
    nodes: Sequence[Node3D],
    config: Config,
) -> List[Node3D]:
    """Load configured table observations before 3D edge inference.

    Configured files are strict inputs by default.  A missing or unsupported
    source therefore fails a geophysical run instead of allowing a hydraulic
    fallback to be labelled as geophysical evidence.  Set
    ``geophysics_strict_inputs=False`` only for an explicitly exploratory run.
    """
    paths = [str(path) for path in getattr(config, "geophysics_ert_profiles", [])]
    bedrock_path = getattr(config, "geophysics_bedrock_elevation_raster", None)
    attached_nodes = list(nodes)
    if bedrock_path:
        bedrock_suffix = Path(str(bedrock_path)).suffix.lower()
        if bedrock_suffix in {".tif", ".tiff"}:
            attached_nodes = attach_bedrock_elevation_raster(
                attached_nodes,
                str(bedrock_path),
                config,
            )
        else:
            paths.append(str(bedrock_path))
    if not paths:
        return attached_nodes

    merged: Dict[str, Dict[str, Any]] = {}
    for path in paths:
        try:
            loaded = load_geophysical_observations(path)
        except FileNotFoundError:
            if getattr(config, "geophysics_strict_inputs", True):
                raise
            logger.warning("Configured geophysical source does not exist: %s", path)
            continue
        for node_id, row in loaded.items():
            merged.setdefault(node_id, {}).update(row)
    return attach_geophysical_observations(attached_nodes, merged)


def _node_depth_below_surface(node: Any, config: Any) -> Optional[float]:
    """Return node depth using the declared positive-down coordinate contract."""
    z = getattr(node, "z", None)
    if z is None:
        return None
    if bool(getattr(config, "geophysics_z_positive_down", True)):
        return float(z)
    elevation = getattr(node, "elevation_m", None)
    if elevation is None:
        return None
    return float(elevation) - float(z)


def _bedrock_depth_below_surface(node: Any, attrs: Mapping[str, Any]) -> Optional[float]:
    """Convert explicit bedrock depth/elevation to positive depth below surface."""
    if attrs.get("bedrock_depth_m") is not None:
        return float(attrs["bedrock_depth_m"])
    bedrock_elevation = attrs.get("bedrock_elevation_m")
    elevation = getattr(node, "elevation_m", None)
    if bedrock_elevation is None or elevation is None:
        return None
    return float(elevation) - float(bedrock_elevation)


def compute_head_probability(
    head_i: float,
    head_j: float,
    sigma_i: float = 0.5,
    sigma_j: float = 0.5,
) -> float:
    """
    Compute probability that flow goes from i to j based on heads.

    Uses normal CDF to model probability of downhill flow.

    Parameters
    ----------
    head_i, head_j : float
        Hydraulic heads at nodes i and j (m)
    sigma_i, sigma_j : float
        Head measurement uncertainties (m)

    Returns
    -------
    float
        Probability in [0, 1]

    Mathematical Implementation
    ---------------------------
    Δh = head_i - head_j
    σ_combined = sqrt(σ_i² + σ_j²)
    z_score = Δh / σ_combined
    P = Φ(z_score)  # Standard normal CDF

    Where:
    - Φ(0) = 0.5: heads equal (uncertain direction)
    - Φ(1) = 0.84: head_i > head_j by 1σ (likely downhill)
    - Φ(3) = 0.999: head_i > head_j by 3σ (very likely downhill)

    Example
    -------
    >>> compute_head_probability(105.0, 100.0, sigma_i=0.5, sigma_j=0.5)
    0.9999...  # 5m drop >> uncertainties, very likely
    """
    delta_h = head_i - head_j
    sigma_combined = math.sqrt(sigma_i**2 + sigma_j**2)

    if sigma_combined < 1e-6:
        # No uncertainty: binary decision
        return 1.0 if delta_h > 0 else 0.0

    # Z-score
    z_score = delta_h / sigma_combined

    # Standard normal CDF approximation
    # Φ(z) ≈ 0.5 * (1 + erf(z/√2))
    # Using a simpler approximation for speed
    if z_score > 6:
        return 1.0
    elif z_score < -6:
        return 0.0
    else:
        # Approximation using tanh
        # Φ(z) ≈ 0.5 * (1 + tanh(0.7978 * z))
        return 0.5 * (1.0 + math.tanh(0.7978 * z_score))


def compute_distance_probability(
    distance_3d: float,
    radius_km: float = 5.0,
) -> float:
    """
    Compute distance decay probability.

    Gaussian decay based on characteristic influence radius.

    Parameters
    ----------
    distance_3d : float
        3D distance between nodes (m)
    radius_km : float
        Characteristic influence radius (km)

    Returns
    -------
    float
        Probability in [0, 1]

    Mathematical Implementation
    ---------------------------
    r = radius_km * 1000  # Convert to meters
    P = exp(-distance_3d² / (2 * r²))

    This is a Gaussian decay:
    - At distance = 0: P = 1.0
    - At distance = r: P = 0.61
    - At distance = 2r: P = 0.14
    - At distance = 3r: P = 0.01

    Example
    -------
    >>> compute_distance_probability(1000, radius_km=5.0)
    0.98  # 1km << 5km, high probability

    >>> compute_distance_probability(10000, radius_km=5.0)
    0.14  # 10km = 2 × 5km, moderate probability
    """
    radius_m = radius_km * 1000.0

    if radius_m < 1e-6:
        # No decay
        return 1.0

    # Gaussian decay
    exponent = -(distance_3d**2) / (2 * radius_m**2)
    return math.exp(exponent)




def geophysical_barrier_check(
    node_i: Any,
    node_j: Any,
    config: Any,
    barriers: Optional[Sequence[Any]] = None,
) -> Tuple[bool, float]:
    """
    Evaluate continuous geophysical barrier / conduit feasibility between two nodes.

    Returns (is_feasible, p_geophys).
    If a resistive dolerite dyke, impermeable fault gouge, or bedrock high cuts across
    the line segment (u, v), applies an exponential barrier penalty:
        p_geophys = exp(- delta_rho_barrier / sigma_struct)
    """
    if not getattr(config, "geophysics_enabled", False):
        return True, 1.0

    p_geophys = 1.0

    # 1. Check explicit structural barrier list
    if barriers is not None:
        for b in barriers:
            if hasattr(b, "crosses") and b.crosses(node_i, node_j):
                penalty = float(getattr(config, "geophysics_barrier_dyke_penalty", 1e-4))
                return False, penalty

    # 2. Check node-level geophysical barrier attributes
    attrs_i = getattr(node_i, "attrs", {}) or {}
    attrs_j = getattr(node_j, "attrs", {}) or {}
    barrier_detected = (
        attrs_i.get("dyke_crosses_to", {}).get(getattr(node_j, "node_id", ""), False)
        or attrs_j.get("dyke_crosses_to", {}).get(getattr(node_i, "node_id", ""), False)
        or attrs_i.get("geophysical_barrier", False)
        or attrs_j.get("geophysical_barrier", False)
    )
    if barrier_detected:
        penalty = float(getattr(config, "geophysics_barrier_dyke_penalty", 1e-4))
        p_geophys *= penalty

    # 3. Check bedrock interface in one explicit coordinate convention. Node3D.z
    # is positive depth below the local surface; a node is below the bedrock
    # interface when its depth is greater than the converted bedrock depth.
    depth_i = _node_depth_below_surface(node_i, config)
    depth_j = _node_depth_below_surface(node_j, config)
    bedrock_depth_i = _bedrock_depth_below_surface(node_i, attrs_i)
    bedrock_depth_j = _bedrock_depth_below_surface(node_j, attrs_j)
    if (
        depth_i is not None
        and bedrock_depth_i is not None
        and depth_i > bedrock_depth_i
    ):
        p_geophys *= 0.1
    if (
        depth_j is not None
        and bedrock_depth_j is not None
        and depth_j > bedrock_depth_j
    ):
        p_geophys *= 0.1

    sigma_struct = float(getattr(config, "geophysics_sigma_struct", 100.0))
    rho_i = float(attrs_i.get("apparent_resistivity_ohm_m", 50.0))
    rho_j = float(attrs_j.get("apparent_resistivity_ohm_m", 50.0))
    delta_rho = abs(rho_i - rho_j)
    if delta_rho > 2.0 * sigma_struct:
        p_geophys *= math.exp(-delta_rho / (3.0 * sigma_struct))

    p_geophys = min(1.0, max(1e-6, p_geophys))
    is_feasible = p_geophys >= float(getattr(config, "edge_p_min", 0.75))
    return is_feasible, p_geophys


def compute_geophysical_transit_time(
    node_i: Any,
    node_j: Any,
    d_3d: float,
    delta_h: Optional[float],
    config: Any,
) -> Tuple[Optional[float], Optional[float]]:
    """
    Calculate geophysical hydraulic conductivity K and advective travel time tau_geophys.

    Governing equations:
        K_sNMR = C_SDR * phi_sNMR^4 * (T2_star)^2
        v_uv = (K / phi_eff) * (delta_h / d_3d)
        tau_uv = d_3d / v_uv
    """
    if (
        not getattr(config, "geophysics_enabled", False)
        or not getattr(config, "geophysics_snmr_k_enabled", False)
    ):
        return None, None
    if delta_h is None or delta_h <= 0.0 or d_3d <= 1e-6:
        return None, None

    attrs_i = getattr(node_i, "attrs", {}) or {}
    phi_sNMR = attrs_i.get("phi_sNMR")
    if phi_sNMR is None:
        phi_sNMR = attrs_i.get("effective_porosity")
    t2_star = attrs_i.get("t2_star_ms")
    c_sdr = float(getattr(config, "geophysics_snmr_csdr", 1.0e-9))

    if t2_star is not None and phi_sNMR is not None:
        k_snmr_m_s = c_sdr * (float(phi_sNMR) ** 4) * (float(t2_star) ** 2)
        k_m_day = k_snmr_m_s * 86400.0
    else:
        # A direct geophysical K observation is allowed, but a hydraulic
        # fallback is not labelled as geophysical transit-time evidence.
        direct_k = attrs_i.get("k_geophys_m_day")
        if direct_k is None:
            return None, None
        k_m_day = float(direct_k)

    # Advective velocity still requires an explicit porosity.  Returning no
    # transit time is preferable to silently substituting a hydraulic default
    # or dividing by an absent sNMR observation.
    if phi_sNMR is None:
        return k_m_day, None
    phi_eff = float(phi_sNMR)
    gradient = float(delta_h) / d_3d
    velocity_m_day = (k_m_day / max(phi_eff, 0.01)) * gradient
    if velocity_m_day <= 1e-9:
        return k_m_day, None

    travel_time_days = d_3d / velocity_m_day
    travel_time_years = travel_time_days / 365.25
    return k_m_day, travel_time_years


def infer_edges_3d_probabilistic(
    nodes: List[Node3D],
    config: Config,
    layer_system: Optional[LayeredAquiferSystem] = None,
    use_haversine: bool = True,
) -> List[Edge3D]:
    """
    Infer probable flow edges in 3D aquifer network.

    Combines hydraulic gradient, distance decay, layer connectivity,
    and screen overlap to compute edge probabilities.

    Parameters
    ----------
    nodes : List[Node3D]
        All nodes with 3D coordinates and hydraulic heads
    config : Config
        Configuration with edge inference parameters:
        - edge_radius_km: maximum search radius
        - edge_p_min: minimum probability threshold
        - edge_max_neighbors: maximum edges per node
        - vertical_anisotropy: anisotropy factor α_v
    layer_system : Optional[LayeredAquiferSystem]
        Multi-layer system definition
    use_haversine : bool
        Use Haversine for geographic coordinates (default: True)

    Returns
    -------
    List[Edge3D]
        Inferred edges with probabilities

    Mathematical Implementation
    ---------------------------
    For each pair (i, j) where i ≠ j:

        1. Compute distances:
           d_xy, d_z, d_3d = compute_3d_distance(i, j, α_v, use_haversine)

        2. Skip if too far:
           if d_3d > config.edge_radius_km * 1000:
               continue

        3. Compute head probability:
           if head data available:
               P_head = compute_head_probability(h_i, h_j, σ_i, σ_j)
               if P_head < config.edge_p_min:
                   continue  # Wrong direction

        4. Compute distance decay:
           P_dist = exp(-d_3d² / (2r²))

        5. Compute layer probability:
           if layer_system:
               P_layer = get_aquitard_probability(system, L_i, L_j)
           else:
               P_layer = 1.0

        6. Compute screen overlap bonus:
           overlap, frac = compute_screen_overlap(i, j)
           P_screen = 0.5 + 0.5 * frac  # Range [0.5, 1.0]

        7. Combined probability:
           P = P_head × P_dist × P_layer × P_screen

        8. Apply threshold and create edge:
           if P >= config.edge_p_min:
               edge = Edge3D(...)
               edges.append(edge)

    9. Keep top-k neighbors per node:
       For each node, keep edges with highest probabilities
       (up to config.edge_max_neighbors)

    Return edges

    Example
    -------
    >>> nodes = [
    ...     Node3D("W1", x=-120.5, y=38.5, z=50, elevation_m=100, hydraulic_head=105),
    ...     Node3D("W2", x=-120.48, y=38.52, z=55, elevation_m=98, hydraulic_head=102),
    ... ]
    >>> edges = infer_edges_3d_probabilistic(nodes, config, layer_system=None)
    >>> len(edges)
    1  # W1 -> W2 (downhill)
    """
    if getattr(config, "geophysics_enabled", False):
        nodes = attach_configured_geophysics(nodes, config)

    edges = []

    # Get config parameters
    radius_km = getattr(config, "edge_radius_km", 5.0)
    p_min = getattr(config, "edge_p_min", 0.75)
    max_neighbors = getattr(config, "edge_max_neighbors", 3)
    
    # Stratigraphic Priority Quotas
    # Default to splitting max_neighbors if specific quotas not set
    max_primary = getattr(config, "edge_max_neighbors_primary", max_neighbors)
    max_secondary = getattr(config, "edge_max_neighbors_secondary", max(1, max_neighbors // 3))
    
    anisotropy = getattr(config, "vertical_anisotropy", 0.1)

    # Track edges per node for max_neighbors constraint
    edges_per_node: Dict[str, List[tuple]] = {node.node_id: [] for node in nodes}

    # Pairwise evaluation
    for i, node_i in enumerate(nodes):
        for j, node_j in enumerate(nodes):
            if i == j:
                continue

            # Compute 3D distance
            d_xy, d_z, d_3d = compute_3d_distance(
                node_i,
                node_j,
                anisotropy_factor=anisotropy,
                use_haversine=use_haversine,
            )

            # Skip if too far
            if d_3d > radius_km * 1000.0:
                continue

            # Head probability
            if node_i.hydraulic_head is not None and node_j.hydraulic_head is not None:
                p_head = compute_head_probability(
                    node_i.hydraulic_head,
                    node_j.hydraulic_head,
                    node_i.head_uncertainty,
                    node_j.head_uncertainty,
                )

                # Skip wrong direction
                if p_head < p_min:
                    continue
            else:
                # No head data: use topographic Bayesian prior (if possible).
                is_ok, p_head, _ = check_hydraulic_feasibility(
                    node_i,
                    node_j,
                    min_gradient=getattr(config, "edge_gradient_min", 1e-4),
                    topo_sigma_depth=getattr(config, "edge_topo_sigma_depth", 5.0),
                    topo_reject_p=getattr(config, "edge_topo_reject_p", 0.1),
                )
                if not is_ok:
                    continue

            # Distance decay
            p_dist = compute_distance_probability(d_3d, radius_km)

            # Layer probability
            if (
                layer_system
                and node_i.aquifer_layer is not None
                and node_j.aquifer_layer is not None
            ):
                p_layer = get_aquitard_probability(
                    layer_system,
                    node_i.aquifer_layer,
                    node_j.aquifer_layer,
                )
            else:
                p_layer = 1.0

            # Geophysical barrier check
            is_geophys_ok, p_geophys = geophysical_barrier_check(node_i, node_j, config)
            if not is_geophys_ok and getattr(config, "geophysics_barrier_enabled", False):
                continue

            # Screen overlap bonus
            overlap, overlap_frac = compute_screen_overlap(node_i, node_j)
            p_screen = 0.5 + 0.5 * overlap_frac  # Range [0.5, 1.0]

            # Combined probability
            p_combined = p_head * p_dist * p_layer * p_screen * p_geophys

            # Apply threshold
            if p_combined < p_min:
                continue

            # Determine layer info
            same_layer = (
                (node_i.aquifer_layer == node_j.aquifer_layer)
                if node_i.aquifer_layer is not None and node_j.aquifer_layer is not None
                else True # Assume same if unknown, maximum entropy
            )

            # Classify edge type
            edge_type = classify_edge_type(d_xy, d_z, same_layer)

            # Compute gradients
            delta_h = (
                float(node_i.hydraulic_head) - float(node_j.hydraulic_head)
                if node_i.hydraulic_head is not None and node_j.hydraulic_head is not None
                else None
            )
            if d_xy > 1e-6:
                h_grad = delta_h / d_xy if delta_h is not None else None
            else:
                h_grad = None

            if d_z > 1e-6:
                v_grad = delta_h / d_z if delta_h is not None else None
            else:
                v_grad = None

            # Geophysical transit time conditioning
            k_geophys, tau_geophys = compute_geophysical_transit_time(
                node_i, node_j, d_3d, delta_h, config
            )

            # Create edge
            edge = Edge3D(
                edge_id=f"{node_i.node_id}->{node_j.node_id}",
                u=node_i.node_id,
                v=node_j.node_id,
                horizontal_distance_m=d_xy,
                vertical_distance_m=d_z,
                distance_3d=d_3d,
                edge_type=edge_type,
                same_layer=same_layer,
                prob_head=p_head,
                prob_distance=p_dist,
                prob_layer=p_layer,
                prob_combined=p_combined,
                horizontal_gradient=h_grad,
                vertical_gradient=v_grad,
                delta_h=delta_h,
                layer_from=node_i.aquifer_layer,
                layer_to=node_j.aquifer_layer,
                prob_geophys=p_geophys,
                tau_geophys_years=tau_geophys,
                k_geophys_m_day=k_geophys,
            )

            # Store edge with probability for later filtering
            # Tuple: (Probability, Edge, IsPrimary)
            is_primary = same_layer
            edges_per_node[node_i.node_id].append((p_combined, edge, is_primary))

    # Keep top-k edges per node using Stratigraphic Priority
    for node_id, edge_list in edges_per_node.items():
        # Sort by probability descending
        edge_list.sort(key=lambda x: x[0], reverse=True)
        
        primary_candidates = [e for p, e, prim in edge_list if prim]
        secondary_candidates = [e for p, e, prim in edge_list if not prim]
        
        # 1. Fill Primary Quota (Same Layer)
        selected = primary_candidates[:max_primary]
        
        # 2. Fill Secondary Quota (Cross Layer / Vertical Leakage)
        selected.extend(secondary_candidates[:max_secondary])
        
        # 3. If total exceeds absolute hard limit, truncate by probability again?
        # Or trust the quotas. Let's trust the quotas but clamp total if needed.
        # usually max_neighbors >= max_primary + max_secondary
        
        for edge in selected:
            edges.append(edge)

    return edges


def build_network_3d(
    samples: List[Dict[str, object]],
    config: Config,
    layer_definition: Optional[Dict[str, object]] = None,
    use_haversine: bool = True,
) -> Network3D:
    """
    Build complete 3D network from sample data.

    Parameters
    ----------
    samples : List[Dict]
        Sample data with required keys:
        - sample_id: str
        - x, y: float (coordinates)
        - z or depth_key: float (depth/elevation)
        - elevation or elevation_key: float (surface elevation)
        - head or head_key: float (hydraulic head)
        - ion concentrations (Ca, Mg, etc.)
    config : Config
        Configuration with network parameters
    layer_definition : Optional[Dict]
        Layer system definition:
        {
            "n_layers": 3,
            "names": ["Shallow", "Intermediate", "Deep"],
            "tops": [0, 30, 100],
            "bottoms": [30, 100, 250],
            "aquitard_p": [0.3, 0.2],
            "anisotropy": [0.2, 0.1, 0.05]
        }
    use_haversine : bool
        Use Haversine distance for geographic coordinates

    Returns
    -------
    Network3D
        Complete 3D network with nodes, edges, and statistics

    Implementation Steps
    --------------------
    1. Convert samples to Node3D objects
    2. Create LayeredAquiferSystem if layer_definition provided
    3. Assign layers to nodes
    4. Infer edges using probabilistic algorithm
    5. Compute summary statistics
    6. Return Network3D

    Example
    -------
    >>> samples = [
    ...     {"sample_id": "W1", "x": -120.5, "y": 38.5, "z": 50, "elevation": 100, "head": 105, "Ca": 2.0, ...},
    ...     {"sample_id": "W2", "x": -120.48, "y": 38.52, "z": 55, "elevation": 98, "head": 102, "Ca": 1.8, ...},
    ... ]
    >>> layer_def = {"n_layers": 2, "names": ["Shallow", "Deep"], "tops": [0, 50], "bottoms": [50, 200], "aquitard_p": [0.3]}
    >>> network = build_network_3d(samples, config, layer_def)
    >>> print(f"{len(network.nodes)} nodes, {len(network.edges)} edges")
    """
    from .layers import create_layer_system_from_dict

    # Get column names from config
    z_key = getattr(config, "z_coordinate_key", "screen_depth")
    elevation_key = getattr(config, "edge_elevation_key", "elevation")
    head_key = getattr(config, "edge_head_key", "head_meas")
    screen_top_key = getattr(config, "screen_top_key", "screen_top")
    screen_bottom_key = getattr(config, "screen_bottom_key", "screen_bottom")
    layer_key = getattr(config, "layer_key", "aquifer_layer")

    # Convert samples to Node3D
    nodes_dict: Dict[str, Node3D] = {}

    for sample in samples:
        node_id = (
            sample.get("site_id") or sample.get("node_id") or sample.get("sample_id")
        )
        if node_id in (None, ""):
            continue
        node_id_str = str(node_id)

        # Extract coordinates
        x_raw = sample.get("x", sample.get("lon", 0.0))
        y_raw = sample.get("y", sample.get("lat", 0.0))
        x = float(x_raw)  # lon or easting
        y = float(y_raw)  # lat or northing
        z = float(sample.get(z_key, sample.get("z", 0.0)))
        elevation = float(sample.get(elevation_key, sample.get("elevation", 0.0)))

        # Hydraulic head
        head = sample.get(head_key, sample.get("head", None))
        if head is not None:
            head = float(head)

        # Screen interval
        screen_top = sample.get(screen_top_key, None)
        screen_bottom = sample.get(screen_bottom_key, None)
        if screen_top is not None:
            screen_top = float(screen_top)
        if screen_bottom is not None:
            screen_bottom = float(screen_bottom)

        # Layer
        aquifer_layer = sample.get(layer_key, None)
        if aquifer_layer is not None:
            aquifer_layer = int(aquifer_layer)

        # Extract concentrations
        ion_order = getattr(
            config,
            "ion_order",
            ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
        )
        concentrations = []
        for ion in ion_order:
            conc = sample.get(ion, 0.0)
            concentrations.append(float(conc))

        # Preserve explicitly supplied metadata (including ERT/TEM/sNMR
        # fields) at the Node3D boundary.  Dropping these columns here would
        # make the later geophysical adapter appear to succeed while the
        # sample-provided observations never reach barrier or transit logic.
        core_keys = {
            "site_id",
            "node_id",
            "sample_id",
            "x",
            "y",
            "lon",
            "lat",
            z_key,
            "z",
            elevation_key,
            "elevation",
            head_key,
            "head",
            screen_top_key,
            "screen_bottom_key",
            layer_key,
            *ion_order,
        }
        node_attrs = {
            str(key): value
            for key, value in sample.items()
            if str(key) not in core_keys
        }

        # Create node
        node = Node3D(
            node_id=node_id_str,
            x=x,
            y=y,
            z=z,
            elevation_m=elevation,
            screen_top=screen_top,
            screen_bottom=screen_bottom,
            hydraulic_head=head,
            head_uncertainty=(
                float(getattr(config, "edge_sigma_meas", 0.5))
                if head is not None
                else float(getattr(config, "edge_sigma_topo", 10.0))
            ),
            aquifer_layer=aquifer_layer,
            concentrations=concentrations,
            attrs=node_attrs,
        )

        nodes_dict[node_id_str] = node

    # Optional Bayesian hierarchical head estimation (fills missing heads using elevation/DTW priors).
    head_inference = getattr(config, "edge_head_inference", "heuristic")
    if head_inference in {"bayesian", "bayesian_mcmc"}:
        dtw_key = getattr(config, "edge_dtw_key", "dtw")
        if head_inference == "bayesian_mcmc":
            posterior = infer_heads_bayesian_mcmc(
                samples,
                node_id_key="sample_id",
                head_key=head_key,
                dtw_key=dtw_key,
                elevation_key=elevation_key,
                sigma_meas=float(getattr(config, "edge_sigma_meas", 0.5)),
                sigma_dtw=float(getattr(config, "edge_sigma_dtw", 1.0)),
                sigma_elev=float(getattr(config, "edge_sigma_elev", 1.0)),
                sigma_topo=float(getattr(config, "edge_sigma_topo", 10.0)),
                dtw_prior_mu=float(getattr(config, "edge_dtw_prior_mu", 5.0)),
                dtw_prior_sigma=float(getattr(config, "edge_dtw_prior_sigma", 5.0)),
                head_prior_mu=float(getattr(config, "edge_head_prior_mu", 0.0)),
                head_prior_sigma=float(
                    getattr(config, "edge_head_prior_sigma", 1000.0)
                ),
                mcmc_draws=int(getattr(config, "bayesian_n_samples", 1000)),
                mcmc_chains=int(getattr(config, "bayesian_n_chains", 2)),
                mcmc_target_accept=float(
                    getattr(config, "bayesian_target_accept", 0.9)
                ),
                mcmc_warmup_fraction=float(
                    getattr(config, "bayesian_warmup_fraction", 0.5)
                ),
            )
        else:
            posterior = infer_heads_bayesian_linear(
                samples,
                node_id_key="sample_id",
                head_key=head_key,
                dtw_key=dtw_key,
                elevation_key=elevation_key,
                sigma_meas=float(getattr(config, "edge_sigma_meas", 0.5)),
                sigma_dtw=float(getattr(config, "edge_sigma_dtw", 1.0)),
                sigma_elev=float(getattr(config, "edge_sigma_elev", 1.0)),
                sigma_topo=float(getattr(config, "edge_sigma_topo", 10.0)),
                dtw_prior_mu=float(getattr(config, "edge_dtw_prior_mu", 5.0)),
                dtw_prior_sigma=float(getattr(config, "edge_dtw_prior_sigma", 5.0)),
                head_prior_mu=float(getattr(config, "edge_head_prior_mu", 0.0)),
                head_prior_sigma=float(
                    getattr(config, "edge_head_prior_sigma", 1000.0)
                ),
            )
        idx_map = posterior.index()
        for node_id, node in nodes_dict.items():
            idx = idx_map.get(node_id)
            if idx is None:
                continue
            node.hydraulic_head = float(posterior.head_mean[idx])
            node.head_uncertainty = float(
                math.sqrt(float(posterior.head_cov[idx, idx]))
            )

    # Create layer system if provided
    layer_system = None
    if layer_definition:
        layer_system = create_layer_system_from_dict(layer_definition)

        # Assign layers to nodes
        depth_key_for_assign = "z"  # Use depth below surface
        assign_layers_to_nodes(nodes_dict, layer_system, depth_key_for_assign)

    # Infer edges
    nodes_list = list(nodes_dict.values())
    if getattr(config, "geophysics_enabled", False):
        attached_nodes = attach_configured_geophysics(nodes_list, config)
        nodes_dict = {node.node_id: node for node in attached_nodes}
        nodes_list = attached_nodes
    edges = infer_edges_3d_probabilistic(
        nodes_list,
        config,
        layer_system,
        use_haversine,
    )

    # Create network
    network = Network3D(
        nodes=nodes_dict,
        edges=edges,
        layer_system=layer_system,
    )

    # Summary statistics computed in __post_init__

    return network
