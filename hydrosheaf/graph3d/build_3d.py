"""
3D network construction and edge inference.
"""

import math
from typing import Dict, List, Optional

from .distance import classify_edge_type, compute_3d_distance, compute_screen_overlap
from .layers import assign_layers_to_nodes, get_aquitard_probability
from .types_3d import Edge3D, LayeredAquiferSystem, Network3D, Node3D


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


def infer_edges_3d_probabilistic(
    nodes: List[Node3D],
    config: "Config",  # type: ignore
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
    edges = []

    # Get config parameters
    radius_km = getattr(config, "edge_radius_km", 5.0)
    p_min = getattr(config, "edge_p_min", 0.75)
    max_neighbors = getattr(config, "edge_max_neighbors", 3)
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
                # No head data: assume bidirectional possibility
                p_head = 0.5

            # Distance decay
            p_dist = compute_distance_probability(d_3d, radius_km)

            # Layer probability
            if layer_system and node_i.aquifer_layer and node_j.aquifer_layer:
                p_layer = get_aquitard_probability(
                    layer_system,
                    node_i.aquifer_layer,
                    node_j.aquifer_layer,
                )
            else:
                p_layer = 1.0

            # Screen overlap bonus
            overlap, overlap_frac = compute_screen_overlap(node_i, node_j)
            p_screen = 0.5 + 0.5 * overlap_frac  # Range [0.5, 1.0]

            # Combined probability
            p_combined = p_head * p_dist * p_layer * p_screen

            # Apply threshold
            if p_combined < p_min:
                continue

            # Determine layer info
            same_layer = (node_i.aquifer_layer == node_j.aquifer_layer) if node_i.aquifer_layer and node_j.aquifer_layer else True

            # Classify edge type
            edge_type = classify_edge_type(d_xy, d_z, same_layer)

            # Compute gradients
            if d_xy > 1e-6:
                h_grad = (node_i.hydraulic_head - node_j.hydraulic_head) / d_xy if node_i.hydraulic_head and node_j.hydraulic_head else None
            else:
                h_grad = None

            if d_z > 1e-6:
                v_grad = (node_i.hydraulic_head - node_j.hydraulic_head) / d_z if node_i.hydraulic_head and node_j.hydraulic_head else None
            else:
                v_grad = None

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
                layer_from=node_i.aquifer_layer,
                layer_to=node_j.aquifer_layer,
            )

            # Store edge with probability for later filtering
            edges_per_node[node_i.node_id].append((p_combined, edge))

    # Keep top-k edges per node
    for node_id, edge_list in edges_per_node.items():
        # Sort by probability descending
        edge_list.sort(key=lambda x: x[0], reverse=True)

        # Keep top max_neighbors
        for _, edge in edge_list[:max_neighbors]:
            edges.append(edge)

    return edges


def build_network_3d(
    samples: List[Dict[str, object]],
    config: "Config",  # type: ignore
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
        sample_id = str(sample.get("sample_id", ""))
        if not sample_id:
            continue

        # Extract coordinates
        x = float(sample.get("x", 0.0))
        y = float(sample.get("y", 0.0))
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
        ion_order = getattr(config, "ion_order", ["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"])
        concentrations = []
        for ion in ion_order:
            conc = sample.get(ion, 0.0)
            concentrations.append(float(conc))

        # Create node
        node = Node3D(
            node_id=sample_id,
            x=x,
            y=y,
            z=z,
            elevation_m=elevation,
            screen_top=screen_top,
            screen_bottom=screen_bottom,
            hydraulic_head=head,
            head_uncertainty=0.5,
            aquifer_layer=aquifer_layer,
            concentrations=concentrations,
        )

        nodes_dict[sample_id] = node

    # Create layer system if provided
    layer_system = None
    if layer_definition:
        layer_system = create_layer_system_from_dict(layer_definition)

        # Assign layers to nodes
        depth_key_for_assign = "z"  # Use depth below surface
        assign_layers_to_nodes(nodes_dict, layer_system, depth_key_for_assign)

    # Infer edges
    nodes_list = list(nodes_dict.values())
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
