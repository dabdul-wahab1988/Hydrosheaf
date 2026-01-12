"""
Layer assignment and connectivity for multi-aquifer systems.
"""

from typing import Dict, List

from .types_3d import LayeredAquiferSystem, Node3D


def compute_layer_probability(
    layer_i: int,
    layer_j: int,
    aquitard_p: float = 0.3,
) -> float:
    """
    Compute probability of flow between aquifer layers.

    The probability decays exponentially with the number of aquitards
    that must be crossed.

    Parameters
    ----------
    layer_i, layer_j : int
        Layer indices (1-indexed)
    aquitard_p : float
        Base probability of crossing one aquitard
        Typical range: 0.1 to 0.5
        - 0.5: leaky aquitard
        - 0.3: moderate confinement
        - 0.1: strong aquitard

    Returns
    -------
    float
        Probability in [0, 1]

    Mathematical Implementation
    ---------------------------
    layer_diff = |layer_i - layer_j|

    If layer_diff == 0:
        P = 1.0  # same layer, no barrier

    If layer_diff == 1:
        P = aquitard_p  # cross one aquitard

    If layer_diff > 1:
        P = aquitard_p^layer_diff  # compound probability

    Physical Interpretation
    -----------------------
    Flow between layers requires water to cross aquitards (low-K units).
    The probability represents the likelihood of sufficient hydraulic
    connection through fractures, windows, or diffuse leakage.

    Crossing multiple aquitards multiplies the resistance:
    - 1 aquitard: P = 0.3
    - 2 aquitards: P = 0.09
    - 3 aquitards: P = 0.027

    Example
    -------
    >>> compute_layer_probability(1, 1, aquitard_p=0.3)
    1.0  # same layer

    >>> compute_layer_probability(1, 2, aquitard_p=0.3)
    0.3  # adjacent layers

    >>> compute_layer_probability(1, 3, aquitard_p=0.3)
    0.09  # skip one layer (0.3²)
    """
    layer_diff = abs(layer_i - layer_j)

    if layer_diff == 0:
        return 1.0  # Same layer, no barrier

    # Exponential decay with number of aquitards
    return aquitard_p**layer_diff


def assign_layers_to_nodes(
    nodes: Dict[str, Node3D],
    layer_system: LayeredAquiferSystem,
    depth_key: str = "z",
) -> None:
    """
    Assign aquifer layers to nodes based on screen depth.

    Modifies nodes in-place by setting aquifer_layer and aquifer_name.

    Parameters
    ----------
    nodes : Dict[str, Node3D]
        Nodes to assign layers to
    layer_system : LayeredAquiferSystem
        Layer definitions
    depth_key : str
        Which depth attribute to use: "z" or "z_mASL"

    Mathematical Implementation
    ---------------------------
    For each node:
        depth = node.z or node.z_mASL (depending on depth_key)
        layer_index = layer_system.assign_layer(depth)
        node.aquifer_layer = layer_index
        node.aquifer_name = layer_system.layer_names[layer_index - 1]

    Example
    -------
    >>> system = LayeredAquiferSystem(
    ...     n_layers=3,
    ...     layer_names=["Shallow", "Intermediate", "Deep"],
    ...     layer_tops=[0, 30, 100],
    ...     layer_bottoms=[30, 100, 300],
    ...     aquitard_p=[0.3, 0.1],
    ...     anisotropy_factors=[0.2, 0.1, 0.05]
    ... )
    >>> nodes = {"W1": Node3D("W1", x=0, y=0, z=50, elevation_m=100)}
    >>> assign_layers_to_nodes(nodes, system)
    >>> nodes["W1"].aquifer_layer
    2  # Intermediate layer (30-100m)
    """
    for node in nodes.values():
        # Get depth
        if depth_key == "z":
            depth = node.z
        elif depth_key == "z_mASL":
            if node.z_mASL is not None:
                # Convert mASL back to depth below surface
                depth = node.elevation_m - node.z_mASL
            else:
                depth = node.z
        else:
            depth = node.z

        # Assign layer
        layer_index = layer_system.assign_layer(depth)
        node.aquifer_layer = layer_index

        # Assign layer name
        if 1 <= layer_index <= len(layer_system.layer_names):
            node.aquifer_name = layer_system.layer_names[layer_index - 1]
        else:
            node.aquifer_name = f"Layer {layer_index}"


def get_anisotropy_for_layer(
    layer_system: LayeredAquiferSystem,
    layer_index: int,
) -> float:
    """
    Get anisotropy factor for a specific layer.

    Parameters
    ----------
    layer_system : LayeredAquiferSystem
        Layer definitions with anisotropy_factors
    layer_index : int
        Layer index (1-indexed)

    Returns
    -------
    float
        Anisotropy factor α_v for the layer

    Notes
    -----
    If layer_index is out of range or anisotropy_factors is shorter
    than expected, returns a default value of 0.1.

    Example
    -------
    >>> system = LayeredAquiferSystem(
    ...     n_layers=3,
    ...     anisotropy_factors=[0.2, 0.1, 0.05],
    ...     ...
    ... )
    >>> get_anisotropy_for_layer(system, 2)
    0.1  # Intermediate layer
    """
    # Convert to 0-indexed
    index = layer_index - 1

    if 0 <= index < len(layer_system.anisotropy_factors):
        return layer_system.anisotropy_factors[index]
    else:
        # Default anisotropy
        return 0.1


def get_aquitard_probability(
    layer_system: LayeredAquiferSystem,
    layer_i: int,
    layer_j: int,
) -> float:
    """
    Get aquitard probability for specific layer pair.

    Allows for heterogeneous aquitard properties (some leakier than others).

    Parameters
    ----------
    layer_system : LayeredAquiferSystem
        Layer definitions with aquitard_p list
    layer_i, layer_j : int
        Layer indices (1-indexed)

    Returns
    -------
    float
        Probability of flow between layers

    Implementation
    --------------
    If same layer:
        return 1.0

    If adjacent layers (diff = 1):
        # Get specific aquitard probability
        aquitard_index = min(layer_i, layer_j) - 1
        return layer_system.aquitard_p[aquitard_index]

    If non-adjacent (diff > 1):
        # Compound probability across multiple aquitards
        P_total = 1.0
        for each aquitard between layer_i and layer_j:
            P_total *= aquitard_p[aquitard_index]
        return P_total

    Example
    -------
    >>> system = LayeredAquiferSystem(
    ...     n_layers=3,
    ...     aquitard_p=[0.3, 0.1],  # [Shallow-Intermediate, Intermediate-Deep]
    ...     ...
    ... )
    >>> get_aquitard_probability(system, 1, 2)
    0.3  # Shallow to Intermediate

    >>> get_aquitard_probability(system, 2, 3)
    0.1  # Intermediate to Deep

    >>> get_aquitard_probability(system, 1, 3)
    0.03  # Shallow to Deep (0.3 × 0.1)
    """
    layer_diff = abs(layer_i - layer_j)

    if layer_diff == 0:
        return 1.0  # Same layer

    # Determine range of aquitards to cross
    layer_min = min(layer_i, layer_j)
    layer_max = max(layer_i, layer_j)

    # Compound probability
    p_total = 1.0
    for layer_from in range(layer_min, layer_max):
        # Aquitard between layer_from and layer_from + 1
        aquitard_index = layer_from - 1

        if 0 <= aquitard_index < len(layer_system.aquitard_p):
            p_total *= layer_system.aquitard_p[aquitard_index]
        else:
            # Default if missing data
            p_total *= 0.3

    return p_total


def create_layer_system_from_dict(layer_def: Dict) -> LayeredAquiferSystem:
    """
    Create LayeredAquiferSystem from dictionary (e.g., loaded from YAML).

    Parameters
    ----------
    layer_def : Dict
        Layer definition with keys:
        - n_layers: int
        - names: List[str]
        - tops: List[float]
        - bottoms: List[float]
        - aquitard_p: List[float]
        - anisotropy: List[float] (optional)

    Returns
    -------
    LayeredAquiferSystem
        Constructed layer system

    Example
    -------
    >>> layer_def = {
    ...     "n_layers": 3,
    ...     "names": ["Shallow", "Intermediate", "Deep"],
    ...     "tops": [0, 30, 100],
    ...     "bottoms": [30, 100, 300],
    ...     "aquitard_p": [0.3, 0.1],
    ...     "anisotropy": [0.2, 0.1, 0.05]
    ... }
    >>> system = create_layer_system_from_dict(layer_def)
    """
    # Default anisotropy if not provided
    if "anisotropy" not in layer_def or not layer_def["anisotropy"]:
        anisotropy = [0.1] * layer_def["n_layers"]
    else:
        anisotropy = layer_def["anisotropy"]

    return LayeredAquiferSystem(
        n_layers=layer_def["n_layers"],
        layer_names=layer_def["names"],
        layer_tops=layer_def["tops"],
        layer_bottoms=layer_def["bottoms"],
        aquitard_p=layer_def["aquitard_p"],
        anisotropy_factors=anisotropy,
    )
