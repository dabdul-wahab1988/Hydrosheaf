"""
Data structures for 3D aquifer networks.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class Node3D:
    """
    A 3D node (well/sample point) in the aquifer network.

    Represents a groundwater sampling location with full 3D coordinates,
    hydraulic data, and chemical measurements.

    Attributes
    ----------
    node_id : str
        Unique identifier for the node
    x : float
        Longitude (degrees) or easting (m) depending on coordinate system
    y : float
        Latitude (degrees) or northing (m) depending on coordinate system
    z : float
        Screen midpoint depth (m below surface, positive down)
    elevation_m : float
        Ground surface elevation (mASL - meters above sea level)
    z_mASL : Optional[float]
        Screen midpoint elevation in mASL (auto-computed if None)
    screen_top : Optional[float]
        Depth to top of screened interval (m below surface)
    screen_bottom : Optional[float]
        Depth to bottom of screened interval (m below surface)
    hydraulic_head : Optional[float]
        Measured hydraulic head (m)
    head_uncertainty : float
        Uncertainty in head measurement (m), default 0.5
    aquifer_layer : Optional[int]
        Layer index (1=shallow, 2=intermediate, 3=deep, etc.)
    aquifer_name : Optional[str]
        Human-readable aquifer name (e.g., "Quaternary", "Cretaceous")
    concentrations : Optional[List[float]]
        Ion concentrations in mmol/L, following config.ion_order
    """

    node_id: str

    # Geographic position
    x: float
    y: float

    # Vertical position
    z: float
    elevation_m: float
    z_mASL: Optional[float] = None

    # Screened interval
    screen_top: Optional[float] = None
    screen_bottom: Optional[float] = None

    # Hydraulic data
    hydraulic_head: Optional[float] = None
    head_uncertainty: float = 0.5

    # Layer information
    aquifer_layer: Optional[int] = None
    aquifer_name: Optional[str] = None

    # Chemical data
    concentrations: Optional[List[float]] = None

    def __post_init__(self):
        """Compute z_mASL if not provided."""
        if self.z_mASL is None and self.elevation_m is not None:
            # z is depth below surface (positive down)
            # z_mASL is elevation above sea level
            self.z_mASL = self.elevation_m - self.z


@dataclass
class Edge3D:
    """
    A 3D edge connecting two nodes in the aquifer network.

    Represents a potential flow path with geometric, probabilistic,
    and hydraulic properties.

    Attributes
    ----------
    edge_id : str
        Unique identifier (typically "u->v")
    u : str
        Upstream node_id
    v : str
        Downstream node_id
    horizontal_distance_m : float
        Horizontal separation between nodes (m)
    vertical_distance_m : float
        Absolute vertical separation |z_u - z_v| (m)
    distance_3d : float
        Anisotropic 3D distance (m)
    edge_type : str
        Classification: "horizontal", "vertical_leakage", or "oblique"
    same_layer : bool
        Whether nodes are in the same aquifer layer
    prob_head : float
        Probability based on hydraulic gradient
    prob_distance : float
        Probability based on distance decay
    prob_layer : float
        Probability based on layer connectivity
    prob_combined : float
        Combined probability P_head × P_dist × P_layer × P_screen
    horizontal_gradient : Optional[float]
        Horizontal hydraulic gradient (dimensionless)
    vertical_gradient : Optional[float]
        Vertical hydraulic gradient (dimensionless)
    layer_from : Optional[int]
        Source layer index
    layer_to : Optional[int]
        Destination layer index
    """

    edge_id: str
    u: str
    v: str

    # Edge geometry
    horizontal_distance_m: float
    vertical_distance_m: float
    distance_3d: float

    # Edge classification
    edge_type: str
    same_layer: bool

    # Probability components
    prob_head: float
    prob_distance: float
    prob_layer: float
    prob_combined: float

    # Hydraulic gradients
    horizontal_gradient: Optional[float] = None
    vertical_gradient: Optional[float] = None

    # Layer transition
    layer_from: Optional[int] = None
    layer_to: Optional[int] = None


@dataclass
class LayeredAquiferSystem:
    """
    Multi-layer aquifer system definition.

    Defines vertical layering structure with aquitards separating
    distinct aquifer units.

    Attributes
    ----------
    n_layers : int
        Number of aquifer layers
    layer_names : List[str]
        Names for each layer (e.g., ["Shallow", "Intermediate", "Deep"])
    layer_tops : List[float]
        Depth to top of each layer (m below surface)
    layer_bottoms : List[float]
        Depth to bottom of each layer (m below surface)
    aquitard_p : List[float]
        Probability of flow across each aquitard
        Length = n_layers - 1 (between adjacent layers)
    anisotropy_factors : List[float]
        Vertical anisotropy α_v for each layer
        Typical range: 0.01 to 0.1 (K_h / K_v ratio indicator)

    Example
    -------
    system = LayeredAquiferSystem(
        n_layers=3,
        layer_names=["Quaternary", "Tertiary", "Cretaceous"],
        layer_tops=[0, 30, 100],
        layer_bottoms=[30, 100, 300],
        aquitard_p=[0.3, 0.1],
        anisotropy_factors=[0.2, 0.1, 0.05]
    )
    """

    n_layers: int
    layer_names: List[str]
    layer_tops: List[float]
    layer_bottoms: List[float]
    aquitard_p: List[float]
    anisotropy_factors: List[float] = field(default_factory=lambda: [0.1])

    def assign_layer(self, depth: float) -> int:
        """
        Assign node to layer based on screen depth.

        Parameters
        ----------
        depth : float
            Screen midpoint depth (m below surface)

        Returns
        -------
        int
            Layer index (1-indexed). Returns deepest layer if beyond range.

        Example
        -------
        >>> system.assign_layer(50)  # depth=50m
        2  # Tertiary layer (30-100m)
        """
        for i, (top, bottom) in enumerate(zip(self.layer_tops, self.layer_bottoms)):
            if top <= depth < bottom:
                return i + 1  # 1-indexed layers
        return self.n_layers  # Default to deepest layer


@dataclass
class Network3D:
    """
    Complete 3D aquifer network.

    Contains all nodes, edges, and layer system definition for a
    3D groundwater monitoring network.

    Attributes
    ----------
    nodes : Dict[str, Node3D]
        All nodes indexed by node_id
    edges : List[Edge3D]
        All inferred flow edges
    layer_system : Optional[LayeredAquiferSystem]
        Multi-layer system definition
    n_horizontal_edges : int
        Count of primarily horizontal edges
    n_vertical_edges : int
        Count of primarily vertical (leakage) edges
    n_cross_layer_edges : int
        Count of edges crossing aquifer layers
    """

    nodes: Dict[str, Node3D]
    edges: List[Edge3D]
    layer_system: Optional[LayeredAquiferSystem] = None

    # Summary statistics
    n_horizontal_edges: int = 0
    n_vertical_edges: int = 0
    n_cross_layer_edges: int = 0

    def __post_init__(self):
        """Compute summary statistics after initialization."""
        self.n_horizontal_edges = sum(1 for e in self.edges if e.edge_type == "horizontal")
        self.n_vertical_edges = sum(1 for e in self.edges if e.edge_type == "vertical_leakage")
        self.n_cross_layer_edges = sum(1 for e in self.edges if not e.same_layer)
