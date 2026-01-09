"""
3D distance calculations for aquifer networks.
"""

import math
from typing import Tuple

from .types_3d import Node3D


def haversine_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Compute great-circle distance between two points on Earth.

    Parameters
    ----------
    lat1, lon1 : float
        First point (degrees)
    lat2, lon2 : float
        Second point (degrees)

    Returns
    -------
    float
        Distance in kilometers

    Mathematical Implementation
    ---------------------------
    Using the haversine formula:

    a = sin²(Δφ/2) + cos(φ1)·cos(φ2)·sin²(Δλ/2)
    c = 2·atan2(√a, √(1-a))
    d = R·c

    Where:
    - φ is latitude, λ is longitude
    - R = 6371 km (Earth's mean radius)
    """
    R = 6371.0  # Earth radius in km

    # Convert to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Differences
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    # Haversine formula
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance_km = R * c
    return distance_km


def compute_3d_distance(
    node_i: Node3D,
    node_j: Node3D,
    anisotropy_factor: float = 0.1,
    use_haversine: bool = True,
) -> Tuple[float, float, float]:
    """
    Compute 3D distance between two nodes with anisotropy.

    Parameters
    ----------
    node_i, node_j : Node3D
        The two nodes
    anisotropy_factor : float
        α_v: vertical distance scaling factor
        Typical range: 0.01 to 0.1
        Represents K_h / K_v ratio indicator
        Small values make vertical flow "harder"
    use_haversine : bool
        If True, use Haversine for geographic coordinates (lat/lon)
        If False, use Euclidean distance for projected coordinates (m)

    Returns
    -------
    Tuple[float, float, float]
        (horizontal_distance_m, vertical_distance_m, total_3d_distance_m)

    Mathematical Implementation
    ---------------------------
    1. Horizontal distance:
       If use_haversine:
           d_xy = haversine(lat_i, lon_i, lat_j, lon_j) × 1000  # km to m
       Else:
           d_xy = sqrt((x_i - x_j)² + (y_i - y_j)²)

    2. Vertical distance:
       Δz = |z_i - z_j|  # absolute depth difference

    3. Anisotropic 3D distance:
       d_3d = sqrt(d_xy² + (Δz / α_v)²)

    Physical Interpretation
    -----------------------
    The anisotropy factor α_v scales vertical distance to account for
    reduced vertical permeability:

    - α_v = 1.0: isotropic (K_h = K_v)
    - α_v = 0.1: K_h = 10 × K_v (typical aquifer)
    - α_v = 0.01: K_h = 100 × K_v (strong aquitard)

    With α_v = 0.1, a 50m vertical separation "feels like" 500m
    horizontally for flow calculations.

    Example
    -------
    >>> node_a = Node3D("A", x=0, y=0, z=0, elevation_m=100)
    >>> node_b = Node3D("B", x=100, y=0, z=50, elevation_m=100)
    >>> d_xy, d_z, d_3d = compute_3d_distance(node_a, node_b, anisotropy_factor=0.1)
    >>> d_xy  # 100.0 m
    >>> d_z   # 50.0 m
    >>> d_3d  # 509.9 m (sqrt(100² + (50/0.1)²))
    """
    # Horizontal distance
    if use_haversine:
        # Assume x, y are lon, lat in degrees
        d_xy = haversine_distance(node_i.y, node_i.x, node_j.y, node_j.x) * 1000.0  # km to m
    else:
        # Euclidean distance for projected coordinates
        d_xy = math.sqrt((node_i.x - node_j.x) ** 2 + (node_i.y - node_j.y) ** 2)

    # Vertical distance (absolute)
    d_z = abs(node_i.z - node_j.z)

    # Anisotropic 3D distance
    # Dividing by anisotropy_factor makes vertical distance "feel" larger
    if anisotropy_factor > 0:
        scaled_vertical = d_z / anisotropy_factor
    else:
        scaled_vertical = d_z  # Fallback if α_v = 0

    d_3d = math.sqrt(d_xy**2 + scaled_vertical**2)

    return d_xy, d_z, d_3d


def compute_screen_overlap(
    node_i: Node3D,
    node_j: Node3D,
) -> Tuple[float, float]:
    """
    Compute screened interval overlap between two wells.

    Wells with overlapping screened intervals are more likely to be
    hydraulically connected.

    Parameters
    ----------
    node_i, node_j : Node3D
        Nodes with screen_top and screen_bottom defined

    Returns
    -------
    Tuple[float, float]
        (overlap_m, overlap_fraction)

        overlap_m : float
            Absolute overlap length in meters
        overlap_fraction : float
            Overlap normalized by shorter screen length
            Range: [0, 1]

    Mathematical Implementation
    ---------------------------
    Given two intervals [top_i, bot_i] and [top_j, bot_j]:

    1. Find overlap bounds:
       overlap_top = max(top_i, top_j)
       overlap_bot = min(bot_i, bot_j)

    2. Compute overlap:
       overlap = max(0, overlap_bot - overlap_top)

    3. Normalize by shorter screen:
       len_i = bot_i - top_i
       len_j = bot_j - top_j
       overlap_frac = overlap / min(len_i, len_j) if min > 0 else 0

    Example
    -------
    >>> node_a = Node3D("A", ..., screen_top=30, screen_bottom=50)
    >>> node_b = Node3D("B", ..., screen_top=40, screen_bottom=60)
    >>> overlap, frac = compute_screen_overlap(node_a, node_b)
    >>> overlap  # 10.0 m (interval [40, 50])
    >>> frac     # 0.5 (10m / 20m, the shorter screen)
    """
    # Check if screen intervals are defined
    if (
        node_i.screen_top is None
        or node_i.screen_bottom is None
        or node_j.screen_top is None
        or node_j.screen_bottom is None
    ):
        return 0.0, 0.0

    # Compute overlap bounds
    overlap_top = max(node_i.screen_top, node_j.screen_top)
    overlap_bottom = min(node_i.screen_bottom, node_j.screen_bottom)

    # Overlap length (non-negative)
    overlap = max(0.0, overlap_bottom - overlap_top)

    # Screen lengths
    len_i = node_i.screen_bottom - node_i.screen_top
    len_j = node_j.screen_bottom - node_j.screen_top

    # Normalize by shorter screen
    min_len = min(len_i, len_j)
    if min_len > 0:
        overlap_frac = overlap / min_len
    else:
        overlap_frac = 0.0

    # Ensure fraction is in [0, 1]
    overlap_frac = min(1.0, max(0.0, overlap_frac))

    return overlap, overlap_frac


def classify_edge_type(
    horizontal_dist: float,
    vertical_dist: float,
    same_layer: bool,
    horizontal_threshold_ratio: float = 10.0,
    vertical_threshold_ratio: float = 0.1,
) -> str:
    """
    Classify edge as horizontal, vertical, or oblique based on geometry.

    Parameters
    ----------
    horizontal_dist : float
        Horizontal distance in meters
    vertical_dist : float
        Vertical distance in meters
    same_layer : bool
        Whether nodes are in the same aquifer layer
    horizontal_threshold_ratio : float
        Ratio d_xy / d_z above which edge is "horizontal"
        Default: 10.0 (horizontal >> vertical)
    vertical_threshold_ratio : float
        Ratio d_xy / d_z below which edge is "vertical"
        Default: 0.1 (vertical >> horizontal)

    Returns
    -------
    str
        One of: "horizontal", "vertical_leakage", "oblique"

    Classification Logic
    --------------------
    ratio = horizontal_dist / (vertical_dist + ε)

    If same_layer and ratio > horizontal_threshold_ratio:
        → "horizontal" (lateral flow within layer)

    If not same_layer and ratio < vertical_threshold_ratio:
        → "vertical_leakage" (aquitard leakage)

    Otherwise:
        → "oblique" (3D flow path)

    Example
    -------
    >>> classify_edge_type(1000, 50, same_layer=True)
    "horizontal"  # ratio = 20 > 10

    >>> classify_edge_type(5, 100, same_layer=False)
    "vertical_leakage"  # ratio = 0.05 < 0.1

    >>> classify_edge_type(100, 50, same_layer=False)
    "oblique"  # ratio = 2.0 (neither extreme)
    """
    # Avoid division by zero
    epsilon = 1e-6
    ratio = horizontal_dist / (vertical_dist + epsilon)

    # Primarily horizontal flow
    if same_layer and ratio > horizontal_threshold_ratio:
        return "horizontal"

    # Primarily vertical leakage
    if not same_layer and ratio < vertical_threshold_ratio:
        return "vertical_leakage"

    # Oblique / 3D flow path
    return "oblique"
