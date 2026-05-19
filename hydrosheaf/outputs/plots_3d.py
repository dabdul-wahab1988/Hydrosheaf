"""
3D Visualization of Flow Networks using PyVista.
"""

from typing import Dict, List, Optional, Union
import logging
from ..graph3d.types_3d import Network3D, Node3D, Edge3D
from ..graph.types import Edge
from .utils import PlotConfig

logger = logging.getLogger(__name__)

def get_layer_color(layer_index: Optional[int]) -> str:
    """Get color for aquifer layer."""
    colors = [
        "white",      # 0/None
        "cyan",       # 1: Shallow
        "dodgerblue", # 2: Intermediate
        "blue",       # 3: Deep
        "navy",       # 4: Very deep
        "purple"      # 5+
    ]
    if layer_index is None or layer_index < 1:
        return "gray"
    if layer_index >= len(colors):
        return colors[-1]
    return colors[layer_index]

def plot_network_3d(
    nodes: Union[Network3D, Dict[str, Node3D], List[Dict[str, object]]],
    edges: Optional[Union[List[Edge3D], List[Edge]]] = None,
    output_path: Optional[str] = None,
    z_exaggeration: float = 10.0,
    show: bool = False,
    convert_latlon: bool = True,
    config: Optional[PlotConfig] = None
):
    """
    Create 3D plot of the flow network.
    
    Parameters
    ----------
    nodes : Network3D or Dict[str, Node3D] or List[sample_dict]
        Network data or nodes.
    edges : List[Edge3D] or List[Edge], optional
        Edges list. If nodes is Network3D, this is ignored.
    output_path : str, optional
        Path to save screenshot (e.g. "network.png").
    z_exaggeration : float
        Vertical exaggeration factor.
    show : bool
        If True, open interactive window.
    convert_latlon : bool
        If True, attempts to convert lat/lon to local meters.
    config: PlotConfig
        Configuration for styling and output.
    """
    try:
        import pyvista as pv
        import numpy as np
    except ImportError:
        logger.warning("PyVista not installed. Skipping 3D plot.")
        print("Error: PyVista is required for 3D plotting. Install with 'pip install pyvista'.")
        return
        
    if config is None:
        config = PlotConfig()

    # Normalize inputs
    node_map: Dict[str, Node3D] = {}
    edge_list: List[object] = []
    
    if hasattr(nodes, "nodes") and hasattr(nodes, "edges"):
        # It's a Network3D object
        node_map = nodes.nodes
        edge_list = nodes.edges
    elif isinstance(nodes, dict):
        node_map = nodes
        if edges:
            edge_list = edges
    elif isinstance(nodes, list):
        # List of dicts (samples), need to convert to basic Node3D-like structs
        # This is a fallback if called from CLI without full 3D build
        for s in nodes:
            nid = str(s.get("site_id") or s.get("sample_id"))
            x = float(s.get("x") or s.get("lon", 0))
            y = float(s.get("y") or s.get("lat", 0))
            z = float(s.get("z") or s.get("screen_depth") or s.get("elevation") or 0)
            layer = int(s.get("aquifer_layer", 1)) if s.get("aquifer_layer") else None
            node_map[nid] = Node3D(
                node_id=nid, x=x, y=y, z=z, elevation_m=0, aquifer_layer=layer
            )
        if edges:
            edge_list = edges

    if not node_map:
        logger.warning("No nodes to plot.")
        return

    # Coordinate conversion (simple local projection)
    lats = [n.y for n in node_map.values()]
    lons = [n.x for n in node_map.values()]
    
    is_latlon = convert_latlon and (
        (min(lats) > -90 and max(lats) < 90) and 
        (min(lons) > -180 and max(lons) < 180) and
        (max(lons) - min(lons) < 10) # range check to avoid confusing small projected coords with latlon
    )
    
    lat0, lon0 = np.mean(lats), np.mean(lons)
    
    def transform(n):
        if is_latlon:
            # Approx conversion
            R = 6371000
            dx = np.radians(n.x - lon0) * R * np.cos(np.radians(lat0))
            dy = np.radians(n.y - lat0) * R
            return dx, dy, n.z
        return n.x, n.y, n.z

    plotter = pv.Plotter(off_screen=bool(output_path and not show))
    
    # Add nodes
    # Calculate sensible radius based on domain size
    coords = np.array([transform(n) for n in node_map.values()])
    if len(coords) > 1:
        span = np.max(coords, axis=0) - np.min(coords, axis=0)
        radius = max(span[0], span[1]) / 100.0
        # Add bounds box if span is significant
        plotter.add_bounding_box(color='black', line_width=1.0)
    else:
        radius = 10.0
    
    radius = max(radius, 1.0) # minimum size

    for node in node_map.values():
        x, y, z = transform(node)
        # Invert Z for depth if positive-down, but multiply by exaggeration
        # Usually z is depth (positive down). Plotting: z_plot = -z * exag
        z_plot = -z * z_exaggeration
        
        sphere = pv.Sphere(radius=radius, center=(x, y, z_plot))
        color = get_layer_color(node.aquifer_layer)
        plotter.add_mesh(sphere, color=color)

    # Add edges
    for edge in edge_list:
        u_id = getattr(edge, "u", None)
        v_id = getattr(edge, "v", None)
        
        if u_id not in node_map or v_id not in node_map:
            continue
            
        u = node_map[u_id]
        v = node_map[v_id]
        
        ux, uy, uz = transform(u)
        vx, vy, vz = transform(v)
        
        p1 = (ux, uy, -uz * z_exaggeration)
        p2 = (vx, vy, -vz * z_exaggeration)
        
        line = pv.Line(p1, p2)
        
        # Color logic
        # If it's Edge3D, use edge_type or same_layer
        # If Edge, assume blue
        color = "gray"
        width = 1
        
        if hasattr(edge, "same_layer"):
            if getattr(edge, "same_layer"):
                color = "skyblue"
                width = 2
            else:
                color = "red" # vertical/cross-layer
                width = 3
        
        plotter.add_mesh(line, color=color, line_width=width)

    # Enhanced Context
    plotter.add_axes()
    plotter.add_text(f"3D Flow Network (Z-Exag: {z_exaggeration}x)", font_size=10)
    plotter.show_grid()  # Show grid lines for spatial reference
    
    # Add scale bar if supported
    if hasattr(plotter, 'add_scale_bar'):
        plotter.add_scale_bar(title="Scale (m)", color="black")
    
    if output_path:
        plotter.screenshot(output_path)
        logger.info(f"Saved 3D plot to {output_path}")
        print(f"Saved 3D plot to {output_path}")
        
    if show:
        plotter.show()
    else:
        plotter.close()
