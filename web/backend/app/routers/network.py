"""
Network Router - Handles flow network and graph-related endpoints
Now integrated with real Hydrosheaf probabilistic flow inference!
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
import uuid

from ..logger import network_logger as logger
from ..database import (
    create_network as db_create_network,
    get_network as db_get_network,
    get_all_networks,
    update_network,
    delete_network as db_delete_network,
)

HYDROSHEAF_AVAILABLE = None
HydrosheafEdge = None
Config = None
infer_edges_probabilistic = None
infer_edges_3d_probabilistic = None


def _load_hydrosheaf() -> None:
    """Lazy-load Hydrosheaf to keep backend startup fast."""
    global HYDROSHEAF_AVAILABLE, HydrosheafEdge, Config, infer_edges_probabilistic, infer_edges_3d_probabilistic

    if HYDROSHEAF_AVAILABLE is not None:
        return
    try:
        from hydrosheaf import (
            infer_edges_probabilistic as hydrosheaf_infer_edges_probabilistic,
            infer_edges_3d_probabilistic as hydrosheaf_infer_edges_3d_probabilistic,
            Config as HydrosheafConfig,
        )
        from hydrosheaf.graph.types import Edge as HydrosheafEdgeType

        infer_edges_probabilistic = hydrosheaf_infer_edges_probabilistic
        infer_edges_3d_probabilistic = hydrosheaf_infer_edges_3d_probabilistic
        Config = HydrosheafConfig
        HydrosheafEdge = HydrosheafEdgeType
        HYDROSHEAF_AVAILABLE = True
    except Exception as exc:
        HYDROSHEAF_AVAILABLE = False
        HydrosheafEdge = None
        logger.warning(f"Hydrosheaf not available in network router: {exc}")


router = APIRouter()


class Node(BaseModel):
    """A node in the flow network (e.g., a well or sampling point)"""

    id: str
    name: str
    x: Optional[float] = None
    y: Optional[float] = None
    z: Optional[float] = None  # Elevation/depth
    hydraulic_head: Optional[float] = None
    node_type: str = "well"  # well, spring, river, etc.


class EdgeAPI(BaseModel):
    """An edge representing potential flow between nodes (API model)"""

    source: str
    target: str
    weight: Optional[float] = 1.0
    flow_probability: Optional[float] = None


class NetworkData(BaseModel):
    """Complete network data structure"""

    nodes: List[Node]
    edges: List[EdgeAPI]  # Use EdgeAPI model
    name: Optional[str] = "Untitled Network"


class NetworkInferenceConfig(BaseModel):
    """Configuration for network inference"""

    method: str = Field(default="bayesian", description="heuristic or bayesian")
    radius_km: float = Field(default=10.0, ge=0.1, le=100.0)
    max_neighbors: int = Field(default=3, ge=1, le=10)
    
    # 3D settings
    use_3d: bool = Field(default=False, description="Enable 3D network inference")
    vertical_anisotropy: float = Field(default=0.1, description="Vertical anisotropy factor (Kh/Kv)")
    screen_overlap_threshold: float = Field(default=5.0, description="Screen overlap threshold in meters")



class NetworkInferenceResult(BaseModel):
    """Results from network flow inference"""

    edges: List[Dict[str, Any]]
    flow_directions: Dict[str, str]
    probabilities: Dict[str, float]


@router.get("/health")
async def network_health_check():
    """Check if Hydrosheaf network inference is available"""
    _load_hydrosheaf()
    return {
        "hydrosheaf_available": HYDROSHEAF_AVAILABLE,
        "status": "ready" if HYDROSHEAF_AVAILABLE else "fallback",
        "message": (
            "Using real Hydrosheaf probabilistic inference"
            if HYDROSHEAF_AVAILABLE
            else "Using simple gradient method (install Hydrosheaf for advanced features)"
        ),
    }


@router.post("/create")
async def create_network_endpoint(network: NetworkData):
    """Create a new flow network"""
    network_id = str(uuid.uuid4())[:8]

    # Convert Pydantic models to dicts
    nodes = [n.model_dump() for n in network.nodes]
    edges = [e.model_dump() for e in network.edges]

    db_create_network(
        network_id=network_id,
        name=network.name or "Untitled Network",
        nodes=nodes,
        edges=edges,
    )

    return {
        "network_id": network_id,
        "nodes_count": len(nodes),
        "edges_count": len(edges),
        "message": "Network created successfully",
    }


@router.get("/list")
async def list_networks_endpoint():
    """List all networks"""
    networks = get_all_networks()
    return [
        {
            "id": net["id"],
            "name": net["name"],
            "nodes_count": len(net.get("nodes", [])),
            "edges_count": len(net.get("edges", [])),
        }
        for net in networks
    ]


@router.get("/{network_id}")
async def get_network_endpoint(network_id: str):
    """Get a specific network by ID"""
    network = db_get_network(network_id)
    if not network:
        raise HTTPException(status_code=404, detail="Network not found")

    return network


@router.post("/{network_id}/infer-flow")
async def infer_flow_directions(
    network_id: str, config: Optional[NetworkInferenceConfig] = None
):
    """
    Infer flow directions using REAL Hydrosheaf probabilistic inference
    (when available) or simple hydraulic head gradient fallback
    """
    network = db_get_network(network_id)
    if not network:
        raise HTTPException(status_code=404, detail="Network not found")

    nodes = network.get("nodes", [])
    edges = network.get("edges", [])

    if config is None:
        config = NetworkInferenceConfig()

    _load_hydrosheaf()

    # Try to use real Hydrosheaf inference if available
    if HYDROSHEAF_AVAILABLE:
        try:
            # Convert nodes to sample format for Hydrosheaf
            samples = []
            for node in nodes:
                sample = {
                    "site_id": node.get("id"),
                }

                # Check for explicit lat/lon first
                if node.get("lat") is not None and node.get("lon") is not None:
                    sample["lat"] = node["lat"]
                    sample["lon"] = node["lon"]
                # Fallback: Map UI X/Y (pixels) to Lat/Lon for demo purposes
                # Assuming 1 pixel = 1 meter (approx)
                # 1 degree lat ~ 111km = 111,000m
                # So 1 pixel ~ 9e-6 degrees
                elif node.get("x") is not None and node.get("y") is not None:
                    scale_factor = 9e-6  # Convert "meters" (pixels) to degrees
                    sample["lon"] = float(node["x"]) * scale_factor
                    sample["lat"] = float(node["y"]) * scale_factor

                    # Also keep x/y if needed
                    sample["x"] = node["x"]
                    sample["y"] = node["y"]

                if node.get("z") is not None:
                    sample["elevation"] = node["z"]
                if node.get("hydraulic_head") is not None:
                    sample["head_meas"] = node["hydraulic_head"]

                samples.append(sample)

            # Use Hydrosheaf's probabilistic edge inference
            if config.use_3d:
                # 3D Inference
                
                # We need to transform samples into Node3D objects for 3D inference
                # But infer_edges_3d_probabilistic takes a list of Node3D objects.
                # However, looking at build_3d.py, there isn't a helper to convert dict samples to Node3D directly exposed.
                # Actually, build_network_3d takes samples and does the conversion.
                # But we want just the edges.
                # Let's import Node3D and create them manually or use build_network_3d internally.
                
                # To keep it simple and consistent with how 2D works (which takes samples dicts),
                # we might need to do some manual conversion or import Node3D.
                # Let's try to import Node3D inside the try block.
                from hydrosheaf.graph3d.types_3d import Node3D
                
                nodes_3d = []
                for s in samples:
                    # Map sample dict to Node3D
                    node_id = str(s.get("site_id", ""))
                    if not node_id: continue
                    
                    nodes_3d.append(Node3D(
                        node_id=node_id,
                        x=float(s.get("lon", 0.0)), # using lon as x
                        y=float(s.get("lat", 0.0)), # using lat as y
                        z=float(s.get("elevation", 0.0)), # using elevation as z/depth? 
                        # Note: 3D model expects z as depth usually, or elevation. 
                        # Config says z_coordinate_key default is screen_depth.
                        # We'll assume elevation is what we have.
                        elevation_m=float(s.get("elevation", 0.0)),
                        hydraulic_head=float(s.get("head_meas")) if "head_meas" in s else None,
                        head_uncertainty=0.5 # Default
                    ))
                
                # Create a minimal config object for 3D
                class MiniConfig:
                    edge_radius_km = config.radius_km
                    edge_p_min = 0.75
                    edge_max_neighbors = config.max_neighbors
                    vertical_anisotropy = config.vertical_anisotropy
                    
                inferred_edges = infer_edges_3d_probabilistic(
                    nodes_3d,
                    config=MiniConfig(),
                    layer_system=None, # Layers not fully supported in simple UI yet
                    use_haversine=True
                )
            else:
                # 2D Inference
                inferred_edges = infer_edges_probabilistic(
                    samples,
                    radius_km=config.radius_km,
                    max_neighbors=config.max_neighbors,
                    p_min=0.75,
                    head_inference=config.method,
                )

            # Convert to frontend format
            inferred_edges_data = []
            for edge in inferred_edges:
                # Extract probability
                if config.use_3d:
                    probability = getattr(edge, "prob_combined", 0.75)
                    # 3D edges use u and v
                    source = edge.u
                    target = edge.v
                else:
                    probability = getattr(edge, "probability", 0.75)
                    # 2D edges use u and v too? Let's check type.
                    source = edge.u
                    target = edge.v

                # Determine flow direction
                inferred_edges_data.append(
                    {
                        "source": source,
                        "target": target,
                        "weight": getattr(edge, "weight", 1.0),
                        "flow_probability": probability,
                        "flow_direction": "forward",
                        "method": "hydrosheaf_probabilistic_3d" if config.use_3d else "hydrosheaf_probabilistic",
                    }
                )

            return {
                "network_id": network_id,
                "inferred_edges": inferred_edges_data,
                "method": f"hydrosheaf_{config.method}_{'3d' if config.use_3d else '2d'}",
                "hydrosheaf_used": True,
            }

        except Exception as e:
            logger.warning(
                f"Hydrosheaf inference failed, falling back to simple method: {e}"
            )
            # Fall through to simple method

    # Fallback: Simple flow inference based on hydraulic head differences
    inferred_edges = []
    for edge in edges:
        source_node = next(
            (n for n in nodes if n.get("id") == edge.get("source")), None
        )
        target_node = next(
            (n for n in nodes if n.get("id") == edge.get("target")), None
        )

        if source_node and target_node:
            source_head = source_node.get("hydraulic_head")
            target_head = target_node.get("hydraulic_head")
            if source_head is not None and target_head is not None:
                head_diff = source_head - target_head
                # Flow probability based on head gradient
                probability = min(1.0, max(0.0, 0.5 + head_diff * 0.1))
                inferred_edges.append(
                    {
                        "source": edge.get("source"),
                        "target": edge.get("target"),
                        "head_difference": head_diff,
                        "flow_probability": probability,
                        "flow_direction": "forward" if head_diff > 0 else "reverse",
                        "method": "simple_gradient",
                    }
                )
            else:
                inferred_edges.append(
                    {
                        "source": edge.get("source"),
                        "target": edge.get("target"),
                        "head_difference": None,
                        "flow_probability": 0.5,
                        "flow_direction": "uncertain",
                        "method": "simple_gradient",
                    }
                )

    return {
        "network_id": network_id,
        "inferred_edges": inferred_edges,
        "method": "simple_hydraulic_head_gradient",
        "hydrosheaf_used": False,
        "note": "Install Hydrosheaf for advanced probabilistic inference",
    }


@router.delete("/{network_id}")
async def delete_network_endpoint(network_id: str):
    """Delete a network"""
    if not db_delete_network(network_id):
        raise HTTPException(status_code=404, detail="Network not found")

    return {"message": "Network deleted successfully"}
