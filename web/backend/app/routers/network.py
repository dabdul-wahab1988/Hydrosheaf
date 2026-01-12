"""
Network Router - Handles flow network and graph-related endpoints
Now integrated with real Hydrosheaf probabilistic flow inference!
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
import sys

HYDROSHEAF_AVAILABLE = None
HydrosheafEdge = None
Config = None
infer_edges_probabilistic = None


def _load_hydrosheaf() -> None:
    """Lazy-load Hydrosheaf to keep backend startup fast."""
    global HYDROSHEAF_AVAILABLE, HydrosheafEdge, Config, infer_edges_probabilistic

    if HYDROSHEAF_AVAILABLE is not None:
        return
    try:
        from hydrosheaf import infer_edges_probabilistic as hydrosheaf_infer_edges_probabilistic, Config as HydrosheafConfig
        from hydrosheaf.graph.types import Edge as HydrosheafEdgeType

        infer_edges_probabilistic = hydrosheaf_infer_edges_probabilistic
        Config = HydrosheafConfig
        HydrosheafEdge = HydrosheafEdgeType
        HYDROSHEAF_AVAILABLE = True
    except Exception as exc:
        HYDROSHEAF_AVAILABLE = False
        HydrosheafEdge = None
        print(f"WARNING: Hydrosheaf not available in network router: {exc}", file=sys.stderr)

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


class NetworkInferenceResult(BaseModel):
    """Results from network flow inference"""
    edges: List[Dict[str, Any]]
    flow_directions: Dict[str, str]
    probabilities: Dict[str, float]


# In-memory storage for networks
networks: Dict[str, NetworkData] = {}


@router.post("/create")
async def create_network(network: NetworkData):
    """Create a new flow network"""
    network_id = network.name.lower().replace(" ", "_")
    networks[network_id] = network

    return {
        "network_id": network_id,
        "nodes_count": len(network.nodes),
        "edges_count": len(network.edges),
        "message": "Network created successfully",
    }


@router.get("/list")
async def list_networks():
    """List all networks"""
    return [
        {
            "id": net_id,
            "name": net.name,
            "nodes_count": len(net.nodes),
            "edges_count": len(net.edges),
        }
        for net_id, net in networks.items()
    ]


@router.get("/{network_id}")
async def get_network(network_id: str):
    """Get a specific network by ID"""
    if network_id not in networks:
        raise HTTPException(status_code=404, detail="Network not found")

    return networks[network_id]


@router.post("/{network_id}/infer-flow")
async def infer_flow_directions(network_id: str, config: Optional[NetworkInferenceConfig] = None):
    """
    Infer flow directions using REAL Hydrosheaf probabilistic inference
    (when available) or simple hydraulic head gradient fallback
    """
    if network_id not in networks:
        raise HTTPException(status_code=404, detail="Network not found")

    network = networks[network_id]

    if config is None:
        config = NetworkInferenceConfig()

    _load_hydrosheaf()

    # Try to use real Hydrosheaf inference if available
    if HYDROSHEAF_AVAILABLE:
        try:
            # Convert nodes to sample format for Hydrosheaf
            samples = []
            for node in network.nodes:
                sample = {
                    'site_id': node.id,
                }
                if node.x is not None and node.y is not None:
                    sample['x'] = node.x
                    sample['y'] = node.y
                if node.z is not None:
                    sample['elevation'] = node.z
                if node.hydraulic_head is not None:
                    sample['head_meas'] = node.hydraulic_head

                samples.append(sample)

            # Use Hydrosheaf's probabilistic edge inference
            hydrosheaf_config = Config(
                edge_radius_km=config.radius_km,
                edge_max_neighbors=config.max_neighbors,
                edge_head_inference=config.method,
            )

            inferred_edges = infer_edges_probabilistic(
                samples,
                hydrosheaf_config,
                radius_km=config.radius_km,
                max_neighbors=config.max_neighbors,
            )

            # Convert to frontend format
            inferred_edges_data = []
            for edge in inferred_edges:
                # Extract probability from edge metadata if available
                probability = getattr(edge, 'probability', 0.75)

                # Determine flow direction based on edge direction
                inferred_edges_data.append({
                    "source": edge.u,
                    "target": edge.v,
                    "weight": edge.weight,
                    "flow_probability": probability,
                    "flow_direction": "forward",
                    "method": "hydrosheaf_probabilistic"
                })

            return {
                "network_id": network_id,
                "inferred_edges": inferred_edges_data,
                "method": f"hydrosheaf_{config.method}",
                "hydrosheaf_used": True,
            }

        except Exception as e:
            print(f"Hydrosheaf inference failed, falling back to simple method: {e}", file=sys.stderr)
            # Fall through to simple method

    # Fallback: Simple flow inference based on hydraulic head differences
    inferred_edges = []
    for edge in network.edges:
        source_node = next((n for n in network.nodes if n.id == edge.source), None)
        target_node = next((n for n in network.nodes if n.id == edge.target), None)

        if source_node and target_node:
            if source_node.hydraulic_head and target_node.hydraulic_head:
                head_diff = source_node.hydraulic_head - target_node.hydraulic_head
                # Flow probability based on head gradient
                probability = min(1.0, max(0.0, 0.5 + head_diff * 0.1))
                inferred_edges.append({
                    "source": edge.source,
                    "target": edge.target,
                    "head_difference": head_diff,
                    "flow_probability": probability,
                    "flow_direction": "forward" if head_diff > 0 else "reverse",
                    "method": "simple_gradient"
                })
            else:
                inferred_edges.append({
                    "source": edge.source,
                    "target": edge.target,
                    "head_difference": None,
                    "flow_probability": 0.5,
                    "flow_direction": "uncertain",
                    "method": "simple_gradient"
                })

    return {
        "network_id": network_id,
        "inferred_edges": inferred_edges,
        "method": "simple_hydraulic_head_gradient",
        "hydrosheaf_used": False,
        "note": "Install Hydrosheaf for advanced probabilistic inference"
    }


@router.delete("/{network_id}")
async def delete_network(network_id: str):
    """Delete a network"""
    if network_id not in networks:
        raise HTTPException(status_code=404, detail="Network not found")

    del networks[network_id]
    return {"message": "Network deleted successfully"}


@router.get("/health")
async def network_health_check():
    """Check if Hydrosheaf network inference is available"""
    _load_hydrosheaf()
    return {
        "hydrosheaf_available": HYDROSHEAF_AVAILABLE,
        "status": "ready" if HYDROSHEAF_AVAILABLE else "fallback",
        "message": "Using real Hydrosheaf probabilistic inference" if HYDROSHEAF_AVAILABLE else "Using simple gradient method (install Hydrosheaf for advanced features)",
    }
