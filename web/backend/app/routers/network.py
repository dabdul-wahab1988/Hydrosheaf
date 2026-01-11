"""
Network Router - Handles flow network and graph-related endpoints
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any

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


class Edge(BaseModel):
    """An edge representing potential flow between nodes"""
    source: str
    target: str
    weight: Optional[float] = 1.0
    flow_probability: Optional[float] = None


class NetworkData(BaseModel):
    """Complete network data structure"""
    nodes: List[Node]
    edges: List[Edge]
    name: Optional[str] = "Untitled Network"


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
async def infer_flow_directions(network_id: str):
    """Infer flow directions based on hydraulic head data"""
    if network_id not in networks:
        raise HTTPException(status_code=404, detail="Network not found")
    
    network = networks[network_id]
    
    # Simple flow inference based on hydraulic head differences
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
                })
            else:
                inferred_edges.append({
                    "source": edge.source,
                    "target": edge.target,
                    "head_difference": None,
                    "flow_probability": 0.5,
                    "flow_direction": "uncertain",
                })
    
    return {
        "network_id": network_id,
        "inferred_edges": inferred_edges,
        "method": "hydraulic_head_gradient",
    }


@router.delete("/{network_id}")
async def delete_network(network_id: str):
    """Delete a network"""
    if network_id not in networks:
        raise HTTPException(status_code=404, detail="Network not found")
    
    del networks[network_id]
    return {"message": "Network deleted successfully"}
