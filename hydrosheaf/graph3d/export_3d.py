"""3D Process Mapping and Export for Hydrosheaf."""

import json
from typing import Any, Dict, List

def export_network_3d_json(
    nodes_info: List[Dict[str, Any]],
    edges_info: List[Dict[str, Any]],
    output_path: str
) -> None:
    """
    Export nodes and edges with 3D coordinates and reaction extents to a JSON format
    suitable for 3D visualization (e.g., in a web-based Three.js or D3-3D viewer).
    """
    data = {
        "nodes": nodes_info,
        "links": edges_info
    }
    
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

def generate_3d_vtp(
    nodes: List[Dict[str, Any]],
    edges: List[Dict[str, Any]],
    output_path: str
) -> None:
    """
    Generate a VTK PolyData (VTP) file for visualization in ParaView.
    """
    try:
        import vtk
    except ImportError:
        print("VTK not installed. Skipping VTP export.")
        return

    # Create points and cells
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    
    # Node attributes
    node_ids = {}
    for i, n in enumerate(nodes):
        x = n.get("x", 0.0)
        y = n.get("y", 0.0)
        z = n.get("z", 0.0)
        points.InsertNextPoint(x, y, z)
        node_ids[n["id"]] = i

    # Edge attributes
    reaction_data = vtk.vtkFloatArray()
    reaction_data.SetName("ReactionNorm")

    for e in edges:
        u_idx = node_ids.get(e["u"])
        v_idx = node_ids.get(e["v"])
        if u_idx is not None and v_idx is not None:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, u_idx)
            line.GetPointIds().SetId(1, v_idx)
            lines.InsertNextCell(line)
            reaction_data.InsertNextValue(e.get("residual_norm", 0.0))

    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points)
    poly_data.SetLines(lines)
    poly_data.GetCellData().AddArray(reaction_data)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    writer.SetInputData(poly_data)
    writer.Write()
