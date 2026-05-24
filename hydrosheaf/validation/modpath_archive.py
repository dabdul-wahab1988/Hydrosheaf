"""Inference pipeline and data standardisation for public MODPATH archives."""

from __future__ import annotations

import os
import tempfile
import yaml
import zipfile
from pathlib import Path
import pandas as pd

from ..physics.modpath import load_modpath_endpoint_records, load_modpath_pathline_points


def scan_modpath_archive(config_path: str) -> dict:
    """Scan and validate the public-archive config."""
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    # load schema
    schema_path = Path(config_path).parent / "archive_schema.yaml"
    if schema_path.exists():
        with open(schema_path, "r", encoding="utf-8") as f:
            schema = yaml.safe_load(f)
        required = schema.get("required_fields", [])
        for field in required:
            if field not in config:
                raise ValueError(f"Config is missing required field: {field}")

    return config


def load_modpath_endpoints(config: dict) -> pd.DataFrame:
    """Load MODPATH endpoints from configured zip file."""
    root = Path(config.get("local_archive_root", "."))
    zip_name = config.get("zip_file", "output.zip")
    zip_path = root / zip_name

    if not zip_path.exists():
        return pd.DataFrame(
            columns=[
                "particle_id",
                "x0",
                "y0",
                "z0",
                "x",
                "y",
                "z",
                "time",
                "initial_cell",
                "final_cell",
                "status",
            ]
        )

    file_in_zip = config.get("endpoint_file_in_zip", "")
    if file_in_zip.endswith(".mplst"):
        return pd.DataFrame(
            columns=[
                "particle_id",
                "x0",
                "y0",
                "z0",
                "x",
                "y",
                "z",
                "time",
                "initial_cell",
                "final_cell",
                "status",
            ]
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        if file_in_zip not in zf.namelist():
            return pd.DataFrame(
                columns=[
                    "particle_id",
                    "x0",
                    "y0",
                    "z0",
                    "x",
                    "y",
                    "z",
                    "time",
                    "initial_cell",
                    "final_cell",
                    "status",
                ]
            )
        with tempfile.TemporaryDirectory() as tmpdir:
            zf.extract(file_in_zip, tmpdir)
            extracted_path = Path(tmpdir) / file_in_zip
            records = load_modpath_endpoint_records(str(extracted_path))

    return pd.DataFrame(
        [
            {
                "particle_id": r.particle_id,
                "x0": r.x0,
                "y0": r.y0,
                "z0": r.z0,
                "x": r.x,
                "y": r.y,
                "z": r.z,
                "time": r.time,
                "initial_cell": r.initial_cell,
                "final_cell": r.final_cell,
                "status": r.status,
            }
            for r in records
        ]
    )


def load_modpath_pathlines(config: dict) -> pd.DataFrame:
    """Load MODPATH pathlines from configured zip file."""
    root = Path(config.get("local_archive_root", "."))
    zip_name = config.get("zip_file", "output.zip")
    zip_path = root / zip_name

    if not zip_path.exists():
        return pd.DataFrame(columns=["particle_id", "x", "y", "z", "time", "cell", "sequence"])

    file_in_zip = config.get("pathline_file_in_zip", "")
    if file_in_zip.endswith(".mplst"):
        return pd.DataFrame(columns=["particle_id", "x", "y", "z", "time", "cell", "sequence"])

    with zipfile.ZipFile(zip_path, "r") as zf:
        if file_in_zip not in zf.namelist():
            return pd.DataFrame(columns=["particle_id", "x", "y", "z", "time", "cell", "sequence"])
        with tempfile.TemporaryDirectory() as tmpdir:
            zf.extract(file_in_zip, tmpdir)
            extracted_path = Path(tmpdir) / file_in_zip
            points = load_modpath_pathline_points(str(extracted_path))

    return pd.DataFrame(
        [
            {
                "particle_id": r.particle_id,
                "x": r.x,
                "y": r.y,
                "z": r.z,
                "time": r.time,
                "cell": r.cell,
                "sequence": r.sequence,
            }
            for r in points
        ]
    )


def standardise_endpoint_columns(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Standardise and validate endpoint columns."""
    cols = ["particle_id", "x0", "y0", "z0", "x", "y", "z", "time", "initial_cell", "final_cell", "status"]
    for col in cols:
        if col not in df.columns:
            df[col] = None
    return df[cols]


def standardise_pathline_columns(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Standardise and validate pathline columns."""
    cols = ["particle_id", "x", "y", "z", "time", "cell", "sequence"]
    for col in cols:
        if col not in df.columns:
            df[col] = None
    return df[cols]


def build_node_mapping(
    endpoints: pd.DataFrame, pathlines: pd.DataFrame, config: dict
) -> pd.DataFrame:
    """Build coordinate details for each cell node found in the files."""
    nodes = {}

    for _, r in endpoints.iterrows():
        if r["initial_cell"] is not None and not pd.isna(r["initial_cell"]):
            c = f"cell_{int(r['initial_cell'])}"
            if c not in nodes:
                nodes[c] = {"node_id": c, "x": r["x0"], "y": r["y0"], "z": r["z0"]}
        if r["final_cell"] is not None and not pd.isna(r["final_cell"]):
            c = f"cell_{int(r['final_cell'])}"
            if c not in nodes:
                nodes[c] = {"node_id": c, "x": r["x"], "y": r["y"], "z": r["z"]}

    for _, r in pathlines.iterrows():
        if r["cell"] is not None and not pd.isna(r["cell"]):
            c = f"cell_{int(r['cell'])}"
            if c not in nodes:
                nodes[c] = {"node_id": c, "x": r["x"], "y": r["y"], "z": r["z"]}

    if not nodes:
        return pd.DataFrame(columns=["node_id", "x", "y", "z"])
    return pd.DataFrame(list(nodes.values()))


def build_reference_edges(
    endpoints: pd.DataFrame, pathlines: pd.DataFrame, node_map: pd.DataFrame, config: dict
) -> pd.DataFrame:
    """Build aggregated directed reference edges from endpoint connectivity.

    Endpoint source/target cells define the reference topology. Pathline records
    are used only to add supporting point counts for the same particles.
    """
    output_cols = [
        "archive_id",
        "edge_id",
        "u",
        "v",
        "source_node",
        "receptor_node",
        "particle_count",
        "endpoint_count",
        "pathline_point_count",
        "support_type",
        "median_travel_time",
        "min_travel_time",
        "max_travel_time",
        "time_unit",
        "travel_time_harmonised",
        "geometry_level",
        "evidence_level",
        "allowed_claim",
        "required_guardrail",
    ]
    pathline_counts = {}
    if not pathlines.empty and "particle_id" in pathlines.columns:
        pathline_counts = pathlines.groupby("particle_id").size().to_dict()

    edges = {}
    for _, r in endpoints.iterrows():
        init = r["initial_cell"]
        final = r["final_cell"]
        if init is not None and final is not None and not pd.isna(init) and not pd.isna(final):
            u = f"cell_{int(init)}"
            v = f"cell_{int(final)}"
            if u == v:
                continue
            key = (u, v)
            edges.setdefault(key, {"particle_ids": [], "times": []})
            if "particle_id" in r and not pd.isna(r["particle_id"]):
                edges[key]["particle_ids"].append(int(r["particle_id"]))
            if "time" in r and not pd.isna(r["time"]):
                edges[key]["times"].append(float(r["time"]))

    if not edges:
        return pd.DataFrame(columns=output_cols)

    rows = []
    for (u, v), data in sorted(edges.items()):
        particle_ids = data["particle_ids"]
        times = sorted(data["times"])
        rows.append(
            {
                "archive_id": config.get("archive_id", config.get("validation_tier", "")),
                "edge_id": f"{u}->{v}",
                "u": u,
                "v": v,
                "source_node": u,
                "receptor_node": v,
                "particle_count": len(particle_ids),
                "endpoint_count": len(particle_ids),
                "pathline_point_count": int(sum(pathline_counts.get(pid, 0) for pid in particle_ids)),
                "support_type": "endpoint_source_receptor_pair",
                "median_travel_time": float(pd.Series(times).median()) if times else None,
                "min_travel_time": times[0] if times else None,
                "max_travel_time": times[-1] if times else None,
                "time_unit": config.get("time_unit", config.get("time_unit_days", "")),
                "travel_time_harmonised": bool(config.get("travel_time_harmonised", False)),
                "geometry_level": "endpoint_projection_only",
                "evidence_level": "model_conditioned_endpoint_topology",
                "allowed_claim": config.get(
                    "allowed_claim",
                    "MODPATH endpoint connectivity provides directed reference topology for graph benchmarking",
                ),
                "required_guardrail": config.get(
                    "required_guardrail",
                    "Endpoint-derived edges are model-conditioned topology evidence, not absolute field truth or full pathline-shape validation",
                ),
            }
        )
    return pd.DataFrame(rows, columns=output_cols)


def summarise_archive_evidence(
    config: dict,
    endpoints: pd.DataFrame,
    pathlines: pd.DataFrame,
    reference_edges: pd.DataFrame,
) -> dict:
    """Summarise archive content metrics."""
    return {
        "archive_name": config.get("archive_name", ""),
        "validation_tier": config.get("validation_tier", ""),
        "n_particles": len(endpoints["particle_id"].unique()) if not endpoints.empty else 0,
        "n_endpoint_records": len(endpoints),
        "n_pathline_points": len(pathlines),
        "n_reference_edges": len(reference_edges),
        "source_doi": config.get("source_doi", ""),
        "source_url": config.get("source_url", ""),
    }
