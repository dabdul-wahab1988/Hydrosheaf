"""Orchestrate validation runs for public MODPATH archives."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import ConvexHull, QhullError

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.physics.modpath import load_modpath_endpoint_records, load_modpath_pathline_points
from hydrosheaf.validation.modpath_archive import scan_modpath_archive

RANDOM_SEED = 20260521
BOOTSTRAP_ITERATIONS = 2000


def _safe_float(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def _bootstrap_ci(
    values: Sequence[float],
    statistic: Callable[[np.ndarray], float],
    *,
    iterations: int = BOOTSTRAP_ITERATIONS,
    random_seed: int = RANDOM_SEED,
) -> Dict[str, float]:
    arr = np.asarray([value for value in values if math.isfinite(float(value))], dtype=float)
    if arr.size == 0:
        return {"estimate": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"), "n": 0.0}
    rng = np.random.default_rng(random_seed)
    samples = np.empty(iterations, dtype=float)
    for idx in range(iterations):
        draw = rng.choice(arr, size=arr.size, replace=True)
        samples[idx] = statistic(draw)
    return {
        "estimate": float(statistic(arr)),
        "ci_low": float(np.nanpercentile(samples, 2.5)),
        "ci_high": float(np.nanpercentile(samples, 97.5)),
        "n": float(arr.size),
    }


def _compressed_cells(cells: Iterable[int]) -> List[int]:
    compressed: List[int] = []
    previous = None
    for cell in cells:
        if cell is None:
            continue
        cell_int = int(cell)
        if previous is None or cell_int != previous:
            compressed.append(cell_int)
        previous = cell_int
    return compressed


def _polygon_area(poly: Sequence[Sequence[float]]) -> float:
    if len(poly) < 3:
        return 0.0
    x = np.asarray([point[0] for point in poly], dtype=float)
    y = np.asarray([point[1] for point in poly], dtype=float)
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2.0)


def _inside(point: Sequence[float], edge_start: Sequence[float], edge_end: Sequence[float]) -> bool:
    return (edge_end[0] - edge_start[0]) * (point[1] - edge_start[1]) >= (
        edge_end[1] - edge_start[1]
    ) * (point[0] - edge_start[0])


def _intersection(
    line_start: Sequence[float],
    line_end: Sequence[float],
    edge_start: Sequence[float],
    edge_end: Sequence[float],
) -> List[float]:
    x1, y1 = line_start
    x2, y2 = line_end
    x3, y3 = edge_start
    x4, y4 = edge_end
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1.0e-12:
        return [float(x2), float(y2)]
    px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom
    py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom
    return [float(px), float(py)]


def _clip_polygon(subject: Sequence[Sequence[float]], clip: Sequence[Sequence[float]]) -> List[List[float]]:
    output = [list(point) for point in subject]
    if len(output) < 3 or len(clip) < 3:
        return []
    for idx, edge_start in enumerate(clip):
        edge_end = clip[(idx + 1) % len(clip)]
        input_list = output
        output = []
        if not input_list:
            break
        previous = input_list[-1]
        for current in input_list:
            if _inside(current, edge_start, edge_end):
                if not _inside(previous, edge_start, edge_end):
                    output.append(_intersection(previous, current, edge_start, edge_end))
                output.append(current)
            elif _inside(previous, edge_start, edge_end):
                output.append(_intersection(previous, current, edge_start, edge_end))
            previous = current
    return output


def _convex_hull(points: Sequence[Sequence[float]]) -> List[List[float]]:
    arr = np.asarray(points, dtype=float)
    arr = arr[np.isfinite(arr).all(axis=1)]
    if arr.shape[0] < 3:
        return []
    unique = np.unique(arr, axis=0)
    if unique.shape[0] < 3:
        return []
    try:
        hull = ConvexHull(unique)
    except QhullError:
        return []
    return unique[hull.vertices].tolist()


def _convex_iou(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> tuple[float, float, float, float]:
    area_a = _polygon_area(a)
    area_b = _polygon_area(b)
    if area_a <= 0 or area_b <= 0:
        return float("nan"), area_a, area_b, float("nan")
    inter = _clip_polygon(a, b)
    area_intersection = _polygon_area(inter)
    union = area_a + area_b - area_intersection
    return (area_intersection / union if union > 0 else float("nan"), area_a, area_b, area_intersection)


def compute_topology_agreement_and_particles(endpoints: pd.DataFrame, pathlines: pd.DataFrame):
    agreement_cols = [
        "edge_id", "source_node", "target_node", "endpoint_particle_count",
        "pathline_particle_count", "reverse_pathline_particle_count",
        "classification", "direction_agrees", "source_receptor_overlap", 
        "travel_time_mean", "travel_time_p10", "travel_time_p90"
    ]
    particles_cols = [
        "particle_id", "source_node", "target_node", "pathline_source_node",
        "pathline_target_node", "endpoint_edge_id", "pathline_edge_id",
        "endpoint_projection_preserved", "n_pathline_points", "pathline_start_time",
        "pathline_end_time", "pathline_elapsed_time", "start_x", "start_y", "end_x",
        "end_y"
    ]
    if endpoints.empty or pathlines.empty:
        return pd.DataFrame(columns=agreement_cols), pd.DataFrame(columns=particles_cols)
        
    path_by_particle = {}
    for r in pathlines.to_dict("records"):
        p_id = int(r["particle_id"])
        path_by_particle.setdefault(p_id, []).append(r)
        
    particle_rows = []
    endpoint_edges: Dict[tuple[str, str], List[int]] = {}
    pathline_edges: Dict[tuple[str, str], List[int]] = {}
    endpoint_times: Dict[tuple[str, str], List[float]] = {}
    endpoint_by_particle = {}
    for r in endpoints.to_dict("records"):
        p_id = int(r["particle_id"])
        endpoint_by_particle[p_id] = r
        
    for p_id, ep in endpoint_by_particle.items():
        init = ep["initial_cell"]
        final = ep["final_cell"]
        if init is None or final is None or pd.isna(init) or pd.isna(final):
            continue
        u = f"cell_{int(init)}"
        v = f"cell_{int(final)}"
        
        pts = path_by_particle.get(p_id, [])
        if not pts:
            continue
        pts_sorted = sorted(pts, key=lambda x: x["sequence"])

        first_cell = pts_sorted[0].get("cell")
        last_cell = pts_sorted[-1].get("cell")
        path_u = f"cell_{int(first_cell)}" if first_cell is not None and not pd.isna(first_cell) else ""
        path_v = f"cell_{int(last_cell)}" if last_cell is not None and not pd.isna(last_cell) else ""
        endpoint_edge = (u, v)
        pathline_edge = (path_u, path_v) if path_u and path_v else None
        endpoint_edges.setdefault(endpoint_edge, []).append(p_id)
        endpoint_times.setdefault(endpoint_edge, []).append(float(ep["time"]))
        if pathline_edge is not None:
            pathline_edges.setdefault(pathline_edge, []).append(p_id)
        
        start_time = float(pts_sorted[0]["time"])
        end_time = float(pts_sorted[-1]["time"])
        elapsed = end_time - start_time
        
        particle_rows.append({
            "particle_id": p_id,
            "source_node": u,
            "target_node": v,
            "pathline_source_node": path_u,
            "pathline_target_node": path_v,
            "endpoint_edge_id": f"{u}->{v}",
            "pathline_edge_id": f"{path_u}->{path_v}" if path_u and path_v else "",
            "endpoint_projection_preserved": bool(pathline_edge == endpoint_edge),
            "n_pathline_points": len(pts_sorted),
            "pathline_start_time": start_time,
            "pathline_end_time": end_time,
            "pathline_elapsed_time": elapsed,
            "start_x": float(ep["x0"]),
            "start_y": float(ep["y0"]),
            "end_x": float(ep["x"]),
            "end_y": float(ep["y"])
        })
        
    df_particles = pd.DataFrame(particle_rows)

    agreement_rows = []
    all_edges = sorted(set(endpoint_edges) | set(pathline_edges))
    for (u, v) in all_edges:
        edge_id = f"{u}->{v}"
        endpoint_pids = set(endpoint_edges.get((u, v), []))
        pathline_pids = set(pathline_edges.get((u, v), []))
        reverse_pathline_count = len(pathline_edges.get((v, u), []))
        if endpoint_pids and pathline_pids:
            classification = "TP"
        elif endpoint_pids:
            classification = "FN"
        else:
            classification = "FP"

        particle_union = endpoint_pids | pathline_pids
        particle_intersection = endpoint_pids & pathline_pids
        overlap = len(particle_intersection) / len(particle_union) if particle_union else float("nan")
        times = endpoint_times.get((u, v), [])
        times_sorted = sorted(times)
        n_times = len(times_sorted)
        tt_mean = sum(times_sorted) / n_times if n_times else float("nan")
        tt_p10 = times_sorted[int(0.1 * (n_times - 1))] if n_times else float("nan")
        tt_p90 = times_sorted[int(0.9 * (n_times - 1))] if n_times else float("nan")
        
        agreement_rows.append({
            "edge_id": edge_id,
            "source_node": u,
            "target_node": v,
            "endpoint_particle_count": len(endpoint_pids),
            "pathline_particle_count": len(pathline_pids),
            "reverse_pathline_particle_count": reverse_pathline_count,
            "classification": classification,
            "direction_agrees": bool(classification == "TP"),
            "source_receptor_overlap": overlap,
            "travel_time_mean": tt_mean,
            "travel_time_p10": tt_p10,
            "travel_time_p90": tt_p90
        })
        
    df_agreement = pd.DataFrame(agreement_rows)
    return df_agreement, df_particles


def compute_graph_priors(endpoints: pd.DataFrame, config: dict):
    priors_cols = [
        "edge_id", "source_node", "target_node", "particle_count", "p_uv", 
        "travel_time_mean", "travel_time_p10", "travel_time_p50", 
        "travel_time_p90", "travel_time_std", "source", "source_file"
    ]
    if endpoints.empty:
        return pd.DataFrame(columns=priors_cols)
        
    edges = {}
    counts_from_u = {}
    for r in endpoints.to_dict("records"):
        init = r["initial_cell"]
        final = r["final_cell"]
        if init is None or final is None or pd.isna(init) or pd.isna(final):
            continue
        u = f"cell_{int(init)}"
        v = f"cell_{int(final)}"
        edges.setdefault((u, v), []).append(float(r["time"]))
        counts_from_u[u] = counts_from_u.get(u, 0) + 1
        
    rows = []
    for (u, v), times in edges.items():
        times_sorted = sorted(times)
        n = len(times_sorted)
        mean = sum(times_sorted) / n
        p10 = times_sorted[int(0.1 * (n - 1))]
        p50 = times_sorted[int(0.5 * (n - 1))]
        p90 = times_sorted[int(0.9 * (n - 1))]
        var = sum((t - mean) ** 2 for t in times_sorted) / max(1, n - 1)
        std = math.sqrt(var)
        p_uv = float(n) / float(counts_from_u[u])
        
        rows.append({
            "edge_id": f"{u}->{v}",
            "source_node": u,
            "target_node": v,
            "particle_count": n,
            "p_uv": p_uv,
            "travel_time_mean": mean,
            "travel_time_p10": p10,
            "travel_time_p50": p50,
            "travel_time_p90": p90,
            "travel_time_std": std,
            "source": f"MODPATH {config.get('modpath_version', 5)} endpoint file",
            "source_file": config.get("endpoint_file_in_zip", "")
        })
        
    return pd.DataFrame(rows)


def compute_topology_summary(df_agreement: pd.DataFrame, df_particles: pd.DataFrame):
    summary_cols = [
        "n_endpoint_edges", "n_pathline_edges", "true_positive_edges", 
        "false_positive_edges", "false_negative_edges", "edge_precision", 
        "edge_recall", "edge_f1", "direction_agreement_rate", 
        "mean_source_receptor_overlap", "pathline_elapsed_time_median", 
        "pathline_elapsed_time_p90", "edge_error_rate"
    ]
    if df_agreement.empty or df_particles.empty:
        return pd.DataFrame(columns=summary_cols)
        
    endpoint_mask = pd.to_numeric(df_agreement["endpoint_particle_count"], errors="coerce").fillna(0).gt(0)
    pathline_mask = pd.to_numeric(df_agreement["pathline_particle_count"], errors="coerce").fillna(0).gt(0)
    tp = int((df_agreement["classification"] == "TP").sum())
    fp = int((df_agreement["classification"] == "FP").sum())
    fn = int((df_agreement["classification"] == "FN").sum())
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    times = df_particles["pathline_elapsed_time"].dropna().tolist()
    pt_med = float(np.median(times)) if times else float("nan")
    pt_p90 = float(np.percentile(times, 90)) if times else float("nan")
    
    row = {
        "n_endpoint_edges": int(endpoint_mask.sum()),
        "n_pathline_edges": int(pathline_mask.sum()),
        "true_positive_edges": tp,
        "false_positive_edges": fp,
        "false_negative_edges": fn,
        "edge_precision": precision,
        "edge_recall": recall,
        "edge_f1": f1,
        "direction_agreement_rate": tp / int(endpoint_mask.sum()) if int(endpoint_mask.sum()) else float("nan"),
        "mean_source_receptor_overlap": float(df_agreement["source_receptor_overlap"].dropna().mean()),
        "pathline_elapsed_time_median": pt_med,
        "pathline_elapsed_time_p90": pt_p90,
        "edge_error_rate": (fp + fn) / len(df_agreement) if len(df_agreement) else float("nan")
    }
    return pd.DataFrame([row])


def build_pathline_structure(endpoints_path: Path, pathline_path: Path) -> pd.DataFrame:
    endpoints = load_modpath_endpoint_records(str(endpoints_path))
    pathline_points = load_modpath_pathline_points(str(pathline_path))
    endpoint_by_particle = {int(record.particle_id): record for record in endpoints}

    by_particle: Dict[int, List[object]] = {}
    for point in pathline_points:
        by_particle.setdefault(int(point.particle_id), []).append(point)

    rows: List[Dict[str, object]] = []
    for particle_id, points in sorted(by_particle.items()):
        points_sorted = sorted(points, key=lambda point: point.sequence)
        cells = _compressed_cells(point.cell for point in points_sorted if point.cell is not None)
        if not cells:
            continue
        endpoint = endpoint_by_particle.get(particle_id)
        endpoint_source = int(endpoint.initial_cell) if endpoint and endpoint.initial_cell is not None else None
        endpoint_target = int(endpoint.final_cell) if endpoint and endpoint.final_cell is not None else None
        start_matches = endpoint_source is not None and cells[0] == endpoint_source
        end_matches = endpoint_target is not None and cells[-1] == endpoint_target
        start_time = _safe_float(points_sorted[0].time)
        end_time = _safe_float(points_sorted[-1].time)
        endpoint_time = _safe_float(endpoint.time) if endpoint else float("nan")
        rows.append(
            {
                "particle_id": particle_id,
                "endpoint_edge_id": f"cell_{endpoint_source}->cell_{endpoint_target}"
                if endpoint_source is not None and endpoint_target is not None
                else "",
                "endpoint_source_node": f"cell_{endpoint_source}" if endpoint_source is not None else "",
                "endpoint_target_node": f"cell_{endpoint_target}" if endpoint_target is not None else "",
                "pathline_start_cell": cells[0],
                "pathline_end_cell": cells[-1],
                "pathline_start_node": f"cell_{cells[0]}",
                "pathline_end_node": f"cell_{cells[-1]}",
                "endpoint_start_x": endpoint.x0 if endpoint else float("nan"),
                "endpoint_start_y": endpoint.y0 if endpoint else float("nan"),
                "endpoint_end_x": endpoint.x if endpoint else float("nan"),
                "endpoint_end_y": endpoint.y if endpoint else float("nan"),
                "pathline_start_x": points_sorted[0].x,
                "pathline_start_y": points_sorted[0].y,
                "pathline_end_x": points_sorted[-1].x,
                "pathline_end_y": points_sorted[-1].y,
                "start_matches_endpoint": start_matches,
                "end_matches_endpoint": end_matches,
                "endpoint_projection_preserved": bool(start_matches and end_matches),
                "n_pathline_points": len(points_sorted),
                "n_compressed_cells": len(cells),
                "n_unique_cells": len(set(cells)),
                "n_cell_transitions": max(len(cells) - 1, 0),
                "has_cell_revisit": len(cells) > len(set(cells)),
                "pathline_start_time": start_time,
                "pathline_end_time": end_time,
                "pathline_elapsed_time": end_time - start_time if math.isfinite(start_time) and math.isfinite(end_time) else float("nan"),
                "endpoint_time": endpoint_time,
            }
        )
    return pd.DataFrame(rows)


def build_capture_envelope_overlap(endpoints_path: Path, pathline_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    endpoints = load_modpath_endpoint_records(str(endpoints_path))
    pathline_points = load_modpath_pathline_points(str(pathline_path))
    first_pathline_by_particle: Dict[int, object] = {}
    for point in sorted(pathline_points, key=lambda item: (int(item.particle_id), int(item.sequence))):
        first_pathline_by_particle.setdefault(int(point.particle_id), point)

    rows: List[Dict[str, object]] = []
    for target_cell in sorted({record.final_cell for record in endpoints if record.final_cell is not None}):
        target_particles = [record for record in endpoints if record.final_cell == target_cell]
        endpoint_points = [[record.x0, record.y0] for record in target_particles]
        pathline_points_xy = []
        endpoint_sources = set()
        pathline_sources = set()
        for record in target_particles:
            endpoint_sources.add(int(record.initial_cell))
            point = first_pathline_by_particle.get(int(record.particle_id))
            if point is not None:
                pathline_points_xy.append([point.x, point.y])
                if point.cell is not None:
                    pathline_sources.add(int(point.cell))
        endpoint_hull = _convex_hull(endpoint_points)
        pathline_hull = _convex_hull(pathline_points_xy)
        iou, endpoint_area, pathline_area, intersection_area = _convex_iou(endpoint_hull, pathline_hull)
        source_union = endpoint_sources | pathline_sources
        source_intersection = endpoint_sources & pathline_sources
        rows.append(
            {
                "target_node": f"cell_{int(target_cell)}",
                "n_particles": len(target_particles),
                "n_endpoint_source_cells": len(endpoint_sources),
                "n_pathline_source_cells": len(pathline_sources),
                "source_cell_jaccard": len(source_intersection) / len(source_union) if source_union else float("nan"),
                "endpoint_hull_area": endpoint_area,
                "pathline_hull_area": pathline_area,
                "intersection_area": intersection_area,
                "capture_envelope_iou": iou,
                "metric_scope": "convex hull of particle release/start points by receptor cell",
                "guardrail": "This is a capture-envelope overlap from MODPATH point clouds, not a model-supplied capture-zone polygon or raster.",
            }
        )
    detail = pd.DataFrame(rows)
    summary = pd.DataFrame(
        [
            {
                "n_receptor_targets": len(detail),
                "particle_weighted_capture_envelope_iou": float(
                    np.average(detail["capture_envelope_iou"], weights=detail["n_particles"])
                )
                if not detail.empty
                else float("nan"),
                "median_capture_envelope_iou": float(detail["capture_envelope_iou"].median()) if not detail.empty else float("nan"),
                "median_source_cell_jaccard": float(detail["source_cell_jaccard"].median()) if not detail.empty else float("nan"),
                "guardrail": "Point-cloud capture envelopes are plotted because no explicit capture-zone polygons/rasters were present in the prepared archive outputs.",
            }
        ]
    )
    return detail, summary


def build_travel_time_rank(agreement: pd.DataFrame, particles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if agreement.empty or particles.empty:
        return pd.DataFrame(), pd.DataFrame()
    particle_edge = particles.copy()
    particle_edge["edge_id"] = particle_edge["source_node"].astype(str) + "->" + particle_edge["target_node"].astype(str)
    grouped = (
        particle_edge.groupby("edge_id", as_index=False)
        .agg(
            pathline_elapsed_time_median=("pathline_elapsed_time", "median"),
            pathline_elapsed_time_p90=("pathline_elapsed_time", lambda values: float(np.nanpercentile(values, 90))),
            pathline_particle_count=("particle_id", "count"),
            median_pathline_points=("n_pathline_points", "median"),
        )
    )
    cols = [
        "edge_id",
        "source_node",
        "target_node",
        "endpoint_particle_count",
        "pathline_particle_count",
        "travel_time_mean",
        "travel_time_p10",
        "travel_time_p90",
        "classification",
        "direction_agrees",
        "source_receptor_overlap",
    ]
    merged = agreement[[col for col in cols if col in agreement.columns]].merge(grouped, on="edge_id", how="inner")
    merged = merged.rename(columns={"travel_time_mean": "endpoint_travel_time_mean"})
    valid = (
        pd.to_numeric(merged["endpoint_travel_time_mean"], errors="coerce").gt(0)
        & pd.to_numeric(merged["pathline_elapsed_time_median"], errors="coerce").gt(0)
    )
    if valid.any():
        endpoint_times = pd.to_numeric(merged.loc[valid, "endpoint_travel_time_mean"], errors="coerce")
        pathline_times = pd.to_numeric(merged.loc[valid, "pathline_elapsed_time_median"], errors="coerce")
        spearman = stats.spearmanr(endpoint_times, pathline_times, nan_policy="omit")
        kendall = stats.kendalltau(endpoint_times, pathline_times, nan_policy="omit")
        summary = pd.DataFrame(
            [
                {
                    "n_edges": int(valid.sum()),
                    "spearman_rho": float(spearman.statistic),
                    "spearman_p": float(spearman.pvalue),
                    "kendall_tau": float(kendall.statistic),
                    "kendall_p": float(kendall.pvalue),
                    "guardrail": "Endpoint and pathline time columns are compared as rank-order evidence unless archive time units are harmonised.",
                }
            ]
        )
    else:
        summary = pd.DataFrame()
    return merged, summary


def build_harmonized_travel_time(structure: pd.DataFrame, graph_priors: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if structure.empty or graph_priors.empty:
        return pd.DataFrame(), pd.DataFrame()
    endpoint = (
        structure.groupby("endpoint_edge_id", as_index=False)
        .agg(
            particle_endpoint_time_median=("endpoint_time", "median"),
            particle_endpoint_time_mean=("endpoint_time", "mean"),
            particle_endpoint_time_p10=("endpoint_time", lambda values: float(np.nanpercentile(values, 10))),
            particle_endpoint_time_p90=("endpoint_time", lambda values: float(np.nanpercentile(values, 90))),
            n_particles=("particle_id", "count"),
        )
        .rename(columns={"endpoint_edge_id": "edge_id"})
    )
    priors = graph_priors.rename(columns={"travel_time_mean": "hydrosheaf_modpath_weight_mean"}).copy()
    merged = endpoint.merge(
        priors[["edge_id", "hydrosheaf_modpath_weight_mean", "travel_time_p10", "travel_time_p90", "particle_count"]],
        on="edge_id",
        how="inner",
    )
    valid = (
        pd.to_numeric(merged["particle_endpoint_time_median"], errors="coerce").gt(0)
        & pd.to_numeric(merged["hydrosheaf_modpath_weight_mean"], errors="coerce").gt(0)
    )
    if valid.any():
        x = pd.to_numeric(merged.loc[valid, "particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(merged.loc[valid, "hydrosheaf_modpath_weight_mean"], errors="coerce")
        spearman = stats.spearmanr(x, y, nan_policy="omit")
        kendall = stats.kendalltau(x, y, nan_policy="omit")
        log_error = np.log10(y) - np.log10(x)
        summary = pd.DataFrame(
            [
                {
                    "n_edges": int(valid.sum()),
                    "spearman_rho": float(spearman.statistic),
                    "spearman_p": float(spearman.pvalue),
                    "kendall_tau": float(kendall.statistic),
                    "kendall_p": float(kendall.pvalue),
                    "median_abs_log10_difference": float(np.nanmedian(np.abs(log_error))),
                    "guardrail": "This compares endpoint total-time particles with MODPATH-derived Hydrosheaf edge weights; it validates prior transfer, not independent travel-time prediction.",
                }
            ]
        )
    else:
        summary = pd.DataFrame()
    return merged, summary


def build_bootstrap_ci(summary: pd.DataFrame, agreement: pd.DataFrame, particles: pd.DataFrame, structure: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    if not summary.empty:
        top = summary.iloc[0]
        for metric in ("edge_precision", "edge_recall", "edge_f1", "direction_agreement_rate", "mean_source_receptor_overlap"):
            rows.append(
                {
                    "metric": metric,
                    "estimate": _safe_float(top.get(metric)),
                    "ci_low": _safe_float(top.get(metric)),
                    "ci_high": _safe_float(top.get(metric)),
                    "n": _safe_float(top.get("n_pathline_edges")),
                    "method": "archive aggregate; no resampling needed for deterministic all-edge summary",
                }
            )

    if not agreement.empty and "source_receptor_overlap" in agreement.columns:
        ci = _bootstrap_ci(pd.to_numeric(agreement["source_receptor_overlap"], errors="coerce").dropna(), np.nanmean)
        rows.append({"metric": "edge_source_receptor_overlap_mean_bootstrap", **ci, "method": "edge bootstrap"})
    if not particles.empty and "pathline_elapsed_time" in particles.columns:
        elapsed = pd.to_numeric(particles["pathline_elapsed_time"], errors="coerce").dropna()
        rows.append({"metric": "particle_elapsed_time_median_bootstrap", **_bootstrap_ci(elapsed, np.nanmedian), "method": "particle bootstrap"})
        rows.append(
            {
                "metric": "particle_elapsed_time_p90_bootstrap",
                **_bootstrap_ci(elapsed, lambda values: float(np.nanpercentile(values, 90))),
                "method": "particle bootstrap",
            }
        )
    if not structure.empty:
        projection = structure["endpoint_projection_preserved"].astype(float)
        rows.append({"metric": "endpoint_projection_preservation_rate_bootstrap", **_bootstrap_ci(projection, np.nanmean), "method": "particle bootstrap"})
        rows.append(
            {
                "metric": "compressed_path_cells_median_bootstrap",
                **_bootstrap_ci(pd.to_numeric(structure["n_compressed_cells"], errors="coerce").dropna(), np.nanmedian),
                "method": "particle bootstrap",
            }
        )
    return pd.DataFrame(rows)


def build_archive_summary(
    config: dict,
    summary: pd.DataFrame,
    agreement: pd.DataFrame,
    particles: pd.DataFrame,
    structure: pd.DataFrame,
    rank_summary: pd.DataFrame,
    capture_summary: pd.DataFrame,
    harmonized_summary: pd.DataFrame,
) -> pd.DataFrame:
    top = summary.iloc[0].to_dict() if not summary.empty else {}
    projection_rate = float(structure["endpoint_projection_preserved"].mean()) if not structure.empty else float("nan")
    median_cells = float(structure["n_compressed_cells"].median()) if not structure.empty else float("nan")
    median_points = float(structure["n_pathline_points"].median()) if not structure.empty else float("nan")
    rank = rank_summary.iloc[0].to_dict() if not rank_summary.empty else {}
    capture = capture_summary.iloc[0].to_dict() if not capture_summary.empty else {}
    harmonized = harmonized_summary.iloc[0].to_dict() if not harmonized_summary.empty else {}
    row = {
        "validation_tier": config.get("validation_tier", ""),
        "archive_name": config.get("archive_name", ""),
        "source_doi": config.get("source_doi", ""),
        "source_url": config.get("source_url", ""),
        "n_endpoint_records": len(endpoints_loaded_global) if "endpoints_loaded_global" in globals() else len(particles),
        "n_pathline_points": len(pathlines_loaded_global) if "pathlines_loaded_global" in globals() else np.nan,
        "n_particles": len(structure),
        "n_endpoint_edges": top.get("n_endpoint_edges", np.nan),
        "n_pathline_edges": top.get("n_pathline_edges", np.nan),
        "true_positive_edges": top.get("true_positive_edges", np.nan),
        "false_positive_edges": top.get("false_positive_edges", np.nan),
        "false_negative_edges": top.get("false_negative_edges", np.nan),
        "edge_precision": top.get("edge_precision", np.nan),
        "edge_recall": top.get("edge_recall", np.nan),
        "edge_f1": top.get("edge_f1", np.nan),
        "direction_agreement_rate": top.get("direction_agreement_rate", np.nan),
        "mean_source_receptor_overlap": top.get("mean_source_receptor_overlap", np.nan),
        "endpoint_projection_preservation_rate": projection_rate,
        "median_compressed_path_cells": median_cells,
        "median_pathline_points_per_particle": median_points,
        "pathline_elapsed_time_median": top.get("pathline_elapsed_time_median", np.nan),
        "pathline_elapsed_time_p90": top.get("pathline_elapsed_time_p90", np.nan),
        "travel_time_spearman_rho": rank.get("spearman_rho", np.nan),
        "travel_time_kendall_tau": rank.get("kendall_tau", np.nan),
        "travel_time_rank_supported": bool(_safe_float(rank.get("spearman_rho", np.nan)) > 0.5),
        "capture_envelope_iou": capture.get("particle_weighted_capture_envelope_iou", np.nan),
        "source_cell_jaccard": capture.get("median_source_cell_jaccard", np.nan),
        "harmonized_travel_time_spearman_rho": harmonized.get("spearman_rho", np.nan),
        "harmonized_travel_time_kendall_tau": harmonized.get("kendall_tau", np.nan),
        "harmonized_travel_time_median_abs_log10_difference": harmonized.get("median_abs_log10_difference", np.nan),
        "travel_time_rank_interpretation": (
            "Travel-time rank is supportive only when endpoint and pathline time definitions are harmonised "
            "and the rank statistic is positive; current value should be treated as a diagnostic limitation."
        ),
        "claim_guardrail": (
            "External archive validates MODPATH endpoint/pathline topology consistency as diagnostic/prior-assisted evidence; "
            "it does not validate field geochemistry, field truth, capture-zone polygons, or independent travel-time prediction, "
            "and does not constitute independent Hydrosheaf validation."
        ),
    }
    return pd.DataFrame([row])


def write_empty_stubs(config: dict, result_dir: Path):
    tier_prefix = config.get("validation_tier", "")
    
    # 1. Edge Agreement
    agreement_cols = [
        "edge_id", "source_node", "target_node", "endpoint_particle_count", 
        "pathline_particle_count", "reverse_pathline_particle_count", 
        "classification", "direction_agrees", "source_receptor_overlap", 
        "travel_time_mean", "travel_time_p10", "travel_time_p90"
    ]
    pd.DataFrame(columns=agreement_cols).to_csv(result_dir / f"{tier_prefix}_edge_agreement.csv", index=False)
    
    # 2. Pathline Particles
    particles_cols = [
        "particle_id", "source_node", "target_node", "n_pathline_points", 
        "pathline_start_time", "pathline_end_time", "pathline_elapsed_time", 
        "start_x", "start_y", "end_x", "end_y"
    ]
    pd.DataFrame(columns=particles_cols).to_csv(result_dir / f"{tier_prefix}_pathline_particles.csv", index=False)
    
    # 3. Pathline Structure
    structure_cols = [
        "particle_id", "endpoint_edge_id", "endpoint_source_node", "endpoint_target_node", 
        "pathline_start_cell", "pathline_end_cell", "pathline_start_node", "pathline_end_node", 
        "endpoint_start_x", "endpoint_start_y", "endpoint_end_x", "endpoint_end_y", 
        "pathline_start_x", "pathline_start_y", "pathline_end_x", "pathline_end_y", 
        "start_matches_endpoint", "end_matches_endpoint", "endpoint_projection_preserved", 
        "n_pathline_points", "n_compressed_cells", "n_unique_cells", "n_cell_transitions", 
        "has_cell_revisit", "pathline_start_time", "pathline_end_time", "pathline_elapsed_time", "endpoint_time"
    ]
    pd.DataFrame(columns=structure_cols).to_csv(result_dir / f"{tier_prefix}_pathline_structure.csv", index=False)
    
    # 4. Capture Envelope Overlap
    capture_cols = [
        "target_node", "n_particles", "n_endpoint_source_cells", "n_pathline_source_cells", 
        "source_cell_jaccard", "endpoint_hull_area", "pathline_hull_area", "intersection_area", 
        "capture_envelope_iou", "metric_scope", "guardrail"
    ]
    pd.DataFrame(columns=capture_cols).to_csv(result_dir / f"{tier_prefix}_capture_envelope_overlap.csv", index=False)
    
    # 5. Capture Envelope Summary
    capture_sum_cols = [
        "n_receptor_targets", "particle_weighted_capture_envelope_iou", "median_capture_envelope_iou", 
        "median_source_cell_jaccard", "guardrail"
    ]
    capture_sum_row = {
        "n_receptor_targets": 0,
        "particle_weighted_capture_envelope_iou": float("nan"),
        "median_capture_envelope_iou": float("nan"),
        "median_source_cell_jaccard": float("nan"),
        "guardrail": "Fallback stub mode."
    }
    pd.DataFrame([capture_sum_row], columns=capture_sum_cols).to_csv(result_dir / f"{tier_prefix}_capture_envelope_summary.csv", index=False)
    
    # 6. Travel Time Rank
    rank_cols = [
        "edge_id", "source_node", "target_node", "endpoint_particle_count", 
        "pathline_particle_count", "endpoint_travel_time_mean", "travel_time_p10", 
        "travel_time_p90", "classification", "direction_agrees", "source_receptor_overlap", 
        "pathline_elapsed_time_median", "pathline_elapsed_time_p90", "median_pathline_points"
    ]
    pd.DataFrame(columns=rank_cols).to_csv(result_dir / f"{tier_prefix}_travel_time_rank.csv", index=False)
    
    # 7. Travel Time Rank Summary
    rank_sum_cols = ["n_edges", "spearman_rho", "spearman_p", "kendall_tau", "kendall_p", "guardrail"]
    rank_sum_row = {
        "n_edges": 0,
        "spearman_rho": float("nan"),
        "spearman_p": float("nan"),
        "kendall_tau": float("nan"),
        "kendall_p": float("nan"),
        "guardrail": "Fallback stub mode."
    }
    pd.DataFrame([rank_sum_row], columns=rank_sum_cols).to_csv(result_dir / f"{tier_prefix}_travel_time_rank_summary.csv", index=False)
    
    # 8. Harmonized Travel Time
    harmonized_cols = [
        "edge_id", "particle_endpoint_time_median", "particle_endpoint_time_mean", 
        "particle_endpoint_time_p10", "particle_endpoint_time_p90", "n_particles", 
        "hydrosheaf_modpath_weight_mean", "travel_time_p10", "travel_time_p90", "particle_count"
    ]
    pd.DataFrame(columns=harmonized_cols).to_csv(result_dir / f"{tier_prefix}_harmonized_travel_time.csv", index=False)
    
    # 9. Harmonized Travel Time Summary
    harmonized_sum_cols = ["n_edges", "spearman_rho", "spearman_p", "kendall_tau", "kendall_p", "median_abs_log10_difference", "guardrail"]
    harmonized_sum_row = {
        "n_edges": 0,
        "spearman_rho": float("nan"),
        "spearman_p": float("nan"),
        "kendall_tau": float("nan"),
        "kendall_p": float("nan"),
        "median_abs_log10_difference": float("nan"),
        "guardrail": "Fallback stub mode."
    }
    pd.DataFrame([harmonized_sum_row], columns=harmonized_sum_cols).to_csv(result_dir / f"{tier_prefix}_harmonized_travel_time_summary.csv", index=False)
    
    # 10. Bootstrap CI
    boot_cols = ["metric", "estimate", "ci_low", "ci_high", "n", "method"]
    boot_rows = []
    for metric in ("edge_precision", "edge_recall", "edge_f1", "direction_agreement_rate", "mean_source_receptor_overlap"):
        boot_rows.append({
            "metric": metric,
            "estimate": float("nan"),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "n": 0,
            "method": "fallback stub"
        })
    pd.DataFrame(boot_rows, columns=boot_cols).to_csv(result_dir / f"{tier_prefix}_bootstrap_ci.csv", index=False)
    
    # 11. Archive Summary
    summary_cols = [
        "validation_tier", "archive_name", "source_doi", "source_url", "n_endpoint_records", 
        "n_pathline_points", "n_particles", "n_endpoint_edges", "n_pathline_edges", 
        "true_positive_edges", "false_positive_edges", "false_negative_edges", "edge_precision", 
        "edge_recall", "edge_f1", "direction_agreement_rate", "mean_source_receptor_overlap", 
        "endpoint_projection_preservation_rate", "median_compressed_path_cells", "median_pathline_points_per_particle", 
        "pathline_elapsed_time_median", "pathline_elapsed_time_p90", "travel_time_spearman_rho", 
        "travel_time_kendall_tau", "travel_time_rank_supported", "capture_envelope_iou", "source_cell_jaccard", 
        "harmonized_travel_time_spearman_rho", "harmonized_travel_time_kendall_tau", 
        "harmonized_travel_time_median_abs_log10_difference", "travel_time_rank_interpretation", "claim_guardrail"
    ]
    summary_row = {
        "validation_tier": config.get("validation_tier", ""),
        "archive_name": config.get("archive_name", ""),
        "source_doi": config.get("source_doi", ""),
        "source_url": config.get("source_url", ""),
        "n_endpoint_records": 0,
        "n_pathline_points": 0,
        "n_particles": 0,
        "n_endpoint_edges": 0,
        "n_pathline_edges": 0,
        "true_positive_edges": 0,
        "false_positive_edges": 0,
        "false_negative_edges": 0,
        "edge_precision": float("nan"),
        "edge_recall": float("nan"),
        "edge_f1": float("nan"),
        "direction_agreement_rate": float("nan"),
        "mean_source_receptor_overlap": float("nan"),
        "endpoint_projection_preservation_rate": float("nan"),
        "median_compressed_path_cells": float("nan"),
        "median_pathline_points_per_particle": float("nan"),
        "pathline_elapsed_time_median": float("nan"),
        "pathline_elapsed_time_p90": float("nan"),
        "travel_time_spearman_rho": float("nan"),
        "travel_time_kendall_tau": float("nan"),
        "travel_time_rank_supported": False,
        "capture_envelope_iou": float("nan"),
        "source_cell_jaccard": float("nan"),
        "harmonized_travel_time_spearman_rho": float("nan"),
        "harmonized_travel_time_kendall_tau": float("nan"),
        "harmonized_travel_time_median_abs_log10_difference": float("nan"),
        "travel_time_rank_interpretation": "No data available in fallback/stub mode.",
        "claim_guardrail": "Fallback stub for scalability/supplementary archive. No active validation results."
    }
    pd.DataFrame([summary_row], columns=summary_cols).to_csv(result_dir / f"{tier_prefix}_archive_summary.csv", index=False)
    print(f"Wrote empty stubs for {tier_prefix}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MODPATH-to-Graph validation pipeline for a public archive.")
    parser.add_argument("--config", type=str, required=True, help="Path to config YAML file.")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path

    # Scan and validate config
    config = scan_modpath_archive(str(config_path))
    tier = config.get("validation_tier", "")
    local_root = PROJECT_ROOT / config.get("local_archive_root", "")
    zip_name = config.get("zip_file", "output.zip")
    zip_path = local_root / zip_name
    endpoint_file = config.get("endpoint_file_in_zip", "")
    pathline_file = config.get("pathline_file_in_zip", "")

    result_dir = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"
    result_dir.mkdir(parents=True, exist_ok=True)

    # Determine if we are running in fallback-stub mode
    is_stub = False
    if tier == "tier_3_long_island" or endpoint_file.endswith(".mplst") or not zip_path.exists():
        is_stub = True

    if is_stub:
        write_empty_stubs(config, result_dir)
        return

    # Extract files and load DataFrames
    with zipfile.ZipFile(zip_path, "r") as zf:
        if endpoint_file not in zf.namelist() or pathline_file not in zf.namelist():
            print(f"[WARNING] Archive {zip_name} does not contain configured files. Falling back to stubs.")
            write_empty_stubs(config, result_dir)
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            zf.extract(endpoint_file, tmpdir)
            zf.extract(pathline_file, tmpdir)
            
            ep_local = Path(tmpdir) / endpoint_file
            pl_local = Path(tmpdir) / pathline_file
            
            # Load DataFrames
            endpoints = load_modpath_endpoint_records(str(ep_local))
            pathlines = load_modpath_pathline_points(str(pl_local))
            
            endpoints_df = pd.DataFrame(
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
                    for r in endpoints
                ]
            )
            pathlines_df = pd.DataFrame(
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
                    for r in pathlines
                ]
            )
            
            # Set globals so build_archive_summary can grab length of lists directly
            global endpoints_loaded_global, pathlines_loaded_global
            endpoints_loaded_global = endpoints
            pathlines_loaded_global = pathlines

            # Compute topology agreement, particles, graph priors
            agreement, particles = compute_topology_agreement_and_particles(endpoints_df, pathlines_df)
            graph_priors = compute_graph_priors(endpoints_df, config)
            summary = compute_topology_summary(agreement, particles)
            
            # Detailed analyses
            structure = build_pathline_structure(ep_local, pl_local)
            capture_detail, capture_summary = build_capture_envelope_overlap(ep_local, pl_local)
            travel_rank, travel_rank_summary = build_travel_time_rank(agreement, particles)
            harmonized_time, harmonized_time_summary = build_harmonized_travel_time(structure, graph_priors)
            bootstrap = build_bootstrap_ci(summary, agreement, particles, structure)
            archive_summary = build_archive_summary(
                config,
                summary,
                agreement,
                particles,
                structure,
                travel_rank_summary,
                capture_summary,
                harmonized_time_summary,
            )

            # Write config-prefixed files
            tier_prefix = config.get("validation_tier", "")
            agreement.to_csv(result_dir / f"{tier_prefix}_edge_agreement.csv", index=False)
            particles.to_csv(result_dir / f"{tier_prefix}_pathline_particles.csv", index=False)
            structure.to_csv(result_dir / f"{tier_prefix}_pathline_structure.csv", index=False)
            capture_detail.to_csv(result_dir / f"{tier_prefix}_capture_envelope_overlap.csv", index=False)
            capture_summary.to_csv(result_dir / f"{tier_prefix}_capture_envelope_summary.csv", index=False)
            travel_rank.to_csv(result_dir / f"{tier_prefix}_travel_time_rank.csv", index=False)
            travel_rank_summary.to_csv(result_dir / f"{tier_prefix}_travel_time_rank_summary.csv", index=False)
            harmonized_time.to_csv(result_dir / f"{tier_prefix}_harmonized_travel_time.csv", index=False)
            harmonized_time_summary.to_csv(result_dir / f"{tier_prefix}_harmonized_travel_time_summary.csv", index=False)
            bootstrap.to_csv(result_dir / f"{tier_prefix}_bootstrap_ci.csv", index=False)
            archive_summary.to_csv(result_dir / f"{tier_prefix}_archive_summary.csv", index=False)
            
            print(f"Successfully processed and wrote outputs for {tier_prefix}: {len(structure)} particles.")


if __name__ == "__main__":
    main()
