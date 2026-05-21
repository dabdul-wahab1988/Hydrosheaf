"""Ingest the M2 external MODPATH archive validation into the M4 evidence package."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import ConvexHull, QhullError


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.physics.modpath import load_modpath_endpoint_records, load_modpath_pathline_points


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
M2_MODPATH_ROOT = REPO_ROOT / "M2" / "m2_benchmark" / "external" / "modpath"
M2_MODPATH_INPUT = M2_MODPATH_ROOT / "input" / "selected_output"
M2_MODPATH_RESULTS = M2_MODPATH_ROOT / "results"
RESULT_DIR = BENCHMARK_ROOT / "results"
RANDOM_SEED = 20260521
BOOTSTRAP_ITERATIONS = 2000


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _edge_id(source: object, target: object) -> str:
    return f"cell_{int(source)}->cell_{int(target)}"


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


def build_pathline_structure(endpoint_file: Path, pathline_file: Path) -> pd.DataFrame:
    endpoints = load_modpath_endpoint_records(str(endpoint_file))
    pathline_points = load_modpath_pathline_points(str(pathline_file))
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
                "endpoint_edge_id": _edge_id(endpoint_source, endpoint_target)
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


def build_capture_envelope_overlap(endpoint_file: Path, pathline_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    endpoints = load_modpath_endpoint_records(str(endpoint_file))
    pathline_points = load_modpath_pathline_points(str(pathline_file))
    endpoint_by_particle = {int(record.particle_id): record for record in endpoints}
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
    manifest: Mapping[str, object],
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
        "validation_tier": "external_modpath_archive",
        "archive_name": "USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5",
        "source_doi": manifest.get("source_doi", "10.5066/F7J102FK"),
        "source_url": manifest.get(
            "source_url",
            "https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute",
        ),
        "n_endpoint_records": manifest.get("n_endpoint_records", len(particles)),
        "n_pathline_points": manifest.get("n_pathline_points", np.nan),
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
            "External archive validates MODPATH endpoint/pathline topology and source-receptor projection; "
            "it does not validate field geochemistry, capture-zone polygons, or travel time until time definitions are harmonised."
        ),
    }
    return pd.DataFrame([row])


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    summary = _read_csv(M2_MODPATH_RESULTS / "modpath_topology_summary.csv")
    agreement = _read_csv(M2_MODPATH_RESULTS / "modpath_topology_agreement.csv")
    particles = _read_csv(M2_MODPATH_RESULTS / "modpath_pathline_particles.csv")
    graph_priors = _read_csv(M2_MODPATH_RESULTS / "modpath_graph_priors.csv")
    manifest_path = M2_MODPATH_RESULTS / "e2_modpath_source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}

    if summary.empty or agreement.empty or particles.empty:
        raise FileNotFoundError(
            "M2 external MODPATH result files are missing. Expected modpath_topology_summary.csv, "
            "modpath_topology_agreement.csv, and modpath_pathline_particles.csv."
        )

    endpoint_file = M2_MODPATH_INPUT / "base-MP.end"
    pathline_file = M2_MODPATH_INPUT / "base-MP.path"
    structure = build_pathline_structure(endpoint_file, pathline_file)
    capture_detail, capture_summary = build_capture_envelope_overlap(endpoint_file, pathline_file)
    travel_rank, travel_rank_summary = build_travel_time_rank(agreement, particles)
    harmonized_time, harmonized_time_summary = build_harmonized_travel_time(structure, graph_priors)
    bootstrap = build_bootstrap_ci(summary, agreement, particles, structure)
    archive_summary = build_archive_summary(
        manifest,
        summary,
        agreement,
        particles,
        structure,
        travel_rank_summary,
        capture_summary,
        harmonized_time_summary,
    )

    agreement.to_csv(RESULT_DIR / "external_modpath_edge_agreement.csv", index=False)
    particles.to_csv(RESULT_DIR / "external_modpath_pathline_particles.csv", index=False)
    structure.to_csv(RESULT_DIR / "external_modpath_pathline_structure.csv", index=False)
    capture_detail.to_csv(RESULT_DIR / "external_modpath_capture_envelope_overlap.csv", index=False)
    capture_summary.to_csv(RESULT_DIR / "external_modpath_capture_envelope_summary.csv", index=False)
    travel_rank.to_csv(RESULT_DIR / "external_modpath_travel_time_rank.csv", index=False)
    travel_rank_summary.to_csv(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv", index=False)
    harmonized_time.to_csv(RESULT_DIR / "external_modpath_harmonized_travel_time.csv", index=False)
    harmonized_time_summary.to_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv", index=False)
    bootstrap.to_csv(RESULT_DIR / "external_modpath_bootstrap_ci.csv", index=False)
    archive_summary.to_csv(RESULT_DIR / "external_modpath_archive_summary.csv", index=False)

    external_manifest = {
        "source": "M2 external MODPATH validation outputs",
        "m2_manifest": manifest,
        "m4_outputs": {
            "archive_summary": str(RESULT_DIR / "external_modpath_archive_summary.csv"),
            "edge_agreement": str(RESULT_DIR / "external_modpath_edge_agreement.csv"),
            "pathline_particles": str(RESULT_DIR / "external_modpath_pathline_particles.csv"),
            "pathline_structure": str(RESULT_DIR / "external_modpath_pathline_structure.csv"),
            "capture_envelope_overlap": str(RESULT_DIR / "external_modpath_capture_envelope_overlap.csv"),
            "capture_envelope_summary": str(RESULT_DIR / "external_modpath_capture_envelope_summary.csv"),
            "travel_time_rank": str(RESULT_DIR / "external_modpath_travel_time_rank.csv"),
            "travel_time_rank_summary": str(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv"),
            "harmonized_travel_time": str(RESULT_DIR / "external_modpath_harmonized_travel_time.csv"),
            "harmonized_travel_time_summary": str(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv"),
            "bootstrap_ci": str(RESULT_DIR / "external_modpath_bootstrap_ci.csv"),
        },
    }
    (RESULT_DIR / "external_modpath_source_manifest.json").write_text(json.dumps(external_manifest, indent=2), encoding="utf-8")
    print(
        "Wrote M4 external MODPATH archive validation: "
        f"{int(archive_summary.iloc[0]['n_particles'])} particles, "
        f"{int(archive_summary.iloc[0]['n_pathline_edges'])} pathline edges."
    )


if __name__ == "__main__":
    main()
