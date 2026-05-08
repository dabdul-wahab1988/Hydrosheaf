"""Run E2 MODPATH topology validation for the M2 benchmark.

This validation uses the public USGS Savage MODFLOW/MODPATH archive already
staged under `external/modpath/input`. Endpoints are converted into Hydrosheaf
graph priors, while compact pathlines provide an independent particle-track
check on source-receptor direction and overlap.
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.physics.modpath import (  # noqa: E402
    ModpathEndpoint,
    ModpathPathlinePoint,
    load_modpath_endpoint_records,
    load_modpath_pathline_points,
)

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = BENCHMARK_ROOT / "external" / "modpath" / "input"
RESULT_DIR = BENCHMARK_ROOT / "external" / "modpath" / "results"
FIGURE_DIR = BENCHMARK_ROOT / "figures"

SOURCE_REPO = "https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute"
SOURCE_DOI = "10.5066/F7J102FK"
ENDPOINT_FILE = INPUT_DIR / "selected_output" / "base-MP.end"
PATHLINE_FILE = INPUT_DIR / "selected_output" / "base-MP.path"


@dataclass(frozen=True)
class EdgeSummary:
    source_node: str
    target_node: str
    particle_count: int
    p_uv: float
    travel_time_mean: float
    travel_time_p10: float
    travel_time_p50: float
    travel_time_p90: float
    travel_time_std: float

    @property
    def edge_id(self) -> str:
        return f"{self.source_node}->{self.target_node}"


def _cell_node(cell: Optional[int], prefix: str = "cell") -> Optional[str]:
    if cell is None:
        return None
    return f"{prefix}_{int(cell)}"


def _quantile(values: Sequence[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.quantile(np.asarray(values, dtype=float), q))


def _endpoint_edge_summaries(records: Iterable[ModpathEndpoint]) -> List[EdgeSummary]:
    times_by_edge: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    count_by_source: Dict[str, int] = defaultdict(int)
    for rec in records:
        source = _cell_node(rec.initial_cell)
        target = _cell_node(rec.final_cell)
        if source is None or target is None or source == target:
            continue
        times_by_edge[(source, target)].append(float(rec.time))
        count_by_source[source] += 1

    summaries: List[EdgeSummary] = []
    for (source, target), times in sorted(times_by_edge.items()):
        arr = np.asarray(times, dtype=float)
        summaries.append(
            EdgeSummary(
                source_node=source,
                target_node=target,
                particle_count=int(len(times)),
                p_uv=float(len(times) / max(1, count_by_source[source])),
                travel_time_mean=float(arr.mean()),
                travel_time_p10=_quantile(times, 0.10),
                travel_time_p50=_quantile(times, 0.50),
                travel_time_p90=_quantile(times, 0.90),
                travel_time_std=float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
            )
        )
    return summaries


def _pathline_endpoint_edges(points: Iterable[ModpathPathlinePoint]) -> Tuple[Dict[Tuple[str, str], int], pd.DataFrame]:
    by_particle: Dict[int, List[ModpathPathlinePoint]] = defaultdict(list)
    for point in points:
        by_particle[int(point.particle_id)].append(point)

    edge_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    rows: List[Dict[str, object]] = []
    for particle_id, particle_points in by_particle.items():
        ordered = sorted(particle_points, key=lambda point: point.sequence)
        first = next((point for point in ordered if point.cell is not None), None)
        last = next((point for point in reversed(ordered) if point.cell is not None), None)
        if first is None or last is None:
            continue
        source = _cell_node(first.cell)
        target = _cell_node(last.cell)
        if source is None or target is None or source == target:
            continue
        edge_counts[(source, target)] += 1
        rows.append(
            {
                "particle_id": particle_id,
                "source_node": source,
                "target_node": target,
                "n_pathline_points": len(ordered),
                "pathline_start_time": first.time,
                "pathline_end_time": last.time,
                "pathline_elapsed_time": last.time - first.time,
                "start_x": first.x,
                "start_y": first.y,
                "end_x": last.x,
                "end_y": last.y,
            }
        )
    return edge_counts, pd.DataFrame(rows)


def _agreement_rows(endpoint_summaries: Sequence[EdgeSummary], pathline_counts: Dict[Tuple[str, str], int]) -> pd.DataFrame:
    endpoint_edges = {(edge.source_node, edge.target_node): edge for edge in endpoint_summaries}
    all_edges = sorted(set(endpoint_edges) | set(pathline_counts))
    rows: List[Dict[str, object]] = []
    for source, target in all_edges:
        endpoint = endpoint_edges.get((source, target))
        path_count = int(pathline_counts.get((source, target), 0))
        reverse_count = int(pathline_counts.get((target, source), 0))
        endpoint_count = int(endpoint.particle_count) if endpoint else 0
        match = endpoint is not None and path_count > 0
        rows.append(
            {
                "edge_id": f"{source}->{target}",
                "source_node": source,
                "target_node": target,
                "endpoint_particle_count": endpoint_count,
                "pathline_particle_count": path_count,
                "reverse_pathline_particle_count": reverse_count,
                "classification": "TP" if match else ("FP" if endpoint is not None else "FN"),
                "direction_agrees": bool(match and reverse_count == 0),
                "source_receptor_overlap": float(min(endpoint_count, path_count) / max(endpoint_count, path_count, 1)),
                "travel_time_mean": endpoint.travel_time_mean if endpoint else float("nan"),
                "travel_time_p10": endpoint.travel_time_p10 if endpoint else float("nan"),
                "travel_time_p90": endpoint.travel_time_p90 if endpoint else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def _summary(agreement: pd.DataFrame, pathline_rows: pd.DataFrame) -> Dict[str, float]:
    tp = int((agreement["classification"] == "TP").sum()) if not agreement.empty else 0
    fp = int((agreement["classification"] == "FP").sum()) if not agreement.empty else 0
    fn = int((agreement["classification"] == "FN").sum()) if not agreement.empty else 0
    denom = max(tp + fp + fn, 1)
    direction = float(agreement.loc[agreement["classification"] == "TP", "direction_agrees"].mean()) if tp else 0.0
    overlap = float(agreement.loc[agreement["classification"] == "TP", "source_receptor_overlap"].mean()) if tp else 0.0
    elapsed = pathline_rows["pathline_elapsed_time"].dropna() if "pathline_elapsed_time" in pathline_rows else pd.Series(dtype=float)
    return {
        "n_endpoint_edges": int((agreement["endpoint_particle_count"] > 0).sum()) if not agreement.empty else 0,
        "n_pathline_edges": int((agreement["pathline_particle_count"] > 0).sum()) if not agreement.empty else 0,
        "true_positive_edges": tp,
        "false_positive_edges": fp,
        "false_negative_edges": fn,
        "edge_precision": float(tp / max(tp + fp, 1)),
        "edge_recall": float(tp / max(tp + fn, 1)),
        "edge_f1": float((2 * tp) / max(2 * tp + fp + fn, 1)),
        "direction_agreement_rate": direction,
        "mean_source_receptor_overlap": overlap,
        "pathline_elapsed_time_median": float(elapsed.median()) if not elapsed.empty else float("nan"),
        "pathline_elapsed_time_p90": float(elapsed.quantile(0.90)) if not elapsed.empty else float("nan"),
        "edge_error_rate": float((fp + fn) / denom),
    }


def _write_figure(pathline_rows: pd.DataFrame, graph_priors: pd.DataFrame, summary: Dict[str, float]) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = FIGURE_DIR / "figure_s2_modpath_to_graph_prior_real_archive.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)

    ax = axes[0]
    sample = pathline_rows
    if len(sample) > 2500:
        sample = sample.sample(2500, random_state=42)
    for _, row in sample.iterrows():
        ax.plot([row["start_x"], row["end_x"]], [row["start_y"], row["end_y"]], color="#8aa6c1", alpha=0.12, linewidth=0.6)
    ax.scatter(sample["start_x"], sample["start_y"], s=6, color="#2f6f4e", alpha=0.5, label="particle starts")
    ax.scatter(sample["end_x"], sample["end_y"], s=6, color="#8f3d38", alpha=0.5, label="particle endpoints")
    ax.set_title("MODPATH particle source-to-receptor tracks")
    ax.set_xlabel("model x")
    ax.set_ylabel("model y")
    ax.legend(loc="best", fontsize=8)

    ax = axes[1]
    top = graph_priors.sort_values("particle_count", ascending=False).head(12).copy()
    labels = [edge.replace("cell_", "") for edge in top["edge_id"]]
    y = np.arange(len(top))
    ax.barh(y, top["particle_count"], color="#456990")
    ax.set_yticks(y, labels=labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("particle count")
    ax.set_title("Hydrosheaf graph priors from endpoints")
    metric = (
        f"TP={summary['true_positive_edges']} FP={summary['false_positive_edges']} FN={summary['false_negative_edges']}\n"
        f"F1={summary['edge_f1']:.2f}; direction={summary['direction_agreement_rate']:.2f}"
    )
    ax.text(0.98, 0.04, metric, transform=ax.transAxes, ha="right", va="bottom", fontsize=9)

    fig.suptitle("Figure S2. MODPATH-to-graph prior validation from USGS Savage archive", fontsize=12)
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)
    return figure_path


def _update_tables(summary: Dict[str, float], report_path: Path) -> None:
    metric = (
        f"TP={int(summary['true_positive_edges'])}; FP={int(summary['false_positive_edges'])}; "
        f"FN={int(summary['false_negative_edges'])}; F1={summary['edge_f1']:.2f}; "
        f"direction agreement={summary['direction_agreement_rate']:.2f}"
    )
    table4 = BENCHMARK_ROOT / "tables" / "table4_validation_design_and_results.csv"
    if table4.exists():
        df = pd.read_csv(table4)
        mask = df["benchmark"] == "MODFLOW/MODPATH topology validation"
        if mask.any():
            df.loc[mask, "performance_metric"] = metric
            df.loc[mask, "m2_status"] = "completed"
            df.loc[mask, "notes"] = f"Report: {report_path.name}"
            df.to_csv(table4, index=False)

    workplan = BENCHMARK_ROOT / "tables" / "external_validation_workplan.csv"
    if workplan.exists():
        df = pd.read_csv(workplan)
        mask = df["validation_tier"] == "E2_modpath_topology"
        if mask.any():
            df.loc[mask, "status"] = "completed"
            df.loc[mask, "notes"] = f"{metric}; report: {report_path.name}"
            df.to_csv(workplan, index=False)

    inventory = BENCHMARK_ROOT / "tables" / "table_s4_benchmark_dataset_inventory.csv"
    if inventory.exists():
        df = pd.read_csv(inventory)
        mask = df["resource"] == "USGS Savage MODFLOW/MODPATH model archive"
        if mask.any():
            df.loc[mask, "m2_status"] = "completed"
            df.loc[mask, "notes"] = metric
            df.to_csv(inventory, index=False)


def main() -> int:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    if not ENDPOINT_FILE.exists() or not PATHLINE_FILE.exists():
        raise FileNotFoundError(
            f"E2 requires {ENDPOINT_FILE} and {PATHLINE_FILE}. Extract the USGS MODPATH archive first."
        )

    endpoints = load_modpath_endpoint_records(str(ENDPOINT_FILE))
    pathline_points = load_modpath_pathline_points(str(PATHLINE_FILE))
    endpoint_summaries = _endpoint_edge_summaries(endpoints)
    pathline_edges, pathline_rows = _pathline_endpoint_edges(pathline_points)
    agreement = _agreement_rows(endpoint_summaries, pathline_edges)
    summary = _summary(agreement, pathline_rows)

    graph_priors = pd.DataFrame(
        [
            {
                "edge_id": edge.edge_id,
                **asdict(edge),
                "source": "USGS Savage MODPATH5 endpoint file",
                "source_file": str(ENDPOINT_FILE),
            }
            for edge in endpoint_summaries
        ]
    )
    if not graph_priors.empty:
        graph_priors = graph_priors.sort_values(["particle_count", "edge_id"], ascending=[False, True])

    graph_priors_path = RESULT_DIR / "modpath_graph_priors.csv"
    agreement_path = RESULT_DIR / "modpath_topology_agreement.csv"
    particles_path = RESULT_DIR / "modpath_pathline_particles.csv"
    summary_path = RESULT_DIR / "modpath_topology_summary.csv"
    report_path = RESULT_DIR / "e2_modpath_topology_report.md"
    manifest_path = RESULT_DIR / "e2_modpath_source_manifest.json"

    graph_priors.to_csv(graph_priors_path, index=False)
    agreement.to_csv(agreement_path, index=False)
    pathline_rows.to_csv(particles_path, index=False)
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    figure_path = _write_figure(pathline_rows, graph_priors, summary)

    report_lines = [
        "# E2 MODPATH Topology Validation Report",
        "",
        f"Run timestamp UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"Source: USGS Savage MODFLOW/MODPATH archive, DOI `{SOURCE_DOI}`.",
        "",
        "## Outputs",
        "",
        f"- Graph priors: `{graph_priors_path}`",
        f"- Topology agreement: `{agreement_path}`",
        f"- Particle pathline summary: `{particles_path}`",
        f"- Figure S2: `{figure_path}`",
        "",
        "## Summary",
        "",
        pd.DataFrame([summary]).to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "Hydrosheaf converts compact MODPATH5 endpoints into directed cell-to-cell graph priors. Compact MODPATH5 pathlines are used as an independent check that endpoint-derived source-receptor edges preserve particle-tracking direction and overlap. This validates the physics-prior ingestion layer; it does not validate geochemical edge inference because the MODPATH archive does not contain paired hydrochemical samples for the same graph nodes.",
        "",
    ]
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "validation_tier": "E2_modpath_topology",
                "source_url": SOURCE_REPO,
                "source_doi": SOURCE_DOI,
                "endpoint_file": str(ENDPOINT_FILE),
                "pathline_file": str(PATHLINE_FILE),
                "n_endpoint_records": len(endpoints),
                "n_pathline_points": len(pathline_points),
                "outputs": {
                    "graph_priors": str(graph_priors_path),
                    "topology_agreement": str(agreement_path),
                    "pathline_particles": str(particles_path),
                    "summary": str(summary_path),
                    "figure_s2": str(figure_path),
                    "report": str(report_path),
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    _update_tables(summary, report_path)
    print(f"E2 completed with {len(graph_priors)} endpoint-derived graph-prior edges.")
    print(f"Topology summary written to {summary_path}")
    print(f"Figure S2 written to {figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
