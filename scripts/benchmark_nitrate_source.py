"""
Benchmark nitrate source inference performance on synthetic data.

This script targets `hydrosheaf.nitrate_source_v2.infer_node_posteriors`
and reports runtime for:
  - hydro_only: hydrochemical-only path
  - with_isotopes: dual-isotope analytical path
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

# Allow running the script without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from hydrosheaf.nitrate_source_v2 import infer_node_posteriors


@dataclass
class EdgeLike:
    u: str
    v: str
    z_labels: List[str]
    z_extents: List[float]
    transport_model: str


def _build_synthetic_nodes(
    node_ids: List[str],
    include_isotopes: bool,
    rng: np.random.Generator,
) -> pd.DataFrame:
    n_nodes = len(node_ids)
    data: Dict[str, Any] = {
        "site_id": node_ids,
        "NO3": rng.lognormal(mean=3.0, sigma=0.6, size=n_nodes),
        "Cl": rng.lognormal(mean=2.0, sigma=0.6, size=n_nodes),
        "K": rng.lognormal(mean=1.4, sigma=0.6, size=n_nodes),
        "PO4": rng.lognormal(mean=0.6, sigma=0.7, size=n_nodes),
        "Fe": rng.lognormal(mean=-0.1, sigma=0.8, size=n_nodes),
        "HCO3": rng.lognormal(mean=3.3, sigma=0.4, size=n_nodes),
        "Ca": rng.lognormal(mean=2.2, sigma=0.5, size=n_nodes),
        "Mg": rng.lognormal(mean=1.6, sigma=0.5, size=n_nodes),
        "Na": rng.lognormal(mean=1.8, sigma=0.5, size=n_nodes),
        "SO4": rng.lognormal(mean=1.7, sigma=0.5, size=n_nodes),
        "d_excess": rng.normal(loc=10.0, scale=4.0, size=n_nodes),
    }
    if include_isotopes:
        data["d15N"] = rng.normal(loc=8.0, scale=4.0, size=n_nodes)
        data["d18O_NO3"] = rng.normal(loc=8.0, scale=5.0, size=n_nodes)
    df = pd.DataFrame(data)
    df.set_index("site_id", inplace=True, drop=False)
    return df


def _build_synthetic_edges(
    node_ids: List[str],
    n_edges: int,
    py_rng: random.Random,
    np_rng: np.random.Generator,
) -> List[EdgeLike]:
    n_nodes = len(node_ids)
    edges: List[EdgeLike] = []
    while len(edges) < n_edges:
        u = node_ids[py_rng.randrange(0, n_nodes)]
        v = node_ids[py_rng.randrange(0, n_nodes)]
        if u == v:
            continue
        edges.append(
            EdgeLike(
                u=u,
                v=v,
                z_labels=["denit"],
                z_extents=[float(np_rng.uniform(-0.2, 1.5))],
                transport_model="evap" if py_rng.random() < 0.08 else "mix",
            )
        )
    return edges


def _time_case(
    nodes_df: pd.DataFrame,
    edges: List[EdgeLike],
    config_overrides: Dict[str, Any],
    warmup: int,
    repeats: int,
) -> Tuple[List[float], int]:
    for _ in range(max(0, warmup)):
        infer_node_posteriors(nodes_df, edges, config_overrides)

    timings: List[float] = []
    result_size = 0
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        result = infer_node_posteriors(nodes_df, edges, config_overrides)
        timings.append(time.perf_counter() - t0)
        result_size = len(result)
    return timings, result_size


def _summarize(timings: List[float], n_nodes: int) -> Dict[str, float]:
    mean_s = statistics.mean(timings)
    min_s = min(timings)
    max_s = max(timings)
    p50_s = statistics.median(timings)
    nodes_per_s = float(n_nodes) / mean_s if mean_s > 0 else 0.0
    return {
        "mean_s": mean_s,
        "min_s": min_s,
        "max_s": max_s,
        "p50_s": p50_s,
        "nodes_per_s": nodes_per_s,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark nitrate source inference.")
    parser.add_argument("--nodes", type=int, default=4000, help="Number of synthetic nodes.")
    parser.add_argument("--edges", type=int, default=10000, help="Number of synthetic edges.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs per mode.")
    parser.add_argument("--repeats", type=int, default=3, help="Timed runs per mode.")
    parser.add_argument(
        "--mode",
        choices=["hydro", "isotopes", "both"],
        default="both",
        help="Benchmark mode.",
    )
    parser.add_argument(
        "--output-json",
        type=str,
        default="",
        help="Optional JSON output path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.nodes < 2:
        raise ValueError("--nodes must be >= 2")
    if args.edges < 1:
        raise ValueError("--edges must be >= 1")

    np_rng = np.random.default_rng(args.seed)
    py_rng = random.Random(args.seed)
    node_ids = [f"S{i}" for i in range(args.nodes)]
    edges = _build_synthetic_edges(node_ids, args.edges, py_rng, np_rng)

    selected_modes: List[str]
    if args.mode == "both":
        selected_modes = ["hydro_only", "with_isotopes"]
    elif args.mode == "hydro":
        selected_modes = ["hydro_only"]
    else:
        selected_modes = ["with_isotopes"]

    all_results: Dict[str, Dict[str, Any]] = {}
    for mode in selected_modes:
        include_isotopes = mode == "with_isotopes"
        nodes_df = _build_synthetic_nodes(node_ids, include_isotopes, np_rng)
        cfg = {"nitrate_isotope_mixing_enabled": include_isotopes, "prior_prob": 0.5}

        timings, n_out = _time_case(nodes_df, edges, cfg, args.warmup, args.repeats)
        summary = _summarize(timings, args.nodes)
        all_results[mode] = {
            "timings_s": timings,
            "n_output_nodes": n_out,
            **summary,
        }

    print("=" * 80)
    print("NITRATE SOURCE BENCHMARK")
    print("=" * 80)
    print(f"nodes={args.nodes} edges={args.edges} repeats={args.repeats} warmup={args.warmup}")
    for mode in selected_modes:
        res = all_results[mode]
        print(
            f"{mode}: mean={res['mean_s']:.3f}s min={res['min_s']:.3f}s "
            f"max={res['max_s']:.3f}s p50={res['p50_s']:.3f}s "
            f"throughput={res['nodes_per_s']:.1f} nodes/s"
        )

    if "hydro_only" in all_results and "with_isotopes" in all_results:
        iso_mean = all_results["with_isotopes"]["mean_s"]
        hydro_mean = all_results["hydro_only"]["mean_s"]
        ratio = iso_mean / hydro_mean if hydro_mean > 0 else float("inf")
        print(f"isotope_to_hydro_ratio={ratio:.3f}x")

    if args.output_json:
        output_path = Path(args.output_json)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "nodes": args.nodes,
            "edges": args.edges,
            "seed": args.seed,
            "warmup": args.warmup,
            "repeats": args.repeats,
            "mode": args.mode,
            "results": all_results,
        }
        output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"wrote_json={output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
