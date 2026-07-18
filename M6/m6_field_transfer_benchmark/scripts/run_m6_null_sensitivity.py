"""Run a separate Hydrosheaf non-flow-null sensitivity analysis for M6.

This extension deliberately leaves the locked M6 field-transfer outputs intact.
It scores the same M6 edge sets with the core Hydrosheaf null model and reports
where chemistry/isotope similarity has a competing no-flow explanation.  It is
an edge-level sensitivity diagnostic, not a new field-validation claim.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import m6_common as m6  # noqa: E402
from hydrosheaf import Config  # noqa: E402
from hydrosheaf.null_models import compute_null_penalty  # noqa: E402


RESULTS = Path(__file__).resolve().parents[1] / "results"
REPORT = Path(__file__).resolve().parents[1] / "docs" / "m6_null_model_sensitivity.md"
ION_ORDER = list(m6.MAJOR_IONS) + ["F"]


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _sample_mapping(row: Mapping[str, object]) -> dict[str, object]:
    """Map the harmonised M6 row to the null-model input contract."""
    out: dict[str, object] = {}
    for ion in ION_ORDER:
        out[ion] = _finite(row.get(ion))
    out["18O"] = _finite(row.get("d18O"))
    out["2H"] = _finite(row.get("d2H"))
    out["lat"] = _finite(row.get("Latitude"))
    out["lon"] = _finite(row.get("Longitude"))
    out["aquifer_unit"] = row.get("Aquifer_Type")
    out["aquifer_layer"] = row.get("Geology_Group")
    out["lithology"] = row.get("Lithology")
    # The chemistry vector is mmol/L.  The anthropogenic screening thresholds
    # are mg/L, so provide explicit companions rather than mixing units.
    for ion in ("NO3", "Cl"):
        value = _finite(row.get(ion))
        if value is not None:
            out[f"{ion}_mg_L"] = value * float(m6.m5.MOLAR_MASS_G_MOL[ion])
    return out


def _null_config() -> Config:
    return Config(
        ion_order=ION_ORDER,
        weights=[1.0] * len(ION_ORDER),
        conservative_weights=[1.0] * len(ION_ORDER),
        isotope_enabled=True,
        null_model_enabled=True,
        evidence_ladder_enabled=True,
        assumption_calibration_enabled=True,
        null_model_weight=0.5,
        null_chemistry_similarity_threshold=0.3,
        null_lithology_weight=0.3,
        null_endmember_weight=0.4,
        null_spatial_weight=0.2,
        null_anthropogenic_weight=0.2,
        # The null-model implementation uses intercept + slope names.
        lmwl_a=7.22,
        lmwl_b=8.66,
    )


def _northern_edge_sets(df: pd.DataFrame) -> dict[str, list[tuple[str, str]]]:
    dry = df[df["season"].eq("Dry")].reset_index(drop=True)
    nodes = {r["site_id"]: r["sample_id"] for r in dry.to_dict("records")}
    graph = pd.read_excel(m6.DATA["northern_ghana"], sheet_name="Graph_Edges")
    provided: list[tuple[str, str]] = []
    for _, row in graph.iterrows():
        source = nodes.get(row.get("Source_Well_ID"))
        target = nodes.get(row.get("Target_Well_ID"))
        if source and target:
            provided.append((source, target))
    rng = np.random.default_rng(m6.SEED)
    return {
        "provided_graph": sorted(set(provided)),
        "chemistry_knn": m6.chemistry_knn_edges(dry, k=3),
        "geographic_nearest": m6.geographic_edges(dry, k=3),
        "random_perturbed": m6.random_edges(dry, n=len(provided), rng=rng),
    }


def _score_edges(dataset: str, frame: pd.DataFrame, edge_sets: dict[str, list[tuple[str, str]]],
                 config: Config, cap: int = 140) -> list[dict[str, object]]:
    lookup = {r["sample_id"]: _sample_mapping(r) for r in frame.to_dict("records")}
    rows: list[dict[str, object]] = []
    for edge_set, edges in edge_sets.items():
        for source, target in list(edges)[:cap]:
            src = lookup.get(str(source))
            tgt = lookup.get(str(target))
            if src is None or tgt is None:
                continue
            score, flags = compute_null_penalty(src, tgt, config)
            score = float(score)
            status = (
                "no_flow_dominant" if score > 0.8
                else "no_flow_competing" if score > 0.5
                else "not_flagged"
            )
            rows.append({
                "dataset": dataset,
                "edge_set": edge_set,
                "source": str(source),
                "target": str(target),
                "null_score": score,
                "null_status": status,
                "null_flags": ";".join(sorted(set(flags))),
            })
    return rows


def _summarise(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, edge_set), group in scores.groupby(["dataset", "edge_set"], sort=True):
        flag_counts = Counter(
            flag for cell in group["null_flags"] for flag in str(cell).split(";") if flag
        )
        rows.append({
            "dataset": dataset,
            "edge_set": edge_set,
            "n_edges": int(len(group)),
            "mean_null_score": float(group["null_score"].mean()),
            "median_null_score": float(group["null_score"].median()),
            "fraction_gt_0_5": float((group["null_score"] > 0.5).mean()),
            "fraction_gt_0_8": float((group["null_score"] > 0.8).mean()),
            "fraction_retained_at_0_5": float((group["null_score"] <= 0.5).mean()),
            "top_null_flags": ";".join(f"{k}:{v}" for k, v in flag_counts.most_common()),
        })
    return pd.DataFrame(rows)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    data = m6.load_all()
    config = _null_config()
    all_rows: list[dict[str, object]] = []

    all_rows.extend(_score_edges("northern_ghana", data["northern_ghana"],
                                 _northern_edge_sets(data["northern_ghana"]), config))
    for dataset in ("talensi", "manu"):
        frame = data[dataset]
        all_rows.extend(_score_edges(
            dataset, frame, {"chemistry_knn": m6.chemistry_knn_edges(frame, k=3)}, config
        ))

    scores = pd.DataFrame(all_rows)
    summary = _summarise(scores)
    score_path = RESULTS / "m6_null_edge_scores.csv"
    summary_path = RESULTS / "m6_null_sensitivity_summary.csv"
    run_path = RESULTS / "m6_null_model_run.json"
    scores.to_csv(score_path, index=False)
    summary.to_csv(summary_path, index=False)

    source_hashes = {name: _sha256(path) for name, path in m6.DATA.items() if path.exists()}
    metadata = {
        "analysis": "M6 null-model sensitivity extension",
        "locked_baseline_unchanged": True,
        "thresholds": {"competing_no_flow": 0.5, "dominant_no_flow": 0.8},
        "config": {
            "null_model_enabled": True,
            "evidence_ladder_enabled": True,
            "assumption_calibration_enabled": True,
            "null_chemistry_similarity_threshold": 0.3,
            "null_lithology_weight": 0.3,
            "null_endmember_weight": 0.4,
            "null_spatial_weight": 0.2,
            "null_anthropogenic_weight": 0.2,
        },
        "source_hashes": source_hashes,
        "outputs": [str(score_path), str(summary_path)],
    }
    run_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M6 Null-Model Sensitivity Extension",
        "",
        "This is a separate edge-level sensitivity analysis. The locked M6 field-transfer",
        "outputs are unchanged. Null scores quantify whether chemistry, recharge, metadata,",
        "or common-source similarity creates a competing no-flow explanation.",
        "",
        "## Outputs",
        "",
        f"- `{score_path.as_posix()}`: edge-level scores and flags.",
        f"- `{summary_path.as_posix()}`: dataset/edge-set summaries.",
        f"- `{run_path.as_posix()}`: configuration and source hashes.",
        "",
        "Scores >0.5 are reported as a competing no-flow explanation; scores >0.8 are",
        "reported as a dominant no-flow explanation. These are screening thresholds from",
        "Hydrosheaf's evidence-classification logic, not calibrated field probabilities.",
        "",
        "## Summary",
        "",
        summary.to_markdown(index=False) if not summary.empty else "No valid edges were scored.",
        "",
        "## Guardrail",
        "",
        "This extension does not establish field process truth. It tests sensitivity to",
        "non-flow explanations and should be interpreted alongside the locked M6 evidence",
        "gate, edge-set sensitivity, and identifiability results.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"\nWrote {score_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
