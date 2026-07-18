"""Run a separate non-flow-null sensitivity analysis for the M7 synthetic twin.

The published M7 integration benchmark remains the locked baseline.  This script
replays the synthetic twin, scores every candidate edge with Hydrosheaf's formal
null model, and reports how a conservative no-flow gate changes the three-stream
joint classification.  Because M7 is synthetic, the isotope, location, and
lithology fields added here are explicitly labelled *constructed metadata*; they
are a mechanism/sensitivity test, not field validation.
"""
from __future__ import annotations

import builtins
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "M6" / "m6_field_transfer_benchmark" / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "M6" / "m6_field_transfer_benchmark" / "scripts"))

from hydrosheaf import Config  # noqa: E402
from hydrosheaf.null_models import compute_null_penalty  # noqa: E402


BASELINE_SCRIPT = SCRIPT_DIR / "run_m7_integration_benchmark.py"
RESULTS = Path(__file__).resolve().parents[1] / "results"
REPORT = Path(__file__).resolve().parents[1] / "docs" / "m7_null_model_sensitivity.md"
ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]
NO3_MG_L = 62.0049
CL_MG_L = 35.45


def _load_baseline_module():
    spec = importlib.util.spec_from_file_location("m7_baseline_for_null", BASELINE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {BASELINE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _stable_hash(value: object) -> int:
    """Deterministic replacement for the baseline generator's built-in hash call."""
    digest = hashlib.sha256(str(value).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _replay_twin(m7):
    """Replay M7 without changing the baseline source or its locked outputs.

    The original generator uses ``hash(node_id)`` to select an archetype. Python
    randomises that hash between processes, so the extension patches it only for
    this in-memory replay and records the stabilisation in its metadata.
    """
    original_hash = builtins.hash
    builtins.hash = _stable_hash
    try:
        return m7.generate(np.random.default_rng(m7.SEED))
    finally:
        builtins.hash = original_hash


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _sample_mapping(node: Mapping[str, object]) -> dict[str, object]:
    """Map a synthetic M7 node to the formal null-model input contract."""
    chem = node.get("chem") or {}
    out = {ion: _finite(chem.get(ion)) for ion in ION_ORDER}

    # Synthetic geographic coordinates derived from x/y metres. They are used
    # only to exercise the spatial-autocorrelation branch of the null model.
    x = _finite(node.get("x")) or 0.0
    y = _finite(node.get("y")) or 0.0
    lat = 40.0 + y / 111_000.0
    lon = -100.0 + x / (111_000.0 * np.cos(np.deg2rad(40.0)))
    out["lat"], out["lon"] = lat, lon

    # Constructed isotope values are internally consistent with the LMWL. A
    # small branch-specific offset creates a controlled shared-recharge signal.
    d18 = -6.0 + 0.4 * (_finite(node.get("s")) or 0.0)
    branch = str(node.get("branch") or "unknown")
    if branch == "fossil":
        d18 -= 0.4
    elif branch == "young_mineralised":
        d18 += 0.2
    out["18O"] = d18
    out["2H"] = 8.66 * d18 + 7.22

    # M7 has no measured lithology/aquifer metadata. Branch labels are a
    # constructed categorical surrogate, not an observed geological fact.
    out["lithology"] = branch
    out["aquifer_unit"] = branch
    out["aquifer_layer"] = branch

    # The core chemistry vector is mmol/L; anthropogenic thresholds are mg/L.
    if out.get("NO3") is not None:
        out["NO3_mg_L"] = out["NO3"] * NO3_MG_L
    if out.get("Cl") is not None:
        out["Cl_mg_L"] = out["Cl"] * CL_MG_L
    return out


def _null_config() -> Config:
    cfg = Config()
    cfg.ion_order = ION_ORDER
    cfg.weights = [1.0] * len(ION_ORDER)
    cfg.conservative_weights = [1.0] * len(ION_ORDER)
    cfg.isotope_enabled = True
    cfg.null_model_enabled = True
    cfg.evidence_ladder_enabled = True
    cfg.assumption_calibration_enabled = True
    cfg.null_model_weight = 0.5
    cfg.null_chemistry_similarity_threshold = 0.3
    cfg.null_lithology_weight = 0.3
    cfg.null_endmember_weight = 0.4
    cfg.null_spatial_weight = 0.2
    cfg.null_anthropogenic_weight = 0.2
    cfg.lmwl_a = 7.22
    cfg.lmwl_b = 8.66
    return cfg


def _prf(frame: pd.DataFrame, column: str) -> dict[str, float]:
    tp = int(((frame[column] == 1) & (frame["is_true"] == 1)).sum())
    fp = int(((frame[column] == 1) & (frame["is_true"] == 0)).sum())
    fn = int(((frame[column] == 0) & (frame["is_true"] == 1)).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    traps = frame[frame["is_true"] == 0]
    trap_accept = float((traps[column] == 1).mean()) if len(traps) else 0.0
    return {"precision": precision, "recall": recall, "f1": f1,
            "trap_accept_rate": trap_accept}


def main() -> None:
    m7 = _load_baseline_module()
    nodes, true_edges, candidates = _replay_twin(m7)
    config = _null_config()
    ids = list(nodes)
    samples = {node_id: _sample_mapping(node) for node_id, node in nodes.items()}

    # Reproduce the baseline three evidence streams on this paired replay.
    inferred = {node_id: m7.single_node_age(nodes[node_id]["3H_obs"]) for node_id in ids}
    audit = m7.audit_graph_age_coherence(
        [(u, v) for u, v, _ in candidates], inferred,
        min_downstream_increase_years=-8.0,
    )
    coherence = {}
    for rec in audit.get("edges", []):
        key = (str(rec.get("upstream")), str(rec.get("downstream")))
        coherence[key] = int(not bool(rec.get("violation", False)))
    true_lengths = [float(np.hypot(nodes[a]["x"] - nodes[b]["x"],
                                   nodes[a]["y"] - nodes[b]["y"]))
                    for a, b in true_edges]
    near_threshold = 1.6 * float(np.median(true_lengths))

    rows = []
    for u, v, label in candidates:
        dist = float(np.hypot(nodes[u]["x"] - nodes[v]["x"],
                              nodes[u]["y"] - nodes[v]["y"]))
        geom_ok = int(nodes[u]["elev"] > nodes[v]["elev"] and dist < near_threshold)
        chem_score = m7.chem_edge_score(nodes[u]["chem"], nodes[v]["chem"])
        chem_ok = int(chem_score >= 1.0)
        age_ok = int(coherence.get((u, v), inferred[v] >= inferred[u]))
        null_score, flags = compute_null_penalty(samples[u], samples[v], config)
        null_score = float(null_score)
        null_ok = int(null_score <= 0.5)
        joint_ok = int(geom_ok and chem_ok and age_ok)
        rows.append({
            "u": u, "v": v, "label": label, "is_true": int(label == "true"),
            "dist": dist, "geom_ok": geom_ok, "chem_r2": chem_score,
            "chem_ok": chem_ok, "age_ok": age_ok, "joint_ok": joint_ok,
            "null_score": null_score, "null_ok_at_0_5": null_ok,
            "null_status": ("no_flow_dominant" if null_score > 0.8 else
                            "no_flow_competing" if null_score > 0.5 else "not_flagged"),
            "null_flags": ";".join(sorted(set(flags))),
            "joint_null_ok": int(joint_ok and null_ok),
        })
    edges = pd.DataFrame(rows)
    edge_path = RESULTS / "m7_null_edge_scores.csv"
    edges.to_csv(edge_path, index=False)

    gains = []
    for stream, col in [("geometry_only", "geom_ok"), ("chemistry_only", "chem_ok"),
                        ("age_only", "age_ok"), ("joint", "joint_ok"),
                        ("joint_plus_null_gate", "joint_null_ok"),
                        ("null_gate_only", "null_ok_at_0_5")]:
        gains.append({"stream": stream, **_prf(edges, col)})
    gain_df = pd.DataFrame(gains)

    trap_rows = []
    for trap_type in ["trapA", "trapB", "spurious"]:
        sub = edges[edges["label"] == trap_type]
        if sub.empty:
            continue
        trap_rows.append({
            "trap_type": trap_type, "n": int(len(sub)),
            "baseline_joint_rejects": float((sub["joint_ok"] == 0).mean()),
            "null_gate_rejects": float((sub["null_ok_at_0_5"] == 0).mean()),
            "joint_plus_null_rejects": float((sub["joint_null_ok"] == 0).mean()),
            "mean_null_score": float(sub["null_score"].mean()),
        })
    trap_df = pd.DataFrame(trap_rows)
    summary_path = RESULTS / "m7_null_sensitivity_summary.csv"
    gain_df.to_csv(summary_path, index=False)
    trap_path = RESULTS / "m7_null_trap_rejection.csv"
    trap_df.to_csv(trap_path, index=False)

    source_hash = hashlib.sha256(BASELINE_SCRIPT.read_bytes()).hexdigest()
    run_path = RESULTS / "m7_null_model_run.json"
    metadata = {
        "analysis": "M7 null-model sensitivity extension",
        "locked_baseline_unchanged": True,
        "synthetic_metadata": True,
        "constructed_fields": ["18O", "2H", "lat", "lon", "lithology", "aquifer_unit"],
        "null_gate": "null_score <= 0.5",
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
        "generator_replay": {
            "seed": m7.SEED,
            "stable_archetype_hash_used": True,
            "reason": "baseline generator uses process-randomised Python hash(node_id)",
        },
        "source_hashes": {"baseline_script": source_hash},
        "outputs": [str(edge_path), str(summary_path), str(trap_path)],
    }
    run_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M7 Null-Model Sensitivity Extension", "",
        "This is a paired, synthetic-twin sensitivity analysis. The locked M7",
        "integration outputs are unchanged. Hydrosheaf's formal null model is",
        "enabled and applied to every candidate edge; `joint_plus_null_gate` retains",
        "an edge only when the original geometry/chemistry/age joint gate passes and",
        "the null score is <=0.5.", "",
        "The isotope, location, and lithology fields are constructed solely to exercise",
        "the null-model branches. They are not observations and do not support a new",
        "field-validation claim.", "",
        "## Outputs", "",
        f"- `{edge_path.as_posix()}`: edge-level scores, flags, and paired gates.",
        f"- `{summary_path.as_posix()}`: precision/recall/F1 sensitivity summary.",
        f"- `{trap_path.as_posix()}`: trap-type rejection comparison.",
        f"- `{run_path.as_posix()}`: configuration, source hash, and replay caveat.", "",
        "## Paired classification summary", "",
        gain_df.to_markdown(index=False), "", "## Trap rejection", "",
        trap_df.to_markdown(index=False) if not trap_df.empty else "No trap rows were generated.",
        "", "## Interpretation guardrail", "",
        "A null score above 0.5 means that similarity has a competing no-flow",
        "explanation under the screening model; it is not a calibrated probability.",
        "Use this extension to qualify the M7 mechanism demonstration, not to replace",
        "the original controlled-twin result.", "",
        "The baseline generator uses Python's process-randomised `hash(node_id)` when",
        "assigning reaction archetypes. This extension uses a deterministic hash only",
        "for its replay and records that choice in the run metadata; the baseline source",
        "and baseline result files were not modified.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(gain_df.to_string(index=False))
    print("\n", trap_df.to_string(index=False))
    print(f"\nWrote {edge_path}\nWrote {summary_path}\nWrote {REPORT}")


if __name__ == "__main__":
    main()
