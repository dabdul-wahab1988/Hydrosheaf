"""Assemble and audit the M2_1 manuscript package.

The script is intentionally small and deterministic. It reads the locked
programme JSON files, copies no hidden state, regenerates the provenance table,
and assembles the authoritative section files into Manuscript-Final.md.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "M2_1"
RUN = ROOT / ".codex_work" / "programme-validation" / "RUN-INTEGRATION-FULL-20260802-15"
SECTION_DIR = PACKAGE / "manuscript" / "sections"
OUTPUT = PACKAGE / "manuscript" / "Manuscript-Final.md"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regenerate_evidence_tables() -> None:
    gates = load_json(RUN / "performance_gates.json")
    execution = load_json(RUN / "execution_gate.json")
    manifest = load_json(RUN / "run_manifest.json")
    artifacts = PACKAGE / "manuscript" / "artifacts"

    package_rows = [
        {"capability": "Topology and candidate connectivity", "package_module": "hydrosheaf.validation.topology; hydrosheaf.inference.network_fit", "evidence": "Component and integrated synthetic audits", "claim": "bounded candidate/topology inference"},
        {"capability": "Age inference", "package_module": "hydrosheaf.validation.age_competent_baseline; calibration; specialist_calibration", "evidence": "Age gate in locked programme run", "claim": "calibrated bounded synthetic age inference"},
        {"capability": "Reaction-family inference", "package_module": "hydrosheaf.validation.reaction_competent_baseline; reaction_rapm", "evidence": "Reaction gate in locked programme run", "claim": "calibrated family-level inference with selective risk"},
        {"capability": "Kinetic prediction and calibration", "package_module": "hydrosheaf.validation.kinetics_specialist_benchmark; hydrosheaf.reactive_transport.kinetic_phreeqc", "evidence": "M8 held-out kinetic gate", "claim": "predictive performance with conditional parameter identification"},
        {"capability": "Model discrepancy", "package_module": "hydrosheaf.validation.discrepancy", "evidence": "Integrated gate requires calibrated discrepancy", "claim": "disagreement is recorded before averaging"},
        {"capability": "Discrete model averaging", "package_module": "hydrosheaf.validation.model_averaging", "evidence": "Integrated gate requires convergence", "claim": "convergent bounded averaging"},
        {"capability": "Prospective measurement selection", "package_module": "hydrosheaf.validation.decision_utility", "evidence": "Six locked truth-blind policy cases", "claim": "bounded cost-adjusted synthetic utility"},
        {"capability": "Claim and failure gates", "package_module": "hydrosheaf.validation.performance_contract; programme_gate", "evidence": "Execution and performance gate records", "claim": "PASS/WEAK/ABSTAIN reporting"},
    ]
    write_csv(artifacts / "package_capability_inventory.csv", list(package_rows[0]), package_rows)

    rows = []
    for key in ("age", "reaction", "kinetics", "integrated"):
        gate = gates["gates"][key]
        observed = gate["observed"]
        if key == "age":
            headline = (
                f"coverage={observed['coverage_including_abstention']:.4f}; "
                f"specialist_MAE={observed['specialist_mae']:.4f}; "
                f"baseline_MAE={observed['baseline_mae']:.4f}"
            )
            calibration = observed["calibrated"]
        elif key == "reaction":
            headline = (
                f"log_loss={observed['multiclass_log_loss']:.4f}; "
                f"baseline_log_loss={observed['baseline_log_loss']:.4f}; "
                f"coverage={observed['coverage']:.4f}"
            )
            calibration = observed["classwise_calibrated"]
        elif key == "kinetics":
            headline = (
                f"RMSE={observed['predictive_rmse']:.4f}; "
                f"interval_coverage={observed['interval_coverage']:.4f}; "
                f"overall_identifiability={observed['identifiability_rate_overall']:.4f}"
            )
            calibration = observed["calibrated"]
        else:
            headline = (
                f"HydroSheaf_utility_per_cost={observed['hydrosheaf_mean_utility_per_cost']:.4f}; "
                f"random={observed['random_mean_utility_per_cost']:.4f}; "
                f"specialist={observed['strongest_specialist_mean_utility_per_cost']:.4f}"
            )
            calibration = observed["discrepancy_calibrated"] and observed["model_averaging_converged"]
        rows.append({
            "component": key,
            "status": gate["status"],
            "claim_scope": gate["requirements"]["claim_scope"],
            "headline_metric": headline,
            "calibration_or_control": str(calibration),
            "field_validation": "not_established",
        })
    write_csv(artifacts / "synthetic_performance_summary.csv", list(rows[0]), rows)

    source_hashes = manifest.get("source_hashes", {})
    hash_rows = [
        {
            "path": path,
            "sha256": value,
            "run_id": manifest.get("run_id", ""),
            "git_revision": manifest.get("git_revision", ""),
            "worktree_dirty": manifest.get("git_worktree_dirty", ""),
        }
        for path, value in sorted(source_hashes.items())
    ]
    write_csv(
        artifacts / "source_hash_inventory.csv",
        ["path", "sha256", "run_id", "git_revision", "worktree_dirty"],
        hash_rows,
    )

    architecture = """# Figure 1 source: HydroSheaf evidence-to-claim flow

```mermaid
flowchart LR
  A[Field or synthetic observations] --> B[Candidate generation]
  B --> C[Topology, age, reaction and kinetic specialists]
  C --> D[Calibration and selective-risk scoring]
  D --> E[Model discrepancy and convergent averaging]
  E --> F[Prospective decision utility]
  F --> G{Claim gate}
  G -->|PASS or bounded PASS| H[Report with provenance]
  G -->|WEAK or ABSTAIN| I[Report limitation and next measurement]
```
"""
    (artifacts / "architecture_flow.md").write_text(architecture, encoding="utf-8")

    if execution.get("status") != "PASS":
        raise RuntimeError("The locked execution gate is not PASS.")


def assemble() -> None:
    ordered = [
        "00-abstract/section.md",
        "01-introduction/section.md",
        "02-methods/section.md",
        "03-results/section.md",
        "04-discussion/section.md",
        "05-conclusion/section.md",
        "06-statements/section.md",
    ]
    contents = []
    for relative in ordered:
        path = SECTION_DIR / relative
        if not path.exists():
            raise FileNotFoundError(path)
        contents.append(path.read_text(encoding="utf-8").strip())
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n\n".join(contents) + "\n", encoding="utf-8")


def prose_word_count(path: Path) -> int:
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"\$\$.*?\$\$", " ", text, flags=re.S)
    text = re.sub(r"^\s*#.*$", " ", text, flags=re.M)
    text = re.sub(r"^\s*\|.*$", " ", text, flags=re.M)
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)", " ", text)
    text = re.sub(r"\[[^\]]*\]\([^)]*\)", " ", text)
    return len(re.findall(r"\b[\w'’-]+\b", text))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--count", action="store_true")
    args = parser.parse_args()
    regenerate_evidence_tables()
    if args.assemble:
        assemble()
    if args.count:
        if not OUTPUT.exists():
            assemble()
        supplement = PACKAGE / "manuscript" / "supplementary" / "Supplementary-Methods.md"
        print(json.dumps({"main": prose_word_count(OUTPUT), "supplementary_methods": prose_word_count(supplement)}, indent=2))


if __name__ == "__main__":
    main()
