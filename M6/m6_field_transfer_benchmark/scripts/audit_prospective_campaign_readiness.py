"""Audit the four Ghana packages against the prospective field protocol.

This is a read-only gate.  It records what can run now (harmonisation,
chemistry transfer, ratio/geology sensitivity) and emits explicit ABSTAIN
statuses for independent age, direct-adjacency, and reaction-truth scoring.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[3]
M7_SCRIPT_DIR = REPO_ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(M7_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(M7_SCRIPT_DIR))

from field_data_inventory import inventory_field_data, render_markdown  # noqa: E402
from hydrosheaf.data.field import load_all_field_datasets  # noqa: E402
from hydrosheaf.validation.integrated_benchmark import (  # noqa: E402
    source_manifest_entry,
    write_json,
)


DEFAULT_CONFIG = REPO_ROOT / "M6" / "m6_field_transfer_benchmark" / "configs" / "ghana_prospective_campaign.json"
DEFAULT_OUTPUT = REPO_ROOT / ".codex_work" / "runs" / "RUN-GHANA-PROSPECTIVE-20260908-01"


def _dataset_statuses(
    inventory: dict[str, Any], canonical_datasets: dict[str, Any]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for dataset in inventory.get("datasets", []):
        name = str(dataset.get("dataset", ""))
        if not dataset.get("exists"):
            rows.append(
                {
                    "dataset": name,
                    "field_transfer": "INVALID",
                    "reaction_screening": "ABSTAIN",
                    "age": "ABSTAIN",
                    "direct_adjacency": "ABSTAIN",
                    "reaction_truth": "ABSTAIN",
                    "reason": "source file is missing",
                }
            )
            continue
        # The core inverse-reaction screen needs a declared, non-empty subset
        # of major ions.  SiO2/Sr/Fe are optional evidence lifts, not a reason
        # to disable all chemistry screening in the sparse external panels.
        normalised_name = "".join(character for character in name.lower() if character.isalnum())
        key = next(
            (
                candidate
                for candidate in canonical_datasets
                if "".join(character for character in candidate if character.isalnum())
                in normalised_name
            ),
            None,
        )
        coverage = canonical_datasets.get(key).coverage() if key else {}
        core_ions = ("Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3")
        core_ion_count = sum(int(coverage.get(ion, 0)) > 0 for ion in core_ions)
        chemistry_ready = core_ion_count >= 4
        rows.append(
            {
                "dataset": name,
                "field_transfer": "RUN",
                "reaction_screening": "RUN" if chemistry_ready else "ABSTAIN",
                "core_ion_count": core_ion_count,
                "age": "ABSTAIN",
                "direct_adjacency": "ABSTAIN",
                "reaction_truth": "ABSTAIN",
                "reason": (
                    "transfer/diagnostic use is allowed; no independent age, flow, "
                    "direct-edge, or reaction-truth labels are supplied"
                ),
            }
        )
    return rows


def _render_readiness(report: dict[str, Any]) -> str:
    lines = [
        "# Ghana prospective campaign readiness",
        "",
        "Read-only audit against `HS-GHANA-PROSPECTIVE-01`.",
        "",
        "| Package | Transfer | Reaction screening | Age | Direct adjacency | Reaction truth |",
        "|---|---|---|---|---|---|",
    ]
    for row in report["dataset_statuses"]:
        lines.append(
            f"| {row['dataset']} | {row['field_transfer']} | {row['reaction_screening']} | "
            f"{row['age']} | {row['direct_adjacency']} | {row['reaction_truth']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The four supplied packages can support harmonisation, chemistry transfer, "
            "ratio diagnostics, mapped-geology sensitivity, and hold-forward tests. "
            "Independent age, direct-adjacency, and reaction-truth scoring remains "
            "ABSTAIN until the protocol's prospective measurements and labels are "
            "collected without leakage.",
            "",
            "## Source inventory",
            "",
            render_markdown(report["field_inventory"]),
        ]
    )
    return "\n".join(lines)


def run(*, repo_root: Path = REPO_ROOT, config: Path = DEFAULT_CONFIG, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    config = config.resolve()
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty run directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    protocol = json.loads(config.read_text(encoding="utf-8"))
    field_inventory = inventory_field_data(repo_root)
    canonical_datasets = load_all_field_datasets(field_root=repo_root / "data" / "FieldData")
    statuses = _dataset_statuses(field_inventory, canonical_datasets)
    sources = []
    for item in protocol.get("current_field_packages", []):
        path = repo_root / str(item["path"])
        sources.append(
            source_manifest_entry(
                path,
                source_id=str(item["package_id"]),
                role=str(item["role"]),
            )
        )
    sources.append(
        source_manifest_entry(
            config,
            source_id="ghana_protocol",
            role="pre-registered prospective campaign contract",
        )
    )
    report: dict[str, Any] = {
        "schema": "ghana-prospective-readiness-v1",
        "run_id": output.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_id": protocol.get("protocol_id"),
        "protocol_status": protocol.get("status"),
        "source_manifest": sources,
        "dataset_statuses": statuses,
        "field_inventory": field_inventory,
        "module_status": {
            "canonical_loader_and_units": "RUN",
            "four_package_transfer": "RUN" if all(row["field_transfer"] == "RUN" for row in statuses) else "PARTIAL",
            "geochemical_ratios": "RUN",
            "mapped_geology_context": "RUN",
            "independent_age_accuracy": "ABSTAIN",
            "direct_adjacency_accuracy": "ABSTAIN",
            "independent_reaction_accuracy": "ABSTAIN",
            "reason_for_abstain": "No independent co-timed labels and no prospective blind-confirmation run are present.",
        },
        "claim_boundary": (
            "Current Ghana packages are field-transfer/context evidence. They do "
            "not validate groundwater age, direct adjacency, or unique reaction "
            "families."
        ),
    }
    write_json(output / "readiness.json", report)
    (output / "readiness.md").write_text(_render_readiness(report), encoding="utf-8")
    write_json(output / "protocol_snapshot.json", protocol)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(repo_root=args.repo_root, config=args.config, output=args.output)
    print(json.dumps(report["module_status"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["DEFAULT_CONFIG", "DEFAULT_OUTPUT", "main", "run"]
