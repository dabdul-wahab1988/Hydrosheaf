"""Build the complete captioned M5 supplementary figure/table source."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "M5" / "m5_inverse_reaction_benchmark"
OUTPUT = ROOT / "M5" / "manuscript" / "supplementary" / "Supplementary-Figures-and-Tables.md"

FIGURES = [
    ("figureS1_r_reaction_dictionary.png", "Reaction dictionary and signed stoichiometric structure."),
    ("figureS2_r_identifiability_diagnostics.png", "Structural identifiability diagnostics by ion panel."),
    ("figureS3_r_panel_method_heatmap.png", "Panel-by-method recovery heatmap and ion-ablation response."),
    ("figureS4_r_reaction_recovery.png", "Reaction-specific recovery under Hydrosheaf-Core evidence gates."),
    ("figureS5_r_measurement_value.png", "Held-out-ion prediction and retrospective next-best-measurement value."),
    ("figureS6_r_regularization_paths.png", "Regularisation paths across the predeclared tuning grid."),
    ("figureS7_r_bootstrap_support.png", "Bootstrap support stability for separated and collinear reactions."),
    ("figureS8_r_core_evidence_gates.png", "Hydrosheaf-Core evidence-gate scores and synthetic truth audit."),
    ("figureS9_r_data_tier_evidence.png", "Reaction evidence across Core, Plus-lite, and Enhanced data tiers."),
    ("figureS10_r_external_field_elri.png", "External-field evidence-lifted resolution audit."),
    ("figureS11_r_phreeqc_archetypes.png", "PHREEQC synthetic archetype behaviour and inverse-model multiplicity."),
]

TABLE_CAPTIONS = [
    "Reaction stoichiometry, families, signs, and saturation-index mappings.",
    "Exact signed reaction-equivalence classes.",
    "Scenario-level PHREEQC design and quality control.",
    "Predeclared hyperparameter grid and selected penalties.",
    "Complete model metrics by method, noise, panel, and archetype.",
    "Reaction-specific support confusion matrices.",
    "Ion-specific held-out error and realised measurement value.",
    "Thermodynamic-bound sensitivity and residual ambiguity.",
    "Scope-specific software environments and reproducibility metadata.",
    "Northern Ghana region-stratified chemistry-only transfer summary.",
    "Strict and relaxed PHREEQC inverse-model feasibility, multiplicity, and recovery.",
    "Hydrosheaf-Core evidence gates and synthetic recovery checks.",
    "Northern Ghana field evidence-gate support by region and reaction.",
    "Core, Plus-lite, and Enhanced measurement-tier recovery.",
    "Reaction-level optional diagnostic evidence by data tier.",
    "Evidence-lifted resolution for ambiguous reaction classes.",
    "External-field evidence-lifted resolution audit.",
]


def esc(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ").strip()


def table(path: Path) -> str:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    width = len(rows[0])
    lines = [
        "| " + " | ".join(esc(x) for x in rows[0]) + " |",
        "| " + " | ".join("---" for _ in range(width)) + " |",
    ]
    for row in rows[1:]:
        lines.append("| " + " | ".join(esc(x) for x in (row + [""] * width)[:width]) + " |")
    return "\n".join(lines)


def main() -> int:
    parts = ["# Supplementary Figures and Tables", ""]
    figure_dir = MODULE / "figures" / "r_publication" / "supplementary"
    for number, (name, caption) in enumerate(FIGURES, start=1):
        # r2m's converter resolves resources from the M5 project root, not from
        # this Markdown file's directory. Keep the resource path in that frame
        # and request a full text-width image in the generated DOCX.
        relative = (Path("m5_inverse_reaction_benchmark") / "figures" / "r_publication" / "supplementary" / name).as_posix()
        parts.extend([f"## Supplementary Figure S{number}", "", f"![Supplementary Figure S{number}. {caption}]({relative}){{width=6.5in}}", ""])
    table_dir = MODULE / "tables"
    for number, caption in enumerate(TABLE_CAPTIONS, start=1):
        path = next(table_dir.glob(f"tableS{number}_*.csv"))
        footnote = " PHREEQC software rows are scoped separately to the locked synthetic benchmark and field saturation-index extension." if number == 9 else ""
        parts.extend([f"## Supplementary Table S{number}", "", f"**Supplementary Table S{number}.** {caption}{footnote}", "", table(path), ""])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
