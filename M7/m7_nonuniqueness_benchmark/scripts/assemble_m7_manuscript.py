"""Assemble the revised M7 manuscript from section and locked artifact sources."""

from __future__ import annotations

import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript"
SECTIONS = MANUSCRIPT / "sections"
TABLES = PROJECT / "tables" / "publication"
OUTPUTS = (
    MANUSCRIPT / "Manuscript-Final-Revised.md",
    MANUSCRIPT / "Manuscript-Final.md",
)

TITLE = (
    "Conditional evidence integration and the incremental contribution of sheaf "
    "structure in controlled-synthetic groundwater benchmarks"
)

SECTION_PATHS = (
    SECTIONS / "00-abstract" / "section.md",
    SECTIONS / "01-introduction" / "section.md",
    SECTIONS / "02-methods" / "section.md",
    SECTIONS / "02-review-robustness" / "section.md",
    SECTIONS / "03-results" / "section.md",
    SECTIONS / "04-discussion" / "section.md",
    SECTIONS / "05-conclusion" / "section.md",
    SECTIONS / "06-availability-statements" / "section.md",
)

FIGURES = {
    "FIG-1": (
        1,
        "figure1_benchmark_and_claim_design.png",
        "Benchmark architecture and claim boundary. The independent synthetic-truth "
        "branch separates official MODFLOW 6/MODPATH 7 flow and pathline generation, "
        "nonlinear chemistry and tracer generation, truth-blind inference and locked "
        "adverse-control scoring from the Northern Ghana branch, which distinguishes "
        "supportable component diagnostics from field interpretations that remain "
        "non-identifiable under the available evidence.",
    ),
    "FIG-2": (
        2,
        "figure2_evidence_integration.png",
        "Evidence integration is conditional. Panel a compares native hydraulic (H), "
        "age (A), chemistry (C), pairwise and fully integrated topology-ranking "
        "performance. Panel b shows case-block incremental PR-AUC contrasts with 95% "
        "bootstrap intervals. Panels c and d show that permuted evidence can reduce "
        "entropy while degrading discrimination and calibration.",
    ),
    "FIG-3": (
        3,
        "figure3_topology_conditions_age.png",
        "Correct topology improves model-conditioned age inference. Panels a and b "
        "show MAE and interval-width changes, panel c reports effective-sample-size "
        "fractions and panel d shows coverage changes. Reversed tritium-only accuracy "
        "contrasts are not interpreted because 8 of 12 cases failed the stability rule.",
    ),
    "FIG-4": (
        4,
        "figure4_reaction_nonuniqueness.png",
        "Recovery among the six planted reaction archetypes is process-dependent. The "
        "core and enhanced panels are compared for modal accuracy, probability assigned "
        "to the planted family, support entropy and effective support. Carbonate "
        "reactions were not recovered under either tested panel.",
    ),
    "FIG-5": (
        5,
        "figure5_ghana_supportability_boundary_m7_only.png",
        "Northern Ghana evidence and claim boundary using evidence from this study only. Panels a and "
        "b map observed evidence to defensible claims; panel c reports 160 wells, 140 "
        "complete seasonal pairs and 320 seasonal observations; panel d reports the "
        "truth-free seasonal hold-forward comparison. No external field-transfer experiment is included.",
    ),
    "FIG-6": (
        6,
        "figure6_sheaf_vs_weighted_graph.png",
        "Incremental contribution of affine sheaf structure. Overall and selected "
        "scenario contrasts compare edge-local, identity, native affine and permuted-map "
        "arms. Error bars are simultaneous 95% intervals controlling all 120 published "
        "representation-benchmark contrasts as one family. Intended printed width is 7.08 in; minimum label "
        "size is 8 pt.",
    ),
    "FIG-7": (
        7,
        "figure7_robust_hybrid_sheaf.png",
        "Local-first/global-fallback estimator diagnostic. The selected nominal hybrid had local weight 1.0 and "
        "is therefore local-first/global-fallback. Error bars are simultaneous 95% "
        "intervals controlling all 560 published estimator-diagnostic contrasts as one family. Intended "
        "printed width is 7.08 in; minimum label size is 8 pt.",
    ),
}

TABLES_MAIN = {
    "TAB-1": (
        1,
        "table1_benchmark_design.md",
        "Design and claim map for the seven audits, distinguishing generator, sample, "
        "claim-bearing comparator, outcomes, lock, multiplicity family and permitted claim.",
    ),
    "TAB-7": (
        2,
        "table2_primary_m7_3_decision.md",
        "Primary process-based integration decision table. Detailed metrics and complete machine-readable "
        "contrasts are reported in Supplementary Tables S1--S6.",
    ),
    "TAB-8": (
        3,
        "table8_sheaf_vs_weighted_graph.md",
        "Locked-test performance of the competence-matched edge-local graph, identity "
        "graph Laplacian, affine sheaf and permuted-map control across 64 cases. "
        "Family-wise contrasts are in Supplementary Table S10.",
    ),
    "TAB-9": (
        4,
        "table9_robust_hybrid_sheaf.md",
        "Locked-test performance of seven estimator-diagnostic arms across 128 fresh cases. The selected "
        "local-first arms had local weight 1.0; family-wise contrasts are in "
        "Supplementary Table S11.",
    ),
}


def _table_body(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip()
    text = re.sub(r"^# Table[^\n]*\n+", "", text)
    return text


def _replace_tokens(text: str) -> str:
    for key, (number, filename, caption) in FIGURES.items():
        text = text.replace(f"[[FIGREF:{key}]]", f"Figure {number}")
        block = (
            f"![](figures/publication/{filename})\n\n"
            f"**Figure {number}.** {caption}"
        )
        text = text.replace(f"[[FIG:{key}]]", block)
    for key, (number, filename, caption) in TABLES_MAIN.items():
        text = text.replace(f"[[TABREF:{key}]]", f"Table {number}")
        block = f"**Table {number}.** {caption}\n\n{_table_body(TABLES / filename)}"
        text = text.replace(f"[[TAB:{key}]]", block)
    unresolved = re.findall(r"\[\[(?:FIG|FIGREF|TAB|TABREF):[^]]+\]\]", text)
    if unresolved:
        raise RuntimeError(f"Unresolved artifact tokens: {sorted(set(unresolved))}")
    return text


def main() -> int:
    body = "\n\n".join(path.read_text(encoding="utf-8").strip() for path in SECTION_PATHS)
    assembled = _replace_tokens(f"# {TITLE}\n\n{body}\n")
    if assembled.count("![](") != 7:
        raise RuntimeError("Expected exactly seven main-text figures")
    for output in OUTPUTS:
        output.write_text(assembled, encoding="utf-8")
        print(f"M7 manuscript -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
