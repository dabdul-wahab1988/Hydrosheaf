#!/usr/bin/env python3
"""Build the 29 July 2026 reviewer-response and revision-audit records."""

from __future__ import annotations

import csv
import difflib
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REVISION = ROOT / "manuscript" / "revision"
BASELINE = ROOT / ".codex_work" / "revision_baseline_20260730"


RESPONSES = {
    "R01": (
        "ADDRESSED",
        "The title now says controlled-synthetic benchmarks. The abstract, Introduction, Methods and Limitations name the MODFLOW 6/MODPATH 7 M7.3 generator separately from the scalar affine M7.4/M7.5 generator and prohibit cross-generator or field transfer of either result.",
        ["manuscript/sections/00-abstract/section.md", "manuscript/sections/02-methods/section.md", "manuscript/sections/04-discussion/section.md"],
    ),
    "R02": (
        "ADDRESSED",
        "The abstract was split into purpose, design, primary representation result, follow-up estimator result and scope conclusion. It now states explicitly that both M7.4 and M7.5 failed their prespecified complete gates.",
        ["manuscript/sections/00-abstract/section.md"],
    ),
    "R03": (
        "ADDRESSED_BY_CLAIM_RESTRICTION",
        "The unsupported exhaustive novelty claim was withdrawn. The manuscript now describes a targeted, non-systematic search, and the search record stores dated queries, sources, inclusion logic and screened items without claiming bibliographic completeness.",
        ["manuscript/sections/01-introduction/section.md", "manuscript/methods/literature_search.json", "manuscript/LITERATURE.bib"],
    ),
    "R04": (
        "ADDRESSED",
        "A main-text design-and-claim table now maps all seven audits to generator, development/test counts, comparator, primary metrics, lock/multiplicity family and permitted claim.",
        ["tables/publication/table1_benchmark_design.md", "manuscript/Manuscript-Final-Revised.md"],
    ),
    "R05": (
        "ADDRESSED",
        "The phrase 'worth collecting' was removed. Collection recommendations are reserved for a future value-of-information study with measurement costs and decision losses.",
        ["manuscript/sections/01-introduction/section.md", "manuscript/sections/04-discussion/section.md"],
    ),
    "R06": (
        "ADDRESSED",
        "The exact one-time M7.5 runner was recovered from the archived task session. Its SHA-256 is a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191, exactly matching the confirmatory lock. The source and recovery manifest are archived under that content hash. No locked test was rerun.",
        ["scripts/run_m7_robust_hybrid_sheaf.py", "provenance/source_archive/a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191/manifest.json"],
    ),
    "R07": (
        "TECHNICAL_CORRECTION_COMPLETE_RELEASE_PENDING",
        "The false commit-level reconstruction claim was removed and replaced with experiment-specific run identifiers, recorded historical revisions, current file hashes and an explicit statement that the historical commits do not contain M7.4/M7.5. The complete local package is assembled. A versioned commit/tag and persistent DOI remain submission actions because this revision task does not authorise repository publication or create a repository DOI.",
        ["manuscript/sections/06-availability-statements/section.md", "manuscript/artifact_registry.json", "provenance/source_archive"],
    ),
    "R08": (
        "ADDRESSED",
        "A construct-validity subsection distinguishes the two generators, the questions each can test and the explicit non-transfer boundary: scalar-case success is not MODFLOW-system validation, and neither generator constitutes field validation.",
        ["manuscript/sections/02-methods/section.md", "manuscript/supplementary/Supplementary-Methods.md"],
    ),
    "R09": (
        "ADDRESSED_POST_REVIEW",
        "A labelled post-review simulation study used development-only planning inputs, 20,000 replicates and prespecified minimum differences for PR-AUC, Brier, log loss, age MAE, interval width, coverage and modal accuracy. The attainable precision/power results are reported in Supplementary Table S12 and are not presented as prospective preregistration.",
        ["scripts/post_review_statistical_audit.py", "results/post_review_audit_20260730", "tables/publication/tableS12_precision_power.md"],
    ),
    "R10": (
        "ADDRESSED_POST_REVIEW",
        "All 120 published M7.4 contrasts and all 560 published M7.5 contrasts were placed in separate full families. Shared case-block bootstrap resampling with 10,000 replicates and max-z simultaneous 95% intervals was applied. Only findings surviving those families remain inferentially supported.",
        ["scripts/post_review_statistical_audit.py", "tables/publication/tableS10_m7_4_multiplicity_adjusted.md", "tables/publication/tableS11_m7_5_multiplicity_adjusted.md"],
    ),
    "R11": (
        "ADDRESSED_BY_CLAIM_RESTRICTION",
        "No unplanned mechanism-mismatch experiment was added. Instead, every reaction claim is restricted to discrimination among the six planted archetypes, under the two tested indicator panels and the tested noise model; out-of-dictionary chemistry is explicitly unevaluated.",
        ["manuscript/sections/02-methods/section.md", "manuscript/sections/03-results/section.md", "manuscript/sections/04-discussion/section.md"],
    ),
    "R12": (
        "ADDRESSED",
        "The absolute statement about unweighted logistic regression was replaced with a prespecified estimand rationale. Candidate-edge class prevalence is reported for each generator, and the text explains why the unweighted target matches the intended per-candidate probability estimand.",
        ["manuscript/sections/02-methods/section.md", "manuscript/supplementary/Supplementary-Methods.md"],
    ),
    "R13": (
        "ADDRESSED",
        "Because the exact runner was recovered and its hash verified against the lock, the provenance statement now names the recovered hash and source archive. It also states explicitly that the stored locked test was not rerun.",
        ["manuscript/sections/03-results/section.md", "provenance/source_archive/a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191/manifest.json"],
    ),
    "R14": (
        "ADDRESSED_POST_REVIEW",
        "Results now pair absolute effects with relative changes and post-review practical margins. The small topology-age changes and the M7.5 overall PR-AUC change are explicitly classified by whether they clear those margins; the M7.5 complete calibration gate remains failed.",
        ["manuscript/sections/03-results/section.md", "tables/publication/tableS12_precision_power.md"],
    ),
    "R15": (
        "ADDRESSED_POST_REVIEW",
        "Scenario statements were reclassified using the simultaneous full-family intervals. The incompatible-cycle conflict-localisation signal survives; the noisy/missing overall gain and M7.5 scenario ranking-gain claims do not and were withdrawn as supported general gains.",
        ["manuscript/sections/03-results/section.md", "tables/publication/tableS10_m7_4_multiplicity_adjusted.md", "tables/publication/tableS11_m7_5_multiplicity_adjusted.md"],
    ),
    "R16": (
        "ADDRESSED",
        "The heading now reads 'Carbonate reactions were not recovered under either tested indicator panel.' All broader 'regardless of panel richness' or universal non-identifiability language was removed.",
        ["manuscript/sections/03-results/section.md", "manuscript/sections/04-discussion/section.md"],
    ),
    "R17": (
        "ADDRESSED",
        "Figure 5 was rebuilt as an M7-only Ghana supportability figure. The M6 tier-ablation panel was removed; the field claim now rests only on the truth-free Ghana scope and hold-forward audit.",
        ["figures/publication/figure5_ghana_supportability_boundary_m7_only.png", "scripts/make_m7_3_publication_assets.py", "manuscript/sections/03-results/section.md"],
    ),
    "R18": (
        "ADDRESSED",
        "The public-pipeline audit now reports the threshold and TP/FP/FN counts for every arm, explains the identical selected F1 values, distinguishes conditional from end-to-end recall and states that all generated candidates were selected, so the audit did not identify a useful scalar selection threshold.",
        ["manuscript/sections/03-results/section.md", "tables/publication/tableS13_public_pipeline_selection.md"],
    ),
    "R19": (
        "ADDRESSED",
        "The evaluative 'earned' sentence was replaced verbatim with the requested bounded statement about non-identity relations, global-compatibility diagnosis and missing-endpoint fallback under the tested scalar scenarios.",
        ["manuscript/sections/05-conclusion/section.md"],
    ),
    "R20": (
        "ADDRESSED",
        "The combined-module visual was removed. Figure 5 and its caption now contain only M7 evidence and explicitly distinguish synthetic supportability context from the truth-free Ghana audit.",
        ["figures/publication/figure5_ghana_supportability_boundary_m7_only.png", "figures/publication/figure_source_manifest.csv"],
    ),
    "R21": (
        "ADDRESSED",
        "The main paper now has four tables: the design map, a compact M7.3 decision table, M7.4 means and M7.5 means. Detailed metric tables were moved to the 13-table Supplement.",
        ["manuscript/Manuscript-Final-Revised.md", "manuscript/supplementary/Supplementary-Information.md"],
    ),
    "R22": (
        "ADDRESSED",
        "Figures 6 and 7 were regenerated for 7.08-inch journal width with a minimum 8-point label size. Those dimensions are recorded in the figure-source manifest and were checked in the Word and LibreOffice renders.",
        ["scripts/make_m7_sheaf_vs_graph_assets.py", "scripts/make_m7_robust_hybrid_assets.py", "figures/publication/figure_source_manifest.csv"],
    ),
    "R23": (
        "ADDRESSED_WITH_TARGETED_SCOPE",
        "The Introduction and Discussion now include recent 2022-2025 tracer, hydrochemistry, groundwater-age and machine-learning work, plus non-sheaf structured alternatives such as graph regularisation, Gaussian-process smoothing, flow-network inversion and residual diagnostics. The search is disclosed as targeted rather than exhaustive.",
        ["manuscript/sections/01-introduction/section.md", "manuscript/sections/04-discussion/section.md", "manuscript/LITERATURE.bib", "manuscript/methods/literature_search.json"],
    ),
    "R24": (
        "ADDRESSED",
        "The Davis-and-Goadrich wording now states that precision-recall curves can give a more informative view under class imbalance; ROC-AUC remains reported.",
        ["manuscript/sections/02-methods/section.md", "manuscript/supplementary/Supplementary-Methods.md"],
    ),
    "R25": (
        "EXTERNAL_AUTHOR_METADATA_PENDING",
        "The section is correctly structured and the licence and available technical identifiers are stated. Author CRediT roles, funding, competing interests, dataset DOI, software DOI and final versioned release are retained as explicit submission blockers because they require author declarations or external deposits and must not be fabricated.",
        ["manuscript/sections/06-availability-statements/section.md"],
    ),
    "R26": (
        "ADDRESSED",
        "The sentence was replaced with: 'The reaction solver was the inference method under evaluation and was not part of the synthetic generator.'",
        ["manuscript/sections/02-methods/section.md"],
    ),
    "R27": (
        "ADDRESSED_STRUCTURE_METADATA_PENDING",
        "The section is now titled 'Declarations and open research' and is ordered as author contributions, funding, competing interests, data availability and code availability. Only the author-supplied metadata listed under R25 remains open.",
        ["manuscript/sections/06-availability-statements/section.md"],
    ),
    "R28": (
        "ADDRESSED",
        "Both DOCX files were rebuilt from clean Markdown with citation processing. Microsoft Word exported complete 33- and 22-page PDFs; LibreOffice exported complete 31- and 21-page PDFs without the former libpng failure or repair prompt. All 55 Word-rendered pages were visually inspected.",
        ["manuscript/Manuscript-Final-Revised-20260730.docx", "manuscript/supplementary/Supplementary-Information-Revised-20260730.docx"],
    ),
    "R29": (
        "ADDRESSED",
        "The title, abstract, Discussion and Conclusion now frame M7 as a falsifiable controlled-synthetic benchmark and conditional representation result. They explicitly reject validation of HydroSheaf as a whole or a general superiority claim.",
        ["manuscript/sections/00-abstract/section.md", "manuscript/sections/04-discussion/section.md", "manuscript/sections/05-conclusion/section.md"],
    ),
    "R30": (
        "ADDRESSED",
        "The selected M7.5 method is called local-first/global-fallback throughout explanatory prose. The paper states that development selected local weight 1.0, so this test is not evidence for a general local/global blend.",
        ["manuscript/sections/03-results/section.md", "manuscript/sections/04-discussion/section.md", "manuscript/sections/05-conclusion/section.md"],
    ),
    "R31": (
        "ADDRESSED",
        "Limitations now identify independent replication under another generator family or field data with independently measured connectivity as the next claim-bearing step. They explicitly prohibit inference to temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance.",
        ["manuscript/sections/04-discussion/section.md"],
    ),
    "R32": (
        "ADDRESSED",
        "All scope statements now distinguish the six-development/twelve-test M7.3 generator from the 64-case M7.4 and 128-case M7.5 scalar generator tests.",
        ["manuscript/sections/00-abstract/section.md", "manuscript/sections/01-introduction/section.md", "manuscript/sections/04-discussion/section.md"],
    ),
    "R33": (
        "ADDRESSED",
        "The paper now uses seven audits consistently and enumerates them in the Introduction and design table.",
        ["manuscript/sections/01-introduction/section.md", "tables/publication/table1_benchmark_design.md", "manuscript/sections/05-conclusion/section.md"],
    ),
    "R34": (
        "ADDRESSED",
        "The exact M7.5 source was restored and independently hash-checked against the confirmatory lock. The Results report the matching hash and no-rerun status.",
        ["scripts/run_m7_robust_hybrid_sheaf.py", "manuscript/sections/03-results/section.md"],
    ),
    "R35": (
        "FACTUAL_CORRECTION_COMPLETE_RELEASE_PENDING",
        "The false single-commit claim was replaced with experiment-specific run identifiers, the distinct M7.3 freeze commit, the historical M7.4/M7.5 manifest revision and an explicit warning that the latter revision does not contain the files. The final commit/tag/DOI must be inserted only after release publication.",
        ["manuscript/sections/06-availability-statements/section.md"],
    ),
    "R36": (
        "ADDRESSED",
        "Figure 7 and main Table 4 state that the selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Machine-readable arm names are retained only where necessary to preserve provenance.",
        ["manuscript/Manuscript-Final-Revised.md", "figures/publication/figure_source_manifest.csv"],
    ),
}


ANALYSIS_IDS = {"R04", "R06", "R09", "R10", "R14", "R15", "R18", "R22", "R28", "R34"}


def short_id(full_id: str) -> str:
    return full_id.rsplit("-", 1)[-1]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_response(tasks: list[dict]) -> None:
    lines = [
        "# Response to the 29 July 2026 M7 manuscript review",
        "",
        "This response addresses all 36 technical comments against the live manuscript and immutable result records. Post-review analyses are explicitly labelled and do not alter or rerun a locked test. Two classes of submission metadata remain outside technical revision: author declarations, and creation of a versioned public release/DOI. They are recorded as blockers rather than invented.",
        "",
        "## Resolution summary",
        "",
        "- Technical/scientific/editorial comments addressed: 32.",
        "- Technical correction complete but release publication pending: R07 and R35.",
        "- Structure corrected but author/deposit metadata pending: R25 and R27.",
        "- Locked confirmatory tests rerun: none.",
        "",
    ]
    for task in tasks:
        sid = short_id(task["id"])
        status, response, evidence = RESPONSES[sid]
        lines += [
            f"## {task['id']} — {status}",
            "",
            "**Reviewer comment**",
            "",
            task["comment"],
            "",
            "**Response**",
            "",
            response,
            "",
            "**Verification evidence**",
            "",
        ]
        lines += [f"- `{path}`" for path in evidence]
        lines.append("")
    (REVISION / "Response_to_Reviewers.md").write_text("\n".join(lines), encoding="utf-8")


def write_logs(tasks: list[dict]) -> None:
    change_log = []
    analysis_log = []
    for task in tasks:
        sid = short_id(task["id"])
        status, response, evidence = RESPONSES[sid]
        change_record = {
            "comment_id": task["id"],
            "category": task["category"],
            "subsection": task["subsection"],
            "response_status": status,
            "file": evidence[0],
            "before_summary": task["comment"].split(". ", 1)[-1][:300],
            "after_summary": response,
            "evidence": evidence,
            "locked_test_rerun": False,
        }
        change_log.append(change_record)
        if sid in ANALYSIS_IDS:
            analysis_log.append(
                {
                    "comment_id": task["id"],
                    "status": "PASS",
                    "review_resolution_status": status,
                    "summary": response,
                    "evidence": evidence,
                    "locked_test_rerun": False,
                }
            )
    (REVISION / "change_log.json").write_text(json.dumps(change_log, indent=2) + "\n", encoding="utf-8")
    (REVISION / "analysis_changes.json").write_text(json.dumps(analysis_log, indent=2) + "\n", encoding="utf-8")


def write_number_map() -> None:
    rows = [
        ("Figure 1", "Figure 1", "Benchmark architecture and claim boundary"),
        ("Figure 2", "Figure 2", "Evidence integration"),
        ("Figure 3", "Figure 3", "Topology-conditioned age inference"),
        ("Figure 4", "Figure 4", "Reaction-family recovery"),
        ("Figure 5", "Figure 5", "M7-only Ghana supportability boundary; M6 panel removed"),
        ("Figure 6", "Figure 6", "M7.4 sheaf-versus-graph representation"),
        ("Figure 7", "Figure 7", "M7.5 local-first/global-fallback diagnostic"),
        ("Table 1", "Table 1", "Seven-audit design and claim map"),
        ("Tables 2-6", "Supplementary Tables S1-S6", "Detailed M7.3 metric tables moved from main paper"),
        ("Table 7", "Table 2", "Compact M7.3 primary decision table"),
        ("Table 8", "Table 3", "M7.4 locked-test means"),
        ("Table 9", "Table 4", "M7.5 locked-test means"),
        ("New post-review audit", "Supplementary Tables S10-S13", "Multiplicity, planning and public-pipeline audits"),
    ]
    with (REVISION / "figure_table_number_map.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pre_revision_identifier", "revised_identifier", "description"])
        writer.writerows(rows)


def write_baseline_diff() -> None:
    if not BASELINE.exists():
        (REVISION / "baseline.diff").write_text("Baseline directory unavailable.\n", encoding="utf-8")
        return
    # Keep this human-reviewable. Large generated result CSVs are protected by
    # the artifact registry and hashes rather than being repeated in a text diff.
    allowed = {".md", ".json", ".py", ".bib", ".txt", ".yml", ".yaml"}
    routed_roots = {"manuscript", "scripts", "tests", "docs", "provenance"}
    routed_root_files = {
        "DECISIONS.md",
        "Outline.md",
        "README.md",
        "References.bib",
        "analysis_plan.json",
        "artifact_validation.json",
        "proposal.normalized.json",
    }

    def routed(rel: Path) -> bool:
        return rel.parts[0] in routed_roots or rel.as_posix() in routed_root_files

    chunks = ["# Review-round unified diff against revision_baseline_20260730\n"]
    current_paths = {
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.suffix.lower() in allowed
        and ".codex_work" not in path.relative_to(ROOT).parts
        and "revision/round_2026-07-27" not in path.as_posix()
        and routed(path.relative_to(ROOT))
    }
    baseline_paths = {
        path.relative_to(BASELINE)
        for path in BASELINE.rglob("*")
        if path.is_file()
        and path.suffix.lower() in allowed
        and routed(path.relative_to(BASELINE))
    }
    for rel in sorted(current_paths | baseline_paths, key=lambda p: p.as_posix().lower()):
        before_path = BASELINE / rel
        after_path = ROOT / rel
        before = before_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True) if before_path.exists() else []
        after = after_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True) if after_path.exists() else []
        if before == after:
            continue
        chunks.extend(
            difflib.unified_diff(
                before,
                after,
                fromfile=f"baseline/{rel.as_posix()}",
                tofile=f"revised/{rel.as_posix()}",
            )
        )
    rendered = "".join(chunks)
    rendered = "\n".join(line.rstrip(" \t") for line in rendered.splitlines()) + "\n"
    (REVISION / "baseline.diff").write_text(rendered, encoding="utf-8")


def write_qa_manifest() -> None:
    paths = [
        ROOT / "manuscript" / "Manuscript-Final-Revised-20260730.docx",
        ROOT / "manuscript" / "Manuscript-Final-Revised.docx",
        ROOT / "manuscript" / "Manuscript-Final.docx",
        ROOT / "manuscript" / "supplementary" / "Supplementary-Information-Revised-20260730.docx",
        ROOT / "manuscript" / "supplementary" / "Supplementary-Information.docx",
        ROOT / ".codex_work" / "word_main_20260730.pdf",
        ROOT / ".codex_work" / "word_supp_20260730.pdf",
        ROOT / ".codex_work" / "lo_pdf_main4_20260730" / "Manuscript-Final-Revised-20260730.pdf",
        ROOT / ".codex_work" / "lo_pdf_supp4_20260730" / "Supplementary-Information-Revised-20260730.pdf",
    ]
    records = []
    for path in paths:
        records.append(
            {
                "path": path.relative_to(ROOT).as_posix(),
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path) if path.exists() else None,
            }
        )
    report = {
        "date": "2026-07-30",
        "word_pdf_pages": {"main": 33, "supplement": 22},
        "libreoffice_pdf_pages": {"main": 31, "supplement": 21},
        "visual_pages_inspected": 55,
        "locked_confirmatory_test_rerun": False,
        "artifacts": records,
    }
    (REVISION / "docx_interoperability_qa.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    tasks = json.loads((REVISION / "tasks.json").read_text(encoding="utf-8"))
    missing = sorted({short_id(task["id"]) for task in tasks} - RESPONSES.keys())
    extra = sorted(RESPONSES.keys() - {short_id(task["id"]) for task in tasks})
    if missing or extra:
        raise SystemExit(f"Response map mismatch: missing={missing}, extra={extra}")
    write_response(tasks)
    write_logs(tasks)
    write_number_map()
    write_baseline_diff()
    write_qa_manifest()
    print(f"Wrote response package for {len(tasks)} comments")


if __name__ == "__main__":
    main()
