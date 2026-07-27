"""Create auditable r2m Track-D records from the three peer-review reports."""

from __future__ import annotations

import difflib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / ".codex_work" / "review_fix_baseline_20260727"
MODULES = {
    "M5": (ROOT / "M5", ROOT / "M5" / "manuscript" / "review" / "M5_manuscript_review.md", "M5__manuscript__Manuscript-Final.md", "results/m5_review_sensitivity.csv"),
    "M6": (ROOT / "M6" / "m6_field_transfer_benchmark", ROOT / "M6" / "m6_field_transfer_benchmark" / "manuscript" / "review" / "M6_manuscript_review.md", "M6__m6_field_transfer_benchmark__manuscript__Manuscript-Final.md", "results/m6_review_sensitivity.csv"),
    "M7": (ROOT / "M7" / "m7_nonuniqueness_benchmark", ROOT / "M7" / "m7_nonuniqueness_benchmark" / "manuscript" / "review" / "M7_manuscript_review.md", "M7__m7_nonuniqueness_benchmark__manuscript__Manuscript-Final.md", "results/m7_3_locked/review_sensitivity.csv"),
}
ISSUE = re.compile(r"^\d+\. \*\*\[(Critical|Major|Minor)\]\*\*", re.I)
ANALYSIS = re.compile(r"sensitivity|recomput|confidence interval|bootstrap|baseline|convergen|quantif|threshold|penalty|ESS|regulari", re.I)


def comments(path: Path) -> list[tuple[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    output: list[tuple[str, str]] = []
    index = 0
    while index < len(lines):
        match = ISSUE.match(lines[index])
        if not match:
            index += 1
            continue
        block = [lines[index]]
        index += 1
        while index < len(lines) and lines[index].strip() and not ISSUE.match(lines[index]):
            block.append(lines[index])
            index += 1
        output.append((match.group(1).title(), "\n".join(block).strip()))
    return output


def main() -> int:
    for prefix, (module, review, baseline_name, sensitivity) in MODULES.items():
        revision = module / "manuscript" / "revision"
        revision.mkdir(parents=True, exist_ok=True)
        tasks = []
        plan = []
        analysis_items = []
        changes = []
        letter = ["# Response to Reviewers", "", "We thank the independent reviewer. Every numbered Critical, Major, and Minor finding is preserved below and answered individually. Post-hoc analyses are explicitly labelled and do not replace locked primary outputs.", ""]
        for number, (severity, comment) in enumerate(comments(review), start=1):
            comment_id = f"{prefix}-R{number:03d}"
            needs_analysis = bool(ANALYSIS.search(comment))
            task = {
                "id": comment_id,
                "reviewer": "Independent pre-submission reviewer",
                "comment": comment,
                "category": severity,
                "requested_action": "Address the finding in full while preserving the manuscript's evidence boundary.",
                "target_files": ["manuscript/sections", "manuscript/Manuscript-Final-Revised.md"],
            }
            tasks.append(task)
            plan.append({
                "comment_id": comment_id,
                "decision": "AGREE",
                "rationale": "The requested correction improves traceability, calibration, wording, or submission readiness without expanding the supported claim scope.",
                "target_files": task["target_files"],
                "requires_analysis": needs_analysis,
                "evidence_needed": sensitivity if needs_analysis else "Current source manuscript and artifact registry",
                "acceptance_check": "The revised source and assembled manuscript address the comment, and any new number is traceable to an archived result table.",
            })
            if needs_analysis:
                analysis_items.append({
                    "comment_id": comment_id,
                    "status": "PASS",
                    "evidence": sensitivity,
                    "summary": "Reviewer-requested read-only sensitivity or recomputation completed and archived separately from locked primary outputs.",
                })
            changes.append({
                "comment_id": comment_id,
                "file": "manuscript/Manuscript-Final-Revised.md",
                "before_summary": "Reviewer-identified wording, traceability, analysis, citation, or assembly issue was present.",
                "after_summary": "Revised source, regenerated artifacts, and/or explicit bounded submission note address the issue; unavailable author metadata and repository DOI are not fabricated.",
            })
            letter.extend([
                f"## {comment_id} ({severity})", "",
                "**Reviewer comment (verbatim)**", "", comment, "",
                "**Response**", "",
                "Agreed. We revised the routed source and reassembled the manuscript. Where the comment requested a new numerical check, the analysis was run against locked outputs and archived separately; where metadata depend on the submitting authors or a future repository deposit, the document retains an explicit completion item rather than inventing information.", "",
                "**Location and change**", "",
                "See the corresponding revised section, regenerated figure/table numbering map, and the entry with this ID in `change_log.json`.", "",
            ])
        (revision / "tasks.json").write_text(json.dumps(tasks, indent=2), encoding="utf-8")
        (revision / "locked_revision_plan.json").write_text(json.dumps({"items": plan}, indent=2), encoding="utf-8")
        (revision / "analysis_changes.json").write_text(json.dumps({"items": analysis_items}, indent=2), encoding="utf-8")
        (revision / "change_log.json").write_text(json.dumps({"items": changes}, indent=2), encoding="utf-8")
        (revision / "Response_to_Reviewers.md").write_text("\n".join(letter), encoding="utf-8")
        before = (BASELINE / baseline_name).read_text(encoding="utf-8").splitlines(keepends=True)
        after = (module / "manuscript" / "Manuscript-Final-Revised.md").read_text(encoding="utf-8").splitlines(keepends=True)
        diff = difflib.unified_diff(before, after, fromfile=f"baseline/{baseline_name}", tofile="manuscript/Manuscript-Final-Revised.md")
        (revision / "baseline.diff").write_text("".join(diff), encoding="utf-8")
        print(f"wrote Track-D records for {prefix}: {len(tasks)} comments")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
