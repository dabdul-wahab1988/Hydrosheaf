"""Normalise the 29 July 2026 M7 review into auditable Track D contracts."""

from __future__ import annotations

import json
import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
REVIEW = PROJECT / "manuscript" / "review" / "M7_MANUSCRIPT_REVIEW_2026-07-29.md"
REVISION = PROJECT / "manuscript" / "revision"


def _target_files(comment: str) -> list[str]:
    lower = comment.lower()
    targets = ["manuscript/sections", "manuscript/Manuscript-Final-Revised.md"]
    if any(word in lower for word in ("supplement", "multiplicity", "sample-size", "precision")):
        targets.append("manuscript/supplementary")
    if any(word in lower for word in ("figure", "table", "printed width")):
        targets.extend(["figures/publication", "tables/publication"])
    if any(word in lower for word in ("runner", "hash", "source lock", "repository commit")):
        targets.extend(["scripts", "provenance"])
    if any(word in lower for word in ("literature", "references", "davis and goadrich")):
        targets.extend(["manuscript/methods/literature_search.json", "manuscript/LITERATURE.bib"])
    return list(dict.fromkeys(targets))


def _requires_analysis(comment: str) -> bool:
    lower = comment.lower()
    return any(
        phrase in lower
        for phrase in (
            "runner does not match",
            "sample-size or precision",
            "multiplicity",
            "simultaneous bootstrap",
            "confusion counts",
            "candidate recall",
        )
    )


def _decision(comment: str) -> str:
    lower = comment.lower()
    if any(
        phrase in lower
        for phrase in (
            "commit the complete",
            "tag a release",
            "deposit",
            "author contributions",
            "competing interests",
            "funding",
            "dataset doi",
            "software doi",
        )
    ):
        return "CLARIFY"
    return "AGREE"


def extract_comments(text: str) -> list[dict]:
    primary = text.split("## 1. Section-by-section critique", 1)[1].split(
        "## 2. Cross-section consistency audit", 1
    )[0]
    consistency = text.split("## 2. Cross-section consistency audit", 1)[1].split(
        "## 3. Research integrity", 1
    )[0]
    records: list[tuple[str, str, str]] = []
    for section_name, block in (("section", primary), ("consistency", consistency)):
        subsection = ""
        for line in block.splitlines():
            if line.startswith("### "):
                subsection = line.removeprefix("### ").strip()
                continue
            match = re.match(
                r"^\d+\. \*\*\[(Critical|Major|Minor)\]\*\*\s+(.*)$", line
            )
            if match:
                records.append((match.group(1), subsection or section_name, line))
    tasks = []
    for index, (category, subsection, comment) in enumerate(records, start=1):
        task_id = f"M7-20260729-R{index:02d}"
        tasks.append(
            {
                "id": task_id,
                "reviewer": "Independent manuscript reviewer, 29 July 2026",
                "comment": comment,
                "category": category,
                "subsection": subsection,
                "requested_action": (
                    "Address the finding in full, preserve the locked evidence boundary, "
                    "and record the exact verification evidence."
                ),
                "target_files": _target_files(comment),
            }
        )
    return tasks


def main() -> int:
    tasks = extract_comments(REVIEW.read_text(encoding="utf-8"))
    if len(tasks) != 36:
        raise RuntimeError(f"Expected 36 review comments, found {len(tasks)}.")
    plan = {
        "schema_version": "1.0",
        "source_review": REVIEW.relative_to(PROJECT).as_posix(),
        "locked_utc_date": "2026-07-30",
        "scope": (
            "All scientific, statistical, computational, provenance, figure, table, "
            "language and interoperability comments. Author identities, funding facts, "
            "persistent DOI minting and publication deposits remain external inputs/actions."
        ),
        "items": [
            {
                "comment_id": task["id"],
                "decision": _decision(task["comment"]),
                "rationale": (
                    "The technical finding is accepted. Where the requested action requires "
                    "author-supplied declarations or an external repository release, the "
                    "manuscript will state the verified current status without inventing it."
                ),
                "target_files": task["target_files"],
                "requires_analysis": _requires_analysis(task["comment"]),
                "evidence_needed": (
                    "A source diff plus a reproducible audit artifact or render/test record "
                    "appropriate to the finding."
                ),
                "acceptance_check": (
                    "The response matrix cites a current file or generated artifact, the "
                    "manuscript statement matches that evidence, and no broader claim remains."
                ),
            }
            for task in tasks
        ],
    }
    REVISION.mkdir(parents=True, exist_ok=True)
    (REVISION / "tasks.json").write_text(
        json.dumps(tasks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (REVISION / "locked_revision_plan.json").write_text(
        json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"comments": len(tasks), "plan_items": len(plan["items"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
