"""Assemble r2m manuscript sections without inline artifact corruption.

Unlike the stock preview assembler, this project-local assembler includes every
CSV row and treats an inline embed token as an in-text reference followed by a
block-level artifact.  This prevents captions and tables from splitting a
sentence and makes the generated Markdown the authoritative DOCX source.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


TOKEN_RE = re.compile(r"\[\[\s*(FIG|TAB|FIGREF|TABREF):([A-Za-z0-9_.-]+)\s*\]\]")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output", default="manuscript/Manuscript-Final.md")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def manuscript_title(root: Path) -> str:
    proposal = root / "proposal.normalized.json"
    if proposal.is_file():
        value = str(read_json(proposal).get("title", "")).strip()
        if value:
            return value
    outline = root / "Outline.md"
    if outline.is_file():
        match = re.search(r"(?m)^#\s+(.+?)\s*$", outline.read_text(encoding="utf-8"))
        if match:
            return match.group(1).strip()
    raise ValueError(f"No manuscript title found under {root}")


def section_text(root: Path) -> str:
    section_root = root / "manuscript" / "sections"
    directories = sorted(
        path for path in section_root.iterdir()
        if path.is_dir() and re.fullmatch(r"\d{2}-.+", path.name)
    )
    if not directories:
        raise ValueError(f"No manuscript sections found under {section_root}")
    parts: list[str] = []
    for directory in directories:
        source = directory / "section.md"
        if not source.is_file():
            raise ValueError(f"Missing section source: {source}")
        parts.append(source.read_text(encoding="utf-8").strip())
        for subsection in sorted(directory.glob("subsection-*.md")):
            parts.append(subsection.read_text(encoding="utf-8").strip())
    return "\n\n".join(part for part in parts if part)


def registry(root: Path) -> dict[str, dict]:
    path = root / "manuscript" / "artifact_registry.json"
    records = read_json(path).get("artifacts", [])
    output: dict[str, dict] = {}
    for record in records:
        artifact_id = str(record.get("id", "")).strip()
        if not artifact_id or artifact_id in output:
            raise ValueError(f"Invalid or duplicate artifact id: {artifact_id!r}")
        artifact_path = (root / str(record["path"])).resolve()
        artifact_path.relative_to(root.resolve())
        if not artifact_path.is_file():
            raise ValueError(f"Missing artifact: {artifact_path}")
        item = dict(record)
        item["_path"] = artifact_path
        output[artifact_id] = item
    return output


def escape_cell(value: str) -> str:
    return value.replace("\\", "\\\\").replace("|", "\\|").replace("\r", " ").replace("\n", " ").strip()


def table_markdown(path: Path) -> str:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return "_(No rows.)_"
    width = len(rows[0])
    output = [
        "| " + " | ".join(escape_cell(cell) for cell in rows[0]) + " |",
        "| " + " | ".join("---" for _ in range(width)) + " |",
    ]
    for row in rows[1:]:
        normalized = (row + [""] * width)[:width]
        output.append("| " + " | ".join(escape_cell(cell) for cell in normalized) + " |")
    return "\n".join(output)


def embed(record: dict, label: str, root: Path) -> str:
    path: Path = record["_path"]
    relative = path.relative_to(root).as_posix()
    caption = str(record["caption"]).strip().rstrip(".")
    if record["kind"] == "figure":
        # Pandoc treats image alternative text as the single visible caption.
        return f"![{label}. {caption}.]({relative})"
    if path.suffix.lower() == ".csv":
        return f"**{label}.** {caption}.\n\n{table_markdown(path)}"
    return f"**{label}.** {caption}.\n\nSource: `{relative}`"


def hydrate(text: str, records: dict[str, dict], root: Path) -> tuple[str, list[dict[str, str]]]:
    matches = list(TOKEN_RE.finditer(text))
    used: list[str] = []
    for match in matches:
        tag, artifact_id = match.groups()
        if artifact_id not in records:
            raise ValueError(f"Unregistered artifact token: {artifact_id}")
        expected = "figure" if tag.startswith("FIG") else "table"
        if records[artifact_id].get("kind") != expected:
            raise ValueError(f"Artifact kind mismatch: {artifact_id}")
        if artifact_id not in used:
            used.append(artifact_id)

    counters = {"figure": 0, "table": 0}
    labels: dict[str, str] = {}
    for artifact_id in used:
        kind = records[artifact_id]["kind"]
        counters[kind] += 1
        labels[artifact_id] = f"{'Figure' if kind == 'figure' else 'Table'} {counters[kind]}"

    embedded: set[str] = set()
    output_parts: list[str] = []
    for paragraph in re.split(r"\n\s*\n", text.strip()):
        queued: list[str] = []

        def replace(match: re.Match[str]) -> str:
            tag, artifact_id = match.groups()
            label = labels[artifact_id]
            if tag.endswith("REF") or artifact_id in embedded:
                return label
            embedded.add(artifact_id)
            block = embed(records[artifact_id], label, root)
            if paragraph.strip() == match.group(0):
                queued.insert(0, block)
                return ""
            queued.append(block)
            return label

        replaced = TOKEN_RE.sub(replace, paragraph).strip()
        if replaced:
            output_parts.append(replaced)
        output_parts.extend(block for block in queued if block)

    missing = [artifact_id for artifact_id in used if artifact_id not in embedded]
    if missing:
        raise ValueError(f"Referenced artifacts have no embed token: {missing}")

    map_rows = [
        {
            "artifact_id": artifact_id,
            "kind": records[artifact_id]["kind"],
            "display_label": labels[artifact_id],
            "artifact_path": records[artifact_id]["_path"].relative_to(root).as_posix(),
        }
        for artifact_id in used
    ]
    return "\n\n".join(output_parts) + "\n", map_rows


def main() -> int:
    args = arguments()
    root = Path(args.project_root).resolve()
    output = (root / args.output).resolve()
    output.relative_to(root)
    records = registry(root)
    body, map_rows = hydrate(section_text(root), records, root)
    output.write_text(f"# {manuscript_title(root)}\n\n{body}", encoding="utf-8")
    revised = root / "manuscript" / "Manuscript-Final-Revised.md"
    revised.write_text(output.read_text(encoding="utf-8"), encoding="utf-8")
    map_path = root / "manuscript" / "revision" / "figure_table_number_map.csv"
    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["artifact_id", "kind", "display_label", "artifact_path"])
        writer.writeheader()
        writer.writerows(map_rows)
    print(f"assembled {output} with {len(map_rows)} artifacts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
