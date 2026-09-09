"""Move main-manuscript display items to the post-reference section.

The authoring source keeps tables and figure captions in explicit Markdown
blocks.  The submission copy requested for this revision places those blocks
after the references, followed by the embedded figure pages.  This helper
performs the same deterministic move in Markdown and in already-built DOCX
files without rewriting table cells, captions, images, or tracked colour
markup.

Examples
--------
python reorder_display_items.py --markdown Manuscript-Final-Revised.md
python reorder_display_items.py --docx Manuscript-Final-Revised.docx
python reorder_display_items.py --docx Manuscript-Colored-Changes.docx
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from docx import Document


def _normalise(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace("\u00a0", " ")).strip()


def _markdown_heading_indices(lines: list[str]) -> dict[str, int]:
    headings: dict[str, int] = {}
    for i, line in enumerate(lines):
        if line.startswith("## "):
            title = line[3:].strip()
            if title in {"Tables", "Figure captions", "5. Discussion", "References"}:
                headings[title] = i
    return headings


def reorder_markdown(path: Path) -> bool:
    """Move Tables and Figure captions after References; return whether changed."""
    original = path.read_text(encoding="utf-8")
    lines = original.splitlines(keepends=True)
    headings = _markdown_heading_indices(lines)
    required = {"Tables", "Figure captions", "5. Discussion", "References"}
    missing = required - headings.keys()
    if missing:
        raise ValueError(f"{path}: missing headings {sorted(missing)}")

    tables_i = headings["Tables"]
    figure_i = headings["Figure captions"]
    discussion_i = headings["5. Discussion"]
    references_i = headings["References"]

    if tables_i > discussion_i:
        # Already in post-Discussion position.  The source is considered
        # canonical only when both display blocks occur after References.
        if tables_i > references_i and figure_i > references_i:
            return False
        raise ValueError(f"{path}: display blocks are not in a movable order")
    if not (tables_i < figure_i < discussion_i):
        raise ValueError(f"{path}: expected Tables, Figure captions, Discussion order")

    display = lines[tables_i:discussion_i]
    retained = lines[:tables_i] + lines[discussion_i:]
    new_references_i = _markdown_heading_indices(retained)["References"]
    before = retained[:new_references_i]
    references_and_tail = retained[new_references_i:]
    prefix = "".join(before).rstrip("\n") + "\n\n"
    middle = "".join(references_and_tail).rstrip("\n")
    display_text = "".join(display).strip("\n")
    updated = prefix + middle + "\n\n" + display_text + "\n"
    if updated == original:
        return False
    path.write_text(updated, encoding="utf-8", newline="")
    return True


def _element_text(element) -> str:
    return _normalise(" ".join(element.itertext()))


def _is_paragraph(element) -> bool:
    return element.tag.rsplit("}", 1)[-1] == "p"


def _has_drawing(element) -> bool:
    return any(child.tag.rsplit("}", 1)[-1] == "drawing" for child in element.iter())


def _find_heading(body, pattern: re.Pattern[str]) -> int:
    for i, element in enumerate(body.iterchildren()):
        if _is_paragraph(element) and pattern.match(_element_text(element)):
            return i
    raise ValueError(f"DOCX: heading matching {pattern.pattern!r} not found")


def reorder_docx(path: Path) -> bool:
    """Move the display-item XML block before the first embedded figure page."""
    doc = Document(path)
    body = doc.element.body
    children = list(body.iterchildren())

    tables_i = _find_heading(body, re.compile(r"^Tables(?: Tables)*$"))
    discussion_i = _find_heading(body, re.compile(r"^5\. Discussion(?: 5\. Discussion)*$"))

    if tables_i > discussion_i:
        # A second invocation is harmless when the display block is already
        # after Discussion.  Validate that it is also after the references.
        references_i = _find_heading(body, re.compile(r"^References(?: References)*$"))
        if tables_i > references_i:
            return False
        raise ValueError(f"{path}: display block is after Discussion but before References")

    # Coloured DOCX files wrap changed blocks in bookmarkStart/bookmarkEnd
    # elements.  Include a directly preceding bookmarkStart so its pair moves
    # with the display block; clean DOCX files simply start at the heading.
    start = tables_i
    if start > 0 and children[start - 1].tag.rsplit("}", 1)[-1] == "bookmarkStart":
        start -= 1
    block = children[start:discussion_i]

    # Figure pages are existing image paragraphs appended by the DOCX builder.
    # Insert captions/tables immediately before the first such paragraph, so
    # the final order is References -> Tables -> Figure captions -> Figures.
    anchor = None
    for element in children:
        if _is_paragraph(element) and _has_drawing(element):
            anchor = element
            break
    if anchor is None:
        # Fallback for a DOCX without embedded figures: insert before sectPr.
        for element in children:
            if element.tag.rsplit("}", 1)[-1] == "sectPr":
                anchor = element
                break
    if anchor is None:
        raise ValueError(f"{path}: cannot find an insertion anchor")

    for element in block:
        body.remove(element)
    anchor_index = body.index(anchor)
    for offset, element in enumerate(block):
        body.insert(anchor_index + offset, element)

    doc.save(path)
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markdown", type=Path, help="Markdown source to reorder in place")
    parser.add_argument("--docx", type=Path, action="append", help="DOCX to reorder in place; repeatable")
    args = parser.parse_args()
    if args.markdown is None and not args.docx:
        parser.error("provide --markdown and/or --docx")

    changed = []
    if args.markdown is not None and reorder_markdown(args.markdown):
        changed.append(str(args.markdown))
    for path in args.docx or []:
        if reorder_docx(path):
            changed.append(str(path))
    print("reordered:" if changed else "already ordered:", ", ".join(changed) or "none")


if __name__ == "__main__":
    main()
