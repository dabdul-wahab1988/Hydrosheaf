"""Remove implicit duplicate figure captions from assembled M7 Markdown."""

from __future__ import annotations

import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
DEFAULTS = (
    PROJECT / "manuscript" / "Manuscript-Final.md",
    PROJECT / "manuscript" / "Manuscript-Final-Revised.md",
)
IMAGE_LINE = re.compile(r"(?m)^!\[[^\n]*\]\(([^\n)]+)\)$")


def main() -> int:
    for path in DEFAULTS:
        text = path.read_text(encoding="utf-8")
        revised, count = IMAGE_LINE.subn(r"![](\1)", text)
        if count != 7:
            raise RuntimeError(f"Expected seven main figures in {path}, found {count}")
        path.write_text(revised, encoding="utf-8")
        print(f"Finalised {path} ({count} figure captions de-duplicated)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
