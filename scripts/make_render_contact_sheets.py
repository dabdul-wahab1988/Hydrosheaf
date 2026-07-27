"""Create contact sheets for complete page-by-page DOCX visual QA."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output")
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--rows", type=int, default=2)
    args = parser.parse_args()
    source = Path(args.source)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    pages = sorted(source.glob("page-*.png"))
    per_sheet = args.columns * args.rows
    for offset in range(0, len(pages), per_sheet):
        batch = pages[offset:offset + per_sheet]
        opened = [Image.open(path).convert("RGB") for path in batch]
        width = max(image.width for image in opened)
        height = max(image.height for image in opened)
        sheet = Image.new("RGB", (args.columns * width, args.rows * (height + 40)), "white")
        draw = ImageDraw.Draw(sheet)
        for index, (path, page) in enumerate(zip(batch, opened)):
            x = (index % args.columns) * width
            y = (index // args.columns) * (height + 40)
            sheet.paste(page, (x, y + 40))
            draw.text((x + 8, y + 8), path.stem, fill="black")
        sheet.save(output / f"sheet-{offset // per_sheet + 1:02d}.jpg", quality=90)
    print(f"{len(pages)} pages -> {len(list(output.glob('sheet-*.jpg')))} sheets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
