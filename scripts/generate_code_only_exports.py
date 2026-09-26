#!/usr/bin/env python3
"""Generate source-only code exports for HydroSheaf.

Produces two clean, reproducible source code text bundles in Code_Only_Exports/:
1. Hydrosheaf_Core_Code.txt:
   All reusable Python source files under hydrosheaf/ (the core package).
2. Hydrosheaf_Analysis_Code.txt:
   All analysis, benchmark, validation, milestone (M1-M9, O3, O4), scripts,
   and test Python files across the repository.

Usage:
    python scripts/generate_code_only_exports.py
    python scripts/generate_code_only_exports.py --check
    python scripts/generate_code_only_exports.py --core-only
    python scripts/generate_code_only_exports.py --analysis-only
"""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "Code_Only_Exports"

# Excluded directory components (caches, build artifacts, environments, ephemeral files)
EXCLUDE_PARTS = {
    ".git",
    "__pycache__",
    "venv",
    ".venv",
    "env",
    ".pytest_cache",
    "pymc_cache",
    "build",
    "dist",
    "hydrosheaf.egg-info",
    ".idea",
    ".vscode",
    "scratch",
    "tmp",
    ".codex_work",
    ".reasonix",
    ".system_generated",
    ".r-lib",
    "draftbychapter",  # Thesis drafting artifacts
}

# Standalone root Python scripts to exclude (one-off scratch tools)
EXCLUDE_ROOT_FILES = {
    "extract_pdf.py",
    "scratch_search_both.py",
    "scratch_search_keys.py",
}

# Analysis / benchmark / test root directories to scan
ANALYSIS_DIR_ROOTS = [
    "M1",
    "M2",
    "M3",
    "M3.1",
    "M4",
    "M5",
    "M6",
    "M7",
    "M8",
    "M9_ttd_graph_paper",
    "O3",
    "O4",
    "scripts",
    "tests",
]

# Curated standalone root Python utility scripts to include in Analysis bundle
ANALYSIS_STANDALONE_FILES = [
    "run_ghana_field_discovery.py",
    "bundle_project_files.py",
    "generate_m2_md_tables.py",
]

DELIMITER_LINE = "=" * 88
SUB_DELIMITER_LINE = "-" * 88


def is_excluded(path: Path) -> bool:
    """Return True if any component of path is in EXCLUDE_PARTS."""
    return any(part in EXCLUDE_PARTS for part in path.parts)


def read_and_normalise(path: Path) -> str:
    """Read file with utf-8 (fallback latin-1) and normalise line endings and trailing whitespace."""
    try:
        raw = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        raw = path.read_text(encoding="latin-1")
    return "\n".join(line.rstrip() for line in raw.splitlines())


def sha256_file(path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def sha256_text(text: str) -> str:
    """Compute SHA-256 hash of a text string."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def collect_core_files(root: Path) -> list[Path]:
    """Collect all live Python source files in the hydrosheaf/ package."""
    core_root = root / "hydrosheaf"
    if not core_root.exists():
        raise FileNotFoundError(f"Core directory not found: {core_root}")

    files = [
        p
        for p in core_root.rglob("*.py")
        if p.is_file() and not is_excluded(p)
    ]
    return sorted(files, key=lambda p: p.relative_to(root).as_posix().lower())


def collect_analysis_files(root: Path) -> list[Path]:
    """Collect all analysis, benchmark, test, and script Python files."""
    files: list[Path] = []

    for rel_dir in ANALYSIS_DIR_ROOTS:
        d = root / rel_dir
        if d.exists():
            for p in d.rglob("*.py"):
                if p.is_file() and not is_excluded(p):
                    files.append(p)

    for rel_file in ANALYSIS_STANDALONE_FILES:
        p = root / rel_file
        if p.exists() and p.is_file() and not is_excluded(p):
            files.append(p)

    # De-duplicate and sort alphabetically
    unique = sorted(set(files), key=lambda p: p.relative_to(root).as_posix().lower())
    return unique


def build_bundle_text(
    title: str,
    files: Sequence[Path],
    root: Path,
    utc_timestamp: str,
) -> tuple[str, list[dict[str, str | int]]]:
    """Assemble the complete bundle string adhering to the canonical Code_Only_Exports format."""
    header_lines = [
        title,
        "Source-only text export; no data, results, figures, tables, reports, or manuscript files are included.",
        f"Generated (UTC): {utc_timestamp}",
        f"Included source files: {len(files)}",
        "Each FILE block contains source text from the live checkout.",
        DELIMITER_LINE,
        "",
    ]
    header = "\n".join(header_lines)

    body_blocks: list[str] = []
    file_manifest: list[dict[str, str | int]] = []

    for p in files:
        rel_path = p.relative_to(root).as_posix()
        content = read_and_normalise(p)
        file_hash = sha256_text(content)

        block = f"FILE: {rel_path}\n{SUB_DELIMITER_LINE}\n{content}\n"
        body_blocks.append(block)

        file_manifest.append({
            "path": rel_path,
            "lines": len(content.splitlines()),
            "sha256": file_hash,
        })

    # Join blocks with delimiter
    separator = f"\n{DELIMITER_LINE}\n\n"
    bundle_text = header + separator.join(body_blocks) + "\n"
    return bundle_text, file_manifest


def export_bundle(
    title: str,
    output_path: Path,
    files: Sequence[Path],
    root: Path,
    utc_timestamp: str,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Generate and write a code bundle, returning summary metadata."""
    print(f"\nProcessing: {output_path.name}")
    print(f"  Target file count: {len(files)}")

    bundle_text, file_manifest = build_bundle_text(title, files, root, utc_timestamp)
    total_bytes = len(bundle_text.encode("utf-8"))
    total_lines = len(bundle_text.splitlines())
    bundle_hash = sha256_text(bundle_text)

    if not dry_run:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(bundle_text, encoding="utf-8")
        print(f"  [SAVED] {output_path}")
    else:
        print(f"  [DRY RUN] Would write to {output_path}")

    print(f"  Size: {total_bytes / (1024 * 1024):.2f} MB ({total_bytes:,} bytes)")
    print(f"  Total lines: {total_lines:,}")
    print(f"  SHA-256: {bundle_hash}")

    return {
        "file_name": output_path.name,
        "path": str(output_path),
        "source_file_count": len(files),
        "total_lines": total_lines,
        "size_bytes": total_bytes,
        "sha256": bundle_hash,
        "files": file_manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate HydroSheaf code-only text exports.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where text bundles will be saved (default: Code_Only_Exports).",
    )
    parser.add_argument(
        "--core-only",
        action="store_true",
        help="Generate only Hydrosheaf_Core_Code.txt.",
    )
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Generate only Hydrosheaf_Analysis_Code.txt.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Dry-run check: verify file lists and report counts without modifying disk.",
    )
    args = parser.parse_args()

    root = REPO_ROOT
    output_dir: Path = args.output_dir
    dry_run: bool = args.check

    utc_now = datetime.now(timezone.utc).isoformat()

    print("=" * 80)
    print("HydroSheaf Code-Only Export Generator")
    print(f"Root: {root}")
    print(f"Output Directory: {output_dir}")
    print(f"Timestamp (UTC): {utc_now}")
    print("=" * 80)

    run_core = not args.analysis_only
    run_analysis = not args.core_only

    manifest: dict[str, Any] = {
        "package": "HydroSheaf Code-Only Exports",
        "generated_utc": utc_now,
        "dry_run": dry_run,
        "bundles": {},
    }

    if run_core:
        core_files = collect_core_files(root)
        core_meta = export_bundle(
            title="HYDROSHEAF CORE REUSABLE IMPLEMENTATION",
            output_path=output_dir / "Hydrosheaf_Core_Code.txt",
            files=core_files,
            root=root,
            utc_timestamp=utc_now,
            dry_run=dry_run,
        )
        manifest["bundles"]["core"] = core_meta

    if run_analysis:
        analysis_files = collect_analysis_files(root)
        analysis_meta = export_bundle(
            title="HYDROSHEAF ANALYSIS AND VALIDATION CODE",
            output_path=output_dir / "Hydrosheaf_Analysis_Code.txt",
            files=analysis_files,
            root=root,
            utc_timestamp=utc_now,
            dry_run=dry_run,
        )
        manifest["bundles"]["analysis"] = analysis_meta

    if not dry_run:
        manifest_path = output_dir / "export_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        print(f"\n[OK] Wrote export manifest: {manifest_path}")

    print("\n[SUCCESS] Code-only export processing completed.")


if __name__ == "__main__":
    main()
