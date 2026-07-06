"""Export a run-specific M5 figure, table, and plotting-code manifest."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re
import sys
from typing import Iterable

import duckdb
import pandas as pd


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = BENCHMARK_DIR / "docs"
TABLES_DIR = BENCHMARK_DIR / "tables"
FIGURES_DIR = BENCHMARK_DIR / "figures"
SCRIPTS_DIR = BENCHMARK_DIR / "scripts"
R_FIGURES_DIR = BENCHMARK_DIR / "r_figures"
DATABASE_PATH = BENCHMARK_DIR / "results" / "m5_results.duckdb"
MANIFEST_CSV = DOCS_DIR / "m5_artifact_manifest.csv"
MANIFEST_MD = DOCS_DIR / "m5_artifact_manifest.md"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_number(path: Path) -> str:
    stem = path.stem
    figure = re.match(r"figure(s?)(\d+)", stem, flags=re.IGNORECASE)
    if figure:
        prefix = "Figure S" if figure.group(1).lower() == "s" else "Figure "
        return f"{prefix}{figure.group(2)}"
    table = re.match(r"table(s?)(\d+)", stem, flags=re.IGNORECASE)
    if table:
        prefix = "Table S" if table.group(1).lower() == "s" else "Table "
        return f"{prefix}{table.group(2)}"
    return ""


def artifact_role(path: Path) -> str:
    rel = path.relative_to(BENCHMARK_DIR)
    parts = rel.parts
    if parts[0] == "tables":
        return "manuscript_table"
    if parts[0] == "docs":
        if path.name == "m5_artifact_manifest.csv":
            return "run_artifact_manifest"
        return "evidence_map"
    if parts[0] == "scripts":
        return "python_workflow_code"
    if parts[0] == "r_figures":
        return "r_plotting_code"
    if parts[0] == "figures":
        if "r_publication" in parts and "supplementary" in parts:
            return "r_supplementary_figure"
        if "r_publication" in parts:
            return "r_main_figure"
        if "publication" in parts and "supplementary" in parts:
            return "database_supplementary_figure"
        if "publication" in parts:
            return "database_main_figure"
        if "supplementary" in parts:
            return "legacy_supplementary_figure"
        return "legacy_main_figure"
    return "artifact"


def iter_artifacts() -> Iterable[Path]:
    patterns = [
        TABLES_DIR.glob("*.csv"),
        FIGURES_DIR.glob("*.png"),
        FIGURES_DIR.glob("*.pdf"),
        (FIGURES_DIR / "supplementary").glob("*.png"),
        (FIGURES_DIR / "supplementary").glob("*.pdf"),
        (FIGURES_DIR / "publication").glob("*.png"),
        (FIGURES_DIR / "publication").glob("*.pdf"),
        (FIGURES_DIR / "publication" / "supplementary").glob("*.png"),
        (FIGURES_DIR / "publication" / "supplementary").glob("*.pdf"),
        (FIGURES_DIR / "r_publication").glob("*.png"),
        (FIGURES_DIR / "r_publication").glob("*.tif"),
        (FIGURES_DIR / "r_publication").glob("*.pdf"),
        (FIGURES_DIR / "r_publication" / "supplementary").glob("*.png"),
        (FIGURES_DIR / "r_publication" / "supplementary").glob("*.tif"),
        (FIGURES_DIR / "r_publication" / "supplementary").glob("*.pdf"),
    ]
    for pattern in patterns:
        yield from pattern
    for path in [
        DOCS_DIR / "02_figures.md",
        DOCS_DIR / "03_tables.md",
        DOCS_DIR / "m5_results_summary.md",
        MANIFEST_CSV,
        SCRIPTS_DIR / "export_m5_results_database.py",
        SCRIPTS_DIR / "make_m5_database_figures.py",
        SCRIPTS_DIR / "export_m5_artifact_manifest.py",
        R_FIGURES_DIR / "README.md",
        R_FIGURES_DIR / "theme_m5.R",
        R_FIGURES_DIR / "plot_m5_publication_figures.R",
    ]:
        if path.exists():
            yield path


def collect_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    generated_utc = datetime.now(timezone.utc).isoformat()
    seen: set[Path] = set()
    for path in sorted(iter_artifacts()):
        resolved = path.resolve()
        if resolved in seen or not path.exists() or not path.is_file():
            continue
        seen.add(resolved)
        rel = path.relative_to(BENCHMARK_DIR).as_posix()
        stat = path.stat()
        rows.append(
            {
                "artifact_role": artifact_role(path),
                "display_number": display_number(path),
                "path": rel,
                "file_stem": path.stem,
                "extension": path.suffix.lower(),
                "size_bytes": stat.st_size,
                "modified_utc": datetime.fromtimestamp(
                    stat.st_mtime, timezone.utc
                ).isoformat(),
                "sha256": file_sha256(path),
                "generated_utc": generated_utc,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values(
            ["artifact_role", "display_number", "path"], kind="stable"
        ).reset_index(drop=True)
    return frame


def write_markdown(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "artifact_role",
        "display_number",
        "path",
        "extension",
        "size_bytes",
    ]
    lines = [
        "# M5 Run Artifact Manifest",
        "",
        "This manifest is generated from the current M5 output folders. It maps",
        "publication figures, supplementary figures, manuscript tables, evidence",
        "maps, and plotting code to actual files from the current run.",
        "",
        "| Role | Display | Path | Type | Bytes |",
        "|---|---:|---|---|---:|",
    ]
    for row in frame[columns].itertuples(index=False):
        lines.append(
            f"| {row.artifact_role} | {row.display_number} | "
            f"`{row.path}` | `{row.extension}` | {row.size_bytes} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def register_in_database(frame: pd.DataFrame, csv_path: Path) -> None:
    if not DATABASE_PATH.exists():
        return
    with duckdb.connect(str(DATABASE_PATH)) as connection:
        connection.register("_m5_artifact_manifest", frame)
        connection.execute(
            "CREATE OR REPLACE TABLE m5_artifact_manifest AS "
            "SELECT * FROM _m5_artifact_manifest"
        )
        connection.unregister("_m5_artifact_manifest")
        existing_tables = {
            row[0] for row in connection.execute("SHOW TABLES").fetchall()
        }
        if "m5_table_catalog" in existing_tables:
            catalog = connection.execute("SELECT * FROM m5_table_catalog").df()
            rel = csv_path.relative_to(BENCHMARK_DIR).as_posix()
            catalog = catalog[catalog["table_name"] != "m5_artifact_manifest"]
            catalog = pd.concat(
                [
                    catalog,
                    pd.DataFrame(
                        [
                            {
                                "table_name": "m5_artifact_manifest",
                                "source_file": rel,
                                "source_kind": "docs",
                                "row_count": len(frame),
                                "column_count": len(frame.columns),
                                "sha256": file_sha256(csv_path),
                                "modified_utc": datetime.fromtimestamp(
                                    csv_path.stat().st_mtime, timezone.utc
                                ).isoformat(),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            connection.register("_m5_table_catalog", catalog)
            connection.execute(
                "CREATE OR REPLACE TABLE m5_table_catalog AS "
                "SELECT * FROM _m5_table_catalog"
            )
            connection.unregister("_m5_table_catalog")
            manifest = pd.DataFrame(
                [
                    {
                        "database": str(
                            DATABASE_PATH.relative_to(BENCHMARK_DIR).as_posix()
                        ),
                        "generated_utc": datetime.now(timezone.utc).isoformat(),
                        "n_tables": int(len(catalog)),
                        "n_rows_total": int(catalog["row_count"].sum()),
                        "tables": catalog["table_name"].tolist(),
                    }
                ]
            )
            connection.register("_m5_database_manifest", manifest)
            connection.execute(
                "CREATE OR REPLACE TABLE m5_database_manifest AS "
                "SELECT * FROM _m5_database_manifest"
            )
            connection.unregister("_m5_database_manifest")


def export_manifest(
    csv_path: Path = MANIFEST_CSV,
    markdown_path: Path = MANIFEST_MD,
    update_database: bool = True,
) -> pd.DataFrame:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    frame = collect_manifest()
    frame.to_csv(csv_path, index=False)
    write_markdown(frame, markdown_path)
    if update_database:
        register_in_database(frame, csv_path)
    return frame


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=MANIFEST_CSV)
    parser.add_argument("--markdown", type=Path, default=MANIFEST_MD)
    parser.add_argument(
        "--no-database",
        action="store_true",
        help="Do not register the manifest in results/m5_results.duckdb.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    frame = export_manifest(args.csv, args.markdown, not args.no_database)
    print(f"Exported M5 artifact manifest with {len(frame)} files: {args.csv}")


if __name__ == "__main__":
    main(sys.argv[1:])
