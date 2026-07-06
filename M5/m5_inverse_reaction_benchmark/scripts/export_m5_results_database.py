"""Export M5 CSV/JSON outputs into a single DuckDB results database."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Iterable

import duckdb
import pandas as pd


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
DATABASE_PATH = RESULTS_DIR / "m5_results.duckdb"


def table_name(path: Path) -> str:
    stem = path.stem.lower()
    name = re.sub(r"[^a-z0-9]+", "_", stem).strip("_")
    if not name:
        name = "table"
    if name[0].isdigit():
        name = f"t_{name}"
    return name


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def csv_files() -> list[Path]:
    return sorted(RESULTS_DIR.glob("*.csv")) + sorted(TABLES_DIR.glob("*.csv"))


def write_dataframe(
    connection: duckdb.DuckDBPyConnection,
    frame: pd.DataFrame,
    name: str,
) -> None:
    connection.register("_m5_frame", frame)
    connection.execute(f'CREATE OR REPLACE TABLE "{name}" AS SELECT * FROM _m5_frame')
    connection.unregister("_m5_frame")


def write_json_table(
    connection: duckdb.DuckDBPyConnection,
    path: Path,
    name: str,
) -> int:
    data = json.loads(path.read_text(encoding="utf-8"))
    frame = pd.json_normalize(data, sep="__")
    write_dataframe(connection, frame, name)
    return len(frame)


def import_csvs(connection: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    used_names: set[str] = set()
    for path in csv_files():
        if not path.exists() or path.stat().st_size == 0:
            continue
        base = table_name(path)
        name = base
        suffix = 2
        while name in used_names:
            name = f"{base}_{suffix}"
            suffix += 1
        used_names.add(name)
        frame = pd.read_csv(path)
        write_dataframe(connection, frame, name)
        rows.append(
            {
                "table_name": name,
                "source_file": str(path.relative_to(BENCHMARK_DIR)),
                "source_kind": path.parent.name,
                "row_count": len(frame),
                "column_count": len(frame.columns),
                "sha256": file_sha256(path),
                "modified_utc": pd.Timestamp.utcfromtimestamp(
                    path.stat().st_mtime
                ).isoformat(),
            }
        )
    return pd.DataFrame(rows)


def import_jsons(connection: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in sorted(RESULTS_DIR.glob("*.json")):
        name = table_name(path)
        row_count = write_json_table(connection, path, name)
        rows.append(
            {
                "table_name": name,
                "source_file": str(path.relative_to(BENCHMARK_DIR)),
                "source_kind": "results",
                "row_count": row_count,
                "column_count": len(
                    pd.json_normalize(json.loads(path.read_text(encoding="utf-8"))).columns
                ),
                "sha256": file_sha256(path),
                "modified_utc": pd.Timestamp.utcfromtimestamp(
                    path.stat().st_mtime
                ).isoformat(),
            }
        )
    return pd.DataFrame(rows)


def create_views(connection: duckdb.DuckDBPyConnection) -> None:
    existing = {
        row[0]
        for row in connection.execute("SHOW TABLES").fetchall()
    }
    if {"benchmark_fits", "mechanism_resolution_scores"}.issubset(existing):
        connection.execute(
            """
            CREATE OR REPLACE VIEW v_primary_full_panel AS
            SELECT *
            FROM benchmark_fits
            WHERE panel = 'full_11' AND noise_level = 0.03
            """
        )
    if {"data_tier_experiment", "data_tier_evidence_lifted_resolution"}.issubset(
        existing
    ):
        connection.execute(
            """
            CREATE OR REPLACE VIEW v_data_tier_summary AS
            SELECT
                data_tier,
                archetype,
                AVG(class_f1) AS mean_class_f1,
                AVG(false_discovery_rate) AS mean_false_discovery_rate,
                AVG(extent_rmse_mmolL) AS mean_extent_rmse_mmolL,
                COUNT(*) AS n
            FROM data_tier_experiment
            GROUP BY data_tier, archetype
            """
        )
    if "external_field_evidence_lifted_resolution" in existing:
        connection.execute(
            """
            CREATE OR REPLACE VIEW v_external_field_elri AS
            SELECT
                dataset,
                data_tier,
                class_id,
                members,
                AVG(evidence_lifted_resolution_index) AS mean_elri,
                MEDIAN(evidence_lifted_resolution_index) AS median_elri,
                AVG(CASE
                    WHEN resolution_status IN (
                        'conditionally_preferred',
                        'evidence_lifted_resolved'
                    )
                    THEN 1.0 ELSE 0.0 END) AS preferred_or_resolved_fraction,
                COUNT(*) AS n
            FROM external_field_evidence_lifted_resolution
            GROUP BY dataset, data_tier, class_id, members
            """
        )


def export_database(database_path: Path = DATABASE_PATH) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if database_path.exists():
        database_path.unlink()
    with duckdb.connect(str(database_path)) as connection:
        csv_catalog = import_csvs(connection)
        json_catalog = import_jsons(connection)
        catalog = pd.concat([csv_catalog, json_catalog], ignore_index=True)
        write_dataframe(connection, catalog, "m5_table_catalog")
        create_views(connection)
        manifest = {
            "database": str(database_path.relative_to(BENCHMARK_DIR)),
            "generated_utc": pd.Timestamp.utcnow().isoformat(),
            "n_tables": int(len(catalog)),
            "n_rows_total": int(catalog["row_count"].sum()) if not catalog.empty else 0,
            "tables": catalog["table_name"].tolist(),
        }
        write_dataframe(connection, pd.json_normalize(manifest), "m5_database_manifest")
    catalog.to_csv(RESULTS_DIR / "m5_results_database_catalog.csv", index=False)
    return database_path


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database",
        type=Path,
        default=DATABASE_PATH,
        help="DuckDB file to write.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    path = export_database(args.database)
    print(f"Exported M5 results database: {path}")


if __name__ == "__main__":
    main(sys.argv[1:])
