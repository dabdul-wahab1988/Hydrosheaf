"""Orchestrate the ingestion of MODPATH public archives (Savage, Great Miami, Long Island) into the M4 benchmark results."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_ROOT / "results"
CONFIG_DIR = BENCHMARK_ROOT / "configs"
SCRIPTS_DIR = BENCHMARK_ROOT / "scripts"

def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Sequential execution of the validation pipeline for the three archives
    configs = ["savage.yaml", "great_miami.yaml", "long_island.yaml"]
    pipeline_script = SCRIPTS_DIR / "run_public_archive_pipeline.py"

    for config in configs:
        config_path = CONFIG_DIR / config
        print(f"Running validation pipeline for {config}...", flush=True)
        subprocess.run(
            [sys.executable, str(pipeline_script), "--config", str(config_path)],
            cwd=REPO_ROOT,
            check=True
        )

    # 2. Concatenate validation summaries into a single archive summary CSV
    print("Concatenating archive summaries...", flush=True)
    summary_files = [
        "tier_1_savage_archive_summary.csv",
        "tier_2_great_miami_archive_summary.csv",
        "tier_3_long_island_archive_summary.csv"
    ]
    
    summaries = []
    for f in summary_files:
        path = RESULT_DIR / f
        if path.exists() and path.stat().st_size > 0:
            summaries.append(pd.read_csv(path))
        else:
            print(f"[WARNING] Summary file {f} was not generated or is empty.", flush=True)

    if summaries:
        combined_summary = pd.concat(summaries, ignore_index=True)
        combined_summary.to_csv(RESULT_DIR / "external_modpath_archive_summary.csv", index=False)
        print(f"Combined summary written to {RESULT_DIR / 'external_modpath_archive_summary.csv'}", flush=True)
    else:
        print("[ERROR] No summaries found to concatenate.", flush=True)

    # 3. Copy Savage detailed validation outputs to main names to avoid breaking downstream scripts
    print("Copying Savage detailed validation outputs to default filenames...", flush=True)
    mappings = {
        "tier_1_savage_edge_agreement.csv": "external_modpath_edge_agreement.csv",
        "tier_1_savage_pathline_particles.csv": "external_modpath_pathline_particles.csv",
        "tier_1_savage_pathline_structure.csv": "external_modpath_pathline_structure.csv",
        "tier_1_savage_capture_envelope_overlap.csv": "external_modpath_capture_envelope_overlap.csv",
        "tier_1_savage_capture_envelope_summary.csv": "external_modpath_capture_envelope_summary.csv",
        "tier_1_savage_travel_time_rank.csv": "external_modpath_travel_time_rank.csv",
        "tier_1_savage_travel_time_rank_summary.csv": "external_modpath_travel_time_rank_summary.csv",
        "tier_1_savage_harmonized_travel_time.csv": "external_modpath_harmonized_travel_time.csv",
        "tier_1_savage_harmonized_travel_time_summary.csv": "external_modpath_harmonized_travel_time_summary.csv",
        "tier_1_savage_bootstrap_ci.csv": "external_modpath_bootstrap_ci.csv"
    }

    for src_name, dst_name in mappings.items():
        src_path = RESULT_DIR / src_name
        dst_path = RESULT_DIR / dst_name
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"Copied {src_name} -> {dst_name}", flush=True)
        else:
            print(f"[WARNING] Savage file {src_name} is missing. Cannot copy.", flush=True)

    # 4. Generate the source manifest documenting the input archives and files
    print("Generating source manifest...", flush=True)
    manifest = {
        "source": "Public MODPATH archive ingestion diagnostics and guarded projection outputs",
        "archives": {
            "tier_1_savage": {
                "name": "USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH 5",
                "doi": "10.5066/F7J102FK",
                "url": "https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute"
            },
            "tier_2_great_miami": {
                "name": "Great Miami River Basin MODFLOW 6/MODPATH 7",
                "doi": "10.5066/P9X4C9R6",
                "url": "https://www.sciencebase.gov/catalog/item/60d0e65fd34e86b19a0a0ff0"
            },
            "tier_3_long_island": {
                "name": "Long Island Nitrogen MODFLOW 6/MODPATH 7",
                "doi": "10.5066/P97VFXZ4",
                "url": "https://www.sciencebase.gov/catalog/item/609da029d34e6a607548b87d"
            }
        },
        "m4_outputs": {
            "archive_summary": str(RESULT_DIR / "external_modpath_archive_summary.csv"),
            "edge_agreement": str(RESULT_DIR / "external_modpath_edge_agreement.csv"),
            "pathline_particles": str(RESULT_DIR / "external_modpath_pathline_particles.csv"),
            "pathline_structure": str(RESULT_DIR / "external_modpath_pathline_structure.csv"),
            "capture_envelope_overlap": str(RESULT_DIR / "external_modpath_capture_envelope_overlap.csv"),
            "capture_envelope_summary": str(RESULT_DIR / "external_modpath_capture_envelope_summary.csv"),
            "travel_time_rank": str(RESULT_DIR / "external_modpath_travel_time_rank.csv"),
            "travel_time_rank_summary": str(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv"),
            "harmonized_travel_time": str(RESULT_DIR / "external_modpath_harmonized_travel_time.csv"),
            "harmonized_travel_time_summary": str(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv"),
            "bootstrap_ci": str(RESULT_DIR / "external_modpath_bootstrap_ci.csv")
        }
    }

    manifest_path = RESULT_DIR / "external_modpath_source_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Manifest written to {manifest_path}", flush=True)
    print("Orchestration pipeline complete.", flush=True)

if __name__ == "__main__":
    main()
