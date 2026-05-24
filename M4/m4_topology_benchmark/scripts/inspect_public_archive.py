"""Scan local public MODPATH archive datasets and inspect their contents."""

from __future__ import annotations

import argparse
import os
import sys
import zipfile
from pathlib import Path

# Add project root to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.validation.modpath_archive import scan_modpath_archive


def inspect_archive(config_path: Path) -> bool:
    print("=" * 60)
    print(f"Inspecting archive config: {config_path.name}")
    print("=" * 60)

    try:
        config = scan_modpath_archive(str(config_path))
    except Exception as e:
        print(f"[ERROR] Config failed schema validation: {e}")
        return False

    archive_name = config.get("archive_name", "Unknown")
    tier = config.get("validation_tier", "unknown")
    modpath_ver = config.get("modpath_version", "unknown")
    local_root = PROJECT_ROOT / config.get("local_archive_root", "")
    zip_name = config.get("zip_file", "output.zip")
    zip_path = local_root / zip_name

    print(f"Archive Name:     {archive_name}")
    print(f"Validation Tier:  {tier}")
    print(f"MODPATH Version:  {modpath_ver}")
    print(f"Archive Root:     {local_root}")
    print(f"Zip File:         {zip_path.name}")

    if not local_root.exists():
        print(f"[ERROR] Local archive root directory does not exist: {local_root}")
        print("\nRecommended Next Actions:")
        print(f"  1. Create the directory: {local_root}")
        print(f"  2. Place public USGS archive zip files inside it.")
        print(f"Processing readiness: False")
        return False

    if not zip_path.exists():
        print(f"[ERROR] Zip file not found at: {zip_path}")
        print("\nRecommended Next Actions:")
        print(f"  1. Download the required archive data from source URL: {config.get('source_url')}")
        print(f"  2. Place the file {zip_name} in {local_root}")
        print(f"Processing readiness: False")
        return False

    print("[SUCCESS] Zip archive file found.")

    endpoint_file = config.get("endpoint_file_in_zip", "")
    pathline_file = config.get("pathline_file_in_zip", "")

    # Open zip and check file entries
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            namelist = zf.namelist()
            print(f"Total files in zip: {len(namelist)}")

            end_found = endpoint_file in namelist
            path_found = pathline_file in namelist

            print(f"Endpoint File in Zip: {endpoint_file} -> {'FOUND' if end_found else 'NOT FOUND'}")
            print(f"Pathline File in Zip: {pathline_file} -> {'FOUND' if path_found else 'NOT FOUND'}")

            if not end_found or not path_found:
                # If we're in Tier 3 (Long Island), it's known that files are missing or .mplst only.
                if tier == "tier_3_long_island":
                    print("[WARNING] Long Island archive does not contain explicit .endpoint7/.pathline7 files.")
                    print("Processing readiness: True (Fallback/Stub mode enabled)")
                    print("\nRecommended Next Actions:")
                    print("  - The pipeline will run in fallback-stub mode using synthetic test/metadata stubs.")
                    return True
                else:
                    print(f"[ERROR] Critical MODPATH outputs are missing from the zip archive.")
                    print("\nRecommended Next Actions:")
                    print(f"  1. Verify if the config contains correct file paths inside output.zip.")
                    print(f"  2. Re-extract or check original model output runs.")
                    print(f"Processing readiness: False")
                    return False

            # Inspect format by peeking
            with zf.open(endpoint_file, "r") as f:
                peek = f.read(500).decode("utf-8", errors="replace")
                print("\nEndpoint Header Peek:")
                print("-" * 40)
                print(peek.strip())
                print("-" * 40)
                if "MODPATH_ENDPOINT_FILE" in peek or "MODPATH7" in peek:
                    print("Detected MODPATH Version format: MODPATH 7")
                elif "compact" in peek or len(peek.splitlines()) > 1:
                    print("Detected MODPATH Version format: MODPATH 5/6/Compact")

    except Exception as e:
        print(f"[ERROR] Failed to open/inspect zip file: {e}")
        return False

    print("\nProcessing readiness: True")
    print("\nRecommended Next Actions:")
    print("  - Run the public archive pipeline script:")
    print(f"    python M4/m4_topology_benchmark/scripts/run_public_archive_pipeline.py --config {config_path}")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect public MODPATH archives.")
    parser.add_argument("--config", type=str, help="Path to archive YAML configuration file.")
    args = parser.parse_args()

    if args.config:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = PROJECT_ROOT / config_path
        inspect_archive(config_path)
    else:
        # Scan all three configurations by default
        config_dir = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "configs"
        configs = ["savage.yaml", "great_miami.yaml", "long_island.yaml"]
        for c in configs:
            inspect_archive(config_dir / c)
            print()


if __name__ == "__main__":
    main()
