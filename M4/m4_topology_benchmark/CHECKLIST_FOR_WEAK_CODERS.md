# Checklist for Weak Coders

If the pipeline fails, check these common issues:

## 1. File Location & Directory Structure
- [ ] Are the raw archive zip files located in `M4/data/Teir_1/`, `M4/data/Teir_2/`, or `M4/data/Teir_3/`? Note the directory typo "Teir" instead of "Tier".
- [ ] Is the output directory `M4/m4_topology_benchmark/results/` writable?
- [ ] Are the config YAML files in `M4/m4_topology_benchmark/configs/`?

## 2. Configuration & Ingestion Schema
- [ ] Does your config YAML file have all required fields?
  - `archive_name` (string)
  - `validation_tier` (string)
  - `modpath_version` (integer)
  - `local_archive_root` (string)
  - `zip_file` (string)
  - `endpoint_file_in_zip` (string)
  - `pathline_file_in_zip` (string)
- [ ] Is `archive_schema.yaml` present in `M4/m4_topology_benchmark/configs/`?

## 3. Execution Failures
- [ ] **FloPy parse errors**: Check if the `.endpoint` or `.pathline` file has headers. If using MODPATH 5/6 compact layouts, ensure the code redirects to `_load_compact_modpath5_endpoint_records` and `_load_compact_modpath5_pathline_points`.
- [ ] **Execution slowness**: Avoid `.iterrows()` loops on large dataframes. Always cast dataframes to records using `.to_dict("records")` before looping.
- [ ] **Long Island / Tier 3**: If running Tier 3 or files are missing, the pipeline expects stub mode. Check that `is_stub` triggers correctly and writes the empty data structures.

## 4. Scientific Rule Compliance
- [ ] **Independent Validation**: Check that candidate edges are inferred completely separate from prior-assisted overlays.
- [ ] **Allowed Claims**: Every output must include guardrail columns. Ensure that IOU is marked as a point-cloud proxy and travel-time is marked as Rank/Diagnostic only.
- [ ] **Output Order**: Do not re-order `external_modpath_archive_summary.csv`! Row 0 MUST be Savage, Row 1 Great Miami, Row 2 Long Island.
