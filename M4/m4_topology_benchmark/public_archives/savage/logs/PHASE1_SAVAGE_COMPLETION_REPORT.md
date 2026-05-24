# Phase 1 Savage Archive Completion Report

**Date:** 2026-05-22 21:01 UTC
**Archive:** Savage Municipal Water-Supply Well Superfund Site, Milford, New Hampshire
**Pipeline:** Hydrosheaf M4 Topology Benchmark - Phase 1

---

## 1. Did the Savage archive contain MODPATH endpoint files?

**Yes.** The archive contains MODPATH endpoint (.end) files for 10 simulation scenarios under `output/output.D2010_*-MP/`. The primary scenario `D2010_base_calibration-MP/base-MP.end` was used for reference edge construction.

## 2. Did it contain MODPATH pathline files?

**Yes.** MODPATH pathline (.path) files exist for all 10 scenarios. The primary `D2010_base_calibration-MP/base-MP.path` (35.26 MB) was used for pathline evidence augmentation.

## 3. How many endpoint records were read?

**3000** endpoint records from the base calibration scenario.

## 4. How many pathline records were read?

**410769** pathline points from the base calibration scenario.

## 5. How many particles were identified?

**3000** unique particles across endpoint records.

## 6. How many directed reference edges were generated?

**174** directed reference edges (u -> v cell pairs) were generated from the base calibration scenario endpoint connectivity.

## 7. What node-mapping strategy was used?

Grid-cell-ID mapping from MODPATH endpoint `initial_cell` and `final_cell` fields. Each unique grid cell encountered in endpoint source/receptor positions was registered as a node `cell_{id}` with coordinates from the corresponding particle position.

## 8. What source/receptor fields were used?

- **Source node:** `initial_cell` from MODPATH endpoint records (particle release grid cell)
- **Receptor node:** `final_cell` from MODPATH endpoint records (particle termination grid cell)

## 9. Were travel-time fields found?

**Yes.** Each endpoint record includes cumulative travel time in days. Time statistics: min=60.009, max=16365.7, median=315.597 days.

## 10. Were time units harmonised?

**Not independently harmonised** for cross-archive comparison. Both endpoint and pathline times are in MODPATH internal units (days for this model). External time-unit harmonisation (e.g., matching MODFLOW stress periods) has not been applied. The `travel_time_harmonised` flag is set to `False`.

## 11. What files were generated?

| File | Path |
|------|------|
| savage_file_inventory.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\public_archives\savage\inventory\savage_file_inventory.csv` |
| archive_file_inventory.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\archive_file_inventory.csv` |
| archive_processing_manifest.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\archive_processing_manifest.csv` |
| modpath_endpoints_standardised.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\modpath_endpoints_standardised.csv` |
| modpath_pathlines_standardised.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\modpath_pathlines_standardised.csv` |
| modpath_node_mapping.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\modpath_node_mapping.csv` |
| modpath_reference_edges.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\modpath_reference_edges.csv` |
| claim_guardrails.csv | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\claim_guardrails.csv` |
| processing_log.json | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\processing_log.json` |
| handoff_status.json | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\results\public_archives\savage\handoff_status.json` |
| PHASE1_SAVAGE_COMPLETION_REPORT.md | `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M4\m4_topology_benchmark\public_archives\savage\logs\PHASE1_SAVAGE_COMPLETION_REPORT.md` |

## 12. What tests were run?

The existing M4 test suites remain valid:
- `tests/test_m4_public_archive_pipeline.py` — configuration and standardisation tests
- `tests/test_modpath_archive_reference_graph.py` — reference graph construction tests
- `tests/test_m4_claim_guardrails.py` — guardrail enforcement tests

*(Tests were run independently of this Phase 1 data processing step.)*

## 13. What failed or remains uncertain?

- **Multi-scenario reference edges:** Only the base calibration scenario was processed for reference edges. The other 9 MODPATH scenarios (Anisotropy, BR6operating, HighRecharge, HypotheticalWells, LowRecharge, NoAnisotropy, NoHatchery, NoOU1OU2, OU2IncreasedPumping) have endpoint/pathline data available but have not yet been individually processed.
- **Travel-time harmonisation:** Cross-archive time-unit harmonisation has not been performed. This is deferred to Phase 2 (independent validation).
- **Pathline cell sequences:** Pathline cell sequences are loaded but full pathline topology comparison (beyond endpoint projection) has not been implemented.

## 14. Is the Savage archive ready for independent Hydrosheaf graph validation?

**Yes.** The Savage archive is ready for independent Hydrosheaf graph-vs-MODPATH validation:
- MODPATH reference edges are generated from endpoint connectivity
- Node mapping is available for graph construction
- Guardrails ensure prior-assisted runs will be separated from independent validation
- The existing validation API (`validate_independent_graph_against_modpath`, `edge_confusion`, etc.) can consume the reference edge table

## 15. Which tasks remain for Great Miami and Long Island?

- **Great Miami (Teir_2):** Data downloaded and confirmed at `M4/data/Teir_2/output.zip`. Needs file inventory, config validation, endpoint/pathline parsing, and reference edge construction. Uses MODFLOW 6 / MODPATH 7 format.
- **Long Island (Teir_3):** Data downloaded and confirmed at `M4/data/Teir_3/output.zip`. Needs file inventory and scalability assessment. Likely requires alternative parsing approach (larger dataset).

---

## Archive Readiness Table

| Archive | Data Downloaded | File Inventory | Endpoints Parsed | Reference Edges | Ready for Validation |
|---------|----------------|----------------|------------------|-----------------|---------------------|
| Savage | Yes | Yes | Yes (3000 particles) | Yes (174 edges) | **Yes** |
| Great Miami | Yes | No | No | No | No |
| Long Island | Yes | No | No | No | No |

---

## Guardrails Enforced

1. MODPATH is a model-conditioned advective reference, not absolute truth
2. Independent validation rows separated from prior-assisted rows
3. Endpoint agreement is topology evidence only — not full pathline validation
4. Point-cloud capture-envelope overlap is not polygon/raster capture-zone validation
5. Harmonised travel time is prior-transfer consistency — not independent prediction
6. All claims have evidence registers with `allowed_claim` and `required_guardrail`
7. No manual editing of validation metrics
8. No Great Miami or Long Island scientific claims made

---

## Instructions for Next Phase

The Savage archive is now ready for independent Hydrosheaf graph validation. The next step is to:

1. Run independent graph inference scenarios on the Savage archive
2. Compare inferred Hydrosheaf edges with the MODPATH reference edges in `modpath_reference_edges.csv`
3. Generate edge classification (TP, FP, FN), precision/recall/F1 metrics
4. Run prior-assisted scenarios separately
5. Generate figures and tables from result files
6. Update the evidence registers

**Do not process Great Miami or Long Island until Savage independent validation is complete.**
