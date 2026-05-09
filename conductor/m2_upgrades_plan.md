# Implementation Plan: M2 Thesis-Level Upgrades

## Objective
Integrate four major thesis-level capabilities into the Hydrosheaf codebase and the M2 manuscript outline. These upgrades bridge the gap from framework validation to advanced inference and discovery:
1.  **Robustness:** Advanced Uncertainty Quantification (Global Sensitivity Analysis).
2.  **Dynamics:** Temporal Sheaf Sections (Seasonal Wet/Dry dynamics).
3.  **Intelligence:** Thermodynamic Logic Gates (Automated Mineral Selection).
4.  **Visual Discovery:** 3D Process Mapping.

## Key Files & Context
-   **Outline:** `papers/M2/M2outline.md`
-   **Intelligence:** `hydrosheaf/models/reactions.py`, `hydrosheaf/phreeqc/constraints.py`, `hydrosheaf/sheaf/directed_section.py`
-   **Dynamics:** `hydrosheaf/temporal/temporal_edge_fit.py`, `hydrosheaf/sheaf/topology_refine.py`
-   **Robustness:** `hydrosheaf/uncertainty/sensitivity.py` (New)
-   **Visual Discovery:** `hydrosheaf/graph3d/export_3d.py` (New), `hydrosheaf/outputs/`
-   **Benchmarks:** `M2/m2_benchmark/scripts/run_e4c_northern_ghana_validation.py`

## Implementation Steps

### Phase 1: M2 Outline Synchronization
-   Modify `papers/M2/M2outline.md` to formally introduce the four upgrades into the manuscript structure.
-   **Sec 2.5.4 & 3.5.1:** Document the preemptive *Thermodynamic Logic Gates* pruning.
-   **Sec 2.6.1 & 3.1.3:** Highlight the *Global Sensitivity Analysis* bridging age variance to reaction stability.
-   **Sec 2.3.4 & 3.6:** Integrate *Temporal Sheaf Sections* explicitly into the Northern Ghana seasonal demonstration.
-   **Sec 2.6.4 & Proposed Figures:** Upgrade Figure 6 to reflect *3D Process Mapping*.

### Phase 2: Intelligence - Thermodynamic Logic Gates
-   **Current State:** The model uses PHREEQC post-optimization to check Saturation Indices (SI).
-   **Upgrade:** Modify `build_reaction_dictionary` in `hydrosheaf/models/reactions.py` and the constraints in `hydrosheaf/phreeqc/constraints.py` to accept pre-computed SI logic masks.
-   **Action:** If a water is heavily supersaturated (e.g., SI_calcite > 0.2), prune calcite dissolution from the candidate dictionary *before* the L1 sparse solver runs. This massively reduces the combinatorial search space and enhances thermodynamic rigor.

### Phase 3: Dynamics - Temporal Sheaf Sections
-   **Current State:** Temporal data is processed, but the core topological graph is largely static or loosely seasonal.
-   **Upgrade:** Enhance `hydrosheaf/temporal/temporal_edge_fit.py` and `hydrosheaf/sheaf/topology_refine.py` to support dynamic, time-aware edge filtering.
-   **Action:** Allow the Sheaf topology to "break" or "form" edges explicitly based on temporal metadata (e.g., separating wet vs. dry season connectivity dynamically) rather than just fitting reactions to static edges at different times.

### Phase 4: Robustness - Global Sensitivity Analysis
-   **Current State:** Bayesian/Bootstrap uncertainty exists for individual modules (like LPM ages).
-   **Upgrade:** Create `hydrosheaf/uncertainty/sensitivity.py`.
-   **Action:** Implement a global sensitivity wrapper that perturbs the upstream age/isotope inputs (e.g., ±10% tritium age) and quantifies how that variance propagates down into the L1 reaction residual and selected mineral phases. This explicitly links temporal uncertainty to geochemical non-uniqueness.

### Phase 5: Visual Discovery - 3D Process Mapping
-   **Current State:** Topologies and networks are represented as 2D spatial graphs.
-   **Upgrade:** Create `hydrosheaf/graph3d/export_3d.py`.
-   **Action:** Implement a VTK (for ParaView) or 3D JSON exporter that takes the solved directed section (nodes with 3D coordinates, edges with reaction labels and extents) and outputs it for 3D visualization.

## Verification & Testing
-   **Tests:** Run `pytest` to ensure no existing workflows break. Write targeted unit tests for `sensitivity.py`, `export_3d.py`, and the updated reaction dictionary pruning.
-   **Benchmark Validation:** Rerun the Northern Ghana benchmark (`run_e4c_northern_ghana_validation.py`) to verify that the Temporal Sheaf Sections and Thermodynamic Logic Gates actively improve (or maintain) the process fits under sparse seasonal data.

## Migration & Rollback
-   All new logic (Logic Gates, Temporal Sheafs, Sensitivity) will be introduced as configurable options (e.g., `use_logic_gates=True`) that default to the new behavior but can be toggled off to rollback to the strict M1/original M2 baseline if necessary.