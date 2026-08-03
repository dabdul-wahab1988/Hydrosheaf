# HydroSheaf code consolidation plan

**Status:** audit draft; no deletion authorised  
**Prepared:** 2026-08-01  
**Purpose:** reduce implementation duplication while preserving locked M2-M8 provenance.

## 1. Current diagnosis

The repository contains one main reusable package, `hydrosheaf/`, surrounded by
seven milestone benchmark packages, manuscript assets, historical bundles, and
provenance snapshots. A file inventory found approximately:

| Area | Python/structured files | Interpretation |
|---|---:|---|
| `hydrosheaf/` | 161 | Main reusable package, including several legacy and optional modules |
| `M2` | 27 | Benchmark and publication package |
| `M3` | 210 | Large age benchmark, diagnostics, manuscript and historical material |
| `M4` | 99 | Topology benchmark and publication material |
| `M5` | 68 | Reaction benchmark and manuscript package |
| `M6` | 68 | Field-transfer benchmark and manuscript package |
| `M7` | 126 | Non-uniqueness and integration benchmark package |
| `M8` | 100 | Calibration and active-learning benchmark package |

The exact duplicate-code audit found no broad duplication of the core scientific
implementation. The exact duplicate beyond empty package initialisers was the
M7 robust-hybrid runner and its content-addressed provenance archive, which is
intentional historical preservation. The apparent redundancy is therefore
mostly duplicated orchestration, overlapping APIs, legacy compatibility, and
multiple result/manuscript copies.

## 2. Target architecture

### Authoritative reusable code

`hydrosheaf/` remains the only authoritative implementation of reusable
inference algorithms:

- `hydrosheaf.api`: public connected-pipeline orchestration;
- `hydrosheaf.inference`: network and edge inference;
- `hydrosheaf.models`: reaction and process models;
- `hydrosheaf.nuclear`: tracer and age models;
- `hydrosheaf.graph` and `hydrosheaf.graph3d`: topology construction;
- `hydrosheaf.uncertainty`: uncertainty, calibration and selective-risk tools;
- `hydrosheaf.validation`: claim, metric, provenance and benchmark guardrails;
- a new programme-level validation module for discrepancy scenarios and common
  `ACTION` / `SET_REPORT` / `ABSTAIN` outputs.

### Thin benchmark packages

M2-M8 should retain their protocol files, input manifests, locked runners,
result files, claim decisions, tests, and manuscript assets. Their Python
scripts should become thin adapters that import the authoritative package and
call shared benchmark utilities. They should not redefine core reaction,
tracer, topology, uncertainty, or provenance algorithms.

### Shared benchmark utilities

Create one versioned benchmark utility layer for:

- run IDs and immutable run ledgers;
- source, input, environment and artifact hashes;
- common random-number handling;
- paired bootstrap contrasts;
- calibration, coverage and selective-risk metrics;
- truth-blind case writing and truth-release checks;
- stable CSV/JSON output; and
- claim-decision records.

Milestone-specific protocols may choose parameters, but the implementation of
these mechanics should not be copied seven times.

### Historical and publication material

`provenance/source_archive/`, old bundles, superseded scripts, rendered
figures, and manuscript snapshots must be clearly labelled as immutable
historical material. They must not be on the active import path. They are not
duplicates to delete casually because they explain how locked results were
produced.

## 3. Likely consolidation targets

### Reaction inference

`hydrosheaf.models.dictionary` is already a compatibility wrapper around
`hydrosheaf.models.reactions`; it should remain a small documented shim rather
than grow new logic. The M5 `m5_common.py` reaction dictionary and metrics look
similar to core functions, but M5 defines a benchmark-specific full reaction
reference and equivalence-class truth. It should not be merged blindly. The
safe solution is to share stoichiometric primitives and metric utilities while
keeping the M5 reference truth explicitly separate.

### Tracer and age inference

`joint_lpm.py`, `multi_tracer.py`, `old_groundwater.py`, and related modules have
different roles, but the public boundary is not yet obvious. The consolidation
task should designate one canonical age-result schema and label the other
functions as forward models, diagnostic utilities, or compatibility adapters.
No function should be removed until all M2/M3/M7 callers and tests use the
canonical schema.

### Topology and network inference

`hydrosheaf.api.fit_network_pipeline` should be the public orchestration entry
point. `inference.network_fit`, `graph.build`, and `graph3d.build_3d` should be
documented as lower-level 2-D, 3-D, and fitting components. Benchmark runners
should not reimplement pipeline sequencing or silently disable stages without
recording the stage status.

### Benchmark runners

Repeated `run_*`, bootstrap, hash, output-writing, and manuscript-export code
across M2-M8 should move behind shared utilities. A runner may define a
scientific protocol and its scenario-specific truth, but it should not own a
second copy of a general-purpose algorithm.

### Data and documentation copies

`M2_ready`, `M3_geochemistry`, old manuscript folders, code bundles, and
rendered exports should be classified as active, historical, or derived. Only
one copy should be authoritative for each current result. Other copies should
carry a manifest pointing to the authority and a checksum, rather than being
edited independently.

## 4. Safe migration sequence

1. Freeze the current M7/M8 dirty worktree boundary; do not mix consolidation
   with manuscript corrections.
2. Generate an import/dependency inventory and classify every candidate as
   `CANONICAL`, `ADAPTER`, `HISTORICAL`, `DERIVED`, or `UNUSED-CANDIDATE`.
3. Add shared interfaces and tests before moving implementation.
4. Move one function family at a time, keeping temporary deprecation shims.
5. Update all callers and benchmark manifests to the new source path.
6. Run core tests, milestone tests, smoke benchmarks, and source-hash checks.
7. Compare locked outputs before and after migration. Any numerical change
   requires a new scientific decision, not silent acceptance.
8. Archive superseded code with its hash and provenance reference. Delete only
   after the inventory, migration, and reproducibility checks pass.

## 5. Deletion rules

Do not delete or rewrite without an explicit archive record:

- locked benchmark runners;
- claim decisions and protocol locks;
- source archives used to explain a historical result;
- input files, result tables, and artifact manifests;
- manuscript sections that cite a specific historical run.

Potentially removable after verification are unreferenced caches, generated
bundles, stale exports, and duplicate presentation files. They still require a
read-only reference check before removal.

## 6. Completion criteria

The cleanup is complete only when:

- the package has one documented public pipeline and one result schema;
- benchmark scripts contain protocol orchestration rather than duplicated core
  algorithms;
- all active imports resolve to canonical package modules;
- historical files are excluded from active imports and clearly labelled;
- the full relevant test suite passes;
- representative locked results are unchanged or scientifically re-decided;
- a clean checkout reproduces the declared smoke and locked outputs; and
- the root README and release metadata describe the consolidated architecture
  without claiming unsupported production capability.

