# Field-data provenance and canonical data contract

**Contract:** `field-data-contract-v1`  
**Checkout audited:** 2026-09-07  
**Scope:** the four source packages under `data/FieldData/`.  The
`derived/` directory is a generated helper and is not a fifth source package.

The active data-in-hand execution scope is summarised in
[`CURRENT_DATA_IN_HAND_SCOPE.md`](CURRENT_DATA_IN_HAND_SCOPE.md).

This document fixes source identity, internal units, join status, and the
conditions under which a field record may enter an analysis.  The hashes below
are the observed SHA-256 values in this checkout; every analytical run must
recompute them and record them in its run manifest.

## Source packages

| Package | Source file(s) and role | Source rows/sheets | SHA-256 (observed) | Status |
|---|---|---:|---|---|
| `NorthenGhana` | `data/FieldData/NorthenGhana/NorthernGhana.xlsx`; canonical Northern Ghana seasonal chemistry | `Dry` 160 + `Wet` 160 = 320 samples | `55ca931f98ba2424d0711697b4b80b18b9a12b97c91f4738aa45f78b318aae0e` | **Canonical M6 input** |
| `Talensi_MiningArea` | `data/FieldData/Talensi_MiningArea/talensi.csv`; sparse mining-area transfer | 63 samples | `1e264b5600d46fc9154c6a2af6a2b7513f274b06b5c9a5074f2247c7d44dfce8` | **Canonical M6 input** |
| `LowerAnayari` | `data/FieldData/LowerAnayari/manu.csv`; sparse external transfer | 41 samples | `9362a641950c22aa6a2e4a7d0b5429da13c6cbe97acc8e03389dcc1d0baa7e05` | **Canonical M6 input** |
| `NorthernGhanaNew` | `data/FieldData/NorthernGhanaNew/compiled UER data_new.xlsx`; alternate compiled workbook (`GW`, `rain`, `monitoring wells`) | `GW`: 240 non-empty rows (237 numbered samples + 3 type-summary rows); `rain`: 61 records; `monitoring wells`: 107 records | `e3ab58394f8d3cc1e62874496cbe526e713e4fbff460957d16dbcd2cdda30f11` | **Available to the new four-package field contract; not merged into locked historical M6 outputs** |

`NorthernGhanaNew` must not be silently merged with `NorthenGhana`:
the workbooks have different schemas, coverage, missingness, and field
semantics.  The refactored `load_all_field_datasets` loader makes all four
packages available as separate, provenance-bearing panels; any combined fit
still requires a separately declared harmonisation and a new run ID.  Locked
historical M6 outputs remain unchanged.

For a new four-package run, load the sources separately, then flatten only
after recording the manifest:

```python
from hydrosheaf import Config, flatten_field_datasets, load_all_field_datasets

datasets = load_all_field_datasets(field_root="data/FieldData")
samples = flatten_field_datasets(datasets)
config = Config(sparse_panel_enabled=True, minimum_observed_ions=4)
```

The sparse-panel switch projects each candidate edge onto ions observed at
both endpoints.  It does not impute missing chemistry; edges with fewer than
the minimum panel are skipped and must be recorded as an abstention by the
caller.

The generated elevation helper is recorded separately for traceability:
`data/FieldData/derived/well_elevations_dem.csv`, 5,504 bytes,
SHA-256 `3400a1fe972a0e13166919a4672cf1a34f6c0f283c748cb1489a00b71ee04663`.
It is a derived DEM comparison, not an observed field-data source and must not
be treated as independent hydraulic-head truth.

### Source-unit observations

The canonical Northern Ghana workbook names most concentration fields with
`_mg_L`, `EC_uS_cm`, and isotope fields with `_per mil`.  The alternate
`NorthernGhanaNew` workbook includes an explicit units row: elevation `m`, EC
`Us/cm`, TDS and dissolved ions `mg/L`, stable isotopes `‰`, and tritium `TU`.
The CSV headers for Talensi and Lower Anayari do not encode units for every
column.  Their interpretation below is the loader contract, not a claim that
the original source metadata was complete.  If source documentation disagrees,
stop the run, record the discrepancy, and issue `ABSTAIN` until harmonisation
is reviewed.

## Canonical units and field mapping

All model-facing records use one explicit unit system:

| Quantity | Canonical unit | Source representations and rule |
|---|---|---|
| Dissolved ions (`Ca`, `Mg`, `Na`, `K`, `HCO3`, `Cl`, `SO4`, `NO3`, `F`, `Fe`, `PO4`) | `mmol/L` | Source dissolved concentrations are interpreted as `mg/L` where the loader declares them; divide by the ion molar mass exactly once. Never infer a missing ion as zero for a quantitative result. |
| Strontium and silica (`Sr_mgL`, `SiO2_mgL`) | `mmol/L` | Source values named `Sr_mg_L`/`SiO2_mg_L` are converted exactly once using their molar masses; the native field and source hash remain in provenance. |
| Charge-balance terms | `meq/L` | Derived from canonical `mmol/L` using absolute ionic charge; charge-balance error is reported in percent. |
| Reaction extents/residuals | `mmol/L` | Computed only after concentration harmonisation. |
| pH | dimensionless | Preserve numeric value; no concentration conversion. |
| Electrical conductivity | `µS/cm` | Source `EC`, `EC_uS_cm`, or `Us/cm`; preserve the source value and normalize the unit label. |
| TDS and salinity | `mg/L` for TDS; source unit for salinity | `Sal` in Talensi is retained as source metadata unless its unit is documented; it is not silently treated as TDS. |
| Temperature | `°C` | Source `Temperature_C`, `Temp`, or `Tempe`. |
| Latitude/longitude | decimal degrees, WGS 84 (`EPSG:4326`) | Northern Ghana `Latitude`/`Longitude`; Lower Anayari `Y coordinate`/`X coordinate`; Talensi latitude plus positive degrees-west longitude converted to signed west longitude by the loader; alternate workbook `Latituted`/`Longitude`. |
| Elevation, borehole depth, static water level, distances | `m` (horizontal distances may be source `km` where named) | Preserve source fields and names; do not compare a source `Distance_*_km` directly with a metre-valued coordinate distance. |
| Stable water isotopes | `‰` (delta notation) | `d18O`/`d2H`, `d18O_permil`/`d2H_permil`, or Greek-symbol columns. These are recharge/source or mixing evidence, not residence-time truth. |
| Tritium | `TU` | Available only in parts of `NorthernGhanaNew`; it is not present in the canonical M6 Northern Ghana workbook. A tritium value is not automatically a groundwater-age posterior. |
| Nitrate isotopes | `‰` | `d15N-NO3 (Air)` and `d18O-NO3 (VSMOW)` in `NorthernGhanaNew`; optional and separate from dissolved `NO3` concentration. |
| Geology geometry | source polygons in `EPSG:32630`; points transformed from `EPSG:4326` | Map metadata only; it is not screened-interval lithology or flow truth. |

The public API default is also `mmol/L` (see
[`docs/INPUTS_REFERENCE.md`](INPUTS_REFERENCE.md)).  A result manifest must
state both each source unit and the canonical unit after conversion; a bare
numeric column is not sufficient.

## Geology join status

The 2026-09-05 metadata join used the full `geology_polygon` layer (1,900
features, source CRS `EPSG:32630`) with points interpreted as WGS 84 decimal
degrees and a primary `within` predicate.  Invalid source geometries were
repaired in memory only.  The aggregate geology-sidecar package hash was
`d7db60451c0ce47137af947dd66ef2a9b5e4fb38affa65742c023739d7c7669e`.

| Source package | Join output | Status counts | Contract interpretation |
|---|---|---|---|
| `NorthenGhana/NorthernGhana.xlsx` | `outputs/2026-09-05_multi_geology_join/NorthernGhanaData_geology_join.csv` | 314 `MATCHED`, 6 `NO_MATCH` of 320 | Matched records may carry mapped-geology metadata; the six unmatched records retain missing geology and cannot be assigned a nearest polygon as if it were a match. |
| `Talensi_MiningArea/talensi.csv` | `.../Talensi_geology_join.csv` | 63 `MATCHED` of 63 | Map metadata available; no depth-resolved lithology or independent flow truth. |
| `LowerAnayari/manu.csv` | `.../LowerAnayari_geology_join.csv` | 41 `MATCHED` of 41 | Map metadata available; no depth-resolved lithology or independent flow truth. |
| `NorthernGhanaNew/compiled UER data_new.xlsx` | `outputs/2026-09-05_northern_ghana_geology_join/NorthernGhanaNew_geology_join.csv` | 236 `MATCHED`, 1 `NO_MATCH` of 237 groundwater samples | Use only this package's own audited sidecar and retain its coordinate-repair flags; never infer a geology assignment from the other Northern Ghana workbook. |

The source geology sidecars used by the join were `geology_polygon.shp`,
`.shx`, `.dbf`, `.prj`, `.sbn`, and `.sbx`; their individual hashes and all
join assumptions are preserved in the three per-dataset join manifests.  The
nearest polygon fields in those outputs are review references, not replacements
for a failed point-in-polygon match.  Mapped surface geology may be used as a
traceable descriptive covariate only after the coordinate and screened-interval
limitations are stated.

When a sidecar supplies `geology_boundary_distance_m`, the optional geology
prior applies the bounded spatial-quality factor
`c = 1 - exp(-distance / geology_boundary_sigma_m)` to the prior strength.
This down-weights points on polygon boundaries; it does not convert map
distance into screened-interval lithology or hydraulic evidence.

## Missingness, quality, and `ABSTAIN`

Missingness is typed, never silently imputed:

- `OBSERVED`: value is present, numeric, unit-declared, and passes the relevant
  range/QC checks;
- `MISSING`: no value was supplied;
- `CENSORED`: a detection-limit value such as `<0.01` was supplied and its
  censoring limit is retained;
- `INVALID`: value or unit cannot be parsed or fails a hard range check;
- `UNMATCHED`: a coordinate/geology join has no polygon match;
- `UNSUPPORTED`: the module requires a field that this package does not carry.

The pipeline may proceed with a lower evidence tier when an optional field is
missing, but it must record the disabled module and reason.  It must emit
`ABSTAIN` for an output that would otherwise imply unsupported evidence.  At a
minimum:

- no age mean/uncertainty, tracer, time axis, or transport distribution →
  `ABSTAIN` for age or direct-adjacency claims;
- no independent head/gradient evidence → do not call elevation a hydraulic
  head measurement; use only the declared elevation proxy, with a limitation
  flag;
- no required ion panel or unresolved units → `ABSTAIN` for the affected
  reaction fit;
- failed charge-balance/QC or invalid/censored values outside the declared
  handling rule → exclude from quantitative scoring and retain a reason code;
- `NO_MATCH`, `BOUNDARY_REVIEW`, or `OVERLAP_REVIEW` geology → do not assign a
  definitive mapped unit; use `UNMATCHED`/review status;
- no independent process, flow-path, or age truth → field outputs are
  screening/transfer diagnostics, not validation metrics.

Zero filling is permitted only where a protocol explicitly defines a measured
zero or a detection-limit model.  It is not a default missing-value policy.
All exclusions, lower-tier activations, and abstentions must appear in a
machine-readable readiness table with counts by package and variable.

## Run IDs and provenance requirements

Every analysis, harmonisation, geology join, or rerun receives a new immutable
run ID, for example:

`RUN-FIELD-DATA-CONTRACT-20260907-01`  
`RUN-M6-GEOLOGY-COVARIATE-20260907-01`

The run directory must be separate from source data and prior locked results.
Its manifest must contain:

1. run ID, UTC timestamp, protocol/contract identifier and commit/ref;
2. absolute source paths, file sizes, SHA-256 hashes, sheet names, row counts,
   and source-unit observations;
3. canonical field mapping, conversion version, missingness/QC rules, and
   package inclusion/exclusion decisions;
4. coordinate columns, source/target CRS, longitude-sign transformations, and
   geology layer/sidecar hashes where applicable;
5. software/runtime/dependency versions and hashes of generated outputs;
6. explicit status for every module (`RUN`, `DISABLED`, `ABSTAIN`, or
   `INVALID`) and a reason for every non-`RUN` status.

Do not overwrite source files, the three 2026-09-05 join outputs, or any locked
M6/M7 artifact.  A changed source hash, unit interpretation, geology join,
missingness rule, or canonical loader requires a new run ID and regenerated
derived outputs.  The field packages provide transfer/application evidence;
they do not supply independent truth for groundwater age, exact flow paths,
direct adjacency, or unique reaction mechanisms.

## Hash check

From the repository root, recompute the source hashes before execution:

```powershell
Get-ChildItem data\FieldData -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Select-Object Path, Hash
```

Record the resulting list rather than relying on this document's snapshot.
