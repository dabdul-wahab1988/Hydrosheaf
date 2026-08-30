# Metric Audit — CAGEO-D-26-00847 (canonical numbers)

All numbers below were regenerated from the current locked pipeline
(git HEAD `463e1ce`, `.venv` Python 3.14.6, `M2/m2_benchmark` and
`M3/m3_age_benchmark` runners). They are the single canonical set that the
revised manuscript text, tables, and figures must quote. Differences from the
submitted manuscript are flagged with the reviewer comment they answer.

## Synthetic benchmark (canonical, 100 realisations)

| Metric | Canonical value | Submitted claim | Comment |
|---|---|---|---|
| Active-reaction R2 (correlation, \|true\|>0.01) | **0.23** | 0.86 (text) / 0.74 (Fig. 3B) | R1-M10; text value had no surviving source; figure value from pre-fix code |
| Active-reaction MAE | 0.37 mmol/L | 0.23 | R1-M10 |
| Active-reaction RMSE | 0.62 mmol/L | 1.38 | R1-M10 |
| False-activation rate (>0.05 mmol/L) | 54.1% | 35.3% | R1-M10, R2-M3 |
| Median chemistry R2 (edge fit) | 1.00 | 1.00 | unchanged (Table 2) |
| Transport model detection | 91.7% (evap 75.5%, mix 99.7%) | not reported | R2-m1 |
| Transport median \|rel. bias\| | 15.8% | — | R2-m1 |
| Isotope pointwise R2 / MAE | -0.03 / 0.58 permil | 0.52 / 0.58 | R2-m1 |
| Isotope edge-mean R2 | 0.99 (noise sigma_delta = 0.71 permil) | 0.99 (0.71) | Unchanged. The noise amplitude is sqrt(2) x the 0.5 permil per-node isotope sigma (make_publication_figures.py:507) and Fig. 3D has always printed 0.71. An earlier version of this row recorded it as 0.07 in error and the manuscript inherited that typo; both corrected 2026-08-20. |
| PHREEQC proxy median RMSE / NSE | 0.02 / 1.00 | 0.05 / 1.00 | R2-m1 |
| PHREEQC feasibility fraction | 0.88 | — | R2-m1 |
| Age: pooled log10 corr-R2 / median AE | 0.98 / 17.86 y | 0.86 / 183.3 y (text) | R2-m1; Table 2 values already pooled |
| Age per class (Table 4) | young R2 0.97/MAE 0.81; intermediate 0.85/15.75; old 0.82/188.49; mixed NA/103.23; fossil NA/1069.81 | same | consistent; kept |
| Age-order consistency index | 0.85 | 0.85 | consistent; kept |

## Topology (canonical)

| Mode | Candidates | TP | FP | FN | P | R | F1 | Comment |
|---|---|---|---|---|---|---|---|---|
| No-prior (head-gradient k2) | 302 | 147 | 155 | 27 | 0.49 | 0.84 | **0.62** | R1-M11/R2-M2; replaces unsourced F1=0.86; matches M4 |
| Prior-assisted | 174 | 174 | 0 | 0 | 1.00 | 1.00 | 1.00 | ingestion fidelity only |

Baselines (results/revision/topology_baselines.csv): all-pairs elevation-drop
F1 = 0.47 (P 0.33, R 0.84 — same recall, 2x FPs vs k2); proximity kNN k2
F1 = 0.00 (306 edges, 0 TP); conservative-tracer ordering not evaluable on the
Savage archive (no hydrochemistry). MODPATH source: USGS Savage Municipal
Water-Supply Well MODFLOW-2005/MODPATH5 archive, DOI 10.5066/F7J102FK.

## Public age (canonical, M3 tracerlpm_parity_agefractions, identifiability-gated)

| Metric | Canonical | Submitted |
|---|---|---|
| Total rows / metric rows | 1272 / **356 identifiable** | 1249 (ungated) |
| Median \|log10 error\| | 0.022 | 0.166 |
| log10 RMSE | 0.235 | 0.799 |
| Within 2x / 10x | 91.9% / 98.6% | 63.0% / 85.3% |
| log10 R2 | 0.960 | 0.628 |
| Non-identifiable rows | 916 (flagged) | not reported |

## Data corrections found by coordinate check (2026-08-20)

**Talensi longitudes had the wrong sign.** Recorded as +0.618 to +0.815, which
places the site at ~0.7 deg E in Togo; Talensi District is at ~0.8 deg W.
Sampling SRTM at the recorded coordinates disagreed with the surveyed
elevations by -80.7 m; at sign-corrected coordinates the disagreement is
-5.5 m. Confirmed as a data-entry error by the corresponding author and fixed
in `data/FieldData/Talensi_MiningArea/talensi.csv` and the working copy
(backups `*.bak_lonsign`). **No reported result changes**: pairwise haversine
distance depends on longitude only through sin^2(dlon/2), which is even, so a
uniform sign flip mirrors the layout without altering any inter-well distance.
Re-running the field pipeline gives an identical edge set and a maximum
difference of 5.7e-08 across all numeric columns (solver noise in
`sheaf_obstruction_energy`). Any geolocated figure built before 2026-08-20
plotted Talensi in the wrong country and must be regenerated.

**PSI was not reproducible.** `analyze_sensitivity_mc` draws from the global
numpy RNG and neither PSI script seeded it, so the reported family distribution
drifted between identical runs (observed: Evaporites 108 -> 111 -> 107). Both
scripts now seed (20260820) and PSI output is byte-reproducible across runs.

## DEM elevation test for Lower Anayari (2026-08-20) -- rejected

Lower Anayari elevations are village constants (41 wells, 6 unique values), so
no within-village pair carries a direction. Per-well elevations were sampled
from SRTM 30 m and ASTER 30 m at the surveyed coordinates
(`data/FieldData/derived/well_elevations_dem.csv`), with Talensi as a control
because it has genuine per-well survey.

| Quantity | Value |
|---|---|
| SRTM error at Talensi (n=63) | bias -0.56 m, **RMSE 7.91 m**, MAE 5.16 m |
| ASTER error at Talensi | bias -6.48 m, RMSE 10.25 m |
| Implied pairwise dz noise (SRTM) | 11.18 m |
| dz needed for p_ij >= 0.75 at configured sigma (1.0 m) | 0.95 m |
| dz needed for p_ij >= 0.75 at measured sigma (7.91 m) | 7.54 m |
| Lower Anayari lateral pairs, DEM \|dz\| | median 6.5 m (IQR 3.0-10.0) |
| "Recovered" at configured sigma | **99 of 100** |
| Recovered at measured sigma | 40 of 100 |
| **\|dz\| exceeding 2x the DEM's own pairwise noise** | **0 of 100** |

Village-mean DEM agrees with the recorded village elevations (r = 0.773, median
\|offset\| 9.2 m), so the recorded values are approximately right, just coarse.
Adopting DEM elevations would have produced 99 confidently-directed edges from
differences that are, measurably, below the DEM's own noise floor. **Rejected**;
surveyed elevations retained and the within-village edges reported as
undirected.

## Field topology corrected (2026-08-20, third-round check)

The field results reported through the second revision were produced by a plain
2-D Euclidean kNN construction (`cKDTree`, k=3 including self) that applied **no
directional test at all**: no p_ij, no `edge_p_min`, no `edge_gradient_min`, no
`edge_radius_km`, and no sheaf refinement (`sheaf_refinement: not_requested`,
and the pipeline returned as many edges as it was given, 208 of 208). That
construction did not match Section 2.3, left 61 edges running uphill and 148 of
208 edges reciprocal (both A->B and B->A), and produced the 7 zero-closure sink
edges. The field tier is now run through the documented pipeline
(`infer_edges(method="probabilistic")` + `refine_edges_with_sheaf` with
`sheaf_cohomology_enabled=True`); the as-run construction is retained for
comparison in `M2/m2_benchmark/scripts/run_m2_field_documented_pipeline.py` and
the previous outputs in `M2/m2_benchmark/results/_as_run_baseline/`.

| Metric | as-run (superseded) | documented (canonical) |
|---|---|---|
| Candidates -> retained | 208 -> 208 (0 rejected) | 572 -> 258 (314 rejected) |
| Uphill edges | 61 | 0 |
| Reciprocal edges | 148 | 56 |
| edge_confidence present | no | yes, all 258 (0.50-1.00) |
| Sheaf refinement | not_requested | completed |
| Median chemistry R2 (Manu / Talensi / all) | 0.60 / 0.77 / 0.71 | 0.53 / 0.82 / 0.70 |
| Median PSI | 1.00 | 0.97 |
| PSI families | Evap 93, Cons 72, Redox 34, Carb 8, Plag 1 | Evap 107, Cons 84, Redox 55, Carb 7, Plag 5 (seeded) |
| Null (zero-closure) edges | 7 | 0 |
| Sheaf cohomology | not computed | Manu H0=0, H1=800, obstruction 0.295, 56 cycles; Talensi H0=10, H1=760, obstruction 0.128, 74 cycles |

## Field topology and canonical field outputs

The documented field construction is the canonical result used by the revised
manuscript: 572 candidates are refined to 258 retained edges, with 121 at
Lower Anayari (Manu) and 137 at Talensi. No retained edge runs against the
elevation proxy, and the canonical output contains zero unresolved null edges.
The previously submitted 208-edge construction is retained only as a
superseded as-run comparison in results/_as_run_baseline/; its metrics must
not be quoted as current field results.

| Metric | Canonical documented output |
|---|---|
| Candidates -> retained | 572 -> 258 (121 Lower Anayari, 137 Talensi) |
| Uphill edges | 0 |
| Reciprocal edges | 56 |
| Median chemistry R2 | 0.70 (Lower Anayari 0.53, Talensi 0.82) |
| Median PSI | 0.97 |
| PSI families | Evaporites 107, Conservative 84, Redox 55, Carbonates 7, Plagioclase 5 |
| Null edges | 0 |
| Sheaf diagnostics | Lower Anayari H0 = 0, H1 = 800, obstruction 0.295; Talensi H0 = 10, H1 = 760, obstruction 0.128 |

The H0 values above are homogeneous nullities. Because both affine
obstruction energies are positive, neither field network has an exact affine
global section at the reported numerical tolerance. The age-overlap audit also
uses explicit interval-overlap and severe flags; it does not apply an
overlap-proportional weight.


## New diagnostics (D3)

- Dictionary: 14 reactions x 11 ions, rank 8/11, condition number inf,
  10 collinear pairs |cos|>0.7 (calcite~anorthite 1.00, CaNa~NaCa 1.00, ...).
- Family-level recovery: R2 0.16, sign-match 36.5%, dominant-family hit 48%.
- Dictionary sensitivity: minus-albite improves active R2 0.09->0.17, MAE
  0.40->0.34; no single removal resolves degeneracy.
- PSI pairs: CaNa/NaCa separated (mean |dPSI| 0.85 Manu / 0.76 Talensi);
  calcite/dolomite both ~0 on field edges.
- Multi-endmember: halite 0.40 -> recovered 0.04; calcite 0.20 -> 0.00;
  chemistry R2 0.999 (mixing absorbs reaction signal).
- Runtime: candidate edges ~3n with kNN pruning; 0.06 s at 320 nodes.

## Consistency rules enforced

1. Text, Table 2, Table 4, Fig. 3B, Fig. 5 quote the canonical values above
   verbatim; all generators read the same result files.
2. Table 5 carries both topology modes; Fig. 2 identical.
3. Table 6 interpretations match the dominant-process column; null edges
   excluded and explained.
4. Section 5.3/5.5 claims reframed: extent-level recovery is limited by
   dictionary degeneracy (quantified); process-family stability (PSI) is the
   supported claim class.
