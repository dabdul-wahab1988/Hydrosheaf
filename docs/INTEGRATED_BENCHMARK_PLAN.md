# Integrated benchmark plan (points 1–3)

This repository now implements three complementary evidence paths. They are
deliberately reported as separate tiers rather than pooled into one apparent
accuracy number. The active execution scope is now data-in-hand analysis;
the detailed scope lock is in
[`docs/CURRENT_DATA_IN_HAND_SCOPE.md`](CURRENT_DATA_IN_HAND_SCOPE.md).

## 1. Published USGS model-reference panels

Run:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_usgs_integrated_reference.py
```

The run ingests the compact tables from the [USGS national groundwater-age
release](https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004)
and the already audited M4 projections of the [Savage MODPATH release](https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute)
and Great Miami MODFLOW/MODPATH release. It writes separate age, edge, and
travel-time tables plus source hashes. Supplying
`--aiken-root data/AikenCounty` additionally parses the local [Aiken integrated
MODFLOW/MODPATH release](https://www.usgs.gov/data/modflow-nwt-and-modpath5-used-evaluate-groundwater-availability-geochemistry-and-flow-pathways):
its compact chemistry/CFC workbooks, well metadata, MODFLOW/MODPATH text
metadata, and explicitly validated MODPATH5 binary endpoint/pathline panels
are written to separate run-scoped CSVs. The multi-gigabyte archives are read
through ZIP members and are not extracted or silently hashed during a normal
run.

These are calibrated/model-derived references used for implementation,
age--transport concordance, pathway extent, and model-reference diagnostics.
Aiken's CFC age intervals and Table 15 pathway summaries remain
model-conditioned; its active/stopped/truncated MODPATH records are explicitly
status-coded. The generated manifest therefore keeps integrated scoring
disabled and preserves the source-kind boundary for every row.

The supplied release was actually replayed in
`.codex_work/runs/RUN-USGS-INTEGRATED-REFERENCE-20260908-07`. The manifest
records 20 wells, 26 field-sample rows, 20 CFC-age rows, 16 recharge-to-well
summaries, 4,785 canonical chemistry records, 61 MODPATH5 endpoints, 5,322
pathline records, and 28 validated binary files. Fifty-three endpoints are
backward-to-recharge and eight are from the separately identified forward run;
35 particles remain active at the stop time and are right-censored. This is a
reproducible model-reference execution, not a new M2 field-accuracy result.

## 2. Independent integrated synthetic truth

Run a non-claim-bearing smoke benchmark with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_integrated_truth_blind.py --quick
```

The independent MODFLOW/MODPATH/chemistry generator creates observations and
sealed edge, age, and process truth. The wrapper audits the persisted blind
observation tables before scoring, records generator/executable hashes, and
never overwrites `m7_3_locked`. A full confirmatory run should use the locked
M7.3 settings and a new run ID; the quick run is only an execution and leakage
smoke test. The locked-settings wrapper is now available by omitting
`--quick`; it completed as
`.codex_work/runs/RUN-INTEGRATED-SYNTHETIC-20260908-05` (18 cases, truth-blind
audit PASS). Its manifest records the full settings and retains the
controlled-synthetic claim boundary.

## 3. Available Ghana data-in-hand panels

The four current packages are audited with:

```powershell
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\audit_prospective_campaign_readiness.py
```

The four current packages are analysed separately under the field-data
contract. They support source harmonisation, missingness/QC accounting,
geochemical-ratio diagnostics, mapped-geology sensitivity, and truth-free
transfer summaries. The package-specific readiness and abstention statuses
are retained in the run manifest; no unsupported topology or reaction label is
created from a missing field.

## Age directness: two-tier implementation

The missing direct-versus-indirect test is implemented as a separate,
run-scoped pair of panels rather than being retrofitted into the locked M7.3
F1/PR-AUC table. Run the full controlled synthetic panel with:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\age_bayes_factor_benchmark.py `
  --output .codex_work\runs\RUN-AGE-BF-CONTROLLED-YYYYMMDD-01
```

It generates independent, truth-sealed direct edges and two-hop reachable
skips; scores `no_age`, endpoint ordering, paired direct/indirect travel-time
Bayes factors, and a permuted-transport control; and writes the complete
truth denominator, case-block bootstrap intervals, and two-tier QA report.
The observed age increment is

\[
\Delta a=a_v-a_u,
\qquad
\sigma_\Delta^2=\sigma_u^2+\sigma_v^2-2\operatorname{Cov}(a_u,a_v)+\sigma_p^2,
\]

where (a_u,a_v) are endpoint ages (years), (sigma_u,sigma_v) are their
one-standard-deviation errors (years), (operatorname{Cov}(a_u,a_v)) is the
endpoint error covariance (years²), and (sigma_p) is process/mixing
uncertainty (years). The direct-versus-indirect score is

\[
\log BF_{D:I}=\log p(\Delta a\mid T_D)-\log p(\Delta a\mid T_I),
\]

with (T_D,T_I) the independently supplied direct and indirect travel-time
means and their RTD/model dispersions included in the two likelihoods. Missing
or overlapping hypotheses produce `ABSTAIN`; they are not converted to
negative edges.

The supplied Aiken release is run separately as a calibrated-model reference:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_aiken_model_conditioned_emulation.py `
  --aiken-root data\AikenCounty `
  --output .codex_work\runs\RUN-AIKEN-MODEL-CONDITIONED-YYYYMMDD-01
```

Aiken contributes declared MODFLOW/MODPATH travel-time distributions, CFC
apparent-age intervals, direction, and censoring diagnostics. Its rows remain
in a separate calibrated-model panel and are not pooled with the synthetic
directness metrics. The combined convenience runner is
`scripts/run_age_adjacency_two_tier.py`.
