# M3 Hydrosheaf-USGS Investigation and Improvement Plan

> **Historical planning document — superseded 2026-07-28.** It predates the
> locked reportability, leakage-control, graph, and withdrawal decisions. Do
> not quote its status or anticipated outcomes as current evidence.

- Generated: 2026-05-15
- WISER source used throughout: `data/wiser_north_america.xlsx`

## Scope

This note answers four questions:

1. Why did Hydrosheaf not improve as much as expected against the USGS benchmark?
2. Which remaining issues are benchmark-harness flaws versus real model limitations?
3. What has already been fixed?
4. What should be done next if Hydrosheaf is to perform better and support stronger M3 claims?

## Files investigated

- `M3/m3_age_benchmark/results/m3_usgs_benchmark_results.csv`
- `M3/m3_age_benchmark/results/m3_usgs_benchmark_run_manifest.json`
- `M3/m3_age_benchmark/results/m3_phase4_screened_full_results.csv`
- `M3/m3_age_benchmark/results/m3_phase4_screened_full_results_manifest.json`
- `M3/m3_age_benchmark/results/m3_phase4_younggas_results.csv`
- `M3/m3_age_benchmark/results/m3_phase4_younggas_results_manifest.json`
- `M3/m3_age_benchmark/results/m3_real_usgs_graph_benchmark.csv`
- `M3/m3_age_benchmark/docs/m3_graph_benchmark_qa.md`
- `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py`
- `M3/m3_age_benchmark/scripts/run_m3_design_matrix.py`
- `M3/m3_age_benchmark/scripts/run_m3_real_usgs_graph_benchmark.py`

## Executive finding

Hydrosheaf did improve, but not in the way earlier summaries might suggest.

- The young-gas screening work is real and useful.
- The full benchmark still does not improve strongly because young-gas problems are only one part of the total error budget.
- A critical benchmark-harness flaw was also identified: recent design-matrix and screened runs were executed on a coarse age grid (`age_steps = 8` or `12`) even though the canonical M3 script default is `35`.

That coarse grid suppresses performance and should not be treated as the final canonical benchmark setting.

## Current benchmark state

### Full screened run currently on disk

`M3/m3_age_benchmark/results/m3_phase4_screened_full_results.csv`

Metrics from the refreshed canonical run:

- rows: `1272`
- finite multi-tracer estimates: `1272`
- median absolute log10 error: `0.382737`
- log10 RMSE: `1.064566`
- within factor 2: `0.409591`
- within factor 10: `0.677673`
- age-grid setting: `age_steps = 35`

This full screened run is now the canonical public-age benchmark, but the result is weaker than earlier provisional notes. It supports screening-level public age agreement only.

### Earlier parity run on disk

`M3/m3_age_benchmark/results/m3_usgs_benchmark_results.csv`

Metrics:

- median absolute log10 error: `0.191433`
- log10 RMSE: `0.867363`
- within factor 2: `0.596244`
- within factor 10: `0.849765`

This earlier parity run used `age_steps = 12`, which is also below the canonical default.

## Critical flaw identified and fixed

### Flaw: benchmark age-grid under-resolution

Evidence:

- `run_m3_usgs_benchmark.py` canonical default: `M3_DEFAULT_AGE_STEPS = 35`
- `design_matrix.yaml` had default `age_steps: 12`
- earlier full screened run manifests recorded `age_steps: 8`
- the earlier 80-row Phase 4 run manifest also recorded `age_steps: 8`

This meant earlier M3 design-matrix outputs were generated at lower resolution than the benchmark code itself defines as canonical. The current `m3_phase4_screened_full_results.csv` has been rerun at `age_steps = 35`.

### Why this matters

A coarse age grid can:

- lock the fit onto poor local age mixtures
- reduce the ability of the screened young-gas path to express better solutions
- make Hydrosheaf look weaker than it is

### Sensitivity check

A controlled rerun of the same screened scenario on the same 80 stratified rows showed:

- `age_steps = 8`:
  - median absolute log10 error: `0.115949`
  - log10 RMSE: `0.459719`
  - within factor 2: `0.7625`
  - within factor 10: `0.9375`

- `age_steps = 35`:
  - median absolute log10 error: `0.105540`
  - log10 RMSE: `0.475314`
  - within factor 2: `0.7875`
  - within factor 10: `0.9500`

Paired comparison on those same 80 rows:

- improved-fraction (`35` vs `8`): `0.4625`
- factor-2 gains/losses: `4 / 2`

Interpretation:

- the coarse grid was suppressing benchmark performance
- raising age-grid resolution does not fix everything
- but it is a real harness flaw and must be corrected before making strong benchmark claims

### Fix applied

The harness has been corrected:

- `M3/m3_age_benchmark/configs/design_matrix.yaml`
  - default age steps changed from `12` to `35`
- `M3/m3_age_benchmark/scripts/run_m3_design_matrix.py`
  - full runs now reject coarse grids below the canonical default unless explicitly overridden with `--allow-coarse-full-grid`
- tests added in `tests/test_m3_usgs_benchmark.py`

## Why Hydrosheaf still does not improve strongly after the gas-screening work

### 1. The young-gas fix only touches part of the benchmark

From the full screened run:

- kept `usgs_dgm`: `1014` rows
- switched to `raw`: `265` rows
- rows with any gas masking: `209`

Masked tracers:

- `SF6`: `191`
- `CFC113`: `14`
- `CFC12`: `6`
- `CFC11`: `5`

Interpretation:

- the screening logic is not affecting the entire benchmark
- it mainly repairs a subset of young-gas pathologies, especially `SF6`
- that is useful, but it cannot transform the whole benchmark by itself

### 2. Many rows are dominated by old-water tracers, not young gases

For old and very-old rows combined:

- `14C_only`: `473`
- `14C+4He`: `165`
- `4He_only`: `62`
- `none`: `21`

Interpretation:

- for a large fraction of the national benchmark, the main decision variables are `14C` and `4He`
- fixing `SF6/CFC` contamination alone cannot move those rows very much

### 3. Improvements are concentrated in younger water

Against the earlier parity file, the screened full run improved only a minority of rows in each age class:

- `modern_le_50`: improved fraction `0.1026`
- `intermediate_50_1000`: improved fraction `0.0945`
- `old_1000_30000`: improved fraction `0.0655`
- `very_old_gt_30000`: improved fraction `0.0448`

Interpretation:

- the recent Hydrosheaf work mostly helps modern and intermediate groundwater
- the remaining benchmark loss is increasingly an old-water problem, not only a gas-correction problem

### 4. The current graph prior is not physically strong enough

The real-USGS graph benchmark is now implemented, but it is a negative result.

Using the full screened pointwise output:

- no candidate graph family improved the pointwise baseline once nonzero regularization was applied
- wrong-direction and randomized controls degraded much more strongly, which is useful as a guardrail

Interpretation:

- the graph benchmark is valid as a falsifiable test
- but `StudyUnit + depth proxy` is not a strong enough hydrogeologic graph to improve the benchmark
- current graph regularization should not be used as an improvement claim

## What this means technically

The remaining issue is now a mix of:

1. benchmark-resolution problems that were suppressing performance
2. tracer-reliability and dissolved-gas-model limitations in Hydrosheaf
3. old-groundwater uncertainty/calibration limits
4. weak national-scale graph edges

In other words, the issue is no longer one simple bug.

## What should be done next

### Priority 1: canonical screened benchmark at proper resolution

This step has now been run at canonical resolution. Re-run it with:

```powershell
python -u M3\m3_age_benchmark\scripts\run_m3_design_matrix.py `
  --scenario screened_dgm_gases `
  --full `
  --output M3\m3_age_benchmark\results\m3_phase4_screened_full_results.csv
```

Then rerun the graph benchmark:

```powershell
python -u M3\m3_age_benchmark\scripts\run_m3_real_usgs_graph_benchmark.py `
  --pointwise-results M3\m3_age_benchmark\results\m3_phase4_screened_full_results.csv `
  --scenario screened_dgm_gases
```

### Priority 2: replace hard gas decisions with soft tracer reliability weighting

Current behavior:

- choose corrected vs raw
- mask obviously bad tracers

What is missing:

- down-weight suspicious tracers instead of treating them as fully trusted until a mask threshold is crossed

Recommended redesign:

- compute tracer-level reliability scores using:
  - historical-atmosphere plausibility
  - disagreement with `3H/3He` or `3H`
  - contamination/supersaturation flags
  - tracer-specific uncertainty inflation
- pass those reliability scores into the joint objective as uncertainty inflation or explicit weights

Expected benefit:

- smoother handling of borderline `SF6/CFC` rows
- fewer brittle young-gas failures
- better modern/intermediate benchmark performance

### Priority 3: add an internal dissolved-gas process model

Hydrosheaf currently depends too much on:

- USGS corrected gas values
- raw fallback
- heuristic screening

What it should estimate directly:

- excess air
- gas loss / degassing
- recharge temperature
- saturation departure
- contamination likelihood

Expected benefit:

- less dependence on external correction choices
- more defensible `SF6/CFC` benchmarking
- stronger scientific claim for Hydrosheaf as a model rather than only a benchmark wrapper

### Priority 4: strengthen the old-groundwater module

The next large gains are likely in old and very-old groundwater.

Needed improvements:

- regional or hydrogeologically grouped `4He` accumulation priors
- stronger probabilistic fusion of `14C` correction candidates
- explicit uncertainty propagation when `14C` and `4He` disagree
- age-class-aware model selection so old-water rows are not influenced too strongly by noisy young tracers that happen to be present

Expected benefit:

- better performance where most remaining national-benchmark difficulty now resides

### Priority 5: treat the graph benchmark conservatively until real edges exist

Do not try to force a positive graph result from `StudyUnit + depth proxy`.

Needed edge information:

- common aquifer / screen interval
- hydraulic head gradient
- known stratigraphic continuity
- basin-specific flow direction

Use current graph benchmark as:

- a guardrail
- a falsification test
- evidence that coarse national proxy edges are not enough

## Recommended claim discipline for M3

What can be claimed now:

- Hydrosheaf now has a defensible screened young-gas pathway with explicit contamination masking
- the benchmark harness now prevents accidental coarse-grid full runs
- the real-USGS graph benchmark exists and is falsifiable
- the full public-age result is screening-level and heterogeneous

What should not be claimed yet:

- that graph regularization improves USGS benchmarking
- that the current `1272`-row screened run is strong TracerLPM equivalence
- that dissolved-gas handling is fully solved

## Practical development sequence

1. use the refreshed full screened benchmark at `35` age steps as the current public-age evidence
2. keep the graph benchmark as neutral/negative evidence under current national proxy edges
3. add soft tracer reliability weighting for young gases
4. rerun the 80-row and full-screened comparisons after weighting changes
5. strengthen old-water probabilistic handling
6. only then revisit model-selection and graph-improvement claims

## Summary

Hydrosheaf did not fail for one single reason.

The investigation shows:

- one real harness flaw was suppressing performance: coarse age-grid settings
- that flaw has now been fixed in the benchmark tooling
- the remaining underperformance is mostly structural:
  - young-gas fixes affect only part of the benchmark
  - old-water rows remain the harder problem
  - the graph prior is too weak and coarse to help at national scale

So the route to better benchmark performance is clear:

- use the corrected benchmark harness
- rerun the canonical screened benchmark at proper resolution
- then improve tracer reliability weighting, dissolved-gas process treatment, and old-groundwater inference
