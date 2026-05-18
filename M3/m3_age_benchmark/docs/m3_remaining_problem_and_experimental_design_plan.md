# M3 Remaining Problem Characterisation and Experimental Design Plan

Date: 2026-05-14

## Purpose

This note characterises the remaining M3 benchmark problem after the recent fixes to:

- calibrated 4He accumulation in the joint LPM path;
- better 14C correction-row selection;
- dissolved-gas-corrected SF6/CFC/3H/3He values from USGS tracer-output tables;
- strict use of `data/wiser_north_america.xlsx` for regional tritium input.

It also gives a detailed plan for turning M3 into a defensible experimental benchmark of Hydrosheaf against USGS groundwater-age estimates.

## Current State

The main M3 objective is to test whether Hydrosheaf can perform multi-tracer groundwater-age inference and whether graph regularization adds value beyond well-wise lumped-parameter models. The current codebase is closer to that goal, but the evidence package is not yet publication-ready.

Current verified facts:

- `load_usgs_national_dataset()` now loads 1,279 unique USGS site rows with no duplicated `site_id` values.
- 680 rows have finite USGS reference ages; 599 rows lack finite reference ages and should not be used for accuracy metrics.
- Age coverage among finite reference-age rows is:
  - modern, <=50 yr: 236 rows;
  - intermediate, 50-1,000 yr: 89 rows;
  - old, 1,000-30,000 yr: 249 rows;
  - very old, >30,000 yr: 106 rows.
- Corrected 14C values are now available for 679 rows, and corrected A0 values for 644 rows.
- USGS dissolved-gas-corrected tracer values are now available in the merged benchmark for many rows:
  - any DGM correction flag: 1,150 rows;
  - DGM tritium: 844 rows;
  - DGM 3H/3He: 383 rows;
  - DGM SF6: 570 rows;
  - DGM CFC-11: 55 rows;
  - DGM CFC-12: 54 rows;
  - DGM CFC-113: 47 rows.
- Calibrated 4He is now usable for 780 rows because the benchmark can combine terrigenic 4He with USGS LPM helium solution-rate information.
- The post-fix smoke run produced finite multi-tracer ages for 11 of 12 diagnostic rows and recovered several old 4He-only samples within a physically plausible range. Example: `CMORPAS1-01` has reference age 741,000 yr and post-fix estimate about 708,333 yr.
- One remaining smoke-row failure is not a general Hydrosheaf failure: the row is modeled as `3He(trit)` without usable tritium, and closed-system 3H/3He apparent age requires both tritium and tritiogenic helium.

The main result CSV at `M3/m3_age_benchmark/results/m3_usgs_benchmark_results.csv` is stale relative to the new fixes. It contains 1,584 rows from an older run where `est_age_multi` is entirely non-finite. It must be regenerated before any manuscript figure, table, or conclusion is trusted.

## Remaining Problem Characterisation

The remaining issue is not simply "the Hydrosheaf model is bad." It is a combination of model physics, benchmark design, and evidence-generation problems.

### 1. Post-fix full benchmark has not yet been frozen

The fixed code has been validated on targeted tests and a 12-row smoke run, but the full USGS benchmark output has not been regenerated and reviewed. The existing result file still reflects the pre-fix failure state. This means all downstream tables and figures that read the stale CSV are invalid until regenerated.

Required decision: treat the smoke run as diagnostic only, not as manuscript evidence.

### 2. Old-groundwater handling is improved but still not complete

The joint LPM path now supports calibrated 4He using sample-specific accumulation rates and background values. This fixes the previous failure where old 4He-rich samples could not be aged realistically.

Remaining gaps:

- 4He accumulation is still represented as a linear rate model. That is acceptable for a benchmark approximation, but not enough for a strong process claim unless uncertainty in rate, lithology, U/Th production, porosity, and crustal flux is propagated.
- 14C selection now prefers geochemically corrected rows, but it is still a ranking heuristic rather than a full geochemical decision tree.
- 14C and 4He are not yet combined in a hierarchical old-water module that treats 14C correction and He accumulation rate as uncertain latent variables.
- Rows with only tritiogenic 3He and no tritium need an explicit "not identifiable" diagnostic rather than a quiet missing estimate.

### 3. Dissolved-gas correction is now used but not independently audited

The benchmark now uses USGS Table 6 corrected tracer outputs where available. That is the right immediate move for fair USGS parity.

Remaining gaps:

- Hydrosheaf does not yet independently reconstruct dissolved-gas corrections from Ne, Ar, N2, temperature, salinity, elevation, and excess-air model choices for all rows.
- Contamination flags for SF6/CFCs are not yet fully integrated into tracer weights or exclusion rules.
- The benchmark needs paired comparisons of raw versus corrected gas tracers to show that dissolved-gas correction improves age inference rather than merely changing inputs.

### 4. Graph regularization evidence is still mostly synthetic

The current `run_m3_age_benchmark.py` graph benchmark is deliberately synthetic and controlled. It is useful for proving the mechanism and negative-control logic, but it cannot support a broad claim that graph regularization improves USGS benchmark performance.

Remaining gaps:

- There is not yet a full real-USGS graph benchmark using aquifer, depth, screen interval, coordinates, and possible flow-direction constraints.
- The current graph results show improvement in only 3 of 9 scenario-graph rows and degradation or no improvement in 6 of 9 rows. That is valuable, but it supports a conditional claim only.
- Randomized graph negative controls exist in the synthetic benchmark, but not yet in the real USGS graph benchmark.

### 5. Publication-ready artefacts are not yet evidence-driven

Several manuscript-ready tables and figures are currently mock, synthetic, or hard-coded:

- Table 3 AICc model-selection values are hard-coded.
- Table 4 regional gradients are hard-coded.
- Supplementary Bayesian diagnostics are hard-coded.
- Figure 2 includes mock tracer convolution.
- Figure 3 uses random synthetic posterior shrinkage.
- Figure 4 creates synthetic distance and single-node ages.
- Figure 7 creates random modern-fraction discordance examples.

These can remain only as conceptual placeholders. They cannot be used as manuscript evidence unless clearly labelled as schematic. The M3 objectives require output tables and figures to be generated from actual benchmark result files.

### 6. Metrics are not yet aligned with the age-scale problem

Raw RMSE in years is dominated by old and very old samples. It is not enough for a tracer-age benchmark spanning 1 yr to more than 700,000 yr.

Required metrics:

- log10 RMSE;
- median absolute log10 error;
- within-factor-2 and within-factor-10 accuracy;
- bias by age class;
- modern/mixed/premodern classification skill;
- uncertainty coverage where intervals are produced;
- performance stratified by tracer set, age class, aquifer, and correction availability.

### 7. M3 scope includes Ar-39 and Kr-85, but USGS support may be limited

The outline names Ar-39 and Kr-85 as part of the M3 tracer library. The current USGS national benchmark is dominated by 3H, 3H/3He, SF6, CFCs, 14C, and 4He. If Ar-39 and Kr-85 are sparse or absent in the USGS tables, M3 must not imply that the USGS benchmark validates those tracers at scale.

Resolution:

- either add a second dataset with Ar-39/Kr-85 coverage;
- or keep Ar-39/Kr-85 as synthetic or method-capability demonstrations, separate from USGS empirical validation.

## Critical Flaws to Fix Before Manuscript Claims

1. Regenerate the full post-fix `m3_usgs_benchmark_results.csv`.

Acceptance criterion: finite `est_age_multi` for a substantial share of rows with finite reference ages, zero duplicated site rows, and explicit failure reasons for non-estimable rows.

2. Replace stale and mock manuscript outputs with evidence-driven outputs.

Acceptance criterion: every manuscript table/figure has a provenance chain to a result CSV, config file, and script. Placeholder schematic figures must be labelled as schematic.

3. Split parity benchmarking from capability benchmarking.

Acceptance criterion: two result modes are reported separately:

- USGS-parity mode: Hydrosheaf fits the USGS-reported LPM family and modeled tracer set.
- Hydrosheaf-selection mode: Hydrosheaf selects among supported TTDs using AICc/BIC or cross-validation.

4. Add explicit non-identifiability diagnostics.

Acceptance criterion: rows such as `3He(trit)` without tritium are reported as "not identifiable from available tracer fields" rather than hidden as generic fit failure.

5. Add full tracer-correction ablation experiments.

Acceptance criterion: results quantify the effect of raw gases versus DGM-corrected gases, raw 14C versus corrected 14C, and default 4He versus calibrated 4He.

6. Build a real-USGS graph benchmark.

Acceptance criterion: graph claims are based on physically motivated USGS edges, randomized negative controls, prior-strength sweeps, and held-out validation by aquifer or region.

7. Add uncertainty calibration.

Acceptance criterion: benchmark outputs include intervals or Monte Carlo ensembles and report coverage, interval width, and whether uncertainty is under- or over-confident.

8. Strengthen old-groundwater physics.

Acceptance criterion: old samples report which old-water constraints were active: corrected 14C, 14C A0 model, calibrated 4He rate, 4He source, and whether old-water age is constrained by one or multiple tracers.

## Experimental Design: Can It Enhance the USGS Benchmark?

Yes. Experimental design is the correct way to make M3 stronger. The main benefit is that it converts the benchmark from a single parity plot into a controlled set of tests that can identify what Hydrosheaf does well, when it fails, and why.

The benchmark should be treated as a factorial experiment with real USGS rows as observational units and controlled modelling choices as experimental factors.

## Proposed Hypotheses

H1. Multi-tracer Hydrosheaf inference outperforms tritium-only inference when at least one old-water tracer or one dissolved-gas-corrected young tracer is available.

H2. Using USGS dissolved-gas-corrected SF6/CFC/3H/3He values improves modern-age performance relative to raw dissolved or uncorrected gas values.

H3. Geochemically corrected 14C improves old-groundwater age agreement relative to raw 14C or a fixed A0 assumption.

H4. Calibrated 4He improves very-old groundwater estimates, especially when 14C is low, absent, or weakly diagnostic.

H5. Graph regularization improves performance only under physically plausible graph topology and appropriate prior strength; randomized or wrong-direction graphs should not improve performance.

H6. Tracer-omission experiments can identify a minimum defensible tracer set for modern, intermediate, old, and very-old groundwater classes.

## Experimental Factors

### Data and correction factors

- Gas correction:
  - raw measured gases;
  - USGS DGM-corrected values;
  - DGM-corrected values with contamination screening.
- 14C correction:
  - raw measured 14C;
  - fixed A0;
  - USGS corrected 14C;
  - selected best correction model;
  - correction-model ensemble.
- 4He handling:
  - disabled;
  - default accumulation rate;
  - USGS calibrated accumulation rate;
  - calibrated rate with uncertainty.

### Model factors

- LPM strategy:
  - USGS-reported LPM family;
  - Hydrosheaf AICc/BIC-selected LPM;
  - fixed baseline families such as PFM, EM, DM, GA.
- Tracer set:
  - 3H only;
  - 3H plus 3H/3He;
  - young gases only;
  - 14C only;
  - 14C plus 4He;
  - full modeled tracer set;
  - leave-one-tracer-out subsets.
- Graph prior:
  - none;
  - weak;
  - medium;
  - strong;
  - randomized graph;
  - wrong-direction graph;
  - aquifer-only graph;
  - aquifer plus depth graph;
  - aquifer plus depth plus hydraulic plausibility graph.

### Stratification factors

- Age class:
  - modern, <=50 yr;
  - intermediate, 50-1,000 yr;
  - old, 1,000-30,000 yr;
  - very old, >30,000 yr.
- Tracer availability class.
- Aquifer or study unit.
- Presence or absence of corrected gases.
- Presence or absence of corrected 14C.
- Presence or absence of calibrated 4He.

## Controls

The benchmark should include controls that can falsify over-claiming:

- Tritium-only baseline.
- Single-well no-graph baseline.
- USGS-reported LPM parity baseline.
- Randomized graph negative control.
- Wrong-direction graph negative control.
- Permuted reference-age negative control.
- Tracer-removal controls.
- Correction-removal controls, such as raw gas instead of DGM-corrected gas.

## Response Metrics

Primary metrics:

- log10 RMSE against USGS reference age;
- median absolute log10 error;
- percentage within factor 2;
- percentage within factor 10;
- Spearman rank correlation;
- age-class classification accuracy.

Secondary metrics:

- raw MAE in years by age class;
- bias by age class;
- number of finite estimates;
- number of non-identifiable rows;
- fit objective and AICc/BIC;
- tracer residuals by tracer;
- graph-age monotonicity violations;
- graph improvement relative to randomized controls.

Uncertainty metrics:

- 50%, 80%, and 95% interval coverage;
- median interval width;
- negative log predictive density if likelihoods are available;
- calibration by age class.

## Statistical Analysis

Use paired comparisons wherever possible because each USGS row can be run under multiple modelling configurations.

Recommended analysis:

- paired bootstrap confidence intervals for metric differences;
- Wilcoxon signed-rank tests on absolute log10 error where distributions are non-normal;
- mixed-effects models with sample nested within aquifer or study unit;
- interaction terms for correction type, tracer set, graph prior, and age class;
- multiple-comparison control when many tracer-removal tests are reported;
- effect-size reporting, not only p-values.

## Detailed Implementation Plan

### Phase 1: Freeze the post-fix USGS parity benchmark

Tasks:

- Run `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py` on the full loaded USGS national dataset.
- Save a run manifest containing date, git status, WISER source path, age-step setting, model strategy, and code version.
- Produce a QA table with row counts, finite estimates, non-identifiable rows, failed rows, and duplicated IDs.
- Compare full-run summaries against the 12-row smoke run.

Deliverables:

- refreshed `results/m3_usgs_benchmark_results.csv`;
- `results/m3_usgs_benchmark_run_manifest.json`;
- `docs/m3_usgs_benchmark_qa.md`.

Acceptance criteria:

- zero duplicated `site_id` rows;
- no silent all-NaN estimate columns;
- explicit failure reason for each non-finite `est_age_multi`;
- metrics reported only for rows with finite reference ages.

### Phase 2: Add a benchmark design matrix

Tasks:

- Create a configuration-driven design matrix for model runs.
- Include factors for gas correction, 14C correction, 4He mode, tracer set, LPM strategy, and graph-prior strategy.
- Add row filters for age class, aquifer, and tracer availability.
- Ensure each run writes one tidy result table with `scenario_id`, `factor_*` columns, inputs used, outputs, diagnostics, and metrics.

Deliverables:

- `M3/m3_age_benchmark/configs/design_matrix.yaml`;
- `M3/m3_age_benchmark/scripts/run_m3_design_matrix.py`;
- `results/m3_design_matrix_results.csv`.

Acceptance criteria:

- each scenario is reproducible from a config row;
- no scenario overwrites another scenario's outputs;
- scenario IDs are stable and manuscript-readable.

### Phase 3: Harden old-groundwater module

Tasks:

- Promote calibrated 4He handling into an explicit old-groundwater helper module.
- Add uncertainty for 4He accumulation rates using Table 4 fields where available.
- Add 14C correction strategy labels and selection diagnostics.
- Add an ensemble mode for 14C correction where multiple plausible correction models exist.
- Add combined 14C plus 4He diagnostics that report agreement, conflict, or single-tracer constraint.

Deliverables:

- `hydrosheaf/nuclear/old_groundwater.py`;
- tests for calibrated 4He, corrected 14C selection, 14C/4He conflict, and non-identifiable cases;
- benchmark output fields for old-water diagnostics.

Acceptance criteria:

- old-water rows identify active constraints;
- 4He-only, 14C-only, and 14C+4He cases are distinguishable in results;
- very-old rows have finite estimates when the relevant tracer information is present.

### Phase 4: Complete dissolved-gas correction experiments

Tasks:

- Keep USGS DGM-corrected values as the parity benchmark input.
- Add raw-versus-corrected gas ablations.
- Add contamination/excess-air flags into tracer weighting or masking.
- Where source dissolved noble-gas fields are sufficient, independently reconstruct a Hydrosheaf DGM correction and compare it against USGS Table 6 values.

Deliverables:

- gas-correction ablation scenarios in the design matrix;
- DGM audit table comparing raw, USGS-corrected, and Hydrosheaf-corrected values;
- tests for DGM value precedence and masking.

Acceptance criteria:

- the effect of dissolved-gas correction is quantified, not assumed;
- contaminated or physically implausible gases are down-weighted or excluded with explicit flags.

### Phase 5: Build the real-USGS graph benchmark

Tasks:

- Construct candidate graph edges from aquifer, study unit, location, screen depth, and available hydraulic/head information.
- Generate graph families: coordinate-only, aquifer-constrained, depth-constrained, and hydraulic-plausibility-constrained.
- Run no-graph and graph-prior scenarios with weak, medium, and strong priors.
- Add randomized and wrong-direction negative controls.
- Evaluate graph regularization on held-out aquifers or study units.

Deliverables:

- `results/m3_real_usgs_graph_edges.csv`;
- `results/m3_real_usgs_graph_benchmark.csv`;
- `docs/m3_graph_benchmark_qa.md`.

Acceptance criteria:

- graph improvement is compared against randomized graph controls;
- strong priors that degrade accuracy are reported, not hidden;
- claims are conditional on graph quality and prior strength.

### Phase 6: Regenerate manuscript tables and figures from evidence

Tasks:

- Replace hard-coded tables with outputs derived from result CSVs.
- Replace mock figures with actual benchmark visualisations.
- Keep conceptual/schematic figures only if labelled as schematic.
- Add figure/table manifest entries stating source data and script for each artefact.

Deliverables:

- refreshed `tables/Manuscript_Ready/*.csv`;
- refreshed `figures/Manuscript_Ready/*.png`;
- updated `MANIFEST.md`;
- updated `docs/m3_results_summary.md`.

Acceptance criteria:

- every numeric claim in the manuscript-ready artefacts can be traced to a result file;
- no random mock data are presented as empirical benchmark results.

### Phase 7: Final M3 evidence gate

Tasks:

- Run targeted tests for joint LPM, USGS loader, dissolved gas, old groundwater, and graph diagnostics.
- Run full M3 design matrix.
- Generate a claim-audit report that maps each manuscript claim to the supporting result.
- Mark unsupported objectives as limitations rather than findings.

Deliverables:

- `docs/m3_claim_audit.md`;
- final benchmark tables and figures;
- test report.

Acceptance criteria:

- all M3 claims are supported by named tables, figures, or diagnostics;
- limitations are explicit for Ar-39/Kr-85, public-data transferability, and graph topology uncertainty;
- no result depends on the full/global WISER workbook.

## Recommended Immediate Next Steps

1. Regenerate the full post-fix USGS benchmark.

2. Add failure-reason fields to the benchmark output before the full run is treated as final.

3. Create the design matrix runner with a small subset first, then scale to the full dataset.

4. Replace manuscript-ready tables and figures that are currently mock or hard-coded.

5. Build the real-USGS graph benchmark only after the point-wise benchmark is stable.

## Claim Guidance for M3

Safe current claim:

Hydrosheaf can now ingest the USGS national groundwater-age benchmark without row multiplication, use North America WISER tritium input, apply USGS dissolved-gas-corrected tracers, select stronger 14C corrections, and use calibrated 4He accumulation for old groundwater cases.

Not yet safe:

Hydrosheaf has not yet demonstrated full-dataset superiority over USGS LPM estimates, because the full post-fix benchmark has not been regenerated and audited.

Not yet safe:

Hydrosheaf graph regularization improves USGS age inference in general. Current graph evidence supports only a conditional claim: graph priors help under physically plausible topology and can degrade performance under wrong topology or excessive prior strength.

Target final claim after the plan is implemented:

Hydrosheaf improves groundwater-age inference under identifiable tracer and graph conditions, especially when corrected multi-tracer inputs and calibrated old-groundwater constraints are available; its benchmark framework also identifies conditions where graph priors or tracer assumptions degrade inference.
