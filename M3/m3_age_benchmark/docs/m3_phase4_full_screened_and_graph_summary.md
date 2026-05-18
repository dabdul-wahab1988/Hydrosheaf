# M3 Phase 4 Full Screened Canonical Benchmark and Real-USGS Graph Benchmark

- Generated: 2026-05-17
- Pointwise full-screened results: `M3/m3_age_benchmark/results/m3_phase4_screened_full_results.csv`
- Pointwise summary: `M3/m3_age_benchmark/results/m3_phase4_screened_full_results_summary.csv`
- Pointwise QA: `M3/m3_age_benchmark/docs/m3_phase4_screened_full_results_qa.md`
- Graph benchmark results: `M3/m3_age_benchmark/results/m3_real_usgs_graph_benchmark.csv`
- Graph benchmark QA: `M3/m3_age_benchmark/docs/m3_graph_benchmark_qa.md`

## Pointwise full-screened canonical result

The full canonical age-grid-35 refresh completed with the `screened_dgm_gases` scenario. Finite multi-tracer estimates were produced for `1272/1272` eligible rows.

Main metrics:

- median absolute log10 error: `0.382737`
- log10 RMSE: `1.064566`
- within factor 2: `0.409591`
- within factor 10: `0.677673`
- calibrated `4He` rows: `838`

This run has only the screened scenario, so paired comparisons against `parity_reported_corrected` are not available in this output.

Interpretation: the full public benchmark is materially weaker than the synthetic age benchmark and must be described as screening-level public tracer-age agreement, not strong TracerLPM equivalence.

## Real-USGS graph benchmark

The real-USGS graph benchmark was run on the full screened pointwise output (`1272` finite nodes, `20` study units, `5155` benchmark edge rows across candidate and control families).

Main result:

- no candidate graph family improved the pointwise baseline once nonzero graph regularization was applied
- best candidate rows at prior weight `0.0` are identical to the pointwise baseline by construction
- weak/medium/strong priors generally degraded RMSE; the smallest nonzero candidate change was nearly neutral but not improved
- wrong-direction and randomized negative controls degraded much more strongly, which is the expected guardrail behavior

Important negative result:

- `parameter_smooth_aicc`, weak prior: `delta_rmse_graph_minus_single = +0.003177`
- `study_unit_coordinate`, weak prior: `+0.159690`
- `coordinate_global`, weak prior: `+0.161825`
- randomized and wrong-direction controls remain negative-control guardrails, not success evidence

Interpretation: the real-USGS graph benchmark is now implemented and useful, but it currently supports a conservative claim: graph regularization is benchmarked and falsifiable, yet it does not improve the full USGS pointwise age fit under the present `StudyUnit + depth-proxy` construction.

## What this fixes

This closes the earlier limitation that Phase 4 had only been demonstrated on a targeted subset. M3 now has:

- explicit contaminated/saturated `SF6/CFC` masking in the screened young-gas path
- a refreshed full canonical screened benchmark on all `1272` eligible USGS rows
- a real-USGS graph benchmark with candidate families and negative controls

## What remains

The main remaining issue is no longer missing infrastructure. It is evidential:

- keep the screened pointwise pathway as the canonical M3 benchmark result
- report the graph benchmark as a negative or neutral finding, not as a success claim
- only revisit graph regularization after adding stronger hydrogeologic edge information than `StudyUnit + depth proxy`
