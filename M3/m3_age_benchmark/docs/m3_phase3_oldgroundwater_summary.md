# M3 Phase 3 Old-Groundwater Summary

Date: 2026-05-15

## What was completed

Phase 3 from the M3 experimental plan is now implemented in code:

- old-groundwater logic was extracted into `hydrosheaf/nuclear/old_groundwater.py`;
- the M3 loader now preserves multiple plausible `14C` correction candidates per sample;
- the benchmark now reports explicit old-groundwater diagnostics:
  `old_groundwater_case`, `old_groundwater_status`, `old_groundwater_active_constraints`,
  apparent `14C` age, apparent `4He` age, and the log10 gap between them;
- the design matrix now includes:
  `oldwater_c14_ensemble`,
  `oldwater_he4_uncertainty`,
  `oldwater_ensemble_he4_uncertainty`.

## Targeted Phase 3 run

Run artefacts:

- `results/m3_phase3_oldwater_results.csv`
- `results/m3_phase3_oldwater_results_summary.csv`
- `results/m3_phase3_oldwater_results_pairwise_deltas.csv`
- `docs/m3_phase3_oldwater_results_qa.md`

Run configuration:

- 80 stratified rows
- 4 scenarios
- 320 output rows

## Main findings

Baseline parity on this 80-row block:

- median absolute log10 error: `0.1664`
- factor-2 accuracy: `0.6875`

Phase 3 scenarios:

- `oldwater_he4_uncertainty`: unchanged from parity on this block
- `oldwater_c14_ensemble`: slightly worse than parity on average
- `oldwater_ensemble_he4_uncertainty`: same as `oldwater_c14_ensemble`

Paired interpretation:

- `14C` ensemble handling changed the effective `14C` source in `17` paired rows;
- rows with more than one plausible `14C` candidate in this block: `16`;
- calibrated `4He` rows in the block: `21`;
- `4He` uncertainty fractions applied:
  - `0.25` in `45` rows
  - `0.50` in `35` rows

Old-groundwater case coverage in the parity scenario:

- `14C_only`: `55`
- `14C+4He`: `16`
- `4He_only`: `5`
- `none`: `4`

Old-groundwater status coverage in the parity scenario:

- `agreement`: `10`
- `tension`: `2`
- `conflict`: `4`
- `single_tracer`: `60`
- `none`: `4`

## Interpretation

Phase 3 fixed an architectural and evidential gap more than a headline-accuracy gap.
Hydrosheaf can now distinguish `14C`-only, `4He`-only, and combined old-water cases,
and it can report when `14C` and `4He` agree or conflict. That was missing before.

However, this targeted run does **not** show a material benchmark improvement from
`14C` ensemble handling or from inflating calibrated `4He` uncertainty. The old-water
module is therefore better characterized and more defensible, but it is not the main
remaining driver of benchmark error.

## What remains next

The next technical bottlenecks are still:

1. young-tracer correction interpretation (`3H`, `3H/3He`, `SF6`);
2. dissolved-gas correction masking/down-weighting for problematic rows;
3. the real-USGS graph benchmark with negative controls;
4. manuscript tables and figures regenerated directly from these run-specific artefacts.
