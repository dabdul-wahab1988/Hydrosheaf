# M3 Phase 4 Young-Gas Screening Summary

Date: 2026-05-15

## What was completed

Phase 4 was advanced from a diagnostic-only state to an implemented benchmark mode:

- `screened_dgm` is now a real gas-correction mode in the M3 benchmark.
- For each row, Hydrosheaf now evaluates:
  - USGS-corrected young gases (`usgs_dgm`)
  - raw young gases (`raw`)
- The screen keeps USGS-corrected inputs by default, but switches to raw when raw gives:
  - materially better young-tracer proxy coherence; and
  - materially better fit support.

This is implemented in:

- `M3/m3_age_benchmark/scripts/run_m3_usgs_benchmark.py`
- `M3/m3_age_benchmark/configs/design_matrix.yaml`

## Targeted Phase 4 run

Run artefacts:

- `results/m3_phase4_younggas_results.csv`
- `results/m3_phase4_younggas_results_summary.csv`
- `results/m3_phase4_younggas_results_pairwise_deltas.csv`
- `docs/m3_phase4_younggas_results_qa.md`

Run configuration:

- 80 stratified rows
- 3 scenarios
- `parity_reported_corrected`
- `ablation_raw_gases`
- `screened_dgm_gases`

## Main results

Scenario metrics:

- parity corrected:
  - median absolute log10 error: `0.1664`
  - factor-2 accuracy: `0.6875`
- raw-gas ablation:
  - median absolute log10 error: `0.1159`
  - factor-2 accuracy: `0.7375`
- screened DGM:
  - median absolute log10 error: `0.1159`
  - factor-2 accuracy: `0.7625`

Paired effect versus parity:

- raw-gas ablation:
  - mean delta log10 error: `-0.0739`
  - gained factor-2 rows: `6`
  - lost factor-2 rows: `2`
- screened DGM:
  - mean delta log10 error: `-0.0828`
  - gained factor-2 rows: `6`
  - lost factor-2 rows: `0`

## Selection behavior

In the screened scenario:

- selected `raw` in `28` rows
- selected `usgs_dgm` in `52` rows

Selection reasons:

- `default_corrected_parity`: `45`
- `corrected_better_young_proxy_coherence`: `7`
- `raw_better_young_proxy_coherence_and_fit`: `14`
- `raw_tiebreak_fit_improvement`: `14`

Age-class pattern:

- modern rows selected `raw` in `17/20`
- intermediate rows selected `raw` in `3/20`
- old rows selected `raw` in `6/20`
- very-old rows selected `raw` in `2/20`

## Interpretation

This is the first M3 result where Hydrosheaf improves on strict USGS-parity mode
without simply forcing all young gases back to raw values.

The new screen captures the strongest modern-row failures from the corrected-gas
parity path, especially the repeated ~`42.65 yr` modern solutions, while preserving
USGS-corrected inputs in most rows.

On this stratified block, `screened_dgm` is better than parity and safer than the
always-raw ablation because it keeps the factor-2 gains while avoiding the factor-2
losses introduced by unconditional raw-gas fallback.

## Remaining next step

The next phase should extend this into a fuller dissolved-gas handling step:

1. propagate screening/masking logic to contaminated or saturated SF6/CFC rows more explicitly;
2. run a broader Phase 4 matrix including `screened_dgm` with larger row coverage;
3. then move to the real-USGS graph benchmark phase.
