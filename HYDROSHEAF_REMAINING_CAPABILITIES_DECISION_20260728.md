# HydroSheaf remaining-capability benchmark decision

Date: 2026-07-28

## Scope

This work implemented the remaining recommended checks after excluding, by
instruction, temporal-series sheaves, three-dimensional graphs, and vadose-zone
capabilities. No substitute public datasets were downloaded for those excluded
modules.

The completed scope was:

1. strict public-pipeline stage completion and failure reporting;
2. public-pipeline/full-sheaf system acceptance using the M7 independent
   MODFLOW 6/MODPATH 7 generator;
3. local-only, global-sheaf, hydraulic-only, and permuted-age ablations;
4. M8 independent numerical transport truth calibrated with the analytical
   production adapter;
5. existing active-learning heuristic versus E-optimal, random, and
   development-oracle controls; and
6. confirm-or-remove decisions for the remaining M8 pilot claims.

All evidence is controlled synthetic. Nothing below is field validation.

## Public pipeline

`fit_network_pipeline` now supports opt-in strict stage completion through:

- `strict_stage_completion`;
- `required_stages`;
- `sheaf_refinement_enabled`; and
- `nuclear_inference_options`.

It reports stage status for latent endmembers, temporal fitting, nuclear age,
sheaf refinement, and network fitting. A requested failed or skipped stage
raises `PipelineStageError` in strict mode. Nuclear ages are inferred before
topology refinement on an edge-free graph to avoid conditioning local ages on
the unknown candidate topology. Caller-owned sample lists are no longer
mutated. Empty posterior results and empty sheaf outputs cannot be marked as
successful.

## M7 strict system acceptance

Authoritative run: `RUN-M7-SYSTEM-20260728-01`

| Check | Result |
|---|---:|
| Fresh independent cases | 6 |
| Mean candidate recall | 0.9815 |
| Requested stages completed | all |
| Metrics finite | yes |
| Generator imports HydroSheaf | no |
| System acceptance | PASS |

The full global sheaf improved topology discrimination relative to local age
alone:

- PR-AUC difference: +0.05865; 95% paired-bootstrap interval +0.03864 to
  +0.07771;
- Brier difference: -0.00849; interval -0.00983 to -0.00723.

However, full-sheaf PR-AUC was 0.3075 versus 0.3272 for hydraulic-only, and
full-sheaf did not beat the permuted-age control. The preregistered incremental
superiority rule therefore failed.

Decision: claim verified strict public-pipeline execution and a conditional
global-section increment over local age alone. Do not claim overall topology
superiority, general HydroSheaf superiority, or field validation.

## M8 independent-model robustness

Authoritative run: `RUN-M8-INDEPENDENT-20260728-01`

The independent implicit finite-volume truth passed the numerical gate: the
240-cell result differed from the 480-cell reference by at most 0.203 declared
observation standard deviations, below the locked 0.25 limit.

Both analytical E-optimality and the 80-replicate development oracle selected
50 d. On 250 locked test replicates:

| Parameter | No new observation | Add 50 d | Paired difference (95% CI) |
|---|---:|---:|---:|
| Dispersivity absolute log10 error | 0.8262 | 0.1674 | -0.6690 (-0.7374, -0.5757) |
| Decay absolute log10 error | 0.1367 | 0.1541 | +0.0210 (+0.0092, +0.0276) |

Decision: 50 d is robustly useful for dispersivity and was the equal-weight
development oracle, but the earlier matched-model claim that it improved both
parameters did not survive independent numerical truth. Report the opposed
parameter effects explicitly.

## Active learning

The existing `rank_next_measurements` routine produced tied categorical
topology recommendations and no explicit transport concentration sampling-time
action. It does not consume candidate-time Jacobians. Applying an ad hoc mapping
would not be a valid benchmark.

Decision: record `NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION`. Retain the
routine as topology/categorical decision support. A transport OED claim requires
a new Jacobian-aware interface and prospective validation.

## M8 disposition

Retain:

- the matched-model transport result, labelled matched-model only;
- the independent-model parameter-specific qualification;
- the PHREEQC kinetic `k*A` structural-confounding result and independent-area
  intervention; and
- the negative active-learning portability result.

Remove from the current contribution set:

- C1 core dissociation;
- C2 coverage under-calibration;
- C3 covariance-fallback detector claim;
- C3b transferable condition-number threshold; and
- C3c coverage-trap headline.

Those five claims remain exploratory pilot observations until a fresh
prospective protocol confirms them. Corresponding pilot Figures 2--5 and 7 and
Tables 2, 3, 5, and 6 are not confirmed manuscript results.

## Reproducibility and validation

- M7 protocol, source inputs, and seeds were hash-locked before the final run.
- M8 protocol, source inputs, numerical grid, development/test split, and
  replicate counts were hash-locked before the final run.
- Independent regeneration produced identical SHA-256 sets for all seven M7
  and all ten M8 registered artifacts.
- Focused final validation: 15 API/M7 tests passed; M8 confirmatory and
  independent tests had also passed in the combined 21-test run; Ruff passed.
- Full repository suite: 660 passed, 2 skipped. Four failures and 13 setup
  errors were confined to legacy `tests/synthetic_data_tests` because
  `data/synthetic/water_chem_full.csv`, `stations_full.csv`, and related
  fixtures are absent. These failures are unrelated to the changes above.
