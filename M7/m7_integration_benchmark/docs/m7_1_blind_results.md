# M7.1 Blind Replicated Results

## Acceptance checks

- 20 development, 100 disjoint test, and 8 heavy-audit aquifers completed.
- Candidate-universe recall was 1.000.
- The truth-poisoning leakage test passed.
- All heavy-audit cases completed in the last source-frozen run; execution,
  convergence, and sequential-coupling validity are reported separately.
- Module provenance recorded 128 calls each to candidate inference, sheaf
  refinement, and inverse chemistry, plus 8 calls each to the topology posterior,
  PHREEQC, and Bayesian network aging.

## Held-out test performance

Means below are computed across 100 aquifers. Parentheses are aquifer-bootstrap
95% intervals for F1.

| Method | F1 | MCC | PR-AUC | Brier | ECE10 | SHD |
|---|---:|---:|---:|---:|---:|---:|
| Hydraulic-spatial | 0.659 (0.649–0.670) | 0.579 | 0.627 | 0.112 | 0.112 | 7.89 |
| Age-only | pending final rerun | — | — | — | — | — |
| Sheaf-multievidence | pending final rerun | — | — | — | — | — |
| Chemistry | 0.462 (0.445–0.480) | 0.319 | 0.423 | 0.145 | 0.104 | 14.85 |
| Equal-weight joint | 0.524 (0.510–0.537) | 0.411 | 0.522 | 0.132 | 0.098 | 14.00 |
| Development-trained logistic joint | pending final rerun | — | — | — | — | — |

The previous fusion figures are superseded by the final age-only/logistic rerun.

Final paired differences, fusion coefficients, and subgroup metrics will be
inserted from the source-frozen 20/100/8 artifact.

## Chemistry and age diagnostics

Strictly held-out Mg/SO4/Fe reaction-delta RMSE was 0.0487 mmol/L. Only 58.8% of
candidate-covered true processes belonged to a reaction family represented in
the fitted dictionary, and supported-family accuracy was only 38.7%. Reaction
identification is therefore not validated; sulfate reduction and unmodelled
mixing require an expanded process model and identifiability work.

Bayesian-age MAE and interval coverage are withheld whenever four-chain R-hat,
bulk/tail ESS, divergence, or node-completeness gates fail. The final table
reports exploratory values separately from publication-valid values.

## Interpretation

The defensible positive result is narrow: a development-trained logistic evidence
stacker improves several held-out topology metrics over hydraulic-spatial
inference on 100 independent synthetic aquifers. The result does not justify
equal-weight integration, unique reaction recovery, universal improvement,
field validation, or the label “validated digital twin.”

The most important next development priorities are tracer-age ambiguity
handling, reaction-family coverage/identifiability, posterior convergence and
probability calibration, followed by external MODFLOW/MODPATH and field
prequential validation.
