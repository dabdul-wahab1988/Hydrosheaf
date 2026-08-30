# M7 public-pipeline system-acceptance protocol

Run ID: `RUN-M7-SYSTEM-20260728-01`

## Question

Does the public `fit_network_pipeline` execute the nuclear-age, sheaf-refinement,
and network-fit stages as one strict system on code-independent synthetic aquifer
cases, and does the resulting evidence improve topology discrimination beyond a
hydraulic-only baseline and an age-permuted adverse control?

This is a controlled-synthetic system-acceptance benchmark. It is not field
validation and cannot support a general-superiority claim.

## Frozen design

- Generator: the existing M7 MODFLOW 6/MODPATH 7 independent generator, whose
  provenance declares `imports_hydrosheaf=false`.
- Fresh locked seeds: 6301--6306.
- Candidate topology: HydroSheaf probabilistic hydraulic candidate generation,
  with truth withheld.
- Evidence conditions:
  - `hydraulic_only`: public pipeline network fit without sheaf refinement;
  - `local_age`: strict nuclear-age, local age-sheaf, and network-fit stages;
  - `full_sheaf`: strict nuclear-age, local plus global chemistry-section sheaf,
    and network-fit stages;
  - `age_permuted`: full-sheaf pipeline after within-case permutation of the
    tritium observations.
- Nuclear inference: edge-free exact-grid local posteriors, 600 draws per chain,
  four chains, and the independently specified synthetic recharge boundary.
- Metrics: candidate recall, PR-AUC, Brier score, log loss, and selected-edge F1.
- Uncertainty: 5,000 paired case bootstraps.

## Gate and claim rule

System acceptance requires all requested stages to report `completed`, finite
metrics, no truth fields in inference rows, generator independence, and mean
candidate recall at least 0.80.

An incremental full-sheaf topology claim is allowed only if all of the following
hold on the locked seeds:

1. the 95% paired-bootstrap interval for `full_sheaf - hydraulic_only` PR-AUC
   is strictly positive;
2. the corresponding Brier-score interval is strictly negative;
3. full-sheaf PR-AUC exceeds age-permuted PR-AUC with a strictly positive paired
   interval; and
4. full-sheaf does not reduce selected-edge F1 relative to hydraulic-only.

The global-section contribution is reported separately as `full_sheaf -
local_age`; it is not inferred from the full-versus-hydraulic contrast.

If the gate passes but the scientific rule fails, the allowed result is limited
to verified system execution and a negative or conditional ablation finding.

## Exclusions

Temporal-series sheaves, three-dimensional graphs, and vadose-zone modules are
out of scope by instruction. No field dataset is used to substitute for them.
