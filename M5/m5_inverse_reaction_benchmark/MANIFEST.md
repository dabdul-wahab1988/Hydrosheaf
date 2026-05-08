# M5 Inverse Reaction Benchmark Manifest

Generated outputs:

- `results/l1_penalty_sensitivity.csv`: sparse reaction fits across L1 penalty
  values.
- `results/missing_ion_sensitivity.csv`: refits after dropping selected ions.
- `results/thermodynamic_bound_violations.csv`: PHREEQC SI-bound compatibility
  checks.
- `tables/table1_inverse_reaction_validation_summary.csv`: manuscript-ready
  summary of the best sparse fit and diagnostic flags.
- `docs/m5_results_summary.md`: wording guardrails and benchmark summary.

Claim rule:

M5 should describe Hydrosheaf as a sparse linear inverse reaction model with
PHREEQC thermodynamic screening and forward-validation diagnostics, not as a
fully coupled nonlinear PHREEQC inverse solver.
