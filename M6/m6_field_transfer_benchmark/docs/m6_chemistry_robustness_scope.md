# Supporting chemistry robustness

The files below were consolidated from the former separate robustness
scaffold:

- `scripts/run_m6_chemistry_robustness.py`
- `scripts/run_m6_chemistry_stress_tests.py`
- `results/m6_regularization_path.csv`
- `results/m6_phase_stability_index.csv`
- `results/m6_loo_structural_uncertainty.csv`

They remain because M2 publication assets consume the regularisation-path and
phase-stability results. They are cross-milestone supporting analyses, not part
of the M6 field-transfer Q1 workflow, and `scripts/run_m6_q1.py` does not invoke
them.

The former robustness figure generator and its three demonstration PNGs were
removed. Those graphics used methodological or randomly generated display data
and are not submission assets.
