# M4 Topology Benchmark Manifest

Generated outputs:

- `results/independent_graph_vs_modpath.csv`: scenario metrics for independent
  Hydrosheaf graph inference against MODPATH reference connectivity.
- `results/edge_classification.csv`: TP/FP/FN/TN edge-level classifications.
- `results/modpath_informed_priors.csv`: MODPATH-informed prior-mode output;
  these rows are explicitly not independent validation.
- `tables/table1_topology_validation_summary.csv`: manuscript-ready summary of
  precision, recall, F1, false-positive rate, and false-negative rate.
- `docs/m4_results_summary.md`: narrative guardrails for M4 claims.

Claim rule:

M4 can claim reduced-order topology reproduction only for independent graph
inference rows. MODPATH-informed priors should be described as a prior-assisted
Hydrosheaf mode, not as independent evidence of topology inference skill.
