# M3 Age Benchmark Manifest

Generated outputs:

- `results/graph_regularization_scenarios.csv`: row-level scenario results.
- `tables/table1_m3_benchmark_summary.csv`: compact table for manuscript use.
- `figures/figure1_graph_regularization_delta_rmse.png`: optional figure created
  when Matplotlib is available.
- `docs/m3_results_summary.md`: narrative summary with claim guardrails.

Claim rule:

Hydrosheaf should only claim graph regularization improves nuclear-tracer age
inference where benchmark rows show lower RMSE than single-node inference and
where negative-control randomized graphs do not show comparable improvement.
