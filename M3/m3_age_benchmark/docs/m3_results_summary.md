# M3 Age Benchmark Results

- Scenario-graph rows tested: 9.
- Rows where graph regularization improved RMSE: 3.
- Rows where graph regularization degraded or failed to improve RMSE: 6.
- Randomized negative-control rows with apparent improvement: 0.

Interpretation guardrail:

Hydrosheaf should claim graph regularization improves age inference only for the
specific tracer/graph/prior settings where RMSE decreases and randomized graphs
do not produce comparable improvement. Degrading scenarios must be reported as
evidence that graph priors can harm inference when topology, recharge structure,
or tracer contamination assumptions are wrong.
