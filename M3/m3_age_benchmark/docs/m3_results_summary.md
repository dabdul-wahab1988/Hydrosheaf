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

## Network-enhanced dating: consistency (real data) + accuracy in the ambiguous regime

The real USGS graph benchmark shows graph priors enforce age-ordering consistency
(severe violations 306 -> ~0) without improving point accuracy, because those tracers
are largely outside the tritium bomb-peak ambiguity zone. A controlled ambiguity
demonstration (`run_m3_network_dating_demo.py`; `m3_network_dating_demo_summary.csv`;
`Suppl_Fig4_Network_Dating_Ambiguity`; QA in `m3_network_dating_demo_qa.md`) establishes
the scoped positive claim: in the tritium-ambiguous regime, flow-ordering raises
within-factor-2 accuracy from 0.63 (single-node) to 0.84 (network). Combined claim:
graph priors deliver flow-consistent age fields everywhere, and improve age accuracy
specifically where the tracer signal is ambiguous.
