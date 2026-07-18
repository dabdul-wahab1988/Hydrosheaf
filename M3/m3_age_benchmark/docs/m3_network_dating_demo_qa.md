# M3 Network-Enhanced Dating — Ambiguity Demonstration (QA)

## Why this exists
The real-data graph benchmark (`m3_real_usgs_graph_benchmark.csv`) shows that graph
priors do **not** improve age accuracy on well-constrained USGS tracers — RMSE and
within-factor-2 are unchanged or slightly worse — and that their only benefit is
enforcing age-ordering consistency (severe violations 306 → ~0). Taken alone that
makes the "network-enhanced dating" contribution look weak. This controlled
demonstration establishes the **condition under which network dating genuinely
improves accuracy**, so the claim is scoped rather than overstated.

## Design
- Controlled twin, 30 nodes in 9 flow chains, known true ages, real Hydrosheaf tritium
  model (`nuclear.lpm.convolve_input`, `nuclear.input_history`), seed 1234.
- Ages deliberately include the **tritium bomb-peak ambiguity zone** (~45–56 yr), where
  a measured value aliases to a much younger solution (e.g. 8.4 TU → ~3 yr *or* ~50 yr),
  bracketed between a clearly-young upstream and clearly-old downstream neighbour.
- Single-node = global best tritium match (degenerate; may pick the wrong alias).
- Network = single-node candidate set filtered by the flow-ordering constraint
  (downstream age ≥ upstream age), selecting the correct alias.
- Capability demonstration on synthetic truth, not field validation.

## Result (`m3_network_dating_demo_summary.csv`, Suppl_Fig4)
| Subset | n | Within-factor-2 single | Within-factor-2 network | RMSE-log single | RMSE-log network |
|---|---:|---:|---:|---:|---:|
| Unambiguous | 11 | 1.00 | 1.00 | 0.04 | 0.04 |
| **Ambiguous** | 19 | **0.63** | **0.84** | **0.51** | **0.27** |
| All | 30 | 0.77 | 0.90 | 0.41 | 0.22 |

Age-ordering violations: single-node 4 → network 1.

## Interpretation (honest, scoped)
- Network-enhanced dating provides **no accuracy benefit where single-node inversion is
  already well-posed** (unambiguous nodes: identical 100%).
- It provides a **substantial accuracy benefit in the tracer-ambiguous (bomb-peak)
  regime** (within-factor-2 0.63 → 0.84), by using graph flow-ordering to select the
  physically consistent age alias.
- This reconciles the real-USGS result: those USGS tracers are largely outside the
  bomb-peak ambiguity zone, so graph priors there deliver **consistency, not accuracy**.

Together the two results give the complete, defensible claim: *graph priors enforce
flow-consistent age fields everywhere, and improve age accuracy specifically where the
tracer signal is ambiguous.*
