# M7 Null-Model Sensitivity Extension

This is a paired, synthetic-twin sensitivity analysis. The locked M7
integration outputs are unchanged. Hydrosheaf's formal null model is
enabled and applied to every candidate edge; `joint_plus_null_gate` retains
an edge only when the original geometry/chemistry/age joint gate passes and
the null score is <=0.5.

The isotope, location, and lithology fields are constructed solely to exercise
the null-model branches. They are not observations and do not support a new
field-validation claim.

## Outputs

- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M7/m7_integration_benchmark/results/m7_null_edge_scores.csv`: edge-level scores, flags, and paired gates.
- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M7/m7_integration_benchmark/results/m7_null_sensitivity_summary.csv`: precision/recall/F1 sensitivity summary.
- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M7/m7_integration_benchmark/results/m7_null_trap_rejection.csv`: trap-type rejection comparison.
- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M7/m7_integration_benchmark/results/m7_null_model_run.json`: configuration, source hash, and replay caveat.

## Paired classification summary

| stream               |   precision |    recall |        f1 |   trap_accept_rate |
|:---------------------|------------:|----------:|----------:|-------------------:|
| geometry_only        |    0.692308 | 0.782609  | 0.734694  |             0.5    |
| chemistry_only       |    0.904762 | 0.826087  | 0.863636  |             0.125  |
| age_only             |    0.75     | 0.782609  | 0.765957  |             0.375  |
| joint                |    1        | 0.565217  | 0.722222  |             0      |
| joint_plus_null_gate |    1        | 0.0434783 | 0.0833333 |             0      |
| null_gate_only       |    0.388889 | 0.304348  | 0.341463  |             0.6875 |

## Trap rejection

| trap_type   |   n |   baseline_joint_rejects |   null_gate_rejects |   joint_plus_null_rejects |   mean_null_score |
|:------------|----:|-------------------------:|--------------------:|--------------------------:|------------------:|
| trapA       |   5 |                        1 |                 0.4 |                         1 |          0.359384 |
| trapB       |   5 |                        1 |                 0   |                         1 |          0.216    |
| spurious    |   6 |                        1 |                 0.5 |                         1 |          0.591737 |

## Interpretation guardrail

A null score above 0.5 means that similarity has a competing no-flow
explanation under the screening model; it is not a calibrated probability.
Use this extension to qualify the M7 mechanism demonstration, not to replace
the original controlled-twin result.

The baseline generator uses Python's process-randomised `hash(node_id)` when
assigning reaction archetypes. This extension uses a deterministic hash only
for its replay and records that choice in the run metadata; the baseline source
and baseline result files were not modified.
