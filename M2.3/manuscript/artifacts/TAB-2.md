| Capability | Package module | Output | Claim scope |
|---|---|---|---|
| Graph construction and candidate edges | hydrosheaf.inference.network_fit | Candidate directed connectivity | Screening-level candidate set |
| Residence-time inference | hydrosheaf.temporal; hydrosheaf.nuclear | Single-node and network-constrained age | Calibrated bounded synthetic inference |
| Reaction candidate generation | hydrosheaf.validation.reaction_competent_baseline | Stoichiometric and thermodynamic families | Family-level candidate set |
| Reaction-family scoring (RAPM) | hydrosheaf.validation.reaction_rapm | Calibrated family probabilities | Calibrated inference with selective risk |
| Kinetic adapter | hydrosheaf.reactive_transport.kinetic_phreeqc | Forward kinetic response | Prediction; conditional parameter identification |
| Interval calibration | hydrosheaf.validation.calibration; specialist_calibration | Coverage, width, selective risk | Frozen on development cases |
| Model discrepancy | hydrosheaf.validation.discrepancy | Disagreement scale | Applied before averaging |
| Discrete model averaging | hydrosheaf.validation.model_averaging | Case-blocked weights | Convergence is a gate, not an assumption |
| Prospective measurement selection | hydrosheaf.validation.decision_utility | Truth-blind cost-adjusted policy | Bounded synthetic utility |
| Claim and failure gates | hydrosheaf.validation.performance_contract; programme_gate | PASS / WEAK / ABSTAIN | Programme-level claim control |
