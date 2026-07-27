# M8 decisions

## 2026-07-27 — Separate exploration from confirmation

The earlier 16-design transport sweep, 30-replicate OED output, pooled
two-parameter errors, and ad hoc parameter-specific correlations were treated
as exploratory. They motivated the confirmatory design but are not inputs to its
reported estimates.

## 2026-07-27 — Use log-scaled, whitened information

Raw `J.T @ J` was rejected as the confirmatory information matrix because the
two parameters have different units and scales. The confirmatory matrix uses
sensitivities to log10 parameters divided by the declared observation standard
deviations.

## 2026-07-27 — Remove start-centred regularisation

The production adapters construct prior means from their base values. When the
base value is also a random starting value, this changes the prior between
replicates. The primary confirmatory experiment therefore retains the production
forward adapter and PESTGLM optimiser but disables these priors.

## 2026-07-27 — Treat kinetics as a structural negative control

The PHREEQC rate laws contain `k*A` multiplicatively. Residence-time sampling
alone cannot separate the two constants. Kinetics is therefore used to test the
limit of sampling design and the effect of an independent surface-area
constraint, not to claim a second transport-like Pareto frontier.

## 2026-07-27 — Presentation-only amendment after B1 review

The immutable confirmatory run passed, but human review rejected two display
choices. The first figure title implied that target-specific criteria selected
different times even though every criterion selected 50 d. The OED table also
showed the first random draw as if it were the strategy's fixed time. The locked
raw outputs and all numerical estimates were preserved. A checksummed
postprocessor produced reviewed artifacts with an explicit common optimum and a
blank candidate time for the varying random strategy.
