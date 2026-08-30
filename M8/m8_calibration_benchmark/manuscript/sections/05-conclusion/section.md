## Conclusions

Calibration is not identification. In controlled synthetic experiments with
production HydroSheaf adapters, we demonstrated that a converged optimisation
with a small residual can coexist with poor parameter recovery, that linearised
95% intervals can undercover even in well-conditioned matched-model designs,
and that the parameter-specific effects of measurement placement do not
survive model-form discrepancy.

The confirmatory results, stated with their negative counterparts as
prominently as their positive ones, are:

1. A 50 d front observation selected by local Fisher criteria improved recovery
   of both dispersivity and decay in the matched analytical transport model
   (dispersivity absolute log10 error 0.2154 to 0.0173; decay 0.0182 to 0.0146),
   but the anticipated split between dispersivity- and decay-targeted next
   observations was not confirmed - every criterion selected the same time.
2. Under independent numerical truth, the same observation robustly improved
   dispersivity (0.8262 to 0.1674) and degraded decay (0.1367 to 0.1541); the
   dual-parameter robustness claim is not supported.
3. In the PHREEQC kinetic adapter, residence-time sampling alone is rank-one
   because the rate law identifies the product of rate constant and reactive
   surface area; an independent surface-area observation restores rank two.
4. The existing categorical active-learning heuristic cannot express an
   actionable transport sampling time and is reported as
   `NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION`.
5. A prospectively specified Bayesian value-of-information workflow for
   topology-measurement design passed every frozen gate: superior to random
   acquisition in Brier score and information per cost, noninferior (not
   superior) to a strong legacy policy within narrow frozen margins, and
   byte-identical on regeneration - with all selected actions being chemistry
   panels, so multimodal usefulness remains undemonstrated.

For modellers, the operational message is that identifiability diagnostics
should accompany every calibration report: the condition number of the
noise-whitened log-parameter Fisher matrix, per-parameter recovery and
coverage jointly with interval width, and explicit model-form sensitivity of
any recommended design. A small residual is a starting point for inference,
not evidence that the parameters have been identified.
