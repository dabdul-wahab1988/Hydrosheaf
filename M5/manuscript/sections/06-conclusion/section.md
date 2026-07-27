# Conclusion

Sparse linear inversion recovered true mineral phases and extents reliably
only where the reaction dictionary was well-conditioned for the available
ion panel; noise, missing ions, and stoichiometric collinearity degraded
exact-phase support well before they degraded concentration reconstruction,
and the L1 penalty traded recall for false discoveries along a path shaped
by the dictionary's own null-space structure rather than by noise alone.
Thermodynamic constraints improved mechanistic recovery only partially: they
reliably removed physically impossible reaction directions but left most
constrained fits below a defensible equivalence-class threshold, and the
conventional PHREEQC inverse baseline's own high model multiplicity
confirmed the same limit from full thermodynamic consistency. The calibrated
Mechanism Resolution Score distinguished identifiable from ambiguous
pathways above chance on an untouched archetype, but only moderately, so it
should inform rather than replace inspection of the underlying diagnostics.
Fluoride carried the greatest realised measurement value in the Northern
Ghana transfer analysis, and next-best-measurement ranking outperformed
fixed and random panels; edge-level chemical confidence was not conditioned
on topology or residence-time uncertainty in this chemistry-only study,
leaving that integration for future work. Identifiability-aware reporting
produced more conservative, better-supported field interpretations than a
single best-fit model on the same Northern Ghana pairs. Reliable inverse
hydrogeochemistry requires evidence for identifiability, not only a low
residual; PHREEQC-style thermodynamic screening is necessary but
insufficient on its own. Hydrosheaf is offered as a sparse linear inverse
reaction model with transparent identifiability diagnostics and PHREEQC
screening and forward-consistency checks, contributing a reusable template
for reproducible groundwater-process inference under sparse, noisy chemical
observation.
