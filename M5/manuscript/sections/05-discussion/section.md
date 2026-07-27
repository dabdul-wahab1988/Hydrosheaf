# Discussion

## Principal finding

Numerical fit, sparsity, and thermodynamic feasibility answered different
questions in this benchmark, and none of them alone certified a uniquely
correct reaction mechanism. Reconstruction accuracy, phase-support accuracy,
and equivalence-class accuracy diverged systematically: 55.0% of the
lowest-residual quartile of fits still had exact-phase F1 below 0.80, the
unpenalised bounded-least-squares comparator achieved the best reconstruction
RMSE of any method while producing the worst phase-support outcome, and the
guarded model's equivalence-class F1 (0.609) consistently exceeded its own
exact-phase F1 (0.563) by a margin traceable to a specific, identified
null-space direction rather than to unexplained noise. Thermodynamic
screening removed a real category of false positives but left 72% of
constrained fits below an equivalence-class F1 of 0.80, and the conventional
PHREEQC inverse baseline's own output, a mean of 185.8 feasible models per
scenario, demonstrated the same underlying multiplicity from a fully
thermodynamically consistent solver. A plausible sparse solution is
therefore not automatically a uniquely supported mechanism; it is, at best,
one member of a structurally determined equivalence class whose size and
composition can themselves be measured and reported. The clearest
illustration was not an exchange-direction ambiguity but an
anorthite/calcite pair: silicate weathering and carbonate dissolution were
exactly indistinguishable under the measured panel because both reactions
affect only calcium and bicarbonate among the measured species, and
additional diagnostics lifted this ambiguity only slightly (Figure 3). This finding extends,
in a sparse-inverse-reaction setting, the broader equifinality result that
many environmental models can be observationally indistinguishable despite
producing very different internal explanations [@beven2001glue;
@beven2006manifesto]: here the indistinguishable alternatives are not
competing parameter sets within one model structure but competing reaction
combinations within one fixed stoichiometric dictionary, and the null-space
analysis makes that indistinguishability exact and quantifiable rather than
only empirically observed across resamples. The moderate transfer accuracy
of the Mechanism Resolution Score (48.9% on the untouched mixed archetype,
against 25% uniform-random and 37.6% majority-class baselines) is best read
alongside, not instead of, the
underlying null-space diagnostics: the score compresses six diagnostics into
a single ordinal prediction for convenience, but the raw rank, coherence,
and stability diagnostics that feed it remain individually inspectable
whenever the compressed score alone is insufficient for a specific
interpretive decision, so a single moderate accuracy number, reported
without the diagnostics behind it, would understate what the framework
actually delivers.

## Advance over existing inverse hydrogeochemical practice

This framework does not replace PHREEQC inverse mole-balance modelling
[@parkhurst2013phreeqc] or NETPATH-style mass-balance reconstruction
[@plummer1994netpath; @plummer1980massbalance], and its lasso, elastic-net,
and stability-selection components are pre-existing statistical tools
[@tibshirani1996lasso; @zou2005elasticnet; @meinshausen2010stability]
combined into a reaction-inference identifiability layer. Conventional
practice, including single-support lasso-type reporting, typically presents
one preferred reaction set per water pair. The results above show that this
presentation style discards information the solver itself already
contains: the same conventional PHREEQC baseline that returned a single
first-minimal model with equivalence-class F1 of 0.512 also generated an
oracle-selectable minimal model reaching 0.586 among its own feasible set,
and the guarded model's structural diagnostics identified the same
non-identifiable directions independently of any specific fit. PHREEQC
already exposes model multiplicity by returning multiple feasible inverse
models. The contribution here is therefore not the discovery of
non-uniqueness, but a shift from merely enumerating feasible models to
hierarchical reporting across four resolution classes, non-identifiable,
partially identifiable, equivalence-class identifiable, and identifiable,
each backed by an explicit diagnostic rather than an implicit confidence
claim. Held-out-ion prediction provided a stronger falsification test than
goodness of fit on the ions used for inversion, because a reaction
combination can minimise the fitting objective while still failing to
predict an independent ion, precisely the situation that conventional
missing-ion sensitivity tests, which merely zero-weight an omitted ion, are
structurally unable to detect. End-member mixing analysis reached a
structurally similar conclusion for a related but distinct problem, showing
that unambiguous identification of source composition from a mixture alone
is impossible without independent constraints [@christophersen1992emma]; the
present framework applies the analogous logic to reaction inference and
makes the resulting ambiguity explicit through equivalence classes rather
than treating non-uniqueness as an unstated caveat.

## Scientific and monitoring implications

Because the diagnostic ions available for a given water pair directly
control which reactions are identifiable, ion selection is itself an
inferential choice with a measurable consequence, not merely a laboratory
convenience. The Mechanism Resolution Score and the reaction-equivalence
classes together prevent unsupported attribution of a specific mechanism,
salinisation, carbonate dissolution, sulphate sourcing, denitrification, or
ion exchange, to a groundwater sample when the available ions cannot
distinguish it from an alternative. For monitoring-programme design, the
retrospective next-best-measurement analysis showed that ranking candidate
ions by expected ambiguity reduction produced larger realised class-F1 gain
than the non-selected candidates in the same scenarios. This is analogous to
information-directed monitoring design, while Sreekanth et al. address
spatial well-network placement rather than geochemical tracer selection
[@sreekanth2017monitoring]. Fluoride's disproportionately high
measurement-value score in the Northern Ghana transfer analysis illustrates
how such a ranking can surface a specific, actionable sampling priority
rather than a generic call for "more data." The practical recommendation
that follows is to report robust reaction classes and named alternative
pathways rather than a single deterministic narrative wherever the
underlying diagnostics indicate persistent ambiguity, and to let the
next-best-measurement ranking, rather than habit, guide which additional
analyte is added to a monitoring panel when resources allow only one. The
gap between the modest gain in equivalence-class F1 across data tiers
(0.606 to 0.614) and the larger gain in conditionally-preferred-or-resolved
fraction (55.1% to 62.8%) has a direct monitoring-design reading: additional
diagnostics of the kind tested here are better justified as evidence for
choosing among already-identified alternatives than as a route to detecting
qualitatively new reactions, so a monitoring programme facing a binding
budget constraint should weigh this distinction explicitly rather than treat
every additional analyte as offering the same kind of benefit.

## Interpretation of the field application

In the Northern Ghana demonstration, robust findings included generally
coherent Hydrosheaf-Core evidence and TDS-consistency scores across four
regions in terrain broadly characterised by crystalline-basement and
Voltaian sedimentary geology in which silicate weathering and cation
exchange are established regional controls [@anku2008ghana;
@banoengyakubo2011ghana], and a uniform classification of all 160 pairs as
partially identifiable rather than fully resolved. Geological plausibility,
independently PHREEQC-derived thermodynamic screening, and the sparse
geochemical evidence gates narrowed the credible reaction space for these
pairs but did not eliminate ambiguity. The single wet-to-dry residual per
well does not identify separate seasonal mechanism supports. The analysis
therefore makes no claim for or against a wet/dry mechanism shift. Held-out-
ion performance in this setting is interpreted strictly as
predictive consistency, not as proof that the reported reaction class is the
true one, because no independent mineralogical or reaction-rate observation
exists to validate it. Northern Ghana therefore demonstrates field transfer
of the chemistry identifiability audit itself, not independent validation of
groundwater connectivity, residence time, or reaction mechanism, and the
external field-transfer datasets' consistently low ELRI values reinforce
that this demonstration operated near the limit of what sparse, real-world
diagnostic coverage currently supports.

## Limitations

The reaction dictionary was fixed and linear, and screening was
equilibrium-based rather than kinetic, so rate-dependent or strongly
disequilibrium processes were not directly represented. The framework
performs fully coupled nonlinear inverse or reactive-transport modelling in
neither its guarded nor its Hydrosheaf-Core form, and thermodynamic-database
uncertainty was addressed only partially, through a single fixed database
with sensitivity analysis treated as a secondary experiment. Transport-model
misspecification remains a possible source of error in the residual defined
by Equation 1, and synthetic-to-field transferability is constrained by the
gap between the controlled PHREEQC benchmark and the heterogeneity of a real
crystalline-basement setting. Field conclusions depend directly on measured
ion coverage, and the Northern Ghana and external field datasets lack
matched age tracers, independently observed flow paths, or direct
mineralogical or reaction-rate observations, so field interpretation in this
study rests on chemistry and thermodynamic plausibility alone rather than on
independently corroborated mechanism. The Mechanism Resolution Score itself
carries a further limitation: its 48.9% four-class transfer accuracy on the
mixed archetype, while above the 25% uniform-random and 37.6% majority-class
baselines, is a moderate
result obtained on one held-out synthetic archetype, and its generalisation
to hydrochemical settings substantially different from the four archetypes
tested here, including many real crystalline-basement and mixed-lithology
aquifers, has not been established. Generated Hydrosheaf graph edges and
heuristic hydraulic residence times used in the Northern Ghana component
carried no populated edge-confidence field and were explicitly excluded
from every age-conditioned analysis in this study, a restriction that
narrowed the field demonstration's scope but was necessary to avoid
implying a validation that the available data could not support.

## Future development

Priority extensions include stability-selected or fully Bayesian structural
ensembles in place of a single calibrated score [@neuman2003bma], kinetic
and affinity-dependent bounds in place of purely equilibrium screening, and
explicit propagation of analytical-measurement and thermodynamic-database
uncertainty through the reported resolution classes. Prospective field
testing of next-best-measurement recommendations, rather than only
retrospective simulation, would test whether the realised value observed
here generalises to genuinely unobserved measurements chosen in advance.
Matched chemistry, age-tracer, hydraulic-head, screen-depth, and
mineralogical observations at the same wells would allow the chemical
identifiability layer developed here to be jointly evaluated with
Hydrosheaf's topology and residence-time layers rather than assessed as a
chemistry-only demonstration, and extension to multi-edge groundwater
networks would test whether edge-level resolution classes remain informative
when propagated across a connected flow system.
