# Introduction

## General studies: groundwater reactions and water-resource interpretation

Identifying which mineral dissolution, precipitation, ion-exchange, redox, and
anthropogenic processes have shaped a groundwater sample is central to
groundwater-resource assessment [@plummer1980massbalance]. Hydrochemical
evolution along a flow path records the joint influence of recharge
composition, lithology, mixing, residence time, climate, and human pressure,
so a defensible reconstruction of that evolution requires separating these
influences rather than fitting a single aggregate signal. Reaction-pathway
interpretation in turn informs contaminant-source attribution, aquifer
vulnerability assessment, treatment-process selection, and long-term
water-quality prediction, all of which depend on knowing which transformation
actually occurred rather than which transformation is merely consistent with
the observed concentrations. Carbonate dissolution, silicate weathering,
cation exchange, evaporite dissolution, redox transformation, and
nitrate-source attribution each imply a different management response, so
misassigning the operative reaction can misdirect remediation spending,
mislabel a naturally elevated solute as contamination, or understate a
genuine anthropogenic source. Inverse hydrogeochemistry has therefore become
the general scientific approach for reconstructing groundwater transformations
that cannot be observed directly. This study develops that approach within
Hydrosheaf, a graph- and sheaf-based framework in which candidate flow
connectivity, residence time, and reaction pathway are separately evidenced
layers attached to a shared directed edge between two sampled points; the
present work is the chemical-reaction layer.

## Current studies: contemporary inverse hydrogeochemical approaches

Contemporary practice draws on PHREEQC inverse mole-balance modelling
[@parkhurst2013phreeqc] and its predecessor NETPATH [@plummer1994netpath],
mineral saturation indices, end-member mixing analysis
[@christophersen1992emma], isotope constraints, sparse regression such as the
lasso and elastic net [@tibshirani1996lasso; @zou2005elasticnet], and Bayesian
model-averaging treatments of structural uncertainty [@neuman2003bma], and
more recent Bayesian mixing frameworks that explicitly propagate source and
tracer uncertainty rather than reporting a single mixing ratio
[@popp2019bayesianmixing; @beria2020hydromix]. Recent hydrochemical
metamodelling has also estimated groundwater-age distributions from chemistry
using symbolic regression and gradient boosting [@tschritter2023age], while
Northern Ghana studies have applied machine learning to groundwater-quality
prediction in the Nabogo Basin and crystalline aquifer terrain
[@apogba2024nabogo; @abu2025ghana]. Those predictive studies target age or
water-quality outcomes; they do not test recovery of a known reaction support
or attach an identifiability class to an inverse reaction estimate.
Thermodynamic screening,
regularisation, and multi-tracer evidence are widely
used to narrow the set of plausible reactions and to make reported pathways
more parsimonious. Recent studies increasingly combine hydrochemistry with
spatial, hydraulic, isotopic, or age information, and graph-based
representations of candidate flow connectivity, developed independently for
transport-scenario screening rather than for chemistry [@moracchini2025graphflow],
offer one route to representing that connectivity explicitly. Even so,
connectivity, travel time, and reaction pathway are typically estimated and
reported as separate outputs rather than as a jointly conditioned inference.
There is a parallel movement toward reproducible workflows, ensemble
interpretation, and data-guided monitoring-network design that explicitly
values the marginal information content of additional measurements
[@sreekanth2017monitoring]. Collectively, these advances improve the
plausibility and parsimony of an inferred reaction set, but plausibility and
parsimony do not by themselves establish that the inferred mechanism is the
unique explanation of the observations. A model that fits well because it is
correct and a model that fits well because the available ions cannot
distinguish it from a competing reaction combination remain, on residual
error alone, indistinguishable.

## Problem statement: numerical fit is not mechanistic identification

The central problem is that multiple stoichiometrically equivalent or
near-equivalent reaction combinations can reproduce the same concentration
difference between two waters. Christophersen and Hooper showed for
end-member mixing that unambiguous identification of source composition from
a mixture alone is impossible without independent constraints
[@christophersen1992emma], and the wider equifinality literature establishes
that many distinct parameter or structural choices can be observationally
indistinguishable in an environmental model [@beven2001glue; @beven2006manifesto].
The same structural indeterminacy applies to sparse inverse reaction
inference: a low concentration residual, a sparse support, or a
thermodynamically feasible solution may still identify the wrong reaction, or
only one member of an unresolved equivalence class. Incomplete ion panels,
analytical noise, uncertain phase dictionaries, transport-model error, mixing,
and correlated stoichiometric vectors all contribute to this non-identifiability.
Conventional missing-ion tests can appear favourable for the wrong reason,
because setting an omitted ion's weight to zero removes it from the objective
function rather than testing whether the fitted reaction predicts it. Without
ground-truth phase-recovery testing, calibrated resolution diagnostics, and
genuine predictive falsification, an inverse hydrogeochemical study can report
an overconfident mechanistic narrative that a numerically equivalent
alternative would equally well support. The practical consequence is not
merely academic: a management decision built on the wrong member of an
equivalence class can misattribute a water-quality trend, target the wrong
recharge source for protection, or miss a genuine early warning of
contamination because the fitted narrative absorbed it into an
indistinguishable natural process.

## Novelty: identifiability-aware and predictively falsifiable reaction inference

This study reframes the inferential target: rather than selecting one
plausible reaction pathway, it determines what level of mechanism the
available observations can actually resolve. Five innovations, delivered and
tested in this manuscript, support this reframing: (1) reaction-equivalence
classes derived from rank, singular values, null-space structure, and mutual
coherence of the reaction dictionary; (2) a live-PHREEQC factorial benchmark
with known phase and extent ground truth; (3) a calibrated Mechanism
Resolution Score; (4) held-out-ion predictive falsification in place of
retained-ion residual comparison; and (5) retrospective next-best-measurement
selection based on expected ambiguity reduction. A sixth position, edge-
conditioned chemical confidence within the wider Hydrosheaf framework, is
introduced conceptually here but is not evaluated in this chemistry-only
study; it is left for the integration with topology and residence-time
evidence described in Future development. Positioned this way, the study
functions as a chemical claim-auditing layer that distinguishes an
identifiable mechanism, a partially resolved reaction class, and a
non-identifiable pathway, rather than as a new reactive-transport solver.

## Aim and objectives

The aim is to develop and validate an identifiability-aware framework that
determines when groundwater reaction mechanisms are uniquely recoverable,
recoverable only as equivalence classes, or unsupported by the available
observations. Five objectives follow. First, construct a diverse
live-PHREEQC benchmark with known reaction phases, directions, and extents.
Second, quantify how stoichiometric structure, analytical noise, missing
ions, transport confounding, and regularisation alter exact-phase and
equivalence-class recovery. Third, develop and independently test a
Mechanism Resolution Score against held-out synthetic archetypes. Fourth,
evaluate held-out-ion prediction and next-best-measurement selection as
falsification and monitoring-design tools. Fifth, demonstrate chemistry-only
identifiability-aware reporting using 320 wet/dry hydrochemical records from
160 Northern Ghana boreholes sampled in March and August 2025, a crystalline
basement setting in which silicate weathering and cation exchange are
established controls on groundwater composition [@anku2008ghana;
@banoengyakubo2011ghana].

## Significance

The evidence base for these claims is primarily the controlled live-PHREEQC
benchmark, where ground truth is known; the Northern Ghana component
demonstrates chemistry-only transfer under realistic data sparsity and does
not itself validate reaction mechanism, so it is not treated as independent
support for the significance claims below. Scientifically, the framework
replaces unqualified reaction-pathway selection
with a measurable resolution hierarchy that separates numerical
reconstruction, thermodynamic feasibility, and mechanistic identification.
Methodologically, it converts equifinality from an implicit caveat into an
explicit, quantified diagnostic attached to every reported reaction.
Practically, it identifies which additional ion or tracer would most improve
a given groundwater interpretation, giving monitoring design an evidence base
rather than a fixed panel. Within Hydrosheaf, it supplies chemically
qualified edge-level evidence without claiming to have independently
validated flow-path topology or groundwater age, preserving the modular
claim boundary between this chemical-reaction layer and Hydrosheaf's
separately reported topology and residence-time layers, which are developed
and validated in companion work. More broadly, it contributes a transferable
template for reproducible, decision-relevant groundwater-process inference
under sparse and noisy chemical observation.
