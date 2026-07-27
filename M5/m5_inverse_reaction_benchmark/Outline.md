# Q1 Journal Manuscript Outline

Development status: revised 24 July 2026. Word budget reset from 5,000 to
**6,500 words of main text (Introduction through Conclusion)** on explicit
user instruction, and journal positioning reconsidered accordingly, since
6,500 words exceeds the de facto main-text ceiling that Nature Portfolio
Article-format journals enforce.

## Recommended Journal Positioning

Primary target: **Water Resources Research** (AGU/Wiley), Research Article.

Rationale: AGU journals permit research articles up to 25 publication units
(1 PU = 500 words or one display item), i.e. approximately 12,500 words, so a
6,500-word main text plus six main figures and one main table (7 PU) totals
approximately 20 PU and fits comfortably without an excess-length fee. WRR is
also the closest topical match for M5's actual contribution: it is a venue
that routinely publishes identifiability, equifinality, and structural/
practical non-uniqueness studies in hydrological and geochemical inverse
problems, so the paper is read against the right comparison literature
instead of being squeezed into a shorter species-of-advance format.

Strong alternative: **Journal of Hydrology** (Elsevier), Research Paper.
Elsevier's guidance allows up to 10,000 words for full-length articles and
explicitly asks for a transferability dimension beyond local documentation,
which matches the Northern Ghana chemistry-only field-transfer component.

Shorter-form fallback if the target is later constrained back to a Nature
Portfolio Article format (approximately 5,000 words): **Communications Earth
& Environment**. Retain the original 5,000-word section allocation (archived
below) as the basis for that condensed variant; do not use it as the working
budget while the target is Water Resources Research or Journal of Hydrology.

Stretch target after substantially stronger field and benchmark validation: **Nature Water**.

## Working Title

**Identifiability-aware reaction inference makes groundwater mechanisms falsifiable**

Alternative titles:

1. **When accurate geochemical fits recover the wrong groundwater reactions**
2. **Groundwater reactions are often identifiable only as stoichiometric equivalence classes**
3. **Predictive falsification exposes non-unique groundwater reaction pathways**

## Role Within the Hydrosheaf PhD

M5 is not a standalone sparse-inversion paper and must not be written as a second framework paper. It answers the third question in the Hydrosheaf inference chain:

1. **M4 - Where can groundwater move?** Candidate directed edges and their evidential support.
2. **M3 - How long has groundwater moved?** Residence-time distributions and age-ordering constraints.
3. **M5 - What happened chemically along a supported edge?** Reaction classes, extents, alternatives, and identifiability.
4. **M6 - How stable is the resulting interpretation?** Regularisation sensitivity, structural uncertainty, and cascading error.

M1 establishes why these questions have historically been treated separately. M2 is the umbrella paper that integrates their outputs into the Hydrosheaf framework and field demonstration.

The M5 inference unit is therefore not an arbitrary sample pair. It is a directed Hydrosheaf edge carrying, where available:

- topology confidence or posterior edge probability from the graph/sheaf layer;
- residence-time estimate and uncertainty from the temporal/tracer layer;
- transport-model parameters;
- sparse reaction extents and alternative supports;
- PHREEQC feasibility diagnostics;
- uncertainty and quality-control flags.

M5 should validate the chemical interpretation placed on this shared edge object. It should not claim to validate edge existence, groundwater age, or the complete Hydrosheaf system.

## Thesis-Level Dependency Contract

**Inputs to M5**

- From M4/code topology layer: directed candidate edges, edge confidence, sheaf consistency, and topology uncertainty.
- From M3/code temporal layer: residence-time estimates or distributions that can constrain time-dependent reaction plausibility.
- From M2 shared design: transport correction, reaction dictionary, field data contract, and common `EdgeResult` reporting.

**Outputs from M5**

- Exact-phase and reaction-equivalence-class support rather than one unqualified mineral set.
- A calibrated Mechanism Resolution Score and identifiability class.
- Reaction-extent estimates with uncertainty and identifiability class.
- PHREEQC feasibility and forward-consistency diagnostics.
- Held-out-ion predictive performance and next-best-measurement priorities.
- Edge-level chemical evidence that M2 can integrate into the process network.
- Failure modes and stability measures that M6 can propagate through the full framework.

**Claim boundary**

M5 validates whether a chemical mechanism can be defensibly assigned to a candidate groundwater connection. It does not prove that the connection is physically real; that is M4's domain. It does not prove the residence-time model; that is M3's domain. It does not establish whole-framework robustness; that is M6's domain.

## Central Scientific Claim

A low concentration residual, sparse support, and thermodynamic feasibility are
insufficient evidence that a groundwater transformation has been assigned the
correct reaction pathway. The defensible inferential target may be an exact
mechanism, a stoichiometric equivalence class, a partially resolved mechanism,
or a non-identifiable pathway. Mechanistic confidence therefore requires
ground-truth recovery tests, null-space diagnostics, regularisation stability,
held-out-ion prediction, PHREEQC forward consistency, and, where available,
conditioning on topology and residence-time evidence.

This claim makes M5 essential to Hydrosheaf rather than adjacent to it. It addresses the framework's chemical trust problem: even when a directed edge is topologically plausible and temporally consistent, equifinality can produce an excellent mass-balance fit while selecting an incorrect or non-unique mineral assemblage.

## Research Questions

1. Under what chemical and stoichiometric conditions can sparse linear inversion recover the true mineral phases and reaction extents used to generate groundwater evolution?
2. How do analytical noise, missing ions, phase collinearity, and the L1 penalty alter reaction support and extent estimates?
3. Do PHREEQC-derived thermodynamic constraints improve mechanistic recovery, or do they only remove physically impossible solutions?
4. Can a calibrated Mechanism Resolution Score distinguish identifiable pathways from ambiguous but numerically well-fitting solutions?
5. Which additional ion or tracer provides the greatest expected reduction in reaction ambiguity, and how do topology and residence-time uncertainty alter edge-level chemical confidence?
6. Does identifiability-aware reporting produce more defensible and reproducible interpretations than a single best-fit model in the Northern Ghana chemistry-only field demonstration?

## Testable Hypotheses

- H1: Concentration-fit accuracy is only weakly related to phase-support accuracy in stoichiometrically ambiguous systems.
- H2: Thermodynamic bounds reduce impossible phase directions but do not eliminate equifinality among feasible reactions.
- H3: Missing diagnostic ions cause abrupt changes in selected reactions even when the observed-ion residual remains small.
- H4: A Mechanism Resolution Score combining rank, nullity, coherence, support stability, held-out-ion error, and thermodynamic feasibility predicts phase-recovery reliability better than residual error alone.
- H5: Retrospective next-best-measurement selection identifies ions that reduce reaction-equivalence-class ambiguity more efficiently than fixed or random analytical panels.
- H6: Identifiability-aware reporting yields more conservative and reproducible field interpretations, while chemical confidence decreases when reaction support is conditioned on uncertain topology or residence time.

## Article Architecture and Word Budget

The target is **6,500 words from Introduction through Conclusion**, excluding the Abstract, references, figure legends, code/data availability statements, and Supplementary Information. Word counts are prose only and exclude headings, equations, and citation markers, consistent with `methods_manifest.json` policy.

| Section | Target words (6,500-word budget) | Archived 5,000-word allocation |
|---|---:|---:|
| 1. Introduction | 1,050 | 800 |
| 2. Methods | 1,550 | 1,200 |
| 3. Results | 2,100 | 1,600 |
| 4. Discussion | 1,550 | 1,200 |
| 5. Conclusion | 250 | 200 |
| **Total** | **6,500** | **5,000** |

The extra 1,500 words are allocated mainly to Results (+500) so that the six
main figures and Table 1 each receive a fully reported quantitative
paragraph, and to Methods and Discussion (+350 each) so the identifiability
diagnostics, the Mechanism Resolution Score, and the equifinality/ELRI
distinction can be explained without compressing them into single sentences.
Introduction and Conclusion grow only in proportion (+250/+50) because their
argumentative structure (six subsections; direct answers to six research
questions) does not require additional length to stay defensible.

Abstract: 180-200 words, unstructured, no references, word-capped
independently of the main-text journal choice.

## Abstract - 180-200 words, excluded from the 6,500-word total

Write only after all analyses are locked.

Required sequence:

1. State the problem: inverse geochemical models can fit concentration changes without uniquely identifying the responsible reactions.
2. State the methodological advance: an identifiability-aware sparse inversion and validation framework integrating reaction-equivalence classes, a Mechanism Resolution Score, held-out-ion falsification, next-best measurement, thermodynamic bounds, and PHREEQC-generated ground truth.
3. State the benchmark scale: number of synthetic scenarios, hydrochemical archetypes, mineral assemblages, noise levels, missing-ion patterns, and the Northern Ghana field demonstration.
4. Report the principal quantitative results: exact-phase and equivalence-class precision/recall, extent error, false discoveries, held-out prediction, resolution-score calibration, and effect of thermodynamic screening.
5. State the significance: the framework prevents overinterpretation of plausible but non-identifiable groundwater reaction pathways.

Avoid claiming a fully coupled nonlinear PHREEQC inverse solver.

# 1. Introduction - 1,050 words

The Introduction must follow the six subsections below in the stated order.
Each subsection should lead logically to the next and avoid presenting study
results before the Methods and Results sections.

## 1.1 General studies: groundwater reactions and water-resource interpretation - 170 words

- Establish the importance of identifying mineral dissolution, precipitation, ion exchange, redox transformation, and anthropogenic inputs for groundwater-resource assessment.
- Explain how hydrochemical evolution records interactions among recharge, lithology, mixing, residence time, climate, and human pressures.
- Show why reaction-pathway interpretation informs contaminant attribution, aquifer vulnerability, treatment decisions, and long-term water-quality prediction.
- Conclude by establishing inverse hydrogeochemistry as the general scientific approach for reconstructing otherwise unobserved groundwater transformations.

## 1.2 Current studies: contemporary inverse hydrogeochemical approaches - 195 words

- Summarise current use of PHREEQC inverse modelling, saturation indices, mass-balance models, endmember mixing analysis, isotope constraints, sparse regression, Bayesian inference, and groundwater-network approaches.
- Explain the value of thermodynamic screening, regularisation, uncertainty analysis, and multi-tracer evidence.
- Note that recent studies increasingly combine hydrochemistry with spatial, hydraulic, isotopic, or age information, but usually report connectivity, travel time, and reactions as separate outputs.
- Identify the current movement towards reproducible workflows, ensemble interpretation, and data-guided monitoring design.
- End by showing that current methods improve plausibility and parsimony but do not automatically establish mechanistic uniqueness.

## 1.3 Problem statement: numerical fit is not mechanistic identification - 210 words

- State the central problem explicitly: multiple stoichiometrically equivalent or near-equivalent reaction combinations can reproduce the same concentration difference.
- Explain why a low residual, sparse support, or thermodynamically feasible solution may still identify the wrong reaction or only one member of an unresolved equivalence class.
- Identify incomplete ion panels, analytical noise, uncertain phase dictionaries, transport-model error, mixing, and correlated stoichiometric vectors as sources of non-identifiability.
- Explain that conventional missing-ion tests can appear favourable when omitted ions are simply removed from the objective rather than predicted.
- State that the absence of ground-truth phase recovery, calibrated resolution diagnostics, and predictive falsification permits overconfident hydrogeochemical narratives.

## 1.4 Novelty: identifiability-aware and predictively falsifiable reaction inference - 180 words

- Introduce the study as a shift from selecting one plausible reaction pathway to determining what level of mechanism is observationally resolvable.
- Define the principal innovations:
  1. reaction-equivalence classes derived from rank, singular values, null-space structure, and mutual coherence;
  2. a live-PHREEQC factorial benchmark with known phase and extent truth;
  3. a calibrated Mechanism Resolution Score;
  4. held-out-ion predictive falsification rather than retained-ion residual comparison;
  5. retrospective next-best-measurement selection based on expected ambiguity reduction;
  6. edge-conditioned chemical confidence within Hydrosheaf.
- Present M5 as the chemical claim-auditing layer that distinguishes an identifiable mechanism, a partially resolved reaction class, and a non-identifiable pathway.

## 1.5 Aim and objectives - 180 words

State the overall aim:

> To develop and validate an identifiability-aware framework that determines
> when groundwater reaction mechanisms are uniquely recoverable, recoverable
> only as equivalence classes, or unsupported by the available observations.

State five objectives:

1. Construct a diverse live-PHREEQC benchmark with known reaction phases, directions, and extents.
2. Quantify the effects of stoichiometric structure, noise, missing ions, transport confounding, and regularisation on exact-phase and equivalence-class recovery.
3. Develop and independently test a Mechanism Resolution Score against held-out synthetic archetypes.
4. Evaluate held-out-ion prediction and next-best-measurement selection as falsification and monitoring-design tools.
5. Demonstrate chemistry-only identifiability-aware reporting using 320 wet/dry records from 160 Northern Ghana boreholes sampled in March and August 2025.

## 1.6 Significance - 115 words

- State the scientific significance: the framework replaces unqualified pathway selection with a measurable resolution hierarchy.
- State the methodological significance: it separates numerical reconstruction, thermodynamic feasibility, and mechanistic identification.
- State the practical significance: it identifies which additional ion or tracer would most improve a groundwater interpretation.
- State the Hydrosheaf significance: it supplies chemically qualified evidence to the shared edge object without claiming independent validation of topology or groundwater age.
- End with the broader value for reproducible and decision-relevant groundwater-process inference.

# 2. Methods - 1,550 words

## 2.1 Conceptual model and mass-balance formulation - 235 words

Define the directed edge as the basic inferential unit. Report its topology evidence and residence-time information before introducing its hydrochemical transformation. Where those upstream quantities are unavailable, label the analysis explicitly as chemistry-only rather than fully integrated Hydrosheaf inference.

Define:

$$
\mathbf{r} = \mathbf{x}_v - \mathbf{x}_{\mathrm{transport}},
$$

and

$$
\hat{\mathbf{z}}
=
\arg\min_{\mathbf{z}}
\left[
\frac{1}{2}
\|\mathbf{W}^{1/2}(\mathbf{r}-\mathbf{S}^{\mathsf T}\mathbf{z})\|_2^2
+\lambda_1\sum_j p_j|z_j|
+\frac{\lambda_2}{2}\|\mathbf{z}\|_2^2
\right],
$$

subject to thermodynamic lower and upper bounds.

Define every symbol, unit, sign convention, weighting rule, and the distinction between dissolution and precipitation. Lock one objective scaling convention and make the coordinate-descent threshold mathematically consistent with it.

## 2.2 Reaction dictionary and solver - 220 words

- Describe construction of the stoichiometric reaction matrix.
- Explain mineral, exchange, nitrate-source, and denitrification vectors.
- Describe coordinate descent, soft thresholding, convergence tolerance, bounds, and penalty scaling.
- State how signed reactions are handled.
- Explain how rank, nullity, singular values, condition number, mutual coherence, and null-space structure diagnose identifiability.
- Define exact and approximate reaction-equivalence classes and explain how class-level recovery differs from exact-phase recovery.
- Report software versions and deterministic settings.

## 2.3 PHREEQC synthetic benchmark design - 340 words

Generate a factorial or space-filling ensemble of PHREEQC scenarios with known ground truth.

Required dimensions:

- aquifer chemistry archetypes: carbonate, crystalline/silicate, evaporitic, and mixed;
- pH, temperature, ionic strength, and initial saturation state;
- one to five simultaneous phases;
- dissolution and precipitation directions;
- reaction extents spanning realistic groundwater ranges;
- mixing and evaporation as transport confounders;
- multiple analytical-noise levels;
- censored and missing-ion patterns;
- alternative thermodynamic databases where feasible.

Use at least several hundred scenarios and repeated random seeds. Preserve PHREEQC input files, true extents, endpoint compositions, saturation indices, and metadata.

## 2.4 Comparator models and ablation design - 235 words

Compare:

- non-negative or bounded least squares;
- LASSO;
- elastic net;
- Hydrosheaf with and without PHREEQC bounds;
- a conventional PHREEQC inverse or combinatorial inverse baseline where technically comparable.

Ablate:

- thermodynamic bounds;
- geologic penalty priors;
- indicator-ion gates;
- L1 and L2 regularisation;
- individual ions and ion groups.

Select regularisation using cross-validation, a one-standard-error rule, stability selection, or a predeclared Pareto criterion. Do not select the minimum-residual model alone.

## 2.5 Evaluation metrics and statistical analysis - 285 words

Primary metrics:

- phase-support precision, recall, F1 score, and false-discovery rate;
- reaction-equivalence-class precision, recall, and coverage;
- signed reaction-direction accuracy;
- reaction-extent RMSE, MAE, relative bias, and interval coverage;
- concentration reconstruction RMSE;
- held-out-ion RMSE and predictive interval coverage;
- thermodynamic-violation rate;
- support stability across bootstrap resamples and penalty values;
- computational time and convergence rate.

Analyse performance using effect sizes and confidence intervals across scenario classes. Evaluate whether residual error predicts support accuracy using regression, calibration curves, and rank correlation. Develop the Mechanism Resolution Score from predeclared diagnostics and calibrate it only on training scenarios. Test calibration, discrimination, and threshold transfer on untouched hydrochemical archetypes. Define identifiable, partially identifiable, equivalence-class identifiable, and non-identifiable reporting classes. Evaluate next-best measurement by the realised reduction in ambiguity and held-out prediction error relative to fixed and random measurement panels.

## 2.6 Northern Ghana chemistry-only field demonstration - 155 words

- Use 320 wet/dry hydrochemical records from 160 boreholes sampled during March and August 2025.
- Describe the four aquifer classes, 12 lithologies, major-ion chemistry, water isotopes, Sr, SiO2, saturation indices, and charge-balance quality classes.
- Use the 294 quantitative inverse-modelling records for primary field analyses and retain the remaining records for sensitivity reporting.
- Stratify support stability and held-out-ion performance by aquifer, lithology, and season.
- Report reaction-equivalence classes and alternative pathways when the measurements do not distinguish individual phases.
- Compare the identifiability-aware result with a conventional single best-fit interpretation.
- State explicitly that generated graph edges, heuristic residence times, and field reaction labels are not independent ground truth; report this component as chemistry-only demonstration and transfer assessment.

## 2.7 Reproducibility and reporting - 80 words

- Archive code, environment lock file, PHREEQC database, generated inputs, outputs, seeds, and figure scripts.
- Provide a one-command workflow.
- Deposit benchmark data and code in a permanent repository.
- Include data, code, and model-card style claim limitations.

# 3. Results - 2,100 words

## 3.1 Benchmark overview and numerical validity - 235 words

- Report scenario counts, convergence, rank/conditioning distribution, and mass-balance checks.
- Verify exact recovery in well-conditioned noise-free cases before presenting stress tests.
- Demonstrate agreement between generated ground truth and extracted PHREEQC endpoint chemistry.

Proposed display: **Figure 1**, workflow and benchmark design.

## 3.2 Good concentration fits can conceal incorrect pathways - 365 words

- Present the central result relating concentration RMSE to phase F1 and extent error.
- Show cases with similarly low residuals but different reaction supports.
- Quantify how often a low residual produces incomplete or false phase recovery.
- Characterise exact and approximate reaction-equivalence classes and null-space alternatives.
- Compare exact-phase recovery with equivalence-class recovery to distinguish solver failure from information-limited ambiguity.

Proposed display: **Figure 2**, residual accuracy versus phase-support accuracy, with representative ambiguous pathways.

## 3.3 Effects of sparsity and stoichiometric identifiability - 325 words

- Show regularisation paths for representative well-conditioned and ambiguous systems.
- Report changes in support, bias, and false discoveries across $\lambda_1$ and $\lambda_2$.
- Relate performance to rank, condition number, and mutual coherence.
- Identify where elastic net improves stability and where it cannot resolve non-identifiability.
- Report Mechanism Resolution Score calibration, discrimination, and performance on held-out aquifer archetypes.

Proposed display: **Figure 3**, regularisation paths, null-space structure, and Mechanism Resolution Score calibration.

## 3.4 Missing ions and minimum analytical information - 340 words

- Report ion-ablation effects on phase recovery and extent estimates.
- Identify reaction-specific diagnostic ions.
- Distinguish low residual on the retained ions from prediction error on withheld ions.
- Propose a minimum analytical panel for defensible recovery of the tested reaction classes.
- Compare expected-information-gain rankings with the realised value of each hidden measurement.
- Report whether next-best measurement outperforms fixed and random analytical panels.

Proposed display: **Figure 4**, ion-ablation and next-best-measurement map showing class recovery, held-out prediction error, and ambiguity reduction.

## 3.5 Value and limits of thermodynamic screening - 340 words

- Compare unconstrained and PHREEQC-bounded fits.
- Quantify reduction in impossible reaction directions and false positives.
- Test whether thermodynamic bounds improve true support recovery.
- Show feasible but non-unique solutions that remain after screening.

Proposed display: **Figure 5**, constrained-versus-unconstrained performance and thermodynamic feasibility.

## 3.6 Northern Ghana chemistry-only demonstration - 310 words

- Present hydrochemical context across aquifer classes, lithologies, and wet/dry seasons.
- Report the best-supported reaction-equivalence classes and their stability.
- Report held-out-ion prediction and retrospective next-best-measurement performance.
- Show which conclusions are robust and which remain ambiguous.
- Compare conventional minimum-residual interpretation with identifiability-aware reporting.
- Separate chemistry-only conclusions from generated topology and unvalidated residence-time outputs.

Proposed display: **Figure 6**, Northern Ghana reaction-class stability, held-out prediction, and alternative plausible pathways.

## 3.7 Computational performance - 185 words

- Report runtime, convergence, and scaling with ions, reactions, and scenario count.
- Demonstrate feasibility for screening large groundwater networks.

Proposed display: **Table 1**, benchmark performance summary across methods and stress conditions.

# 4. Discussion - 1,550 words

## 4.1 Principal finding - 285 words

- State directly that numerical fit, sparsity, and thermodynamic feasibility answer different questions.
- Explain why a plausible sparse solution is not automatically a uniquely supported mechanism.
- Tie every statement to quantitative results.

## 4.2 Advance over existing inverse hydrogeochemical practice - 310 words

- Contrast the framework with conventional PHREEQC inverse applications that identify plausible reaction assemblages.
- Explain that the contribution is not replacement of PHREEQC, but an independent identifiability and stress-testing layer.
- Contrast single-support LASSO reporting with structural-uncertainty-aware reaction inference.
- Explain the conceptual advance from exact-phase selection to hierarchical reporting of exact mechanisms, equivalence classes, partial resolution, and non-identifiability.
- Discuss why held-out-ion prediction provides a stronger falsification test than goodness of fit on the ions used for inversion.

## 4.3 Scientific and monitoring implications - 325 words

- Explain how ion selection controls inferential credibility.
- Discuss implications for designing groundwater monitoring programmes.
- Show how the framework prevents unsupported attribution of salinisation, carbonate dissolution, sulphate sources, denitrification, or ion exchange.
- Recommend reporting robust reaction classes and alternative pathways rather than one deterministic story.
- Explain how next-best-measurement rankings translate uncertainty diagnostics into an actionable sampling decision.

## 4.4 Interpretation of the field application - 235 words

- Distinguish robust field conclusions from unresolved alternatives.
- Explain how geological evidence and thermodynamics narrow, but may not eliminate, ambiguity.
- Compare wet/dry and aquifer-stratified reaction-class stability without interpreting seasonal association as direct flow-path evolution.
- Interpret held-out-ion performance as predictive consistency rather than proof of reaction truth.
- Avoid claiming validation where no independent ground truth exists.
- State that Northern Ghana demonstrates field transfer of the chemistry audit, not independent validation of connectivity, groundwater age, or reaction mechanism.

## 4.5 Limitations - 245 words

Address explicitly:

- fixed linear stoichiometric dictionary;
- equilibrium-based thermodynamic screening;
- absence of fully coupled nonlinear kinetic inversion;
- uncertainty in thermodynamic databases and phase selection;
- possible transport-model misspecification;
- synthetic-to-field transferability;
- dependence on measured-ion coverage;
- absence of matched age tracers and independent field flow paths;
- field interpretation without direct reaction-rate or mineralogical-abundance observations.

## 4.6 Future development - 150 words

- stability-selected or Bayesian structural ensembles;
- kinetic and affinity-dependent bounds;
- uncertainty propagation from analytical measurements and thermodynamic databases;
- prospective field testing of the next-best-measurement recommendations;
- matched chemistry, age-tracer, hydraulic-head, screen-depth, and mineralogical observations;
- multi-edge groundwater-network applications.

# 5. Conclusion - 250 words

- Answer each research question directly.
- State the practical message: reliable inverse hydrogeochemistry requires evidence for identifiability, not only a low residual.
- Summarise the role of PHREEQC screening as necessary but insufficient.
- Present Hydrosheaf as a sparse linear inverse reaction model with transparent diagnostics and PHREEQC screening/forward checks.
- End with the broader contribution to reproducible groundwater-process inference.
- Introduce no new evidence or citations.

# Display-Item Plan

## Main Figures

1. Identifiability-aware inverse hydrogeochemical workflow.
2. Concentration fit versus exact-phase and equivalence-class recovery.
3. Regularisation paths, null-space structure, and Mechanism Resolution Score.
4. Missing-ion falsification and next-best-measurement value.
5. Effect and limitations of PHREEQC thermodynamic constraints.
6. Northern Ghana chemistry-only demonstration with alternative plausible pathways.

## Main Table

1. Comparative performance of inverse methods across benchmark regimes.

## Supplementary Figures

- S1. Full reaction dictionary and stoichiometric correlations.
- S2. Singular-value spectra, null-space bases, and reaction-equivalence classes.
- S3. Scenario distributions and PHREEQC quality control.
- S4. Noise sensitivity by aquifer archetype.
- S5. Extent bias for each phase.
- S6. Full regularisation paths.
- S7. Mechanism Resolution Score calibration and held-out-archetype testing.
- S8. Complete ion-ablation and held-out-ion results.
- S9. Next-best-measurement rankings and realised information gain.
- S10. Thermodynamic-database sensitivity.
- S11. Bootstrap support frequencies.
- S12. Runtime and convergence diagnostics.
- S13. Northern Ghana hydrochemical quality-control plots.
- S14. Northern Ghana reaction-equivalence-class ensemble.

## Supplementary Tables

- S1. Reaction stoichiometry, sign convention, and diagnostic ions.
- S2. Reaction-equivalence-class definitions and null-space coefficients.
- S3. PHREEQC scenario parameters.
- S4. Hyperparameter grids and selection rules.
- S5. Complete model-comparison metrics.
- S6. Confusion matrices by exact phase and equivalence class.
- S7. Missing-ion, held-out prediction, and measurement-value performance.
- S8. Thermodynamic-bound definitions.
- S9. Software and computational environment.
- S10. Northern Ghana data, 2025 sampling design, and quality-control summary.

# Appendix or Supplementary Methods

## Supplementary Methods word budget

Neither Water Resources Research nor Journal of Hydrology imposes a strict
word limit on Supporting Information, so the Supplementary Methods length is
set by content need rather than a journal cap: it must independently
reproduce every procedure that Methods (1,550 words) can only reference. The
eight subsections below total approximately **4,250 words**, sized so that
each subsection can carry a full derivation, algorithm, or protocol rather
than a compressed summary. Allocate no fewer than 300 words to any
subsection; subsections 3-6 carry the load-bearing technical content
(rank/null-space diagnostics, the Mechanism Resolution Score, held-out-ion
falsification, and the PHREEQC generator) and are weighted accordingly.

| Subsection | Target words |
|---|---:|
| 1. Derivation of the coordinate-descent update | 450 |
| 2. Bound handling and sign conventions | 350 |
| 3. Rank, null-space, singular-value, condition-number, mutual-coherence, and reaction-equivalence-class diagnostics | 700 |
| 4. Stability-selection and Mechanism Resolution Score algorithms | 750 |
| 5. Held-out-ion falsification and next-best-measurement algorithm | 650 |
| 6. PHREEQC synthetic-data generator and species extraction | 650 |
| 7. Northern Ghana preprocessing, unit conversion, and 2025 sampling-date verification | 400 |
| 8. Complete reproducibility protocol | 300 |
| **Total** | **4,250** |

1. Derivation of the coordinate-descent update.
2. Proof or explanation of bound handling and sign conventions.
3. Rank, null-space, singular-value, condition-number, mutual-coherence, and reaction-equivalence-class diagnostics.
4. Stability-selection and Mechanism Resolution Score algorithms.
5. Held-out-ion falsification and next-best-measurement algorithm.
6. PHREEQC synthetic-data generator and species extraction.
7. Northern Ghana preprocessing, unit conversion, and 2025 sampling-date verification.
8. Complete reproducibility protocol.

# Data and Analysis Status After M5 Implementation

The executable M5 package now contains the live-PHREEQC benchmark,
identifiability diagnostics, held-out prediction, measurement-value analysis,
calibrated Mechanism Resolution Score, Northern Ghana transfer analysis, and
all planned main and supplementary displays.

## Implemented evidence

- 240 live-PHREEQC scenarios spanning carbonate, crystalline, evaporitic, and
  mixed hydrochemical archetypes, with known reaction supports, directions,
  and extents.
- 18,000 factorial inverse fits across five comparator methods, three
  analytical-noise levels, five ion panels, and transport-confounding levels.
- A 16-reaction, 11-ion dictionary with exact signed equivalence classes,
  panel-specific rank, nullity, condition, coherence, and singular-value
  diagnostics.
- Exact-phase, equivalence-class, direction, extent, reconstruction,
  held-out-ion, thermodynamic, stability, convergence, and runtime metrics.
- Mechanism Resolution Score calibration on carbonate, crystalline, and
  evaporitic scenarios and untouched transfer testing on the mixed archetype.
- A conventional PHREEQC `INVERSE_MODELING` baseline with strict 5%
  uncertainty and documented relaxed 20% fallback, preserving all scenario
  input/output files and parsed model multiplicity.
- Retrospective next-best-measurement and bootstrap support-stability analyses.
- 160 wet-to-dry Northern Ghana pairs (all sampled boreholes) representing
  four administrative regions; no independent aquifer-type or lithology
  classification is available for these boreholes. All source sampling dates
  are verified as March and August 2025.
- Six data-derived main figures, 14 supplementary figures, one main table, and
  11 supplementary tables in PNG/PDF and CSV formats.

## Principal locked results

- At 3% analytical noise with the full 11-ion panel, the training-calibrated
  Hydrosheaf guarded model achieved mean exact-phase F1 of 0.563 and mean
  equivalence-class F1 of 0.609 while reducing false-discovery rate from
  0.465 for the legacy thermodynamic elastic-net model to 0.361.
- The conventional PHREEQC inverse baseline was feasible for 45.8% of
  scenarios under strict 5% uncertainty and 99.6% after relaxed 20% fallback,
  but its first-minimal equivalence-class F1 was 0.512 and it generated a mean
  185.8 feasible inverse models per scenario.
- Among the lowest concentration-residual quartile, 55.0% of guarded fits still had
  exact-phase F1 below 0.80, directly supporting the paper's central claim.
- Thermodynamic bounds eliminated incompatible fitted directions but 72% of
  constrained fits remained below equivalence-class F1 of 0.80.
- The calibrated MRS achieved 48.9% four-class accuracy on the untouched mixed
  archetype. This is useful but moderate discrimination and must be reported
  without inflation.
- Fluoride had the highest retrospective measurement-value score in the
  Northern Ghana transfer analysis.
- All 160 Northern Ghana wet-to-dry pairs were classified as partially
  identifiable by the transferred calibration, supporting conservative
  class-level rather than unique-phase reporting.

## Remaining work before final Nature-Portfolio submission

1. Add alternative thermodynamic-database sensitivity if the necessary
   databases and phase mappings can be locked.
2. Obtain prospective or independent field evidence, preferably matched
   mineralogy, age tracers, hydraulic heads, screen depths, and repeat
   chemistry.
3. Externally validate or refine the MRS because held-out four-class accuracy
   is currently moderate.
4. Lock the public Northern Ghana dataset citation and permanent code/data DOI.

## Field evidence boundary

Northern Ghana supports field transfer of chemistry-only reaction-class
stability, held-out-ion prediction, and measurement-value analysis. It does
not independently validate groundwater connectivity, residence time, or
reaction truth because matched age tracers, independent flow paths, and direct
mineral-reaction observations remain unavailable.

# Claim Guardrails

Use:

> Hydrosheaf is an identifiability-aware sparse linear inverse reaction model that uses PHREEQC-derived thermodynamic screening and independent forward-validation diagnostics.

Do not use:

> Hydrosheaf is a fully coupled nonlinear PHREEQC inverse or reactive-transport solver.

Do not equate:

- low residual with correct mechanism;
- sparse support with unique support;
- thermodynamic feasibility with pathway validation;
- a field-consistent interpretation with ground-truth recovery.
