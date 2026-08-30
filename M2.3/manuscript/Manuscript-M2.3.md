# HydroSheaf: a claim-gated evidence-integration framework for groundwater connectivity, residence-time and reaction inference in data-limited aquifers

# Abstract

Groundwater interpretation in data-limited aquifers combines hydraulic, tracer and hydrochemical evidence that each admit several compatible explanations, so an integrated framework can become more confident without becoming more correct. HydroSheaf v0.5.1 was implemented as a modular Python framework that retains competing explanations until an explicit claim gate, calibrates intervals on frozen development cases, estimates model discrepancy before averaging, and returns PASS, WEAK or ABSTAIN rather than an unconditional answer. Its central result is that inference targets divide by identifiability, and that the division can be demonstrated rather than assumed. Point reaction extents are structurally non-identifiable from major-ion chemistry: the stoichiometry matrix of the declared reaction dictionary has rank 8 against 14 unknown terms, leaving a six-dimensional null space of extent combinations that produce identical chemistry. Empirically the inversion carried no explanatory power for extent magnitude and assigned non-trivial extents to 54.1% of truth-zero terms while detecting only 42.1% of active ones. Recast as calibrated reaction-family probabilities with abstention, the same evidence passed a predeclared gate with a multiclass log loss of 0.896 against 2.852 for a competence-matched baseline and a false commitment rate of 0.038. Kinetic prediction passed while parameter identification abstained in 83.3% of cases under rate-constant and surface-area collinearity. Network-constrained residence-time inference raised log-space R² from 0.922 to 0.983 relative to single-node inversion while halving interval width. Held-out uncalibrated agreement with published lumped-parameter ages was screening-level (R² 0.482), and fitting a correction raised R² to 0.962 while measuring emulation rather than agreement. Application to 264 samples from northern Ghana returned an evidence ceiling rather than a validation, and identified strontium and silica as reducing non-identifiable wells from 51.9% to 0.6%.

**Keywords:** groundwater inverse modelling; equifinality; uncertainty calibration; abstention; hydrogeochemistry; reproducible research software

# 1. Introduction

Groundwater interpretation is an inverse problem in which the quantity of interest is never observed directly. Hydraulic heads constrain potential gradients, environmental tracers constrain mixtures of recharge histories, and dissolved chemistry records the accumulated effect of mixing, mineral dissolution and precipitation, ion exchange and redox transformation. In well-instrumented settings these streams are combined within a calibrated numerical model. In the data-limited aquifers that supply much of rural sub-Saharan Africa they are not: boreholes are sampled once, screen intervals are unrecorded, heads are not repeated, and environmental age tracers are absent altogether [@macdonald2012africa; @banoengyakubo2011ghana]. The analyst is left to reason from a single hydrochemical snapshot.

The difficulty is not primarily the absence of a sophisticated solver. It is that several structurally different explanations remain compatible with the same limited observations. Groundwater age distributions inferred from environmental tracers are model-dependent, and mean age can be non-unique when tracer information is sparse or when alternative lumped-parameter models reproduce the same concentrations [@jurgens2012tracerlpm; @mccallum2015limitations]. Additional tracers tighten the constraint but leave the estimate conditional on assumed recharge histories, mixing structure and the treatment of measurement and model error [@visser2013multitracer]. Inverse hydrogeochemical calculation has the same character: PHREEQC and comparable tools identify reaction sets that close an observed mass balance, but a closing set is not thereby a unique reaction history [@parkhurst2013phreeqc]. Environmental modelling has described this equifinality for three decades [@beven2001glue; @beven2006manifesto], and the accompanying epistemological point is older still — confirmation of a natural-system model from a finite set of observations is necessarily partial [@oreskes1994validation].

Equifinality has a specific and under-appreciated consequence for integrated frameworks. Combining evidence streams that each admit multiple explanations does not automatically narrow the explanation set; it can instead produce a single confident answer whose confidence is an artefact of the combination rule. A framework can become more precise without becoming more correct. Model discrepancy, if unrepresented, biases parameter estimates and understates their uncertainty [@kennedy2001calibration; @brynjarsdottir2014discrepancy]. Correlated evidence counted twice narrows an interval without adding information. A candidate flow path chosen by proximity and then confirmed by a chemistry fit against that same path is circular.

Existing tools address parts of this problem well and were not designed to address the rest. MODFLOW 6 solves groundwater flow on a discretised domain and MODPATH 7 traces particles through the resulting field [@langevin2017modflow6; @pollock2016modpath]; both are conditional on a conceptual model, boundary conditions and a parameterisation that data-limited settings rarely support. TracerLPM fits lumped age-distribution models to environmental tracers [@jurgens2012tracerlpm], which the datasets considered here do not contain. PHREEQC performs speciation, reaction-path and inverse calculations for a nominated pair of waters along an assumed path [@parkhurst2013phreeqc]. Graph-based representations of hydrogeological systems have been proposed as a way to make the path assumption itself testable [@moracchini2025graphflow], and sheaf-theoretic language has been used to describe the consistency of overlapping observations on a network [@robinson2017sheaves; @hansen2019spectral]. What has been missing is an implementation in which the resulting claim is gated: where the software states, from evidence rather than from the analyst's judgement, which of its own outputs are supported.

This paper documents such an implementation. HydroSheaf represents groundwater observations as nodes and candidate directed edges, attaches hydraulic, isotopic, chemical, thermodynamic and kinetic evidence to those candidates, and — this is the central design commitment — retains the candidate set rather than collapsing it to a single path before the evidence has been evaluated. Calibration is fitted on development cases and frozen; disagreement between components is estimated and reported before any averaging; predictions that are not supported by the observation record are recorded as abstentions rather than imputed. A programme-level contract then returns PASS when every predeclared criterion is met, WEAK when a finite result falls below a performance floor, and ABSTAIN when a result is withheld as unsupported, non-identifiable or not comparable.

The contribution is therefore an executable evidence-and-claim architecture rather than a new forward model, and the paper's evidence is organised to test that architecture rather than to advertise it. The organising question is:

> Which groundwater inference targets are recoverable from incomplete hydraulic, tracer and hydrochemical evidence, and can a framework determine that boundary from the evidence itself rather than by assumption?

Six objectives follow. The implemented package and its evidence-to-claim pipeline are documented (OBJ-1). The three Ghanaian datasets to which the framework is applied are characterised in full, so that a reader can independently judge what they support (OBJ-2). Component recovery is quantified against a 100-realisation synthetic benchmark with known truth, separating recoverable from non-recoverable targets (OBJ-3). Calibrated age, reaction-family, kinetic and integrated decision performance are scored against thresholds fixed before scoring (OBJ-4). Framework output is compared with evidence generated outside the study — published lumped-parameter age estimates and a particle-tracking directed topology (OBJ-5). Finally, the framework is applied to the measured field data and the resulting evidence ceiling is reported (OBJ-6).

The expected contribution is conditional, and deliberately so. A positive result on independent synthetic cases supports a bounded statement about the tested generators and observation mechanisms; it does not establish field effectiveness or superiority over specialist tools. A negative result is equally informative, because it identifies an inference target that the available evidence cannot support. Both outcomes occur in what follows.

# 2. Data

## 2.1 Regional setting

All field measurements came from northern Ghana, within the semi-arid Sudano-Sahelian belt between approximately 8.8 and 11.1 degrees north. Groundwater in this region occurs principally in the weathered regolith and fractured bedrock of the Precambrian Birimian and Voltaian formations, where storage is limited, transmissivity is heterogeneous over short distances, and boreholes are typically the sole dry-season water source [@macdonald2012africa; @banoengyakubo2011ghana]. Recharge is concentrated in a single wet season. These are precisely the conditions under which a calibrated numerical flow model is unavailable and interpretation must proceed from hydrochemistry, and they are therefore the conditions the framework is designed for.

Three datasets were used Table 1, comprising 264 water samples from 264 distinct sampling points distributed across the Northern, North East, Upper East and Upper West regions Fig. 2. Every coordinate was validated by a point-in-polygon test against the national and regional boundaries rather than by inspection of the plotted map; this test located a sign error in one source file, described below. They differ deliberately in completeness and analytical quality, because the purpose of the field application was to test how the framework behaves as the evidence base degrades, not to assemble the best available dataset. Administrative boundaries used for mapping were taken from the geoBoundaries open database and are distributed under CC BY 4.0.

## 2.2 Dataset descriptions

The **Northern Ghana panel** is an author-compiled regional dataset covering 160 boreholes across 31 districts in four administrative regions, with 29 recorded variables. It is the most complete of the three: physical setting (ground elevation 141 to 340 m, borehole depth 20.0 to 99.9 m, static water level 7.8 to 35.0 m below ground), field determinations (pH 6.51 to 7.89, electrical conductivity 427 to 1,192 µS/cm, water temperature 27.0 to 31.5 °C), the eight major ions, fluoride, strontium, dissolved silica, and both stable water isotopes. Nitrate ranges from 0.57 to 160.5 mg/L and fluoride from 0.04 to 3.88 mg/L, so a subset of wells exceeds the World Health Organization guideline values for both. Strontium and silica are carried by this panel alone, which is material because both act as corroborating constraints on silicate-weathering interpretations.

One property of this panel requires explicit statement. The workbook presents its 160 wells as two sheets labelled Dry and Wet, giving 320 rows. The seasonal separation is a reconstruction, not an independent resampling: static water level, ground elevation, borehole depth and distance-to-river are identical in 160 of 160 wells across the two sheets Fig. 3, and a genuine wet-season resampling could not return the same water table at every well. The chemistry is measured; the seasonal attribute is constructed. Accordingly, the Dry sheet is treated throughout as the primary measured panel, the Wet sheet is excluded from every inferential result, the measured unit is the well rather than the seasonal sample, and no seasonal-change or repeated-head result is reported anywhere in this paper. All Northern Ghana statistics below and in the Results refer to the 160-well primary panel.

**Lower Anayari** contains 41 samples with 22 variables from a published catchment study [@zakaria2021anayari]. It carries the major ions, fluoride, iron and both stable isotopes, but no strontium, silica or borehole depth. It is the most dilute of the three (total dissolved solids 115 to 368 mg/L) and the most sodic: sodium reaches 168 mg/L and bicarbonate 512 mg/L. One fluoride value is left-censored below the reporting limit and was retained as censored rather than substituted.

**Talensi** contains 63 samples with 23 variables from a published study of an artisanal mining area [@chegbeleh2020talensi]. One correction was applied. All 63 longitudes are stored as positive values, which places every sample outside Ghana, whereas Talensi District lies near 0.8 degrees west in the Upper East Region. The sign was therefore negated, after which all 63 samples fall inside the Upper East Region polygon. The correction was applied in the analysis code rather than to the source file, and it does not affect any reported result: reflection about the prime meridian preserves pairwise distance exactly, and the candidate edges are generated from inter-sample distances. It is the only dataset carrying redox potential, recorded for 54 of 63 samples and spanning −205.7 to +22.4 mV, which indicates strongly reducing conditions in part of the area. It carries neither fluoride, strontium nor silica. Its pH range extends to 10.6 and its bicarbonate to 434 mg/L, with one reported zero.

## 2.3 Analytical quality

Charge-balance error was computed independently from the reported concentrations rather than taken from any source document. Concentrations were converted to milliequivalents per litre using equivalent weights, and the error was expressed as the difference between summed cations and summed anions as a percentage of their total. Only the ions each dataset actually reports enter its own balance. Samples were then tiered as quantitative (absolute error at most 5%), screening (5 to 10%) or exploratory (above 10%), following conventional practice for hydrochemical data screening [@fritz1994charge].

The three datasets separate sharply Fig. 2. The Northern Ghana panel has a median absolute error of 1.54%, with 151 of 160 wells quantitative, 7 screening and 2 exploratory. Lower Anayari has a median of 3.13%, with 36 quantitative and 5 screening. Talensi has a median absolute error of 29.9%, with no sample meeting the quantitative tier, 5 screening and 58 exploratory.

The Talensi imbalance is systematic rather than random, and is an anion excess: median summed cations are 2.66 meq/L against 4.49 meq/L of anions, with bicarbonate alone contributing 3.38 meq/L. It is not a reporting-unit artefact. Re-interpreting the bicarbonate column as alkalinity expressed as calcium carbonate worsens the median error to −36.0%, and interpreting it as already in milliequivalents worsens it to −97.5%. The most consistent explanation is unmeasured cationic species — the dataset reports iron only at trace level and no aluminium, ammonium or manganese, all plausible in a reducing, mining-affected setting — or an alkalinity determination not matched to the cation suite. Whatever the cause, the consequence for this study is fixed: Talensi supports screening and failure-mode diagnostics, and no quantitative inverse mass-balance claim.

## 2.4 Hydrochemical and isotopic character

Dominant-ion facies were assigned from the milliequivalent composition, with a cation or anion group labelled dominant when it exceeded half of its total and labelled mixed otherwise Fig. 2. The three datasets are compositionally distinct. The Northern Ghana panel is dominated by mixed-cation waters, 64.4% mixed-mixed and 34.4% mixed-bicarbonate, consistent with silicate weathering under limited residence time without a single controlling mineral phase. Lower Anayari is overwhelmingly sodium-bicarbonate (90.2%), the signature of advanced silicate weathering or cation exchange. Talensi is mostly mixed-bicarbonate (65.1%) with a sodium-bicarbonate subset (28.6%) and small calcium-bicarbonate and mixed-chloride groups.

Stable isotopes were interpreted only as evidence of recharge source, evaporation and mixing. All three local regressions fall well below the global meteoric water line slope of 8 [@craig1961meteoric]: 5.33 for Northern Ghana (R² = 0.50), 3.69 for Lower Anayari (R² = 0.75) and 3.09 for Talensi (R² = 0.80). Median deuterium excess is 5.18‰, 8.78‰ and 8.50‰ respectively. Slopes between 3 and 5.5 are diagnostic of evaporative enrichment during recharge, which is expected where rainfall infiltrates through a hot unsaturated zone or via surface ponding. The Lower Anayari and Talensi datasets each contain a small number of strongly enriched samples, reaching δ¹⁸O values of +0.64‰ and +4.48‰, that are consistent with an evaporated surface-water contribution.

These isotope data support statements about recharge and mixing. They do not constrain groundwater age. δ¹⁸O and δ²H are not environmental age tracers, and no interpretation in this paper treats them as such.

## 2.5 What the data do not contain

The framework's field claim is bounded by the variables that are absent, and these are listed explicitly Fig. 3. No dataset contains environmental age tracers of any kind — no tritium, radiocarbon, chlorofluorocarbon, sulfur hexafluoride or noble-gas determination. None records screen-top and screen-bottom intervals, so vertical exchange cannot be resolved and a sample cannot be attributed to a specific depth horizon. None provides a repeated hydraulic-head time series. None supplies an independently verified directed connectivity structure against which an inferred network could be scored. None provides mineralogical, petrographic or isotopic evidence of reaction mechanism against which an inferred reaction could be checked. No laboratory covariance metadata or replicate analyses are available, so measurement uncertainty had to be assigned from declared analytical conventions rather than estimated from the data.

Consequently the field application in this paper can demonstrate ingestion, analytical screening, candidate generation, reaction-family plausibility, measurement-value diagnostics and non-identifiability reporting. It cannot validate groundwater age, directed connectivity or reaction mechanism, and no such validation is claimed.

## 2.6 Reference data generated outside this study

Two external sources were used for comparison. Published lumped-parameter groundwater-age estimates for United States public-supply aquifers provided a reference against which the residence-time component could be scored on data neither generated nor curated by this study, across 20 study units. A published MODFLOW and MODPATH archive for a municipal supply well provided a particle-tracking directed topology of 174 reference edges [@langevin2017modflow6; @pollock2016modpath]. Neither reference is treated as physical truth. A fitted lumped-parameter output is itself a model interpretation of tracer data [@jurgens2012tracerlpm], and a particle path is conditional on the flow model that produced it; both are used as independent points of comparison, which is a weaker and more defensible role.

# 3. Software and methods

## 3.1 Architecture and evidence representation

HydroSheaf version 0.5.1 was treated as the production software surface. It is a modular Python package exposing an end-to-end network-inference entry point together with a validation namespace containing candidate generators, specialist baselines, calibration, discrepancy estimation, model averaging, kinetic benchmarking, decision utility and claim gates Table 2. The repository was separated analytically into reusable package code, benchmark control scripts, historical analysis records and machine-readable run artefacts; only the first was treated as a package capability.

The runtime represents observations as nodes and hypothesised flow connections as candidate directed edges Fig. 1. A candidate edge is a hypothesis generated from a declared candidate universe, never an observed flow path. Hydraulic direction, tracer compatibility, chemical residuals, thermodynamic limits and kinetic response contribute separate evidence records to that hypothesis. The design commitment that distinguishes the framework is that the candidate set and the compatible outcome set are retained until the claim gate, so that a conflict can be reported as model disagreement, numerical non-convergence, non-identifiability or abstention rather than silently resolved by the combination rule.

Two competence-matched baselines carry the comparator role and are described here because a weak comparator would make any contrast meaningless. The age baseline is an independent lumped-parameter inversion that generates its own candidate universe over the same mean-age, distribution-shape and mixing grid, receives the same observed tracers, uncertainty conventions and missingness pattern, and returns an interval or an abstention under the same acceptance rule. The reaction baseline is a fixed-comparator family classifier that receives the same chemistry and the same family dictionary and returns a probability vector without the on-off evidence layer or the fitted calibration. Neither baseline receives sealed truth, framework-selected candidates or any post-hoc correction, and neither is penalised for declining to infer a quantity absent from its own observation contract. Both are defined in the machine-readable baseline registry accompanying the locked run.

Evidence channels return distributions or compatible sets rather than point values wherever the observation record allows. Missing variables are never imputed as positive evidence. The source of every candidate universe is recorded, so that an apparent performance advantage cannot be attributed to a hidden candidate-set advantage — a control that matters because an integrated method can otherwise appear to improve simply by receiving a larger or better-curated set of options than its comparator.

## 3.2 Component inference

Residence-time inference proceeds in two modes. A single-node lumped-parameter inversion estimates a mean residence time from the tracer record at one location using piston-flow, exponential and dispersion transit-time distributions. A network-constrained Bayesian mode additionally conditions on the candidate graph, requiring that inferred ages be consistent with the directed structure along retained edges. The contrast between the two isolates the contribution of the network constraint itself.

Reaction inference uses two independent layers. A stoichiometric and thermodynamic candidate generator enumerates reaction-family explanations from the observed chemistry, including an explicit null option representing no supported family. A sparse inverse component then fits reaction extents under a bounded objective combining an L1 penalty for parsimony with an L2 term for stability, subject to non-negativity and thermodynamic direction constraints where the declared chemistry supports them [@tibshirani1996lasso; @zou2005elasticnet; @parkhurst2013phreeqc]. Because point extents proved not to be identifiable from this evidence (Section 4.2), a second representation was scored: a regularised adjusted plus-minus (RAPM) layer that treats evidence channels as on-off contributions to a reaction-family score and returns a calibrated probability over families. RAPM is project-specific notation for this scoring layer and is not an established external groundwater method.

The kinetic adapter constructs PHREEQC-compatible kinetic blocks and evaluates a forward response over elapsed time. Predictive concentrations and parameter identification are scored separately, because the rate law depends on the product of the rate constant and the reactive surface area. When only that product is constrained, the two parameters are structurally confounded and additional residence-time observations cannot separate them; the adapter reports a prediction and abstains from the parameter values unless an independent surface-area measurement is supplied.

## 3.3 Calibration, discrepancy and abstention

Interval calibration is fitted on development cases and frozen before locked scoring. Coverage is reported including abstentions in the denominator, alongside the acceptance rate and the selective risk computed over committed outputs only. Reporting these together is deliberate: a method can otherwise obtain an excellent selective score by declining every difficult case, and a narrow interval with poor coverage is not a successful age engine. The distinction between marginal coverage and coverage after selection follows selective and conformal prediction, though the intervals here are calibrated only within the declared programme [@angelopoulos2023conformal].

Where components disagree about the same target, a discrepancy calibrator estimates the scale required to cover development residuals and applies it before averaging, rather than treating disagreement as independent additional evidence [@kennedy2001calibration; @brynjarsdottir2014discrepancy]. Discrete model weights are fitted on development observations only by case-blocked log score, so that a case with many candidate edges cannot dominate the fit through row count alone [@neuman2003bma]. The averaging report retains the component distributions, the pairwise total-variation disagreement and the compatible outcome set; a frozen threshold of 0.25 marks material disagreement, above which the mixture is retained for audit but not promoted to a single reportable interpretation. Convergence is assessed by the constrained simplex Karush–Kuhn–Tucker residual and objective stopping criteria. The unconstrained gradient norm is retained as a diagnostic but is not the correct stationarity measure on a simplex, and the two are reported separately rather than interchangeably.

Unsupported likelihoods, missing required fields and incompatible outputs are written as auditable failure or abstention records before any interval arithmetic or policy ranking. An abstention is never counted as a correct prediction. False commitment denotes a committed interpretation incompatible with the sealed outcome under the declared outcome set.

## 3.4 Validation design

Three tiers of evidence were used, in ascending order of independence from the framework.

The **recovery benchmark** comprises 100 realisations generated from a single locked ground-truth specification of nodes, edges, transport processes and reaction extents. Its replication structure must be stated precisely, because it governs how the results may be read. The truth is fixed: the same ten candidate edges, six transport parameter values, fourteen site residence times and fifteen distinct reaction-extent values recur in every realisation, and no truth value varies across realisations. What varies is measurement noise, drawn from a seeded generator per realisation with a major-ion relative standard deviation of 0.035. The 100 realisations are therefore noise replicates of one scenario, not 100 independent aquifers. Scored rows accordingly overstate the units of independent information by two orders of magnitude, and every recovery statistic in Section 4.2 is reported both at row level and aggregated within distinct truth units, with the clustered value taken as the honest one. This tier establishes whether an inference target is recoverable when the generating process is known exactly; it cannot establish transfer to other geometries.

A structural analysis accompanies the empirical one, because a failed inversion by one solver does not by itself establish that a target is non-identifiable. The stoichiometry matrix $A$ was assembled from the declared reaction dictionary, with one row per ion and one column per reaction term, and its rank, null-space dimension and condition number were computed by singular value decomposition. A non-trivial null space means that distinct extent vectors produce identical predicted chemistry, so no estimator can separate them from chemistry alone regardless of penalty, solver or sample size.

The **locked programme** is a controlled-synthetic run in which three independent generator families produce observations from sealed states withheld until scoring: an analytic lattice, a multilayer mixing generator, and a MODFLOW 6 and MODPATH 7 generator coupled to nonlinear chemistry and tracer equations. The generators do not import HydroSheaf. Twelve cases were used for development and six were held out. Observation stresses included structured missingness, left censoring, combined stress, and permutation of the tracer, chemistry and head associations; permutation preserves marginal values while destroying case-specific association, so a method that still appears to perform is drawing on leakage rather than signal. Calibration, thresholds and policy tuning were restricted to development cases. Age, reaction, kinetic and integrated gates were evaluated against thresholds fixed before scoring, and are reported beside those thresholds in Section 4.3. A versioned analysis plan archived these rules before assembly; it was not an external preregistration.

The **external comparison** scores framework output against evidence generated outside the study. For residence time, inferred ages were compared with published lumped-parameter estimates under leave-one-study-unit-out folds, so that no site contributed to any quantity used in its own prediction. Two variants are reported and must not be conflated. The uncalibrated strict-parity variant compares the framework's own estimate with the published value directly and is the primary result. The calibrated variant applies a ridge correction fitted on the training folds; because that correction is fitted to reproduce the reference outputs, it measures how well the framework can emulate the reference, not whether it independently agrees with it. For topology, inferred directed edges were compared with a particle-tracking reference in two modes, and these are likewise never pooled: a no-prior mode inferring edges from elevation and head gradient alone, and a prior-assisted mode that ingests the reference structure and therefore measures ingestion fidelity rather than inference.

Prospective decision utility was evaluated by presenting each policy with candidate measurement actions, costs and pre-measurement distributions. Policies selected an action or abstained before any hidden state was released; only then was the outcome attached and utility minus cost computed, with paired bootstrap contrasts against random and specialist policies [@chaloner1995design; @sreekanth2017monitoring]. Because the resampling unit is the six locked cases, the resulting intervals are within-run descriptive summaries and not confidence intervals over aquifers or generator families.

## 3.5 Computation, provenance and reproducibility

Computation is owned by Python; every reported statistic was recomputed from primary run records and primary data files by scripts under `analysis/python/`, which write tidy read-only CSV exports. All figures were produced in R from those exports and recompute no reported statistic. Figure colours use a colour-vision-deficiency-safe categorical palette verified for lightness, chroma and adjacent-pair separation, with identity always carried by a legend or direct label rather than colour alone.

No value in this paper was inherited from an earlier internal report. Where recomputation disagreed with an earlier internal value, the recomputed value stands and the disagreement is recorded Table 4.

Two provenance limitations are stated rather than resolved. The locked programme run was generated before all source changes were committed: its manifest records the generating revision, cryptographic hashes for 31 source files and a flag indicating an uncommitted working tree. The hashes therefore provide traceability, but they do not substitute for regeneration from a committed tree, which remains an outstanding action. Second, the Northern Ghana seasonal attribute is a reconstruction, as set out in Section 2.2, and is excluded from all inference.

# 4. Results

## 4.1 Implemented package and execution

The installed package loaded as HydroSheaf 0.5.1, exposing the capabilities and claim scopes listed in Table 2. The locked programme execution gate returned PASS: all stages completed without runtime error, the independent generator-quality audit passed with no blocking or major finding across 18 generator records, stress-test records were complete, and the external MODFLOW 6 and MODPATH 7 execution gate passed. The run comprised 12 development and six held-out cases from three independent generator families. The provenance record retained source hashes and the generating revision, and recorded that source changes had not been committed at generation time.

## 4.2 What the recovery benchmark recovers, and what it does not

The 100-realisation benchmark separates the framework's inference targets sharply Fig. 4.

Because the truth is fixed across realisations, statistics are reported at both levels: over scored rows, and aggregated within the distinct truth units that carry the independent information (six transport parameters, fourteen site residence times, twenty-one active edge-by-reaction combinations). The clustered figure is the one to read.

**Transport parameters are recoverable.** Across 550 scored parameter rows the median absolute error was 0.058 in fraction units (interquartile range 0.030 to 0.084); aggregating within the six distinct parameter units gives 0.064. The correct transport model — evaporation against mixing — was selected in 91.7% of cases. Fifty rows returned no recovered value and were recorded as such rather than dropped.

**Residence time is recoverable, and the network constraint measurably helps.** Over 1,400 site-realisations spanning fourteen distinct site residence times, the single-node lumped inversion achieved a log-space R² of 0.922 with a mean absolute error of 477 years. The network-constrained Bayesian mode improved both, to a log-space R² of 0.983 and a mean absolute error of 183 years, while halving the mean credible-interval width from 1,714 to 823 years. The improvement holds across every residence-time class: median absolute relative error fell from 0.177 to 0.081 for young waters, from 0.260 to 0.136 for intermediate waters, from 0.433 to 0.196 for old waters and from 0.450 to 0.238 for fossil waters. The largest relative errors in both modes occur in the mixed class (0.574 and 0.344). This is the configuration a lumped model represents least well, and since binary mixtures of young and old water are common in fractured basement settings, it is also the case where field application would be least reliable.

**Point reaction extents are not identifiable, and this is structural rather than a solver failure.** The stoichiometry matrix assembled from the declared reaction dictionary has 14 columns, one per reaction term, but only 8 of the 11 declared ions participate in any reaction, and its rank is 8. The null space therefore has dimension 6 Table 5. Six independent directions in extent space produce exactly no change in predicted chemistry, with maximum residual ion change at machine precision (of order 10⁻¹⁶ mmol/L). These directions are chemically interpretable rather than numerical artefacts: the largest pairs calcium-sodium against sodium-calcium exchange, another trades calcite against albite and dolomite, and a third trades sulfate reduction against calcite and gypsum. Any estimator working from major-ion chemistry alone must choose among these indistinguishable combinations by penalty rather than by evidence. No regularisation weight, solver or sample size can remove a null space; it can only pick a point within it.

The empirical result follows the structural one. Over the 2,100 scored active reaction terms the coefficient of determination was −0.023, and aggregating within the twenty-one distinct active units gives 0.011 — indistinguishable from zero either way, meaning the inversion carries essentially no information about extent magnitude beyond the mean. The complementary failure is sharper still. Of the 4,900 truth-zero terms, 54.1% were assigned a recovered magnitude above 0.05 mmol/L, while only 42.1% of genuinely active terms were detected at the same threshold. The inversion asserts reactions that did not occur more often than it detects reactions that did, and this ordering holds at every threshold tested from 0.01 to 0.20 mmol/L Table 5. The 0.05 mmol/L threshold is not arbitrary: propagating the benchmark's 3.5% major-ion analytical noise through typical concentrations gives a detection floor of approximately 0.035 mmol/L, so the threshold sits just above the level at which noise alone could masquerade as a reaction.

This is the paper's pivotal result. It is reported in preference to the values circulated in earlier internal reporting of the same benchmark, which are not reproducible from the shipped record under the benchmark's own definitions Table 4.

## 4.3 The same evidence, gated

The locked programme scores a different representation of the same chemical evidence: not point extents, but calibrated probabilities over reaction families with an explicit abstention option. All four component gates met their predeclared thresholds Fig. 5, Table 3.

For **reaction families**, multiclass log loss over 142 locked records was 0.896 against 2.852 for the competence-matched baseline, with coverage 0.739 against a floor of 0.25, maximum classwise expected calibration error 0.271 against a ceiling of 0.35, selective risk 0.219 against a ceiling of 0.4, and a false commitment rate of 0.038 against a ceiling of 0.10. The contrast with Section 4.2 is the substantive finding: the evidence that cannot identify how much of a reaction occurred can support a calibrated statement about which family of reactions is compatible, provided the method is permitted to abstain.

For **age**, coverage including abstention was 0.966 against a target of 0.95, with a relative interval width of 0.544 against a ceiling of 1.5 and selective risk of 2.10 years against a ceiling of 12. The specialist mean absolute error was 4.24 years against 7.65 for the competence-matched baseline, at an acceptance rate of 0.966.

For **kinetics**, predictive root-mean-square error was 0.119 against a ceiling of 0.25 and interval coverage was 1.0. Parameter identification behaved exactly as the structural analysis predicts: conditional identification was 1.0 among cases supplied with an independent surface-area observation, but overall identification was 0.167 and the parameter abstention rate 0.833. The adapter predicted the response and withheld the parameters when only their product was constrained.

For the **integrated policy**, mean cost-adjusted utility per unit cost was 1.478, against 0.246 for a uniform random policy and −0.004 for the strongest specialist comparator, with zero observed false commitment across the six locked prospective cases. Paired differences were 1.23 (95% interval 0.99 to 1.45) against random and 1.48 (1.22 to 1.75) against the specialist. These intervals resample the six locked cases and are therefore within-run descriptive summaries, not population bounds. No general-superiority gate was defined or passed, and none is claimed.

## 4.4 Comparison with externally generated evidence

Against **published lumped-parameter age estimates**, held-out uncalibrated agreement was screening-level: over 675 sites with finite paired values across 20 study units, log-space R² was 0.482, median absolute log₁₀ error 0.228, log₁₀ root-mean-square error 1.014, and 56.4% of estimates fell within a factor of two of the published value, with 81.6% within a factor of ten Fig. 6. Bias was negligible (−0.036 in log₁₀ units), so the disagreement is dispersion rather than systematic offset.

The 675 paired sites are 53% of the 1,272 reference records, and the attrition requires comment because it is not random. All 597 unpaired records were lost because the published reference carried no total age; the framework returned an estimate for every site, so none was lost on the estimate side. The attrition is nevertheless uneven across aquifer groups, from 9.1% for glacial and 20.7% for carbonate aquifers to 53.3% for sandstone and 75.8% for western unconsolidated aquifers, and unpaired sites are deeper on average (median 148 m against 122 m). The comparison is therefore representative of the aquifer types for which published lumped-parameter ages exist, and is weighted toward carbonate and glacial settings. It should not be read as uniform coverage of the reference release.

Applying a ridge correction fitted on the training folds raised R² to 0.962, reduced the log₁₀ root-mean-square error to 0.275 and brought 76.9% of sites within a factor of two. That improvement is real but it measures a different quantity. The correction is fitted to reproduce the published outputs, so the calibrated result quantifies how well the framework can emulate the reference, not whether it independently agrees with it. Reporting only the calibrated number would overstate independent agreement by a factor of two in explained variance. The honest headline is the uncalibrated one.

Against the **particle-tracking directed topology**, no-prior inference from elevation and head gradient recovered 147 of 174 reference edges, giving a recall of 0.845, but proposed 302 edges in total, giving a precision of 0.487 and a false discovery rate of 0.513. The resulting F1 was 0.618. The framework therefore finds most of the reference structure while proposing slightly more false connections than true ones, which is the expected behaviour of a screening step and is reported as such. The prior-assisted mode reproduces the reference exactly (precision, recall and F1 of 1.0); this measures ingestion fidelity, contains no independent information, and is reported separately for that reason. No topology-superiority claim is made.

## 4.5 Field application and evidence ceiling

Applied to the two measured field datasets, the chemistry workflow generated 208 candidate edges with an overall median in-sample chemistry closure of 0.713 Fig. 7. Of these, 161 were fitted with a mixing transport model and carry a site-identifying local end-member, 95 at Talensi and 66 at Lower Anayari; the remaining 47 were fitted with an evaporation model, which references no local end-member and therefore carries no site label. The exclusion is thus by transport model rather than by quality, and the excluded edges are not systematically worse (median closure 0.658 against 0.713 overall). Per-site medians are reported over the 161 site-assigned edges and are 0.795 at Talensi and 0.599 at Lower Anayari; all other field statistics use the full 208. These are closures of a candidate edge against the same chemistry used to fit it. They are not held-out predictions and they are not evidence of physical connectivity. That eight of the 161 site-assigned edges close worse than the sample mean, with coefficients as low as −0.52, is itself useful: the workflow surfaces candidate pairings that its own chemistry cannot support.

The number of retained reaction alternatives per candidate edge is the more informative output. The modal edge retains three reaction terms above 0.05 mmol/L and the distribution extends to eight, so the workflow reports a set of compatible explanations rather than a single mechanism. Given the extent-recovery result of Section 4.2, this is the correct behaviour: a single reported mechanism at these data densities would be an artefact of the fitting procedure.

A cumulative measurement-tier ablation over the 160 Northern Ghana wells quantifies which determinands actually buy identifiability Fig. 7. With major ions alone, 60.0% of wells are non-identifiable at the reaction-family level. Adding stable isotopes reduces this to 51.9%, and adding fluoride leaves it unchanged at 51.9%. Adding strontium and dissolved silica reduces it to 0.6%, and full metadata to zero. The mean resolution score moves only from 68.4 to 71.0 across the entire sequence, and the mean number of retained reaction alternatives stays between 5.4 and 6.2. Identifiability is therefore gained almost entirely from two determinands, and it is gained without narrowing the compatible explanation set — the wells become classifiable, not resolved to a single mechanism. For survey design in this setting the implication is specific: strontium and silica are worth their analytical cost, and additional stable-isotope or fluoride determinations are not, if the objective is reaction-family identifiability.

The Talensi result also illustrates the value of independent analytical screening. Because 58 of its 63 samples fall in the exploratory charge-balance tier, every downstream Talensi product is confined to screening use by the framework itself, without analyst intervention.

The overall field outcome is therefore an evidence ceiling rather than a validation. The datasets support ingestion, quality tiering, candidate generation, reaction-family plausibility, measurement-value diagnostics and explicit non-identifiability reporting. They do not support age, connectivity or mechanism claims, and the framework returns abstentions for all three.

# 5. Discussion

## 5.1 Identifiability, not accuracy, is the organising variable

The results do not divide into successes and failures of the software. They divide by whether the inference target is identifiable from the evidence supplied.

Two of the three failures reported here are structural, and can be demonstrated rather than merely observed. For reaction extents, the stoichiometry matrix has rank 8 against 14 unknown terms, leaving a six-dimensional null space of extent combinations that produce identical chemistry. For kinetics, the rate law depends on the product of the rate constant and the reactive surface area, so the sensitivity columns are collinear and the numerical rank falls to one. In both cases the non-identifiability is a property of the forward map, not of the estimator: a different solver, a better penalty or a hundred times more data would leave it untouched. The kinetic case shows what actually breaks such a degeneracy, since conditional identification rises to 1.0 when an independent surface-area measurement is supplied. The remedy for structural non-identifiability is a different kind of measurement, not more of the same kind.

This distinction matters practically because the two failure modes look identical in a results table. An inversion that returns a poor coefficient of determination might be badly tuned, or it might be estimating something the data cannot distinguish. Reporting the rank of the forward map alongside the fit statistic separates them, and it is a cheap check that inverse hydrogeochemical studies could adopt generally.

This reframes what an integrated framework is for. Its value is not that it produces an answer for every target, but that it distinguishes the targets its evidence can reach from those it cannot, and reports the distinction rather than absorbing it. The reaction result is the clearest demonstration. The same chemical observations that leave a six-dimensional null space in extent space, and that consequently yield no explanatory power for point extents, support a calibrated family-level statement with a log loss of 0.896 against a baseline of 2.852 and a false commitment rate of 0.038. Nothing was added between those two results. What changed was the estimand and the availability of an abstention. Reaction families are identifiable precisely because they are unions of the directions the null space mixes together, so a statement at family resolution does not require separating what the chemistry cannot separate. An inference target that is unreachable at one resolution can be reachable at a coarser one, and locating that resolution is a scientific act rather than a presentational one.

This is a familiar move in data-limited hydrology, where transferability and degradation under missing information are treated as scientific questions in their own right rather than as obstacles to be assumed away [@sivapalan2003pub; @hrachowitz2013pub]. The contribution here is to make the choice of resolution an output of the evidence rather than a decision the analyst takes in advance.

The practical consequence for data-limited hydrogeology is direct. A reported reaction extent in millimoles per litre, derived from a single hydrochemical snapshot along an assumed flow path, should be treated as a fitted quantity and not as a measured one, irrespective of how well the mass balance closes. The 208 field candidate edges close with a median coefficient of determination of 0.713 while retaining a modal three compatible reaction terms; those two facts are consistent, and only the second is informative about the subsurface.

## 5.2 Relationship to existing tools

HydroSheaf does not replace the specialist forward models and is not positioned against them. MODFLOW 6 and MODPATH 7 solve flow and trace particles given a conceptual model and parameterisation [@langevin2017modflow6; @pollock2016modpath]; TracerLPM fits lumped age distributions to environmental tracers [@jurgens2012tracerlpm]; PHREEQC performs speciation, reaction-path and inverse calculations [@parkhurst2013phreeqc]. Each is stronger than the corresponding HydroSheaf component within its own domain and given its own data requirements. The framework's contribution lies where those requirements are unmet: retaining candidate explanations, comparing component outputs, quantifying disagreement, and deciding when not to collapse alternatives. Graph representations have been advanced for related reasons [@moracchini2025graphflow], and sheaf-theoretic language describes the consistency of overlapping observations on a network [@robinson2017sheaves; @hansen2019spectral]; the usage here is sheaf-inspired, and no formal sheaf theorem is claimed or demonstrated.

The topology comparison locates the contribution precisely. Inferring direction from elevation and head gradient alone recovered 84.5% of the reference edges but produced more false connections than true ones (precision 0.487). As a competitor to particle tracking this is a poor result. As a screening step that narrows an unconstrained set of well pairs to a candidate set for further evidence, with its false discovery rate stated, it is a useful one. The distinction matters because the prior-assisted mode reproduces the reference exactly, and reporting that number as an inference result — rather than as the ingestion-fidelity check it is — would misrepresent the method's capability entirely. Keeping the two modes separate is a reporting requirement, not a stylistic preference.

## 5.3 Calibrating to a reference is not agreeing with it

The external age comparison carries a methodological warning that generalises well beyond this framework. Held-out uncalibrated agreement with published lumped-parameter estimates gave a log-space R² of 0.482. Fitting a ridge correction on the training folds raised it to 0.962. Both numbers are correct, and they answer different questions. The first asks whether the framework's independent estimate agrees with an external interpretation; the second asks whether the framework carries enough signal to be mapped onto that interpretation once a correction is fitted. Only the first is evidence of independent agreement, and reporting the second alone would double the apparent explained variance.

This distinction is easy to lose in a benchmarking pipeline, because the calibrated variant is legitimately useful — a well-calibrated emulator is a real capability — and because it produces the better-looking figure. It is worth stating plainly that a comparison is only independent to the extent that no quantity used in the prediction was fitted to the target, and that leave-one-unit-out folding protects against site-level leakage but not against the fitting of the correction itself. Neither reference used here is physical truth in any case: a fitted lumped-parameter output is a model interpretation of tracer data [@jurgens2012tracerlpm], and a particle path is conditional on the flow model that generated it [@oreskes1994validation].

## 5.4 Abstention as a reported quantity

Treating abstention as an output rather than a failure changes what the software can be trusted to do. An abstention costs coverage: the age component's coverage of 0.966 is reported with abstentions in the denominator and beside its acceptance rate, precisely so that a method cannot obtain a strong selective score by declining every difficult case. Under this accounting the kinetic component's overall identification rate of 0.167 is not a poor result but an accurate one, since 83.3% of locked cases genuinely lacked the independent surface-area information required to separate the confounded parameters.

The same logic applies to model discrepancy. Where components disagreed beyond the frozen total-variation threshold, the mixture was retained for audit rather than promoted to a single interpretation, following the argument that unrepresented model inadequacy biases estimates and understates uncertainty [@kennedy2001calibration; @brynjarsdottir2014discrepancy]. This is the mechanism by which an integrated framework avoids becoming more confident merely by combining more evidence, and it is why equifinality is handled by reporting the compatible set rather than by selecting within it [@beven2001glue; @beven2006manifesto].

## 5.5 Which measurements buy identifiability

The most directly transferable result in this study is also the cheapest to act on. Across the 160 Northern Ghana wells, reaction-family identifiability is bought by two determinands. Major ions alone leave 60.0% of wells non-identifiable. Adding stable isotopes reduces that to 51.9%; adding fluoride changes it not at all. Adding strontium and dissolved silica reduces it to 0.6%.

The structural analysis of Section 4.2 explains why this should be expected rather than surprising. Strontium and silica participate in weathering reactions that the major-ion suite constrains only through the same linear combinations that span the null space, so they contribute rows to the stoichiometry matrix that are not already in its row space. Stable isotopes and fluoride largely re-express information the major ions already carry, about mixing and about a single geogenic source respectively. A determinand improves identifiability only to the extent that it adds an independent constraint, and which determinands do so is computable in advance from the reaction dictionary rather than discoverable only after sampling.

Two qualifications bound the result. The mean number of retained reaction alternatives stays between 5.4 and 6.2 across the whole tier sequence and the mean resolution score moves only from 68.4 to 71.0, so wells become classifiable rather than resolved to a single mechanism. And the finding is conditional on the declared dictionary. Neither qualification weakens the practical implication: for survey design in comparable crystalline basement terrain, strontium and silica are worth their analytical cost and further isotope or fluoride determinations are not, if the objective is reaction-family identifiability. That is a budgeting statement a framework reporting its own evidence boundary can produce and one that always returns an answer cannot.

## 5.6 Limitations

Several limitations bound every claim above.

The **recovery benchmark is one scenario, replicated**. Its 100 realisations share a single fixed truth, so its independent information is carried by six transport parameters, fourteen site residence times and twenty-one active reaction terms. The structural identifiability result does not depend on this, since it is a property of the reaction dictionary rather than of the data, but every empirical recovery statistic does. Whether the recovery of transport parameters and residence times transfers to other geometries is untested.

The locked programme rests on **six held-out cases** from three generator families. Its uncertainty is dominated by generator-domain coverage rather than by sampling, and the paired intervals in Section 4.3 are within-run descriptive summaries rather than population bounds. Broader geometries, heterogeneity regimes, tracer histories and reaction mechanisms are required before any statement about transfer to arbitrary aquifers. Six cases are also too few to discriminate a passing method from a failing one with any power, so the programme is better read as a working test of the claim architecture than as a performance evaluation.

The **null-space result is conditional on the declared reaction dictionary**. A different dictionary, or the addition of determinands that participate in reactions the present ions do not constrain, would change the rank and could reduce the degeneracy. That is a testable prediction rather than a caveat, and it is the mechanism behind the strontium and silica result reported below.

The **field application validates nothing**, by construction. No dataset contains environmental age tracers, screen intervals, repeated heads, independent connectivity truth or reaction-mechanism truth (Section 2.5). Field results are ingestion, screening and non-identifiability diagnostics.

The **Northern Ghana seasonal attribute is reconstructed**, and no seasonal or repeated-head result is reported. The measured unit is the well.

**Talensi is analytically limited**, with a systematic anion excess placing 58 of 63 samples in the exploratory tier. Its use is confined to screening and failure-mode diagnostics.

The **reaction dictionary is fixed**, so a process outside the declared families cannot be represented, and the family-level result of Section 4.3 is conditional on that dictionary. Perturbing the dictionary itself, and adding held-out mechanisms, would test the component's scientific scope rather than its coefficients.

The **locked run was generated from an uncommitted working tree**. Its hashes give traceability but not clean-tree regeneration, which remains outstanding.

**No general-superiority gate exists** for any component. Where a comparator contrast is positive it is reported as a bounded within-programme result.

## 5.7 Way forward

Four steps follow directly. First, regenerate the locked programme from a committed tree and compare result hashes. Second, give topology a dedicated contract with independent hydraulic-gradient and particle-tracking baselines and explicit no-prior controls, so that the screening role established here can be quantified rather than inferred. Third, expand the locked case set and declare a superiority or non-inferiority rule before scoring. Fourth — and this is the step the field data themselves identify — a prospective sampling campaign carrying environmental age tracers, screened intervals, repeated heads and mineralogical or redox corroboration. The measurement-value diagnostics reported here exist to inform the design of exactly that campaign, and they already give it a concrete constraint: strontium and dissolved silica move the non-identifiable fraction from 51.9% to 0.6%, whereas further stable-isotope or fluoride determinations move it not at all. Until that stage, the supported claim is bounded controlled-synthetic performance together with field readiness, and not field validation.

# 6. Conclusions

HydroSheaf was implemented as a claim-gated evidence-integration framework that links candidate connectivity, residence-time, reaction and kinetic inference to calibration, model-discrepancy reporting, abstention and cost-aware measurement selection, and it was evaluated across three tiers of increasing independence.

The evaluation divides by identifiability rather than by component, and in two cases the division can be demonstrated from the forward map rather than inferred from a fit. Point reaction extents are structurally non-identifiable from major-ion chemistry: the stoichiometry matrix of the declared dictionary has rank 8 against 14 unknown terms, leaving a six-dimensional null space of extent combinations that produce identical predicted chemistry. The empirical result follows, with no explanatory power for extent magnitude and non-trivial extents asserted on 54.1% of truth-zero terms against only 42.1% detection of active ones. Kinetic parameters are non-identifiable for the same kind of reason, the rate constant and reactive surface area entering only as a product, and identification recovers fully when an independent surface-area measurement breaks the collinearity. Where the forward map does constrain the target, the framework performs: transport parameters were recovered to a median absolute error of 0.058 with 91.7% correct model selection, and network-constrained residence-time inference improved on single-node inversion in every age class, raising log-space R² from 0.922 to 0.983 while halving interval width. Recast at a resolution the evidence can support, calibrated reaction-family probabilities with an abstention option passed a predeclared gate with a log loss of 0.896 against 2.852 and a false commitment rate of 0.038.

Against externally generated evidence the framework performs at screening level: held-out uncalibrated agreement with published lumped-parameter ages gave a log-space R² of 0.482 with 56.4% of estimates within a factor of two, and no-prior directed-topology inference recovered 84.5% of reference edges at a false discovery rate of 0.513. Fitting a correction to the published ages raised R² to 0.962, but that quantity measures emulation of the reference and not independent agreement, and the two are reported separately for that reason.

The contribution is therefore an architecture that determines, from evidence rather than assumption, which of its own outputs are supported, together with a demonstration that the boundary falls in a specific and computable place. Reporting the rank of the forward map beside the fit statistic distinguishes an estimator that is badly tuned from one that is estimating something the data cannot separate, and that check is inexpensive enough for inverse hydrogeochemical studies to adopt generally. The field application in northern Ghana returns an evidence ceiling rather than a validation: the datasets contain no environmental age tracers, screen intervals, repeated heads, connectivity truth or mechanism truth, and the framework abstains on all of them while identifying strontium and silica as the determinands that would most improve identifiability. The next stage is a clean-tree regeneration of the locked programme with wider generator coverage, followed by prospective field sampling carrying the independent age, head, screen and process evidence that the present data lack. Until then the supported claim is bounded controlled-synthetic performance together with field readiness, and not field validation or superiority over specialist forward models.

# 7. Statements

## Code and data availability

HydroSheaf source code, package metadata, tests, benchmark protocols, the analysis
scripts that generate every reported value, and the R scripts that generate every
figure are archived at `https://doi.org/10.5281/zenodo.PLACEHOLDER-DOI` and
developed at `https://github.com/PLACEHOLDER-ORG/hydrosheaf`. The repository
carries an open-source licence, an English README with installation and usage
instructions, a declared dependency lock, worked tutorials, a user guide
describing inputs, outputs and options, and synthetic test cases sufficient to
reproduce the main results. Files are provided uncompressed.

The controlled-synthetic evidence is linked to a machine-readable provenance
manifest recording the generating source revision and cryptographic hashes for 31
source files. That run was generated before all source changes were committed;
the manifest records this, and regeneration from a committed tree is an
outstanding action rather than a completed one.

Three field datasets were used. Lower Anayari and Talensi derive from the
published studies cited in Section 2 and are redistributed subject to those
studies' terms. The Northern Ghana panel is an author-compiled regional dataset:
its chemistry is measured, but the dry and wet sheets do not represent
independent seasonal resampling, and the seasonal attribute is a reconstruction.
Static water level is identical across the two sheets in all 160 wells. That
panel is therefore used only at the well level, its seasonal dimension is
excluded from every inferential result, and no seasonal-change or repeated-head
finding is reported. No new field samples were collected for this study.

## Author contributions

To be completed with named CRediT roles before submission: conceptualisation,
methodology, software, validation, formal analysis, data curation, writing —
original draft, writing — review and editing, and supervision.

## Ethics, competing interests and funding

The study used existing measurements and controlled synthetic data, and involved
no new human or animal subject procedures. The authors declare no competing
interests. Funding is to be declared in the final submission metadata.


# Tables


**Table 1.** Inventory of the three measured field datasets. Charge-balance error was recomputed independently from the reported ions. Northern Ghana is summarised as its 160-well primary measured panel; its second seasonal panel is reconstructed and excluded from inference.

| Dataset | Samples | Sites | Variables | Median &#124;CBE&#124; (%) | Quant. / screen. / expl. | Sr, SiO<sub>2</sub> | Age tracers | Supported use |
|---|---:|---:|---:|---|---|---|---|---|
| Northern Ghana | 160 | 160 | 29 | 1.54 | 151 / 7 / 2 | yes | no | Regional chemistry, quality tiering, measurement-value ablation |
| Lower Anayari | 41 | 41 | 22 | 3.13 | 36 / 5 / 0 | no | no | Sparse chemistry and candidate screening |
| Talensi | 63 | 63 | 23 | 29.85 | 0 / 5 / 58 | no | no | Failure-mode and screening-level transfer |


**Table 2.** Implemented package capabilities, the module providing each, and the claim scope each supports.

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


**Table 3.** Locked controlled-synthetic gate results against thresholds fixed before scoring. Age selective risk is in years; kinetic predictive RMSE is in mmol/L. The kinetic identification rate is conditional on an independent surface-area measurement; overall identification was 0.167 with a parameter abstention rate of 0.833.

| Component | Metric | Observed | Predeclared threshold | Status |
|---|---|---:|---|---|
| Age | Coverage including abstention | 0.9655 | at least 0.95 | met |
| Age | Relative interval width | 0.5445 | at most 1.5 | met |
| Age | Selective risk | 2.0967 | at most 12 | met |
| Age | Acceptance rate | 0.9655 | at least 0.5 | met |
| Reaction family | Coverage | 0.7394 | at least 0.25 | met |
| Reaction family | Maximum classwise calibration error | 0.2710 | at most 0.35 | met |
| Reaction family | Selective risk | 0.2190 | at most 0.4 | met |
| Reaction family | False commitment rate | 0.0385 | at most 0.1 | met |
| Kinetics | Predictive RMSE | 0.1191 | at most 0.25 | met |
| Kinetics | Interval coverage | 1.0000 | at least 0.9 | met |
| Kinetics | Conditional identification rate | 1.0000 | at least 0.8 | met |
| Integrated decision | Locked prospective cases | 6 | at least 6 | met |
| Integrated decision | Observed false commitment | 0.0000 | at most 0.1 | met |
| Age | Specialist MAE vs baseline 7.6548 years | 4.2396 | lower than baseline | met |
| Reaction family | Log loss vs baseline 2.8520 | 0.8964 | lower than baseline | met |


**Table 4.** Quantities where recomputation from primary records disagreed with earlier internal reporting. In every case the recomputed value stands.

| Quantity | Earlier internal value | Recomputed value | Resolution |
|---|---|---|---|
| Reaction extent recovery, R2 (active terms) | 0.57 (M2 abstract) | -0.023 | Recomputed from reaction_recovery.csv using the benchmark's own active definition (&#124;truth&#124; > 0.01 mmol/L, n = 2100). The earlier value is not reproducible from the shipped record under this or any tested subset definition. |
| Reaction extent recovery, MAE (active terms) | 0.33 mmol/L (M2 abstract) | 0.374 mmol/L | Recomputed under the same active definition. |
| Inactive-term leakage above 0.05 mmol/L | 32.7% (M2 abstract) | 54.1% | Recomputed over all inactive truth terms (n = 4900). |
| Public age benchmark, n | 1,249 (M2 abstract) / 1,272 (M2 summary) | 675 finite held-out pairs | Recomputed from the parity record; earlier counts mixed row totals with finite-pair counts. |
| Public age benchmark, R2 | 0.71 (M2 abstract) | 0.482 uncalibrated held-out | Earlier abstract reported a value reproducible from neither the uncalibrated nor the calibrated record. |
| Public age benchmark, median &#124;log10 error&#124; | 0.17 (M2 abstract) / 0.383 (M2 summary) | 0.228 uncalibrated held-out | 0.17 corresponds to the calibrated emulation, which is not an independent comparison (D4). |
| Northern Ghana records | 320 seasonal field records | 160 wells; second seasonal panel reconstructed | Author confirmed the seasonal split was reconstructed (D2). |
| Talensi sample longitudes | Positive (east) in the source file and in M2/M2_1 | Negated; Talensi District lies near 0.8 deg west | As delivered all 63 samples plot outside Ghana. After negation all 63 fall inside the Upper East Region polygon. Because reflection about the prime meridian is an isometry, pairwise sample distances are unchanged, so the distance-based candidate-edge results are unaffected; only absolute position and the study-area map change. |
| Field candidate-edge chemistry closure | 0.711 (M2) / 0.713 (M2_1), described as overall R2 | median 0.713 over 208 edges; mean 0.611 | The two earlier values are the median and a differently pooled summary of the same quantity. The median is reported and named as such. |


**Table 5.** Structural identifiability of the reaction-extent inverse problem, and detection behaviour as a function of threshold. The stoichiometry matrix is rank deficient, so a six-dimensional family of extent vectors produces identical predicted chemistry and no estimator can separate them from major-ion data alone. At every threshold tested the inversion asserts extents on truth-zero terms more often than it detects genuinely active ones. Propagated analytical noise gives a detection floor near 0.035 mmol/L.

| Quantity | Value | Note |
|---|---|---|
| Reaction dictionary terms | 14 | Unknown extents to be recovered |
| Ions in the declared order | 11 | Maximum number of mass-balance constraints |
| Ions appearing in at least one reaction | 8 | Constraints that actually bind: Ca, Mg, Na, HCO3, Cl, SO4, NO3, Fe |
| Rank of the stoichiometry matrix | 8 | Independent constraints available |
| Rank over binding ions only | 8 | Unchanged by dropping non-participating ions |
| Null-space dimension | 6 | Directions in extent space that leave the chemistry unchanged |
| Condition number of the stoichiometry matrix | 13.65 | Ratio of largest to smallest non-zero singular value |
| Detection threshold 0.01 mmol/L | 64.1% / 49.7% | False assertion on truth-zero terms / detection of active terms |
| Detection threshold 0.02 mmol/L | 61.3% / 47.8% | False assertion on truth-zero terms / detection of active terms |
| Detection threshold 0.05 mmol/L | 54.1% / 42.1% | False assertion on truth-zero terms / detection of active terms |
| Detection threshold 0.1 mmol/L | 44.1% / 30.5% | False assertion on truth-zero terms / detection of active terms |
| Detection threshold 0.2 mmol/L | 29.7% / 19.5% | False assertion on truth-zero terms / detection of active terms |
| Detection threshold 0.5 mmol/L | 15.5% / 1.4% | False assertion on truth-zero terms / detection of active terms |


# Figure captions


**Fig. 1.** Evidence-to-claim architecture of the HydroSheaf framework. Observations enter as nodes and candidate directed edges; candidate sets are retained through specialist scoring, calibration and discrepancy-aware averaging, and are resolved only at the claim gate, which returns PASS, WEAK or ABSTAIN. Calibration is fitted on development cases and frozen; locked cases are scored after prediction. (`artifacts/figures/FIG-1_*.pdf`)


**Fig. 2.** The three measured field datasets. (a) Sampling locations across the four northern administrative regions, with neighbouring countries in grey and an inset locating the study frame within Ghana; administrative boundaries from the geoBoundaries database (CC BY 4.0). (b) Independently recomputed charge-balance error, with bands marking the 5% and 10% acceptance tiers; Talensi shows a systematic anion excess. (c) Dominant-ion facies as a share of each dataset. (d) Stable isotope composition with dataset-specific local regressions against the global meteoric water line; all slopes fall below 8, indicating evaporative enrichment during recharge. Northern Ghana is shown as its 160-well primary measured panel. (`artifacts/figures/FIG-2_*.pdf`)


**Fig. 3.** Variable availability and the resulting evidence ceiling. (a) Presence of each canonical determinand in each dataset; numbers give percentage completeness where a carried variable is incomplete. (b) What the available variables can and cannot support. No dataset carries environmental age tracers, screen intervals, repeated heads, independent connectivity truth or reaction-mechanism truth. (`artifacts/figures/FIG-3_*.pdf`)


**Fig. 4.** Component recovery across 100 synthetic realisations with known truth. (a) Transport parameters are recovered. (b) Residence time is recovered, and network constraints improve on single-node inversion in every age class. (c) Point reaction extents are not recovered: the active-term coefficient of determination is negative and over half of truth-zero terms are assigned a non-trivial extent. Dashed lines are 1:1; dotted lines in (c) mark the 0.05 mmol/L leakage threshold. (`artifacts/figures/FIG-4_*.pdf`)


**Fig. 5.** Locked controlled-synthetic programme results. (a) Observed component metrics against thresholds fixed before scoring. (b) Cost-adjusted utility per unit cost for each policy over six locked truth-blind cases. (c) Paired differences with 95% intervals; because the resampling unit is the six locked cases these are within-run descriptive intervals and not population bounds. (`artifacts/figures/FIG-5_*.pdf`)


**Fig. 6.** Comparison against externally generated evidence. (a) Agreement with published lumped-parameter age estimates under leave-one-study-unit-out folds. The calibrated panel applies a correction fitted to the reference and therefore measures emulation, not independent agreement. Dotted lines bound a factor of two. (b) Directed-topology comparison against a particle-tracking reference; the prior-assisted mode ingests the reference and measures ingestion fidelity only. (`artifacts/figures/FIG-6_*.pdf`)


**Fig. 7.** Field application and its ceiling. (a) In-sample chemistry closure per candidate edge; these are closures against the chemistry used to fit them, not held-out predictions. (b) Number of reaction terms retained per candidate edge, showing that a set of compatible explanations is reported rather than a single mechanism. (c) Cumulative measurement-tier ablation over the 160 Northern Ghana wells: strontium and silica reduce the non-identifiable fraction from 51.9% to 0.6%, whereas isotopes and fluoride do not. (`artifacts/figures/FIG-7_*.pdf`)
