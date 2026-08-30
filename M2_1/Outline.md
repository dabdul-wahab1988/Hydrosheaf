# M2_1 manuscript outline and word contract

## Positioning

M2_1 is a software and methods research article, not a second field-validation
paper and not a universal benchmark claim. It documents the construction of
HydroSheaf, reports what the current package does under a locked synthetic
programme, and demonstrates how its claim gates behave when applied to real,
data-limited Ghanaian measurements. A target journal has not yet been locked, so
generic research-article structure and APA 7 citation keys are used.

## Working title

**HydroSheaf: an auditable evidence-integration framework for groundwater age,
connectivity and reaction inference with controlled-synthetic validation**

## Central question

> Can a modular evidence-integration framework combine groundwater connectivity,
> residence-time, hydrochemical reaction and kinetic evidence while quantifying
> disagreement, calibrating uncertainty, abstaining when evidence is insufficient,
> and making only claims supported by independent synthetic tests and available
> field measurements?

## Novel contribution to test

The novelty is the implemented, auditable validation architecture: independent
candidate generation and specialist comparators; calibrated age and reaction
outputs; RAPM/on-off reaction scoring; explicit model discrepancy; convergent
discrete model averaging; kinetic identifiability reporting; truth-blind
cost-adjusted prospective measurement selection; and a programme-level claim gate
that can issue PASS, WEAK or ABSTAIN. The paper must not describe this as a new
universal groundwater solver or as formal sheaf theory unless the mathematical
requirements of a formal sheaf are separately demonstrated. The preferred wording
is “sheaf-inspired” or “sheaf-style evidence integration”.

## Main-manuscript budget: 6,000 words

The count includes the abstract and short statements section. It excludes title,
keywords, headings, references, tables, figure captions, equations, code/data
availability links and Supplementary Methods. Counts are prose-only.

| Section | Target words | Purpose |
|---|---:|---|
| Abstract | 250 | Problem, implemented contribution, bounded synthetic result, field boundary |
| 1. Introduction | 900 | Scientific gap, groundwater context, novelty and objectives |
| 2. Methods | 1,400 | Architecture, evidence streams, validation programme, field data and reproducibility |
| 3. Results | 1,600 | Package capabilities, synthetic gates, field application and audit outcomes |
| 4. Discussion | 1,400 | Interpretation, novelty, comparator meaning, failures, limitations and way forward |
| 5. Conclusion | 250 | Exact claim and next validation stage |
| 6. Statements | 200 | Data/code availability, ethics, competing interests, funding and author contributions |
| **Total** | **6,000** | |

## Main sections

### Abstract, 250 words

1. State the problem: coupled groundwater inference can become more confident
   without becoming more correct.
2. State the implemented contribution: modular HydroSheaf software with
   calibrated specialist comparators, discrepancy handling, model averaging,
   abstention and decision-utility gates.
3. State the synthetic design: three independent generator families, development
   and locked test cases, truth-blind scoring, negative controls and prospective
   policy evaluation.
4. Report the current bounded results: age coverage 0.9655 with specialist MAE
   4.24 years versus baseline MAE 7.65 years; reaction log loss 0.896 versus
   2.852, coverage 0.739 and false commitment 0.038; kinetics predictive RMSE
   0.119 with conditional identifiability and explicit `k×A` confounding; and
   integrated mean utility per cost 1.478 versus 0.246 for random, with positive
   paired lower bounds.
5. State the field boundary: Ghana data demonstrate workflow operation and
   chemistry-level limitations, not field validation.

### 1. Introduction, 900 words

1. Groundwater interpretation from incomplete hydraulic, tracer and chemical
   evidence.
2. Why ordinary best-fit integration is vulnerable to equifinality, structural
   mismatch and false confidence.
3. Position HydroSheaf relative to MODFLOW/MODPATH, PHREEQC, specialist age
   models, sparse reaction inversion and graph/sheaf-style sensor integration.
4. The gap: previous M2 described a framework, whereas M2_1 must document the
   implemented claim-gated software and its current evidence.
5. Aim and five objectives: architecture, specialist competence, disagreement
   handling, prospective decisions and field scope.
6. Explicit hypotheses and claim boundary.

### 2. Methods, 1,400 words

1. Software architecture and package boundary: core inference, validation,
   comparators, calibration, discrepancy, decision utility and provenance.
2. Evidence streams: topology, age, reaction, kinetics and measurement actions;
   evidence is retained as candidate sets and distributions rather than silently
   collapsed to one path.
3. Specialist comparators: independent age candidate generation with atmospheric
   history and excess-air variants; reaction candidate generation with
   stoichiometric/thermodynamic and RAPM on/off layers.
4. Uncertainty and abstention: calibration, selective risk, classwise reliability,
   unsupported-likelihood abstention, false-commitment control and failure
   records.
5. Model discrepancy and averaging: disagreement is estimated and reported before
   discrete model weights are averaged; convergence is a gate, not an assumption.
6. M8 kinetic adapter and OED: transport/kinetic forward models, predictive
   intervals, parameter recovery and explicit `k×A` confounding.
7. Independent synthetic programme: three generator families, development/locked
   split, held-out mechanisms, missingness/censoring/permutation controls,
   truth-blind selectors and cost-adjusted policy scoring.
8. Field-data application: Northern Ghana, Talensi and Lower Anayari inventory,
   charge-balance/variable-availability audit, field inference units and the
   precise non-validation boundary.
9. Reproducibility and statistical reporting: run ID, source hashes, environment,
   paired bootstrap intervals, no post-hoc gate changes, and no commit claim where
   the manifest did not record one.

### 3. Results, 1,600 words

1. What is now in the package, with a compact module-to-capability table.
2. Synthetic generator and programme audit: independent generators, truth-blind
   rows, source independence, held-out test structure and execution PASS.
3. Age: calibrated coverage, width, selective risk, specialist versus baseline,
   and residual weaknesses.
4. Reaction: RAPM/on-off calibration, log loss, coverage, ECE, selective risk,
   and why this is bounded family inference rather than unique mineral proof.
5. Kinetics: predictive recovery and intervals, conditional identification, and
   the `k×A` structural confounding result.
6. Integrated HydroSheaf: discrepancy calibration, model-averaging convergence,
   false-commitment control, prospective six-case utility and comparator contrasts.
7. Real field data: 320 Northern Ghana seasonal records, 41 Lower Anayari and 63
   Talensi records; missingness/QC and current field results; what the data allow
   and what they do not.
8. Reproducibility result and an explicit list of claims that remain ABSTAIN.

### 4. Discussion, 1,400 words

1. What the synthetic PASS means and does not mean.
2. Why the contribution is a validation architecture and software capability,
   not a claim of universal specialist superiority.
3. Age and reaction: how calibration and abstention convert overconfidence into a
   measurable failure signal; why coverage, selective risk and comparator quality
   matter.
4. Kinetics and model discrepancy: why structural non-identifiability cannot be
   repaired by more time points alone.
5. Integrated decision utility: why prospective, cost-adjusted, truth-blind
   comparisons are stronger than post-hoc uncertainty reduction.
6. Field implications for Ghana: workflow readiness and measurement value, not age
   or connectivity validation.
7. Remaining weaknesses: six locked cases, generator domain, calibration transfer,
   no field data yet, no formal sheaf theorem, no universal superiority, and the
   missing Git commit in the locked run manifest.
8. Way forward: rerun after final code freeze, add more held-out field-like
   generators, external retrospective reference data, prospective field sampling,
   stable release and independent reproduction.

### 5. Conclusion, 250 words

Answer the central question in three sentences: HydroSheaf now has a coherent
implemented software and validation contract; it passes bounded controlled-
synthetic gates with calibrated failure behaviour; and field work remains a
separate validation stage. State the exact safe claim and the next experiment.

### 6. Statements, 200 words

Data/code availability, third-party data provenance, ethics statement, competing
interests, funding, author contributions and the absence of new human sampling.

## Supplementary Methods budget: 3,000 words

| Section | Target words | Content |
|---|---:|---|
| S1. Provenance and data contract | 300 | Paths, versions, source hashes and field-data inventory |
| S2. Evidence representation | 350 | Candidate edges, distributions, units and missingness |
| S3. Age specialist and calibration | 450 | Candidate histories, likelihoods, intervals, abstention |
| S4. Reaction specialist and RAPM | 450 | Templates, features, on/off scores, calibration and ECE |
| S5. Discrepancy and model averaging | 350 | Discrepancy scale, weights, convergence and failure records |
| S6. Kinetics and OED | 350 | Forward model, interval scoring, parameter recovery, `k×A` |
| S7. Synthetic generators and controls | 450 | Independent generators, held-out cases, stresses and policy scoring |
| S8. Field-data processing | 250 | Harmonisation, QC, field inference limits |
| S9. Reproduction and claim gates | 50 | Regeneration and PASS/ABSTAIN rules |
| **Total** | **3,000** | |

## Planned display items

1. Main Figure 1: HydroSheaf architecture and evidence-to-claim pipeline.
2. Main Figure 2: synthetic validation design and truth-blind data flow.
3. Main Figure 3: specialist gate results for age, reaction and kinetics.
4. Main Figure 4: discrepancy/model averaging and prospective utility comparison.
5. Main Figure 5: field-data availability and claim boundary.
6. Main Table 1: package capability inventory and source module.
7. Main Table 2: locked synthetic performance gates.
8. Main Table 3: field-data inventory and supportable claims.
9. Supplementary tables: full thresholds, source hashes, generator families,
   field variable crosswalk and reviewer issue-resolution ledger.

## Claims that must not appear

- “HydroSheaf is the best groundwater inference engine.”
- “HydroSheaf is validated in the field” based on the current FieldData.
- “The field graph is the true flow network.”
- “Stable isotopes in the field data establish groundwater age.”
- “RAPM identifies a unique mineral reaction.”
- “The M8 `k` parameter is identified without independent surface-area evidence.”
- “Sheaf theory has been formally proven by the current implementation.”
- “PASS on six locked synthetic cases establishes universal superiority.”
