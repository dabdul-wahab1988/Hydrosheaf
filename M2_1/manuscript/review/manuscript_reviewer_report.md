# M2_1 manuscript reviewer report

## Review basis and recommendation

This report was prepared after reading the assembled main manuscript, the Supplementary Methods, the claim ledger, the analysis plan, the field-data inventory, the locked synthetic JSON artefacts and the source-hash manifest. It follows the manuscript-reviewer workflow, including section-level critique, cross-section consistency audit, research-integrity and adversarial verification, reporting-standard assessment, reproducibility assessment, data-availability assessment, citation spot checks and rubric scoring. The report records both the first-pass issues and the resolution audit because the manuscript was revised during review.

The manuscript is a software and methods paper. Its defensible contribution is an implemented evidence-to-claim architecture that keeps specialist candidate generation, calibration, discrepancy, selective risk, abstention, model averaging and measurement selection auditable. It is not a claim that HydroSheaf is a universal groundwater-inference engine. After the resolution pass, no critical scientific or integrity issue remained. The recommendation is **Return for minor revisions** before external submission, limited to a clean-commit reproducibility rerun, a public release location, rendered figure delivery and completion of author-specific contribution metadata. Those tasks do not require changing the current scientific claim boundary.

## 1. Title and Abstract

1. **[Minor] The title was initially too easy to read as a general groundwater-inference claim.** The revised title now says “with controlled-synthetic validation” in `Manuscript-Final.md`, line 1, which correctly signals the evidence domain. The abstract also states that no general topology-superiority gate was passed and that the Ghana application was not field validation, line 3. This resolved the main scope risk. Revision instruction: retain the scope-qualified title in the submission metadata and do not shorten it to a title that implies field validation or universal comparison.

2. **[Minor] The abstract remains information-dense.** The single paragraph at line 3 reports age, reaction, kinetics, integrated utility, topology status and three field datasets. The numbers are useful, but the reader must process several metrics before reaching the final claim boundary. The current Table 1 provides a better place for exact precision. Revision instruction: for a journal submission, preserve the current headline values but consider moving one secondary metric, such as the reaction false-commitment rate, to the table or abstract keywords if the target journal imposes a tighter abstract limit.

## 2. Introduction

1. **[Major] The paper needed to distinguish an implemented architectural contribution from a claim of a new inference algorithm.** The initial component list could have been read as a collection of existing methods. The revised text now says that HydroSheaf is “an executable evidence and claim architecture” and explicitly states that it is not a replacement for MODFLOW, MODPATH, TracerLPM or PHREEQC, line 13. The comparator-fairness paragraph and central question at lines 15 to 23 make the contribution testable. Revision instruction: keep this architecture-level contribution statement and do not reintroduce algorithm-superiority language during journal adaptation.

2. **[Major] The age discussion needed a direct boundary between environmental-tracer age inference and stable-isotope source or mixing evidence.** The manuscript now states that a stable isotope can support recharge or mixing plausibility without establishing age, line 7, and later states that the Ghana columns do not contain environmental age tracers, line 23. This is consistent with the USGS TracerLPM description, which treats age distributions as model-dependent interpretations of environmental tracer data and notes that mean age can be non-unique when information is limited. Revision instruction: retain the explicit tracer list and do not describe the Ghana stable-isotope columns as age validation.

3. **[Minor] The literature argument needed a current uncertainty and selective-prediction reference.** The manuscript cited the foundational calibration and discrepancy work, and the methods now cite a modern conformal-prediction review when distinguishing coverage, selection and abstention, line 43. The remaining task is editorial rather than scientific. Revision instruction: preserve the primary sources in `LITERATURE.bib` and add any target-journal-required recent groundwater uncertainty references without replacing the foundational citations.

## 3. Methodology

1. **[Major] Model discrepancy was central to the stated contribution and initially lacked a sufficiently explicit operational formula.** The revised Supplementary Methods now gives the dilation rule `[e-s(e-L), e+s(U-e)]`, identifies the case-blocked order statistic used to fit the scale, and states that the scale is fitted on development data only, Supplementary Methods, lines 61 to 67. The model-averaging objective is also given as a case-equal negative log score under a simplex weight floor. Revision instruction: retain these equations and ensure that any rendered supplement includes the same definitions and the corresponding `discrepancy_calibration.json` and `model_averaging.json` artefacts.

2. **[Major] The first draft transcribed several programme thresholds incorrectly.** The locked contract uses age target coverage 0.95, relative-width ceiling 1.5, reaction coverage floor 0.25 and selective-risk ceiling 0.4. The corrected Methods text reports those values and the conditional kinetic identification rule at line 59. It also states that thresholds are safety floors rather than evidence of optimality. Revision instruction: regenerate all tables and manuscript text from the locked contract after any future threshold change; do not hand-edit these values in a journal draft.

3. **[Major] Comparator reproducibility depended on information outside the main text.** The Methods described independent candidate generation, common observation contracts and truth-blind scoring, but exact baseline identifiers and tuning values were not all visible in the manuscript. This is now addressed by explicit pointers to `baseline_registry.json`, `specialist_comparators.json` and `case_metrics.json` in Supplementary Methods, lines 87 and 101. Revision instruction: include these machine-readable artefacts in the public release and preserve their SHA-256 fingerprints beside the manuscript version.

4. **[Major] The locked run was generated from a dirty worktree and therefore cannot support a clean-rebuild claim.** The manuscript reports the run revision, source hashes and `git_worktree_dirty=true` at line 69, and it does not claim that the current checkout regenerated the result. This is transparent, but it remains a reproducibility limitation. Revision instruction: perform one clean-commit rerun, compare result hashes with `RUN-INTEGRATION-FULL-20260802-15`, and replace the current run identifier in the final submission only if the comparison is documented.

5. **[Minor] The programme plan was local rather than externally preregistered.** The Methods now names `M2_1-PLAN-20260802` and explicitly says that it was not an external preregistration, line 55. This is adequate for a software validation package, but it limits protection against hindsight in a formal journal setting. Revision instruction: archive the analysis plan and its hash in the public release, and use a timestamped repository release or preregistration for the next locked benchmark.

## 4. Results and Discussion

1. **[Major] Six locked cases cannot support a population-level performance claim.** The integrated result reports utility values and paired lower bounds at line 106, but the manuscript now states that the six-case bootstrap intervals are within-run descriptive summaries rather than confidence intervals over all aquifers or generator families. The Discussion also states that generator-domain coverage dominates uncertainty, line 162. Revision instruction: retain the six-case limitation in the abstract-adjacent claim boundary and expand the locked case set before any universal or population-level comparison is attempted.

2. **[Major] Connectivity was present in the software inventory but did not receive a general-superiority result.** The Results explicitly state that no topological superiority result was included and that the prior-assisted diagnostic was not pooled with independent comparisons, line 79. The abstract repeats this boundary at line 3. Revision instruction: keep topology in the manuscript as a capability and component diagnostic until a separate no-prior hydraulic and particle-tracking contract has passed.

3. **[Major] The field chemistry results could be mistaken for connectivity or mechanism validation.** The manuscript now calls the 208 edges candidate-screening constructions, reports their in-sample nature and denies physical-connectivity interpretation, line 114. It also states that the seasonal result was not a chronological forecast and that R² was a fit diagnostic, lines 124 and 126. Revision instruction: retain “candidate edge” and “in-sample closure” in every field-results caption, table and downstream manuscript version.

4. **[Major] The reaction result could be overstated because the pooled gate is easier to pass than a complete mechanism test.** The reported log loss, coverage and false-commitment values are traceable at line 98, while the Discussion now calls the output an improved calibrated family diagnostic, not a complete reaction-inference engine, and identifies the fixed dictionary and remaining calibration limitations. The row-level records and per-family diagnostics are retained in the locked artefacts. Revision instruction: report per-family confusion, calibration and held-out-mechanism results in a future expanded benchmark before describing reaction inference as complete or causal.

5. **[Minor] The model-averaging result needed a clear distinction between constrained convergence and the raw gradient diagnostic.** The Results report KKT and objective criteria as the constrained convergence record while retaining the raw gradient norm of 0.786 against the nominal raw-gradient tolerance, line 108. The Discussion explains why those diagnostics are not interchangeable. Revision instruction: preserve both diagnostics in future reports and never describe the current averaging fit as perfectly optimised.

## 5. Tables and Figures

1. **[Minor] The first draft had no publication-facing result table.** The revised manuscript now contains Table 1 with component status, locked row or case counts, headline metrics and qualification, lines 81 to 87. Revision instruction: render this table in the target journal style and keep the exact values linked to `synthetic_performance_summary.csv` and the locked performance gate.

2. **[Minor, resolved in final package] The architecture figure was initially a source diagram rather than a rendered asset.** Figure 1 is identified at line 90 and its Mermaid source is stored in `manuscript/artifacts/architecture_flow.md`. The final package now includes the self-contained SVG `manuscript/artifacts/architecture_flow.svg` alongside the source file.

## 6. References

1. **[Minor] The bibliography needed a type and metadata audit.** The USGS TracerLPM entry is now typed as a technical report, field-study citations were added to the Methods, and the Sreekanth monitoring-design DOI was added. The citation-key check found no cited key missing from `LITERATURE.bib`. Revision instruction: run the final bibliography through the target journal's CSL or BibTeX processor and verify report titles, DOI casing and access dates.

2. **[Minor] The manuscript needed primary or official support for its comparator descriptions.** The current bibliography includes USGS TracerLPM, MODFLOW 6, MODPATH 7 and PHREEQC documentation, as well as the primary discrepancy and experimental-design sources. Revision instruction: retain these primary sources and avoid replacing them with secondary software summaries when the manuscript is adapted.

3. **[Minor] RAPM could be mistaken for an established external groundwater method.** The Methods now states that RAPM is project-specific notation for the reaction scoring layer, line 37. Revision instruction: preserve that sentence and do not cite or describe RAPM as an established literature method unless a separate method paper supports that claim.

## 7. Editorial and Language Quality

1. **[Minor] Several paragraphs carry many clauses and identifiers.** For example, the locked-run paragraph at line 55 includes generator families, case splits, stress transformations, calibration restrictions and truth release in one unit. It is readable but dense. Revision instruction: split this paragraph into design, leakage control and scoring paragraphs if the target journal permits more space.

2. **[Minor] Technical notation required consistency checking.** The draft previously alternated between `(kA)` and `k×A`, and between `(R^2)` and R². The current assembled manuscript uses `k×A` and R² consistently in the prose and tables. Revision instruction: preserve this notation and apply the same convention to figure labels, supplementary equations and source-code documentation.

3. **[Minor] Author-contribution language is still provisional.** The Statements section says that the final submission should replace the provisional contribution statement with author-specific roles, line 186. Revision instruction: insert named CRediT roles before external submission; this is a metadata task, not a scientific revision.

## 8. Originality, Contribution and Significance

1. **[Major] The paper's contribution is integration and claim control, not a new formal sheaf theorem or a single new groundwater solver.** The manuscript is clear about this at lines 13 and 31, and it supports the claim with executable gates, independent candidate generation, model discrepancy, abstention and prospective utility. Revision instruction: frame the paper for a software, methods or reproducibility venue that values auditable validation, and do not submit it as a universal benchmark of MODFLOW, MODPATH, TracerLPM or PHREEQC.

2. **[Major] The contribution would be marginal if the package merely repeated existing component diagnostics.** The locked programme adds a distinct result by combining competence-matched specialist candidates, calibration, disagreement control, truth-blind policy selection and explicit claim withholding. The Discussion states that this architecture is the contribution and identifies the missing topology and field stages, lines 130 to 162. Revision instruction: keep the architecture diagram, claim ledger and failure-signal examples in the main narrative so the reader can see what is added beyond a collection of established solvers.

3. **[Minor] Significance remains conditional on wider held-out domains.** The current six-case programme supports competence under declared generators, not transfer to arbitrary aquifers. Revision instruction: state the generator-domain boundary in the abstract, conclusion and software metadata until the expanded locked benchmark is complete.

## Cross-Section Consistency Audit

The audit was repeated after the resolution edits. The abstract, Table 1, Results and claim ledger agree on age coverage 0.9655, specialist MAE 4.24 years, baseline MAE 7.65 years, reaction log loss 0.896 versus 2.852, reaction coverage 0.739, kinetic RMSE 0.119, kinetic overall identification 0.167, integrated utility 1.478, random utility 0.246 and specialist utility approximately -0.0039. The exact unrounded values are retained in the locked JSON artefacts.

The field counts agree across the abstract, Methods, Results, Table 2 and field inventory: 320 Northern Ghana records from 160 wells, 41 Lower Anayari records and 63 Talensi records. The Northern Ghana quality categories sum to 320. The 208 current M2 candidate edges are consistently described as 82 Lower Anayari plus 126 Talensi in-sample chemistry closures.

The status distinction is consistent: the component and integrated programme gate is `PASS` within its bounded synthetic contract, while the explicit universal-superiority status is `ABSTAIN_NO_SUPERIORITY_GATE`, and field validation is deferred. The dirty-worktree limitation is stated in Methods, Results, README and the claim ledger. The notation audit found no remaining `(kA)`, `k/A` or `R^2` inconsistency in the assembled manuscript after regeneration.

## Research Integrity Red Flag Scan and adversarial verification

No evidence of fabrication, falsification, image manipulation, citation fabrication or selective deletion of contradictory headline results was found. The dirty-worktree flag is a reproducibility weakness, not evidence of fraud. The manuscript reports negative or limiting outcomes, including no topology-superiority gate, 0.167 overall kinetic parameter identification, the raw-gradient qualification, the weak seasonal advantage over the expanding mean-delta baseline and the absence of field truth.

The adversarial verification ledger was as follows. The age headline values were traced to `performance_gates.json` and reproduced to the reported rounding. The reaction values, including coverage and false commitment, were traced to the same locked gate and the 142-record result. The kinetic RMSE, interval coverage, conditional identification and abstention rate were traced to the M8 benchmark and performance gate. The integrated utility values and lower bounds were traced to `prospective_decision_benchmark.json`; the six-case resampling limitation was checked against its paired case count. Field row counts and missingness were checked against the live CSV and workbook inventory. All five headline claim paths were verified as traceable, while universal superiority and field validation were verified as withheld.

The cooked-work check found no suspiciously perfect field results or unexplained removal of difficult cases. The synthetic results are not all perfect, and the manuscript exposes the dirty snapshot and raw optimisation diagnostic. The appropriate integrity conclusion is **no integrity concern**, with a required clean-commit rerun for reproducibility.

## Reporting Standards Compliance

No single clinical, observational or laboratory reporting standard applies because the manuscript is a computational software-validation study with a secondary audit of existing field records. STROBE is not an appropriate primary standard because no new field sampling or epidemiological inference was performed. The package follows the relevant computational expectations by stating the software version, generator families, development and locked split, truth-blind policy selection, observation stresses, metrics, costs, source hashes and field claim boundary. A public release should add a repository DOI, environment manifest and rendered figure before submission.

## Reproducibility Assessment

Software identity and package boundary were **Present**. The manuscript identifies HydroSheaf 0.5.1, the public pipeline, validation modules and the separation between reusable package code, benchmark scripts, milestone archives and generated artefacts.

Data and execution provenance were **Partially present**. The run ID, source revision, SHA-256 inventory, generator metadata, executable metadata and dirty-worktree flag are recorded. A clean-commit rerun and public release URL remain outstanding.

Validation design was **Present** for the bounded claim. The package records three generator families, 12 development cases, six locked cases, stress controls, held-out scoring, calibration scope, truth release timing and claim gates. Generalisation beyond those generator families is not claimed.

Statistical reporting was **Partially present**. Metrics, calibration, selective risk, false commitment and paired bootstrap settings are recorded, and the six-case resampling limitation is now explicit. Future releases should provide per-generator numerical tables and a larger locked case set.

Execution instructions were **Present** in the Supplementary Methods and generation script, but external execution remains **Partially present** until the clean snapshot is rerun and the results are compared by hash.

## Data Availability and Transparency Assessment

Raw field data availability is **Partially adequate**. The working repository contains the original CSV and XLSX files, their inventory and variable crosswalk. A public access URL, licence or source-study permission record has not yet been assigned.

Analysis-code availability is **Partially adequate**. The package, validation modules, scripts, lock file and M2_1 generation script are present locally. A public versioned release and DOI are still required for independent readers.

Supplementary materials are **Adequate for the current working package**. The supplement includes equations, gate thresholds, discrepancy calibration, model averaging, kinetic identifiability, generator controls, field limits, comparator artefact paths and regeneration instructions.

Author contributions are **Partially adequate** because the Statements section still contains a provisional CRediT placeholder. Named roles must be inserted before submission.

## Citation Integrity Check

Citation keys were checked mechanically and every citation in the assembled main manuscript and supplement resolves to an entry in `LITERATURE.bib`. The bibliography contains no citation that was invented for the report. The USGS TracerLPM record was spot-checked against the official USGS record, which describes multiple lumped-parameter models, environmental tracer interpretation and non-uniqueness under limited information. Kennedy and O'Hagan was checked against the journal record, which describes calibration uncertainty and model discrepancy. Oreskes, Shrader-Frechette and Belitz was checked against the Science abstract, which states that confirmation of natural-system models is partial and model results are non-unique. MODFLOW 6, PHREEQC and Bayesian experimental design records were also checked against USGS, official documentation or the original journal record.

The literature sources are used for methodological context rather than as evidence that HydroSheaf itself has been validated. No citation cartel, excessive self-citation, source misrepresentation or salami-slicing signal was found. The only remaining citation task is journal-format conversion and final access-date checking.

## Scoring

| Section | Score | Maximum | Basis |
|---|---:|---:|---|
| Title and Abstract | 4 | 5 | Scope is now explicit; the abstract remains dense. |
| Introduction | 8 | 10 | Clear problem and boundary; contribution is architectural and literature-supported. |
| Methodology | 21 | 25 | Strong contract and artefact trail; clean rerun and public release remain. |
| Results and Discussion | 25 | 30 | Traceable and candid; six locked cases and absent topology gate limit reach. |
| Tables and Figures | 4 | 5 | Tables and figure source are present; rendered figure remains. |
| References | 4 | 5 | Primary sources and metadata are good; final journal formatting remains. |
| Editorial and Language Quality | 4 | 5 | Clear UK English with some dense paragraphs. |
| Originality, Contribution and Significance | 12 | 15 | A real architecture-level contribution, with conditional transfer and no formal theorem claim. |
| **Total** | **82** | **100** | **Minor revisions** under the rubric. |

## Critical Flaw Assessment

No critical flaw was identified. The methodology is appropriate for the bounded question because the synthetic truth is independent, locked scoring is separated from development calibration, and field records are not treated as truth. The six-case design and dirty-worktree provenance limit the strength of generalisation and reproduction, but they do not invalidate the reported controlled-synthetic results and are addressable by a clean rerun and expanded cases. The conclusions follow the evidence because the manuscript withholds topology superiority, field validation, unique reactions and universal engine superiority. The contribution is meaningful as an implemented and auditable validation architecture, even though it is not a formal sheaf theorem or a replacement for specialist forward models. No fabrication or citation-integrity override was triggered.

## Decision

The re-reviewed package is scientifically suitable for submission after minor revisions. The remaining actions are a clean-commit regeneration and hash comparison, public release location and access metadata, and named author contributions. Figure 1 has now been rendered as a self-contained SVG in the final package. No redesign of the current experiments or weakening of the bounded claim is required.

## Structured recommendation answers

1. **Recommendation:** Return for minor revisions. The core manuscript and claim boundary are acceptable; release and metadata tasks remain.

2. **Study design appropriateness:** Yes for the declared controlled-synthetic question. The design is not a field-validation design, and the manuscript correctly says so.

3. **Methods reproducibility:** No, but the remaining clean-rerun, public-release and rendering points can be addressed without changing the study design.

4. **Statistics and uncertainty treatment:** Yes, appropriate for a bounded software benchmark, provided the six-case intervals remain labelled descriptive rather than population-level.

5. **Guidance on overstated claims:** Yes. The report and manuscript specify the wording for topology, field validation, reaction mechanism, kinetic parameters and universal superiority.

6. **Presentation clarity:** Yes. The main narrative is clear, with a few dense paragraphs that merit light editing.

7. **Editorial and language quality:** Generally well-written but needs a light editing pass before publication.

8. **Originality and contribution:** Adequate contribution with minor gaps. The architecture-level contribution is real, while broader transfer and topology evidence remain pending.

9. **Reproducibility:** Partially reproducible, with the dirty-worktree snapshot and public release still to be resolved.

10. **Research integrity:** No concerns. The adversarial scan found traceable metrics, candid limitations and no integrity red flags.

11. **Data availability and transparency:** Partially transparent with notable gaps in public release, access conditions and author metadata.

12. **Citation integrity:** No concerns. Citation keys resolve and spot checks supported the cited methodological claims.

## Summary letter to the authors

M2_1 now reads as a software and methods manuscript rather than as an inflated field-validation paper. Its strongest feature is the explicit separation between what the current synthetic programme demonstrates and what the Ghana measurements cannot yet identify. The age, reaction-family, kinetic and integrated results are tied to a locked run, specialist baselines, calibration, abstention and machine-readable evidence. The manuscript also reports limiting outcomes instead of hiding them, including conditional kinetic identification, the absence of a topology-superiority gate, the six-case decision benchmark and the dirty-worktree provenance limitation.

The main scientific contribution is the implemented claim-gated evidence architecture. That is a legitimate contribution, but it must remain the contribution. The paper should not imply that RAPM is an established external method, that reaction-family probabilities identify unique mineral processes, that candidate chemistry edges are measured flow paths, or that positive six-case utility establishes a universal groundwater engine. The current text handles these boundaries well.

Before external submission, complete the clean-commit rerun and compare hashes, publish the code and synthetic artefacts at a stable URL or DOI, insert author-specific CRediT roles and perform the target-journal formatting pass. The architecture figure is already rendered in the final package. With those steps completed, M2_1 will be a defensible account of how HydroSheaf was built and what its current software can support on controlled synthetic data, while preserving the correct decision to defer field validation.

## Final reader-facing terminology pass

The first review did not separately test whether a reader without access to the project workspace could interpret internal labels. A final terminology audit found internal milestone labels, repository paths, a hidden run identifier, raw revision hashes, machine status strings and several unexplained abbreviations in the assembled text. These were removed or rewritten. PASS, WEAK and ABSTAIN are now defined in the Introduction; RAPM, KKT, R², RMSE, MAE and k×A are defined at first use or written in full; and the field application is described as an archive-based data audit rather than by internal directory names. A pattern search over the final manuscript and supplement found no remaining Codex paths, milestone labels, raw run identifiers, internal status strings or repository-specific source paths. This issue is resolved.
