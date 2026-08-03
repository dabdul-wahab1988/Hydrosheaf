# Re-Review Report: Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks

## Part A: Comment-by-comment assessment

### M7-20260729-R01

**Original comment.**

1. **[Major]** The phrase “single-generator groundwater benchmark” is not accurate for the assembled paper. The M7.3 experiments used the MODFLOW 6, MODPATH 7, tracer and chemistry generator described at main-text lines 221 to 242, whereas M7.4 used a separate scalar graph-case generator that imported no HydroSheaf code, as stated at lines 348 to 358 and Supplementary lines 352 to 361. The title, Introduction line 203 and Limitations lines 1134 to 1136 therefore collapse two distinct generator systems and three different locked-test sizes into one. This matters because the claimed scope and independence of the evidence depend on which generator produced which result. Revise the title to refer to “controlled-synthetic benchmarks” or “two controlled-synthetic generator systems,” then identify M7.3 and M7.4/M7.5 separately in the abstract and limitations.

**Authors' claimed response.**

The title now says controlled-synthetic benchmarks. The abstract, Introduction, Methods and Limitations name the MODFLOW 6/MODPATH 7 M7.3 generator separately from the scalar affine M7.4/M7.5 generator and prohibit cross-generator or field transfer of either result.

**What was found in the revised manuscript.**

The title now reads 'Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks.' The aim names 'two independent controlled-synthetic generator systems,' and the Methods has a dedicated non-transfer subsection.

**Assessment: ADEQUATELY ADDRESSED.**

The revision fixes both the title-level error and the deeper scope problem. A reader can now tell which generator supports each claim.

### M7-20260729-R02

**Original comment.**

2. **[Minor]** The abstract reports the important positive, null and adverse findings, but it compresses seven audits and several estimands into one dense paragraph. At lines 8 to 26, the evidence-panel results, age results, reaction results, M7.4 representation test and M7.5 estimator test arrive without an explicit primary-versus-secondary hierarchy. A reader can miss that M7.4 rejected overall benefit against the claim-bearing comparator and that M7.5 also failed its complete three-outcome gate. Divide the abstract into purpose, design, primary representation result, follow-up result and scope conclusion, and use the words “failed the prespecified complete gate” for both M7.4 and M7.5.

**Authors' claimed response.**

The abstract was split into purpose, design, primary representation result, follow-up estimator result and scope conclusion. It now states explicitly that both M7.4 and M7.5 failed their prespecified complete gates.

**What was found in the revised manuscript.**

The abstract uses five labelled paragraphs. The primary paragraph states that M7.4 'failed the prespecified complete gate,' and the follow-up paragraph uses the same wording for M7.5.

**Assessment: ADEQUATELY ADDRESSED.**

The hierarchy and gate failures are now unmistakable in the abstract.

### M7-20260729-R03

**Original comment.**

1. **[Major]** The originality claim at main-text lines 102 to 109 is not supported by a reproducible literature search. `manuscript/methods/literature_search.json` records general web queries but no bibliographic databases, search dates by database, returned counts, deduplication, inclusion criteria or screening record for the claim that no groundwater benchmark combined an external generator, a locked split and adverse controls. The claim may be correct, but the current evidence cannot establish its completeness. Either replace it with a narrower statement that the authors did not identify such a benchmark in a targeted, non-systematic search, or archive a reproducible search across at least Web of Science or Scopus, Crossref, Google Scholar and a preprint index, with queries, dates, counts and screened records.

**Authors' claimed response.**

The unsupported exhaustive novelty claim was withdrawn. The manuscript now describes a targeted, non-systematic search, and the search record stores dated queries, sources, inclusion logic and screened items without claiming bibliographic completeness.

**What was found in the revised manuscript.**

The Introduction now says 'In a targeted, non-systematic search' and names the searched source types. The versioned search JSON records the 30 July 2026 targeted queries and expressly disclaims an exhaustive absence claim.

**Assessment: ADEQUATELY ADDRESSED.**

The authors selected the reviewer-approved claim-restriction route and no longer imply a systematic or exhaustive review.

### M7-20260729-R04

**Original comment.**

2. **[Major]** The objective structure obscures which experiments were confirmatory. Main-text lines 159 to 190 list five objectives, but objective five contains M7.4 and a second fresh-seed diagnostic, while the manuscript later describes seven linked audits. M7.3, the public-pipeline run, M7.4 and M7.5 also used different generators, case counts, primary outcomes and claim gates. This makes it difficult to distinguish the original scientific questions from later estimator diagnosis. Add a one-page design table at the end of the Introduction that maps each audit to its generator, development and test counts, primary comparator, primary metrics, lock point, multiplicity family and permitted claim.

**Authors' claimed response.**

A main-text design-and-claim table now maps all seven audits to generator, development/test counts, comparator, primary metrics, lock/multiplicity family and permitted claim.

**What was found in the revised manuscript.**

Main Table 1 contains seven audit rows and the requested columns for generator, development/test counts, comparator, metrics, lock or multiplicity family, and permitted claim.

**Assessment: ADEQUATELY ADDRESSED.**

The new table contains every requested design field and resolves the confirmatory-versus-diagnostic ambiguity.

### M7-20260729-R05

**Original comment.**

3. **[Minor]** The practical significance is framed more strongly than the decision analysis supports. Lines 198 to 203 acknowledge that translation to a field decision protocol remains future work, yet the surrounding text says the benchmark identifies which evidence was “worth collecting.” No value-of-information threshold, cost function or field decision loss was defined. Replace “worth collecting” with “improved the prespecified synthetic-benchmark outcomes,” and reserve collection recommendations for a future study that defines measurement cost and decision consequences.

**Authors' claimed response.**

The phrase 'worth collecting' was removed. Collection recommendations are reserved for a future value-of-information study with measurement costs and decision losses.

**What was found in the revised manuscript.**

The phrase 'worth collecting' is absent. The revised text says the experiments 'do not establish a value-of-information threshold or a measurement-cost recommendation.'

**Assessment: ADEQUATELY ADDRESSED.**

The practical claim is now aligned with the absence of measurement costs and field decision loss.

### M7-20260729-R06

**Original comment.**

1. **[Critical]** The current M7.5 runner does not match the confirmatory lock. `m7_robust_hybrid_confirmatory.lock.json`, line 6, requires SHA-256 `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`; the current `scripts/run_m7_robust_hybrid_sheaf.py` hashes to `994e36954775c0577dd6f7a7655e9cf6d90d4cec2463a2cf4fc19891d8af7c12`. Direct execution stops at source line 150 with `Confirmatory-bound artifact changed: runner_sha256`. Stored output hashes remain intact, but the exact code that passed the one-time lock is not recoverable from the present checkout. Recover the exact locked runner from backup or an earlier task snapshot and archive it under its content hash. If it cannot be recovered, state that the stored outputs are an unverifiable locked-run record, create a new protocol and confirmatory lock, and run one new test set with previously unused seeds.

**Authors' claimed response.**

The exact one-time M7.5 runner was recovered from the archived task session. Its SHA-256 is a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191, exactly matching the confirmatory lock. The source and recovery manifest are archived under that content hash. No locked test was rerun.

**What was found in the revised manuscript.**

The active runner and its archived copy both hash to a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191. This equals the confirmatory-lock value. The recovery manifest identifies the archived task-session source and records that no locked test was rerun.

**Assessment: ADEQUATELY ADDRESSED.**

The exact locked source is recoverable and content-addressed, while the locked outputs remain untouched.

### M7-20260729-R07

**Original comment.**

2. **[Critical]** The software availability statement is factually wrong for M7.4 and M7.5. Main-text lines 1249 to 1259 say the later protocols, runners, tests and result generators are held in the repository at protocol-freeze commit `d336e87`. Both M7.4 and M7.5 manifests instead record Git revision `53beb46034d5230c1a061341a5cf2175d9af858e`, and the M7.4/M7.5 protocol, runner and result files are absent even from that commit and from the current remote `main`. The recorded revision therefore cannot recreate either later experiment. Commit the complete M7.4/M7.5 protocol, code, tests, manifests and immutable outputs, tag a release, deposit that release in a persistent repository, and cite separate M7.3 and M7.4/M7.5 freeze identifiers in the manuscript.

**Authors' claimed response.**

The false commit-level reconstruction claim was removed and replaced with experiment-specific run identifiers, recorded historical revisions, current file hashes and an explicit statement that the historical commits do not contain M7.4/M7.5. The complete local package is assembled. A versioned commit/tag and persistent DOI remain submission actions because this revision task does not authorise repository publication or create a repository DOI.

**What was found in the revised manuscript.**

The availability statement now correctly says that commit 53beb460... does not contain the later files and is not their release identifier. The complete files are present locally with run IDs and hashes, but the text also states that the new commit, release, and persistent DOI remain required.

**Assessment: PARTIALLY ADDRESSED.**

The factual availability statement and local technical package are corrected, but the comment also required a committed, tagged, persistently deposited release. That external release has not yet been created, so a reader cannot cite or retrieve an immutable public M7.4/M7.5 package.

**Guidance for further revision.**

Create one versioned repository commit containing the complete M7.4/M7.5 protocols, runners, tests, manifests, and immutable outputs. Tag the release, deposit it in a persistent repository, obtain the software and data identifiers, and replace the pending language with those exact identifiers before submission.

### M7-20260729-R08

**Original comment.**

3. **[Major]** The paper needs to separate the construct validity of its two generator systems. M7.3 used a flow, tracer and chemistry generator, whereas M7.4/M7.5 used scalar affine graph cases in four planted scenarios. Both are legitimate capability tests, but success in the scalar generator does not validate performance in the MODFLOW-based generator, and neither validates a field aquifer. Add a methods subsection that names the two generator systems, states which objectives each can test, and defines the non-transfer boundary between them.

**Authors' claimed response.**

A construct-validity subsection distinguishes the two generators, the questions each can test and the explicit non-transfer boundary: scalar-case success is not MODFLOW-system validation, and neither generator constitutes field validation.

**What was found in the revised manuscript.**

The Methods subsection 'Generator construct validity and non-transfer boundary' states what M7.3 can test, what the scalar M7.4/M7.5 generator can test, and that neither transfers to the other or to field validity.

**Assessment: ADEQUATELY ADDRESSED.**

The non-transfer boundary is explicit and repeated where the results are interpreted.

### M7-20260729-R09

**Original comment.**

4. **[Major]** No prospective sample-size or precision rationale is reported for six development cases and twelve locked M7.3 cases, 64 M7.4 test cases or 128 M7.5 test cases. Bootstrap intervals quantify the realized uncertainty but do not show why these sample sizes were adequate to detect a scientifically relevant difference. This is especially important because several findings rest near zero and the M7.3 exact tests are discrete at 12 blocks. Add a simulation-based precision analysis using development-only effect distributions, state the minimum difference of interest for each primary metric, and report the probability that each locked design would exclude the prespecified null or non-inferiority margin.

**Authors' claimed response.**

A labelled post-review simulation study used development-only planning inputs, 20,000 replicates and prespecified minimum differences for PR-AUC, Brier, log loss, age MAE, interval width, coverage and modal accuracy. The attainable precision/power results are reported in Supplementary Table S12 and are not presented as prospective preregistration.

**What was found in the revised manuscript.**

The revision openly states that no practical margins or prospective power analysis were prespecified. It then labels a 20,000-replicate simulation as post-review planning and reports probabilities of excluding zero and clearing the planning margins in the main Results and Table S12.

**Assessment: ADEQUATELY ADDRESSED.**

The requested planning analysis is quantitative and is correctly labelled as post-review rather than retrospective preregistration.

### M7-20260729-R10

**Original comment.**

5. **[Major]** M7.3 received an explicit 24-test Benjamini-Hochberg analysis, but the M7.4 and M7.5 scenario claims rely on unadjusted percentile intervals drawn from much larger contrast matrices. Supplementary Table S7 contains 120 rows and Table S9 contains 560 rows. Prespecified directions reduce researcher discretion but do not remove family-wise selection risk when multiple scenarios and metrics are interpreted. Define the exact confirmatory family for M7.4 and M7.5, then provide simultaneous bootstrap intervals or adjusted permutation tests for the scenario statements that are retained in the abstract, results and conclusion.

**Authors' claimed response.**

All 120 published M7.4 contrasts and all 560 published M7.5 contrasts were placed in separate full families. Shared case-block bootstrap resampling with 10,000 replicates and max-z simultaneous 95% intervals was applied. Only findings surviving those families remain inferentially supported.

**What was found in the revised manuscript.**

The Methods defines separate complete families of 120 M7.4 and 560 M7.5 contrasts, with 10,000 shared case-block resamples and two-sided max-z simultaneous 95% intervals. Tables S10 and S11 provide the adjusted results.

**Assessment: ADEQUATELY ADDRESSED.**

The multiplicity problem is resolved for the entire published contrast families, not only a selected subset.

### M7-20260729-R11

**Original comment.**

6. **[Major]** The reaction-family experiment uses the same six-family vocabulary for generation and scoring. Supplementary lines 172 to 186 show that predicted reactions were mapped to the six planted process families, while the generator planted those same archetypes. Code independence prevents leakage, but a shared ontology can inflate recovery for well-separated families and does not test out-of-dictionary chemistry. Add a mechanism-mismatch sensitivity that perturbs stoichiometry, mineral assemblage and reaction combinations outside the scoring dictionary, or restrict the claim to discrimination among the six planted archetypes.

**Authors' claimed response.**

No unplanned mechanism-mismatch experiment was added. Instead, every reaction claim is restricted to discrimination among the six planted archetypes, under the two tested indicator panels and the tested noise model; out-of-dictionary chemistry is explicitly unevaluated.

**What was found in the revised manuscript.**

The paper now defines the reaction estimand as discrimination among the six planted archetypes under the specified stoichiometry, mineral assemblage, noise model, and two panels. It explicitly says that out-of-dictionary reactions were not tested.

**Assessment: ADEQUATELY ADDRESSED.**

The narrower claim is scientifically defensible and directly follows the alternative offered by the reviewer.

### M7-20260729-R12

**Original comment.**

7. **[Minor]** The Supplementary Methods describe unweighted logistic regression as “the only choice consistent with reporting calibration honestly.” That statement is too absolute. Weighting changes the fitted target and can harm raw calibration, but weighted estimators can still be assessed and recalibrated on an untouched development fold. Rewrite this as a prespecified design choice, report the class prevalence in each generator, and explain why the unweighted target matched the intended probability estimand.

**Authors' claimed response.**

The absolute statement about unweighted logistic regression was replaced with a prespecified estimand rationale. Candidate-edge class prevalence is reported for each generator, and the text explains why the unweighted target matches the intended per-candidate probability estimand.

**What was found in the revised manuscript.**

The Methods reports M7.3 prevalence as 54/414 in development and 108/827 in locked test, and one-third for the scalar generator. It describes unweighted fitting as a prespecified choice for the generated per-candidate probability estimand, not as the only honest method.

**Assessment: ADEQUATELY ADDRESSED.**

The estimand, prevalence, and design rationale are now precise without making an absolute methodological claim.

### M7-20260729-R13

**Original comment.**

1. **[Major]** The statement at main-text lines 856 to 865 that M7.5 passed its execution and provenance gate cannot remain unqualified while the current runner fails the confirmatory hash. The stored case tables support the numerical contrasts, but the present archive does not support the stronger provenance statement. Suspend that sentence until the exact runner is restored, or revise it to state that the stored manifest reports a pass but the current source archive does not verify the runner hash.

**Authors' claimed response.**

Because the exact runner was recovered and its hash verified against the lock, the provenance statement now names the recovered hash and source archive. It also states explicitly that the stored locked test was not rerun.

**What was found in the revised manuscript.**

The Results reports the restored runner hash, states that it matches the confirmatory lock, and says the stored locked test was not rerun. This is supported by the content-addressed archive manifest.

**Assessment: ADEQUATELY ADDRESSED.**

The provenance sentence is now independently verifiable from the restored bytes and lock file.

### M7-20260729-R14

**Original comment.**

2. **[Major]** Statistical precision is not paired with practical magnitude. Complete topology reduced mean absolute age error by 0.062 years relative to a 2.764-year baseline and by 0.164 years relative to a 4.750-year baseline, approximately 2.2% and 3.5%. Interval-width reductions were also small relative to the baseline widths. Similarly, the M7.5 PR-AUC difference of 0.0200 is about 3.1% of the edge-local mean, while the complete calibration gate failed. Add relative effects, baseline dispersion and a prespecified minimum practically important difference to the Results, then state whether each result cleared that threshold.

**Authors' claimed response.**

Results now pair absolute effects with relative changes and post-review practical margins. The small topology-age changes and the M7.5 overall PR-AUC change are explicitly classified by whether they clear those margins; the M7.5 complete calibration gate remains failed.

**What was found in the revised manuscript.**

The Results now reports relative effects and margins. For example, the M7.5 PR-AUC mean difference is given as 3.09% of the edge-local mean and is said not to exceed the post-review 0.02 planning margin. The topology-age effects are treated the same way.

**Assessment: ADEQUATELY ADDRESSED.**

Statistical intervals are now paired with relative magnitude and explicit, honestly post-review planning margins.

### M7-20260729-R15

**Original comment.**

3. **[Major]** The incompatible-cycle and noisy/missing findings at lines 892 to 900 are presented as supported scenario gains without multiplicity control across the full scenario-metric matrix. These are useful diagnostic observations, but their confirmatory status is weaker than the prose implies. Apply the family correction requested under Methodology, or relabel these findings as prespecified exploratory diagnostics and remove them from the abstract until corrected intervals are available.

**Authors' claimed response.**

Scenario statements were reclassified using the simultaneous full-family intervals. The incompatible-cycle conflict-localisation signal survives; the noisy/missing overall gain and M7.5 scenario ranking-gain claims do not and were withdrawn as supported general gains.

**What was found in the revised manuscript.**

After the full-family correction, the paper retains the incompatible-cycle conflict-localisation result but withdraws support for noisy/missing overall gains and M7.5 scenario ranking gains. The abstract and conclusion follow those corrected decisions.

**Assessment: ADEQUATELY ADDRESSED.**

The manuscript's supported scenario claims now match the simultaneous intervals.

### M7-20260729-R16

**Original comment.**

4. **[Major]** The heading “Carbonate reactions remain non-identifiable regardless of panel richness” and related discussion extend beyond the tested panels. The study compared one core and one enhanced indicator set under one reaction dictionary and one noise model. The zero recovery is a valid finding for those conditions, but it does not establish non-identifiability regardless of isotopic, mineralogical or thermodynamic measurements. Rename the subsection “Carbonate reactions were not recovered under either tested indicator panel” and replace every broader formulation with the same boundary.

**Authors' claimed response.**

The heading now reads 'Carbonate reactions were not recovered under either tested indicator panel.' All broader 'regardless of panel richness' or universal non-identifiability language was removed.

**What was found in the revised manuscript.**

The Results and Discussion heading is now exactly 'Carbonate reactions were not recovered under either tested indicator panel.' The supporting paragraph confines the result to six planted archetypes, 3% noise, and the two tested panels.

**Assessment: ADEQUATELY ADDRESSED.**

The universal non-identifiability claim has been fully withdrawn.

### M7-20260729-R17

**Original comment.**

5. **[Major]** Figure 5 mixes the M7 Ghana audit with a companion M6 field-transfer experiment. The caption discloses this, but placing the M6 tier ablation inside an M7 result figure invites readers to treat all four panels as evidence from the same protocol. Move the M6 panel to a clearly labelled contextual figure in the Supplement, or present it in a boxed comparison with its own protocol, sample, outcome and citation. The M7 field claim must rest only on the truth-free Ghana scope audit and hold-forward analysis.

**Authors' claimed response.**

Figure 5 was rebuilt as an M7-only Ghana supportability figure. The M6 tier-ablation panel was removed; the field claim now rests only on the truth-free Ghana scope and hold-forward audit.

**What was found in the revised manuscript.**

Figure 5 was regenerated as 'Northern Ghana evidence and claim boundary using M7 evidence only.' The M6 panel is absent; the four panels now cover evidence availability, defensible claims, workbook coverage, and the truth-free hold-forward audit.

**Assessment: ADEQUATELY ADDRESSED.**

The M7 field figure no longer imports evidence from another module.

### M7-20260729-R18

**Original comment.**

6. **[Minor]** The public-pipeline table reports the same selected F1 value, 0.4222, for all four arms, despite large changes in probability scores and log loss. Candidate recall is also 0.9815 rather than complete. These observations suggest that the selected threshold was insensitive in six cases and that the system result is conditional on candidate generation. Add thresholds and confusion counts for each arm, explain the identical F1 result, and state explicitly that the system audit tested scoring only on the candidates that were recovered.

**Authors' claimed response.**

The public-pipeline audit now reports the threshold and TP/FP/FN counts for every arm, explains the identical selected F1 values, distinguishes conditional from end-to-end recall and states that all generated candidates were selected, so the audit did not identify a useful scalar selection threshold.

**What was found in the revised manuscript.**

The Results says all 198 generated candidates were retained in every arm, with 53 TP, 145 FP, and 0 conditional FN. It reports conditional F1 0.4223, end-to-end F1 0.4206, one pre-scoring missed truth edge, and recall 0.9815. Table S13 gives the arm-level thresholds and counts.

**Assessment: ADEQUATELY ADDRESSED.**

The identical F1 values and candidate-generation conditioning are now fully explained.

### M7-20260729-R19

**Original comment.**

7. **[Minor]** The conclusion says the sheaf layer “has therefore earned a bounded scientific-workflow role” at lines 1203 to 1205. “Earned” is evaluative language rather than a test result. Replace it with: “The experiments support use of the sheaf layer as a prespecified model of non-identity relations and as a global-compatibility diagnostic under the tested scalar scenarios, with global fallback when endpoint evidence is missing.”

**Authors' claimed response.**

The evaluative 'earned' sentence was replaced verbatim with the requested bounded statement about non-identity relations, global-compatibility diagnosis and missing-endpoint fallback under the tested scalar scenarios.

**What was found in the revised manuscript.**

The Conclusion contains the requested sentence stating that the experiments support the sheaf layer as a prespecified model of non-identity relations and a global-compatibility diagnostic under the tested scalar scenarios, with global fallback for missing endpoints.

**Assessment: ADEQUATELY ADDRESSED.**

The evaluative language was replaced with a test-bounded conclusion.

### M7-20260729-R20

**Original comment.**

1. **[Major]** Figure 5 combines data and claims from two modules, as described above. Its provenance is disclosed in the caption and artifact registry, but the visual grouping still implies a common experiment. Separate the M6 panel or mark it inside the panel title and legend as external companion evidence, including the independent sample size and protocol identifier.

**Authors' claimed response.**

The combined-module visual was removed. Figure 5 and its caption now contain only M7 evidence and explicitly distinguish synthetic supportability context from the truth-free Ghana audit.

**What was found in the revised manuscript.**

The revised Figure 5 is M7-only and the figure-source manifest identifies it as an M7 field-supportability boundary. There is no M6 evidence inside the result figure.

**Assessment: ADEQUATELY ADDRESSED.**

The visual and its provenance now describe one M7 audit only.

### M7-20260729-R21

**Original comment.**

2. **[Minor]** The Word render is clean, but the main paper contains seven figures and nine tables in 30 pages. Table 7 spans three pages and interrupts the transition to the representation results. Retain the design table, the primary M7.3 decision table, Table 8 and Table 9 in the main paper, and move detailed metric tables to the Supplement, where the complete CSVs are already cited.

**Authors' claimed response.**

The main paper now has four tables: the design map, a compact M7.3 decision table, M7.4 means and M7.5 means. Detailed metric tables were moved to the 13-table Supplement.

**What was found in the revised manuscript.**

The revised main manuscript contains seven figures and four tables. The design map, compact M7.3 decision table, M7.4 means, and M7.5 means remain in the main text; detailed results are in thirteen supplementary tables.

**Assessment: ADEQUATELY ADDRESSED.**

The main-text table burden is reduced exactly as requested without removing complete supplementary results.

### M7-20260729-R22

**Original comment.**

3. **[Minor]** Figures 6 and 7 are legible at full-page width but their four-panel labels and several axis labels become small in the 30-page Word layout. Increase label and tick sizes for final journal dimensions, and include the intended printed width in the figure-generation manifest so that legibility can be checked automatically.

**Authors' claimed response.**

Figures 6 and 7 were regenerated for 7.08-inch journal width with a minimum 8-point label size. Those dimensions are recorded in the figure-source manifest and were checked in the Word and LibreOffice renders.

**What was found in the revised manuscript.**

The revised figure-source manifest records 7.08 inches and a minimum 8 point label size for Figures 6 and 7. Both figures are legible in the Word and LibreOffice renderings.

**Assessment: ADEQUATELY ADDRESSED.**

Final-size legibility is now measurable and documented.

### M7-20260729-R23

**Original comment.**

1. **[Major]** Eighteen unique references are too few to support the paper’s broad positioning across groundwater joint inversion, age dating, reactive transport, probabilistic calibration, graph inference and sheaf methods. The present search record is strongest for software documentation and selected method citations but weak for competing groundwater graph or structured-residual approaches. Expand the literature review using the reproducible search requested above, and add the closest non-sheaf alternatives rather than only mathematical sheaf sources and one recent neural-network preprint.

**Authors' claimed response.**

The Introduction and Discussion now include recent 2022-2025 tracer, hydrochemistry, groundwater-age and machine-learning work, plus non-sheaf structured alternatives such as graph regularisation, Gaussian-process smoothing, flow-network inversion and residual diagnostics. The search is disclosed as targeted rather than exhaustive.

**What was found in the revised manuscript.**

The Introduction and Discussion add recent groundwater tracer, hydrochemistry, age-model, and machine-learning sources from 2022 to 2025. They also discuss non-sheaf alternatives, including graph regularisation, Gaussian-process smoothing, flow-network inversion, and residual diagnostics.

**Assessment: ADEQUATELY ADDRESSED.**

The positioning is materially broader and more current while remaining honest about the targeted search scope.

### M7-20260729-R24

**Original comment.**

2. **[Minor]** Davis and Goadrich are cited at main-text lines 416 to 420 and Supplementary lines 231 to 236 as showing that PR-AUC is more informative under imbalance. The paper’s primary statement concerns precision-recall curves under highly skewed data and the relation between ROC and PR spaces, not an unconditional ranking of scalar AUC summaries. Change the wording to “precision-recall curves can give a more informative view under class imbalance,” and retain ROC-AUC as the manuscript already does.

**Authors' claimed response.**

The Davis-and-Goadrich wording now states that precision-recall curves can give a more informative view under class imbalance; ROC-AUC remains reported.

**What was found in the revised manuscript.**

The Methods and Supplement now say that precision-recall curves 'can give a more informative view under class imbalance' and explicitly retain ROC-AUC rather than claiming unconditional scalar superiority.

**Assessment: ADEQUATELY ADDRESSED.**

The citation now supports the precise claim made.

### M7-20260729-R25

**Original comment.**

1. **[Major]** The declarations are submission blockers. Main-text lines 1213 to 1222 contain placeholders for author contributions and competing interests; lines 1233 to 1240 and 1257 to 1260 contain data and software DOI placeholders. There is no funding statement. Replace every placeholder with final contributor roles, competing-interest and funding declarations, dataset DOI, software DOI, licence and versioned release before submission.

**Authors' claimed response.**

The section is correctly structured and the licence and available technical identifiers are stated. Author CRediT roles, funding, competing interests, dataset DOI, software DOI and final versioned release are retained as explicit submission blockers because they require author declarations or external deposits and must not be fabricated.

**What was found in the revised manuscript.**

The declarations are correctly ordered and the MIT licence is stated, but author roles, funding, competing interests, dataset DOI, software DOI, and final release remain explicit pre-submission requirements. These details are not yet supplied.

**Assessment: PARTIALLY ADDRESSED.**

The revision appropriately refuses to invent declarations, but the original submission blockers remain unresolved. Explicit notes that metadata are missing are accurate interim text, not final declarations.

**Guidance for further revision.**

The named authors must supply a final CRediT statement, funding declaration, and competing-interest declaration. After the repository deposits are made, insert the dataset DOI, software DOI, licence, and versioned release, then remove the interim instructions to authors.

### M7-20260729-R26

**Original comment.**

2. **[Minor]** Main-text line 318 reads, “The reaction solver was evaluated inference, not part of the synthetic generator.” The sentence is ungrammatical and obscures an important independence claim. Replace it with: “The reaction solver was the inference method under evaluation and was not part of the synthetic generator.”

**Authors' claimed response.**

The sentence was replaced with: 'The reaction solver was the inference method under evaluation and was not part of the synthetic generator.'

**What was found in the revised manuscript.**

The Methods now says: 'The reaction solver was the inference method under evaluation and was not part of the synthetic generator.'

**Assessment: ADEQUATELY ADDRESSED.**

The grammatical error and the inference-generator independence statement are both corrected.

### M7-20260729-R27

**Original comment.**

3. **[Minor]** “Open Research” at line 1211 contains authorship, interests, data and code statements rather than only open-research material. Rename the section “Declarations and open research,” and order it as author contributions, funding, competing interests, data availability and code availability.

**Authors' claimed response.**

The section is now titled 'Declarations and open research' and is ordered as author contributions, funding, competing interests, data availability and code availability. Only the author-supplied metadata listed under R25 remains open.

**What was found in the revised manuscript.**

The section is titled 'Declarations and open research' and appears in the requested order: author contributions, funding, competing interests, data availability, and code availability.

**Assessment: ADEQUATELY ADDRESSED.**

The requested title and ordering are present. Missing declaration content is assessed separately under R25.

### M7-20260729-R28

**Original comment.**

4. **[Major]** The current `Manuscript-Final.docx` fails LibreOffice headless conversion with `libpng error: Write Error`, although all seven embedded PNGs can be read and Microsoft Word exports a complete 30-page PDF. The supplement exports to 17 pages in Word. This is an interoperability defect in the submission package, not a scientific flaw. Rebuild both DOCX files from clean source, test them in Word and LibreOffice, and keep the package only when both applications produce complete PDFs without repair prompts or conversion errors.

**Authors' claimed response.**

Both DOCX files were rebuilt from clean Markdown with citation processing. Microsoft Word exported complete 33- and 22-page PDFs; LibreOffice exported complete 31- and 21-page PDFs without the former libpng failure or repair prompt. All 55 Word-rendered pages were visually inspected.

**What was found in the revised manuscript.**

Clean DOCX files were rebuilt. Word produced complete 33-page and 22-page PDFs, while LibreOffice produced complete 31-page and 21-page PDFs without the former libpng failure. All 55 Word-rendered pages were inspected, including the final references and Table S13.

**Assessment: ADEQUATELY ADDRESSED.**

The former interoperability defect is absent in two independent rendering applications.

### M7-20260729-R29

**Original comment.**

1. **[Major]** The contribution is real but narrower than the manuscript’s software-level framing sometimes suggests. M7.4 shows exact identity-limit nesting, information in native affine maps relative to permuted maps, and better planted-conflict localisation. It does not show overall benefit against the edge-local comparator. M7.5 shows ranking gains in two planted scenarios but fails the complete calibration gate. Frame the contribution as a falsifiable benchmark and a conditional representation result, not as validation of HydroSheaf as a whole.

**Authors' claimed response.**

The title, abstract, Discussion and Conclusion now frame M7 as a falsifiable controlled-synthetic benchmark and conditional representation result. They explicitly reject validation of HydroSheaf as a whole or a general superiority claim.

**What was found in the revised manuscript.**

The abstract states that the experiments do not establish general predictive superiority or field validity. The Conclusion says: 'They do not validate HydroSheaf as a whole or establish general predictive superiority over weighted graphs.'

**Assessment: ADEQUATELY ADDRESSED.**

The contribution is now a conditional representation and diagnostic result, not a framework-level validation claim.

### M7-20260729-R30

**Original comment.**

2. **[Major]** Development selected a local weight of 1.0 for both M7.5 hybrid arms, as reported at lines 856 to 860. The selected method therefore did not combine local and global residuals when both endpoints were observed. Its added capability was global fallback for missing endpoints plus map-sensitive conflict information. Use “local-first/global-fallback estimator” consistently and avoid presenting it as evidence that a general local/global blend is superior.

**Authors' claimed response.**

The selected M7.5 method is called local-first/global-fallback throughout explanatory prose. The paper states that development selected local weight 1.0, so this test is not evidence for a general local/global blend.

**What was found in the revised manuscript.**

The Results, Discussion, Conclusion, Figure 7 caption, and Table 4 caption call the selected method local-first/global-fallback and report selected local weight 1.0. The paper says this is not evidence that a general local/global blend is superior.

**Assessment: ADEQUATELY ADDRESSED.**

The nomenclature now matches what the selected weight actually did.

### M7-20260729-R31

**Original comment.**

3. **[Major]** The study has no independent cross-generator replication of the M7.4/M7.5 representation result and no field truth for topology, age or reaction mechanisms. This does not erase the contribution, but it limits its likely impact to a method and capability paper. State that the next claim-bearing step is replication under a different generator family or a field dataset with independently measured connectivity, and do not infer temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance from M7.

**Authors' claimed response.**

Limitations now identify independent replication under another generator family or field data with independently measured connectivity as the next claim-bearing step. They explicitly prohibit inference to temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance.

**What was found in the revised manuscript.**

The Limitations identify a different generator family or field data with independently measured connectivity as the next claim-bearing step. They explicitly exclude temporal, three-dimensional, vadose-zone, vector-stalk, and active-learning performance from the M7 evidence.

**Assessment: ADEQUATELY ADDRESSED.**

The limitation and next validation step are concrete, and unsupported capability extrapolation is excluded.

### M7-20260729-R32

**Original comment.**

1. **[Major]** The title, Introduction line 203 and Limitations lines 1134 to 1136 say or imply one generator and a fixed twelve-case locked test, while Methods lines 348 to 390 describe a separate scalar generator with 64 and 128 locked cases. Revise every scope statement to distinguish M7.3 from M7.4/M7.5.

**Authors' claimed response.**

All scope statements now distinguish the six-development/twelve-test M7.3 generator from the 64-case M7.4 and 128-case M7.5 scalar generator tests.

**What was found in the revised manuscript.**

The abstract, Methods, design table, and Limitations consistently distinguish the 6/12-case M7.3 system from the 32/64-case M7.4 and 64/128-case M7.5 scalar designs.

**Assessment: ADEQUATELY ADDRESSED.**

The previously conflicting generator and test-size descriptions are consistent.

### M7-20260729-R33

**Original comment.**

2. **[Major]** Introduction line 207 says “seven linked audits,” whereas Conclusion line 1180 says “six linked audits.” There are seven if the four original experiments, public-pipeline audit, M7.4 and M7.5 are counted. Use seven throughout and enumerate them once.

**Authors' claimed response.**

The paper now uses seven audits consistently and enumerates them in the Introduction and design table.

**What was found in the revised manuscript.**

The Introduction enumerates seven linked audits, Table 1 contains seven audit rows, and the Conclusion begins 'Across seven linked audits.' No conflicting six-audit statement remains.

**Assessment: ADEQUATELY ADDRESSED.**

The audit count is now internally consistent.

### M7-20260729-R34

**Original comment.**

3. **[Critical]** Results line 856 reports that the M7.5 provenance gate passed, but the current runner fails the confirmatory lock. Restore the exact source or qualify the statement until a new locked run is completed.

**Authors' claimed response.**

The exact M7.5 source was restored and independently hash-checked against the confirmatory lock. The Results report the matching hash and no-rerun status.

**What was found in the revised manuscript.**

The exact runner source is restored, the current SHA-256 matches the lock, and the Results states this verification. Unit tests exercise the recovered runner without executing the locked test set.

**Assessment: ADEQUATELY ADDRESSED.**

The critical source-lock inconsistency has been resolved without a new confirmatory run.

### M7-20260729-R35

**Original comment.**

4. **[Critical]** Code availability cites commit `d336e87` for M7.4/M7.5, while both later manifests report `53beb460...` and neither commit contains the later files. Replace the single commit claim with accurate, experiment-specific release identifiers.

**Authors' claimed response.**

The false single-commit claim was replaced with experiment-specific run identifiers, the distinct M7.3 freeze commit, the historical M7.4/M7.5 manifest revision and an explicit warning that the latter revision does not contain the files. The final commit/tag/DOI must be inserted only after release publication.

**What was found in the revised manuscript.**

The manuscript no longer presents d336e87 or 53beb460... as a valid release identifier for all experiments. It gives the M7.3 freeze separately and accurately describes the local M7.4/M7.5 run records, but no new versioned release identifier exists yet.

**Assessment: PARTIALLY ADDRESSED.**

The false commit claim is gone, which is an important correction, but the requested accurate experiment-specific release identifier does not yet exist for M7.4/M7.5. A local run ID and a file hash do not replace a retrievable versioned release.

**Guidance for further revision.**

After publishing the complete M7.4/M7.5 package, cite its commit, tag, and persistent identifier separately from the M7.3 d336e87 freeze. Verify that a clean checkout of that release can regenerate or validate every stated artifact.

### M7-20260729-R36

**Original comment.**

5. **[Minor]** The manuscript calls the M7.5 arm a hybrid while reporting a selected local weight of 1.0. The methods explain the missing-endpoint fallback, but figure and table labels do not. Rename the selected arm in explanatory prose or add “selected local weight = 1.0” to Figure 7 and Table 9 captions.

**Authors' claimed response.**

Figure 7 and main Table 4 state that the selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Machine-readable arm names are retained only where necessary to preserve provenance.

**What was found in the revised manuscript.**

Figure 7 says the selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Table 4 gives the same qualification. Internal machine arm names remain only in provenance-oriented supplementary material.

**Assessment: ADEQUATELY ADDRESSED.**

The explanatory labels prevent the reader from mistaking weight 1.0 for a fitted two-source blend.

## Part B: Current assessment summary

The authors engaged seriously with the technical review. The revised paper is not a cosmetic rewrite. It recovers and archives the exact M7.5 runner, adds full-family simultaneous inference for M7.4 and M7.5, supplies an honestly labelled post-review precision analysis, reduces the main table burden, rebuilds the Ghana figure without M6 evidence, and narrows every major scientific claim to what the controlled generators can establish. Most importantly, the revised paper now answers the ordinary-weighted-graph question directly. The answer is conditional: the sheaf layer encodes non-identity maps and provides global compatibility and conflict localisation, but it does not outperform the strong edge-local comparator overall under the complete gate.

The statistical interpretation is much stronger. The full 120- and 560-contrast families prevent selective scenario reporting, and the adverse heterogeneous-affine and null overall estimator results remain visible. The planning margins are not presented as if they had been prespecified. The distinction between the MODFLOW/MODPATH generator and the scalar affine generator is also clear throughout, which prevents capability tests from being mistaken for field validation or cross-generator replication.

The remaining issues are operational but publication-critical. M7.4 and M7.5 still lack a retrievable versioned repository release and persistent identifier. The declarations also still require author-supplied roles, funding, competing interests, and final data and software identifiers. These omissions do not invalidate the numerical findings, but they prevent the package from meeting the reproducibility and submission requirements set by the original critical comments.

The trajectory is strongly positive. Once the release and declaration work is completed and checked from a clean checkout, no further scientific reanalysis appears necessary for this review round.

## Part C: New issues

The revisions did not introduce a new scientific or statistical contradiction. The added precision and practical-magnitude analysis is clearly labelled as post-review planning, and the manuscript does not use it to recast the locked tests as adequately powered. The different Word and LibreOffice page counts reflect layout pagination; both render complete text, figures, tables, and references.

## Part D: Structured recommendation

**Adequacy summary.** 33 of 36 comments (91.7%) were adequately addressed. 3 were partially addressed and none were inadequately addressed. Two overlapping critical reproducibility comments remain partially open because the versioned M7.4/M7.5 release does not yet exist; one major declarations comment also remains partially open.

**1. Recommendation: Return for major revisions.** The scientific and statistical revisions are sufficient, but the unresolved critical release requirements trigger one more major round under the review decision framework. The required work is bounded: publish and verify the complete versioned M7.4/M7.5 package, mint the persistent identifiers, and replace the interim declaration text with final author-supplied statements.

**2. Study design and evidential support: Yes.** The design is appropriate for the bounded controlled-synthetic questions, the adverse controls are suitable, and the conclusions now remain within the evidence supported by each generator.

**3. Methods reproducibility: No, but it can be addressed with revision.** The methods and local artifacts are sufficiently detailed, and the exact locked runner is recovered. Repetition by an independent reader is still impeded by the absence of a complete versioned public release for M7.4/M7.5.

**4. Statistics and uncertainty: Yes.** The revision uses case-block resampling, full-family simultaneous intervals, transparent post-review planning, and practical-magnitude reporting. The retained claims match those analyses.

**5. Guidance on overstated claims: This was not needed further.** The earlier overstatements have been rewritten. The current conclusion gives a suitably bounded conditional representation claim.

**6. Presentation clarity: Yes.** The paper is dense but coherent, the main table count is manageable, the figures are readable at the recorded dimensions, and both DOCX files render completely in Word and LibreOffice.
