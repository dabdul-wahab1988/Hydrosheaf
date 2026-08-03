# Response to the 29 July 2026 M7 manuscript review

This response addresses all 36 technical comments against the live manuscript and immutable result records. Post-review analyses are explicitly labelled and do not alter or rerun a locked test. Two classes of submission metadata remain outside technical revision: author declarations, and creation of a versioned public release/DOI. They are recorded as blockers rather than invented.

## Resolution summary

- Technical/scientific/editorial comments addressed: 32.
- Technical correction complete but release publication pending: R07 and R35.
- Structure corrected but author/deposit metadata pending: R25 and R27.
- Locked confirmatory tests rerun: none.

## M7-20260729-R01 — ADDRESSED

**Reviewer comment**

1. **[Major]** The phrase “single-generator groundwater benchmark” is not accurate for the assembled paper. The M7.3 experiments used the MODFLOW 6, MODPATH 7, tracer and chemistry generator described at main-text lines 221 to 242, whereas M7.4 used a separate scalar graph-case generator that imported no HydroSheaf code, as stated at lines 348 to 358 and Supplementary lines 352 to 361. The title, Introduction line 203 and Limitations lines 1134 to 1136 therefore collapse two distinct generator systems and three different locked-test sizes into one. This matters because the claimed scope and independence of the evidence depend on which generator produced which result. Revise the title to refer to “controlled-synthetic benchmarks” or “two controlled-synthetic generator systems,” then identify M7.3 and M7.4/M7.5 separately in the abstract and limitations.

**Response**

The title now says controlled-synthetic benchmarks. The abstract, Introduction, Methods and Limitations name the MODFLOW 6/MODPATH 7 M7.3 generator separately from the scalar affine M7.4/M7.5 generator and prohibit cross-generator or field transfer of either result.

**Verification evidence**

- `manuscript/sections/00-abstract/section.md`
- `manuscript/sections/02-methods/section.md`
- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R02 — ADDRESSED

**Reviewer comment**

2. **[Minor]** The abstract reports the important positive, null and adverse findings, but it compresses seven audits and several estimands into one dense paragraph. At lines 8 to 26, the evidence-panel results, age results, reaction results, M7.4 representation test and M7.5 estimator test arrive without an explicit primary-versus-secondary hierarchy. A reader can miss that M7.4 rejected overall benefit against the claim-bearing comparator and that M7.5 also failed its complete three-outcome gate. Divide the abstract into purpose, design, primary representation result, follow-up result and scope conclusion, and use the words “failed the prespecified complete gate” for both M7.4 and M7.5.

**Response**

The abstract was split into purpose, design, primary representation result, follow-up estimator result and scope conclusion. It now states explicitly that both M7.4 and M7.5 failed their prespecified complete gates.

**Verification evidence**

- `manuscript/sections/00-abstract/section.md`

## M7-20260729-R03 — ADDRESSED_BY_CLAIM_RESTRICTION

**Reviewer comment**

1. **[Major]** The originality claim at main-text lines 102 to 109 is not supported by a reproducible literature search. `manuscript/methods/literature_search.json` records general web queries but no bibliographic databases, search dates by database, returned counts, deduplication, inclusion criteria or screening record for the claim that no groundwater benchmark combined an external generator, a locked split and adverse controls. The claim may be correct, but the current evidence cannot establish its completeness. Either replace it with a narrower statement that the authors did not identify such a benchmark in a targeted, non-systematic search, or archive a reproducible search across at least Web of Science or Scopus, Crossref, Google Scholar and a preprint index, with queries, dates, counts and screened records.

**Response**

The unsupported exhaustive novelty claim was withdrawn. The manuscript now describes a targeted, non-systematic search, and the search record stores dated queries, sources, inclusion logic and screened items without claiming bibliographic completeness.

**Verification evidence**

- `manuscript/sections/01-introduction/section.md`
- `manuscript/methods/literature_search.json`
- `manuscript/LITERATURE.bib`

## M7-20260729-R04 — ADDRESSED

**Reviewer comment**

2. **[Major]** The objective structure obscures which experiments were confirmatory. Main-text lines 159 to 190 list five objectives, but objective five contains M7.4 and a second fresh-seed diagnostic, while the manuscript later describes seven linked audits. M7.3, the public-pipeline run, M7.4 and M7.5 also used different generators, case counts, primary outcomes and claim gates. This makes it difficult to distinguish the original scientific questions from later estimator diagnosis. Add a one-page design table at the end of the Introduction that maps each audit to its generator, development and test counts, primary comparator, primary metrics, lock point, multiplicity family and permitted claim.

**Response**

A main-text design-and-claim table now maps all seven audits to generator, development/test counts, comparator, primary metrics, lock/multiplicity family and permitted claim.

**Verification evidence**

- `tables/publication/table1_benchmark_design.md`
- `manuscript/Manuscript-Final-Revised.md`

## M7-20260729-R05 — ADDRESSED

**Reviewer comment**

3. **[Minor]** The practical significance is framed more strongly than the decision analysis supports. Lines 198 to 203 acknowledge that translation to a field decision protocol remains future work, yet the surrounding text says the benchmark identifies which evidence was “worth collecting.” No value-of-information threshold, cost function or field decision loss was defined. Replace “worth collecting” with “improved the prespecified synthetic-benchmark outcomes,” and reserve collection recommendations for a future study that defines measurement cost and decision consequences.

**Response**

The phrase 'worth collecting' was removed. Collection recommendations are reserved for a future value-of-information study with measurement costs and decision losses.

**Verification evidence**

- `manuscript/sections/01-introduction/section.md`
- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R06 — ADDRESSED

**Reviewer comment**

1. **[Critical]** The current M7.5 runner does not match the confirmatory lock. `m7_robust_hybrid_confirmatory.lock.json`, line 6, requires SHA-256 `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`; the current `scripts/run_m7_robust_hybrid_sheaf.py` hashes to `994e36954775c0577dd6f7a7655e9cf6d90d4cec2463a2cf4fc19891d8af7c12`. Direct execution stops at source line 150 with `Confirmatory-bound artifact changed: runner_sha256`. Stored output hashes remain intact, but the exact code that passed the one-time lock is not recoverable from the present checkout. Recover the exact locked runner from backup or an earlier task snapshot and archive it under its content hash. If it cannot be recovered, state that the stored outputs are an unverifiable locked-run record, create a new protocol and confirmatory lock, and run one new test set with previously unused seeds.

**Response**

The exact one-time M7.5 runner was recovered from the archived task session. Its SHA-256 is a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191, exactly matching the confirmatory lock. The source and recovery manifest are archived under that content hash. No locked test was rerun.

**Verification evidence**

- `scripts/run_m7_robust_hybrid_sheaf.py`
- `provenance/source_archive/a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191/manifest.json`

## M7-20260729-R07 — TECHNICAL_CORRECTION_COMPLETE_RELEASE_PENDING

**Reviewer comment**

2. **[Critical]** The software availability statement is factually wrong for M7.4 and M7.5. Main-text lines 1249 to 1259 say the later protocols, runners, tests and result generators are held in the repository at protocol-freeze commit `d336e87`. Both M7.4 and M7.5 manifests instead record Git revision `53beb46034d5230c1a061341a5cf2175d9af858e`, and the M7.4/M7.5 protocol, runner and result files are absent even from that commit and from the current remote `main`. The recorded revision therefore cannot recreate either later experiment. Commit the complete M7.4/M7.5 protocol, code, tests, manifests and immutable outputs, tag a release, deposit that release in a persistent repository, and cite separate M7.3 and M7.4/M7.5 freeze identifiers in the manuscript.

**Response**

The false commit-level reconstruction claim was removed and replaced with experiment-specific run identifiers, recorded historical revisions, current file hashes and an explicit statement that the historical commits do not contain M7.4/M7.5. The complete local package is assembled. A versioned commit/tag and persistent DOI remain submission actions because this revision task does not authorise repository publication or create a repository DOI.

**Verification evidence**

- `manuscript/sections/06-availability-statements/section.md`
- `manuscript/artifact_registry.json`
- `provenance/source_archive`

## M7-20260729-R08 — ADDRESSED

**Reviewer comment**

3. **[Major]** The paper needs to separate the construct validity of its two generator systems. M7.3 used a flow, tracer and chemistry generator, whereas M7.4/M7.5 used scalar affine graph cases in four planted scenarios. Both are legitimate capability tests, but success in the scalar generator does not validate performance in the MODFLOW-based generator, and neither validates a field aquifer. Add a methods subsection that names the two generator systems, states which objectives each can test, and defines the non-transfer boundary between them.

**Response**

A construct-validity subsection distinguishes the two generators, the questions each can test and the explicit non-transfer boundary: scalar-case success is not MODFLOW-system validation, and neither generator constitutes field validation.

**Verification evidence**

- `manuscript/sections/02-methods/section.md`
- `manuscript/supplementary/Supplementary-Methods.md`

## M7-20260729-R09 — ADDRESSED_POST_REVIEW

**Reviewer comment**

4. **[Major]** No prospective sample-size or precision rationale is reported for six development cases and twelve locked M7.3 cases, 64 M7.4 test cases or 128 M7.5 test cases. Bootstrap intervals quantify the realized uncertainty but do not show why these sample sizes were adequate to detect a scientifically relevant difference. This is especially important because several findings rest near zero and the M7.3 exact tests are discrete at 12 blocks. Add a simulation-based precision analysis using development-only effect distributions, state the minimum difference of interest for each primary metric, and report the probability that each locked design would exclude the prespecified null or non-inferiority margin.

**Response**

A labelled post-review simulation study used development-only planning inputs, 20,000 replicates and prespecified minimum differences for PR-AUC, Brier, log loss, age MAE, interval width, coverage and modal accuracy. The attainable precision/power results are reported in Supplementary Table S12 and are not presented as prospective preregistration.

**Verification evidence**

- `scripts/post_review_statistical_audit.py`
- `results/post_review_audit_20260730`
- `tables/publication/tableS12_precision_power.md`

## M7-20260729-R10 — ADDRESSED_POST_REVIEW

**Reviewer comment**

5. **[Major]** M7.3 received an explicit 24-test Benjamini-Hochberg analysis, but the M7.4 and M7.5 scenario claims rely on unadjusted percentile intervals drawn from much larger contrast matrices. Supplementary Table S7 contains 120 rows and Table S9 contains 560 rows. Prespecified directions reduce researcher discretion but do not remove family-wise selection risk when multiple scenarios and metrics are interpreted. Define the exact confirmatory family for M7.4 and M7.5, then provide simultaneous bootstrap intervals or adjusted permutation tests for the scenario statements that are retained in the abstract, results and conclusion.

**Response**

All 120 published M7.4 contrasts and all 560 published M7.5 contrasts were placed in separate full families. Shared case-block bootstrap resampling with 10,000 replicates and max-z simultaneous 95% intervals was applied. Only findings surviving those families remain inferentially supported.

**Verification evidence**

- `scripts/post_review_statistical_audit.py`
- `tables/publication/tableS10_m7_4_multiplicity_adjusted.md`
- `tables/publication/tableS11_m7_5_multiplicity_adjusted.md`

## M7-20260729-R11 — ADDRESSED_BY_CLAIM_RESTRICTION

**Reviewer comment**

6. **[Major]** The reaction-family experiment uses the same six-family vocabulary for generation and scoring. Supplementary lines 172 to 186 show that predicted reactions were mapped to the six planted process families, while the generator planted those same archetypes. Code independence prevents leakage, but a shared ontology can inflate recovery for well-separated families and does not test out-of-dictionary chemistry. Add a mechanism-mismatch sensitivity that perturbs stoichiometry, mineral assemblage and reaction combinations outside the scoring dictionary, or restrict the claim to discrimination among the six planted archetypes.

**Response**

No unplanned mechanism-mismatch experiment was added. Instead, every reaction claim is restricted to discrimination among the six planted archetypes, under the two tested indicator panels and the tested noise model; out-of-dictionary chemistry is explicitly unevaluated.

**Verification evidence**

- `manuscript/sections/02-methods/section.md`
- `manuscript/sections/03-results/section.md`
- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R12 — ADDRESSED

**Reviewer comment**

7. **[Minor]** The Supplementary Methods describe unweighted logistic regression as “the only choice consistent with reporting calibration honestly.” That statement is too absolute. Weighting changes the fitted target and can harm raw calibration, but weighted estimators can still be assessed and recalibrated on an untouched development fold. Rewrite this as a prespecified design choice, report the class prevalence in each generator, and explain why the unweighted target matched the intended probability estimand.

**Response**

The absolute statement about unweighted logistic regression was replaced with a prespecified estimand rationale. Candidate-edge class prevalence is reported for each generator, and the text explains why the unweighted target matches the intended per-candidate probability estimand.

**Verification evidence**

- `manuscript/sections/02-methods/section.md`
- `manuscript/supplementary/Supplementary-Methods.md`

## M7-20260729-R13 — ADDRESSED

**Reviewer comment**

1. **[Major]** The statement at main-text lines 856 to 865 that M7.5 passed its execution and provenance gate cannot remain unqualified while the current runner fails the confirmatory hash. The stored case tables support the numerical contrasts, but the present archive does not support the stronger provenance statement. Suspend that sentence until the exact runner is restored, or revise it to state that the stored manifest reports a pass but the current source archive does not verify the runner hash.

**Response**

Because the exact runner was recovered and its hash verified against the lock, the provenance statement now names the recovered hash and source archive. It also states explicitly that the stored locked test was not rerun.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `provenance/source_archive/a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191/manifest.json`

## M7-20260729-R14 — ADDRESSED_POST_REVIEW

**Reviewer comment**

2. **[Major]** Statistical precision is not paired with practical magnitude. Complete topology reduced mean absolute age error by 0.062 years relative to a 2.764-year baseline and by 0.164 years relative to a 4.750-year baseline, approximately 2.2% and 3.5%. Interval-width reductions were also small relative to the baseline widths. Similarly, the M7.5 PR-AUC difference of 0.0200 is about 3.1% of the edge-local mean, while the complete calibration gate failed. Add relative effects, baseline dispersion and a prespecified minimum practically important difference to the Results, then state whether each result cleared that threshold.

**Response**

Results now pair absolute effects with relative changes and post-review practical margins. The small topology-age changes and the M7.5 overall PR-AUC change are explicitly classified by whether they clear those margins; the M7.5 complete calibration gate remains failed.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `tables/publication/tableS12_precision_power.md`

## M7-20260729-R15 — ADDRESSED_POST_REVIEW

**Reviewer comment**

3. **[Major]** The incompatible-cycle and noisy/missing findings at lines 892 to 900 are presented as supported scenario gains without multiplicity control across the full scenario-metric matrix. These are useful diagnostic observations, but their confirmatory status is weaker than the prose implies. Apply the family correction requested under Methodology, or relabel these findings as prespecified exploratory diagnostics and remove them from the abstract until corrected intervals are available.

**Response**

Scenario statements were reclassified using the simultaneous full-family intervals. The incompatible-cycle conflict-localisation signal survives; the noisy/missing overall gain and M7.5 scenario ranking-gain claims do not and were withdrawn as supported general gains.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `tables/publication/tableS10_m7_4_multiplicity_adjusted.md`
- `tables/publication/tableS11_m7_5_multiplicity_adjusted.md`

## M7-20260729-R16 — ADDRESSED

**Reviewer comment**

4. **[Major]** The heading “Carbonate reactions remain non-identifiable regardless of panel richness” and related discussion extend beyond the tested panels. The study compared one core and one enhanced indicator set under one reaction dictionary and one noise model. The zero recovery is a valid finding for those conditions, but it does not establish non-identifiability regardless of isotopic, mineralogical or thermodynamic measurements. Rename the subsection “Carbonate reactions were not recovered under either tested indicator panel” and replace every broader formulation with the same boundary.

**Response**

The heading now reads 'Carbonate reactions were not recovered under either tested indicator panel.' All broader 'regardless of panel richness' or universal non-identifiability language was removed.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R17 — ADDRESSED

**Reviewer comment**

5. **[Major]** Figure 5 mixes the M7 Ghana audit with a companion M6 field-transfer experiment. The caption discloses this, but placing the M6 tier ablation inside an M7 result figure invites readers to treat all four panels as evidence from the same protocol. Move the M6 panel to a clearly labelled contextual figure in the Supplement, or present it in a boxed comparison with its own protocol, sample, outcome and citation. The M7 field claim must rest only on the truth-free Ghana scope audit and hold-forward analysis.

**Response**

Figure 5 was rebuilt as an M7-only Ghana supportability figure. The M6 tier-ablation panel was removed; the field claim now rests only on the truth-free Ghana scope and hold-forward audit.

**Verification evidence**

- `figures/publication/figure5_ghana_supportability_boundary_m7_only.png`
- `scripts/make_m7_3_publication_assets.py`
- `manuscript/sections/03-results/section.md`

## M7-20260729-R18 — ADDRESSED

**Reviewer comment**

6. **[Minor]** The public-pipeline table reports the same selected F1 value, 0.4222, for all four arms, despite large changes in probability scores and log loss. Candidate recall is also 0.9815 rather than complete. These observations suggest that the selected threshold was insensitive in six cases and that the system result is conditional on candidate generation. Add thresholds and confusion counts for each arm, explain the identical F1 result, and state explicitly that the system audit tested scoring only on the candidates that were recovered.

**Response**

The public-pipeline audit now reports the threshold and TP/FP/FN counts for every arm, explains the identical selected F1 values, distinguishes conditional from end-to-end recall and states that all generated candidates were selected, so the audit did not identify a useful scalar selection threshold.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `tables/publication/tableS13_public_pipeline_selection.md`

## M7-20260729-R19 — ADDRESSED

**Reviewer comment**

7. **[Minor]** The conclusion says the sheaf layer “has therefore earned a bounded scientific-workflow role” at lines 1203 to 1205. “Earned” is evaluative language rather than a test result. Replace it with: “The experiments support use of the sheaf layer as a prespecified model of non-identity relations and as a global-compatibility diagnostic under the tested scalar scenarios, with global fallback when endpoint evidence is missing.”

**Response**

The evaluative 'earned' sentence was replaced verbatim with the requested bounded statement about non-identity relations, global-compatibility diagnosis and missing-endpoint fallback under the tested scalar scenarios.

**Verification evidence**

- `manuscript/sections/05-conclusion/section.md`

## M7-20260729-R20 — ADDRESSED

**Reviewer comment**

1. **[Major]** Figure 5 combines data and claims from two modules, as described above. Its provenance is disclosed in the caption and artifact registry, but the visual grouping still implies a common experiment. Separate the M6 panel or mark it inside the panel title and legend as external companion evidence, including the independent sample size and protocol identifier.

**Response**

The combined-module visual was removed. Figure 5 and its caption now contain only M7 evidence and explicitly distinguish synthetic supportability context from the truth-free Ghana audit.

**Verification evidence**

- `figures/publication/figure5_ghana_supportability_boundary_m7_only.png`
- `figures/publication/figure_source_manifest.csv`

## M7-20260729-R21 — ADDRESSED

**Reviewer comment**

2. **[Minor]** The Word render is clean, but the main paper contains seven figures and nine tables in 30 pages. Table 7 spans three pages and interrupts the transition to the representation results. Retain the design table, the primary M7.3 decision table, Table 8 and Table 9 in the main paper, and move detailed metric tables to the Supplement, where the complete CSVs are already cited.

**Response**

The main paper now has four tables: the design map, a compact M7.3 decision table, M7.4 means and M7.5 means. Detailed metric tables were moved to the 13-table Supplement.

**Verification evidence**

- `manuscript/Manuscript-Final-Revised.md`
- `manuscript/supplementary/Supplementary-Information.md`

## M7-20260729-R22 — ADDRESSED

**Reviewer comment**

3. **[Minor]** Figures 6 and 7 are legible at full-page width but their four-panel labels and several axis labels become small in the 30-page Word layout. Increase label and tick sizes for final journal dimensions, and include the intended printed width in the figure-generation manifest so that legibility can be checked automatically.

**Response**

Figures 6 and 7 were regenerated for 7.08-inch journal width with a minimum 8-point label size. Those dimensions are recorded in the figure-source manifest and were checked in the Word and LibreOffice renders.

**Verification evidence**

- `scripts/make_m7_sheaf_vs_graph_assets.py`
- `scripts/make_m7_robust_hybrid_assets.py`
- `figures/publication/figure_source_manifest.csv`

## M7-20260729-R23 — ADDRESSED_WITH_TARGETED_SCOPE

**Reviewer comment**

1. **[Major]** Eighteen unique references are too few to support the paper’s broad positioning across groundwater joint inversion, age dating, reactive transport, probabilistic calibration, graph inference and sheaf methods. The present search record is strongest for software documentation and selected method citations but weak for competing groundwater graph or structured-residual approaches. Expand the literature review using the reproducible search requested above, and add the closest non-sheaf alternatives rather than only mathematical sheaf sources and one recent neural-network preprint.

**Response**

The Introduction and Discussion now include recent 2022-2025 tracer, hydrochemistry, groundwater-age and machine-learning work, plus non-sheaf structured alternatives such as graph regularisation, Gaussian-process smoothing, flow-network inversion and residual diagnostics. The search is disclosed as targeted rather than exhaustive.

**Verification evidence**

- `manuscript/sections/01-introduction/section.md`
- `manuscript/sections/04-discussion/section.md`
- `manuscript/LITERATURE.bib`
- `manuscript/methods/literature_search.json`

## M7-20260729-R24 — ADDRESSED

**Reviewer comment**

2. **[Minor]** Davis and Goadrich are cited at main-text lines 416 to 420 and Supplementary lines 231 to 236 as showing that PR-AUC is more informative under imbalance. The paper’s primary statement concerns precision-recall curves under highly skewed data and the relation between ROC and PR spaces, not an unconditional ranking of scalar AUC summaries. Change the wording to “precision-recall curves can give a more informative view under class imbalance,” and retain ROC-AUC as the manuscript already does.

**Response**

The Davis-and-Goadrich wording now states that precision-recall curves can give a more informative view under class imbalance; ROC-AUC remains reported.

**Verification evidence**

- `manuscript/sections/02-methods/section.md`
- `manuscript/supplementary/Supplementary-Methods.md`

## M7-20260729-R25 — EXTERNAL_AUTHOR_METADATA_PENDING

**Reviewer comment**

1. **[Major]** The declarations are submission blockers. Main-text lines 1213 to 1222 contain placeholders for author contributions and competing interests; lines 1233 to 1240 and 1257 to 1260 contain data and software DOI placeholders. There is no funding statement. Replace every placeholder with final contributor roles, competing-interest and funding declarations, dataset DOI, software DOI, licence and versioned release before submission.

**Response**

The section is correctly structured and the licence and available technical identifiers are stated. Author CRediT roles, funding, competing interests, dataset DOI, software DOI and final versioned release are retained as explicit submission blockers because they require author declarations or external deposits and must not be fabricated.

**Verification evidence**

- `manuscript/sections/06-availability-statements/section.md`

## M7-20260729-R26 — ADDRESSED

**Reviewer comment**

2. **[Minor]** Main-text line 318 reads, “The reaction solver was evaluated inference, not part of the synthetic generator.” The sentence is ungrammatical and obscures an important independence claim. Replace it with: “The reaction solver was the inference method under evaluation and was not part of the synthetic generator.”

**Response**

The sentence was replaced with: 'The reaction solver was the inference method under evaluation and was not part of the synthetic generator.'

**Verification evidence**

- `manuscript/sections/02-methods/section.md`

## M7-20260729-R27 — ADDRESSED_STRUCTURE_METADATA_PENDING

**Reviewer comment**

3. **[Minor]** “Open Research” at line 1211 contains authorship, interests, data and code statements rather than only open-research material. Rename the section “Declarations and open research,” and order it as author contributions, funding, competing interests, data availability and code availability.

**Response**

The section is now titled 'Declarations and open research' and is ordered as author contributions, funding, competing interests, data availability and code availability. Only the author-supplied metadata listed under R25 remains open.

**Verification evidence**

- `manuscript/sections/06-availability-statements/section.md`

## M7-20260729-R28 — ADDRESSED

**Reviewer comment**

4. **[Major]** The current `Manuscript-Final.docx` fails LibreOffice headless conversion with `libpng error: Write Error`, although all seven embedded PNGs can be read and Microsoft Word exports a complete 30-page PDF. The supplement exports to 17 pages in Word. This is an interoperability defect in the submission package, not a scientific flaw. Rebuild both DOCX files from clean source, test them in Word and LibreOffice, and keep the package only when both applications produce complete PDFs without repair prompts or conversion errors.

**Response**

Both DOCX files were rebuilt from clean Markdown with citation processing. Microsoft Word exported complete 33- and 22-page PDFs; LibreOffice exported complete 31- and 21-page PDFs without the former libpng failure or repair prompt. All 55 Word-rendered pages were visually inspected.

**Verification evidence**

- `manuscript/Manuscript-Final-Revised-20260730.docx`
- `manuscript/supplementary/Supplementary-Information-Revised-20260730.docx`

## M7-20260729-R29 — ADDRESSED

**Reviewer comment**

1. **[Major]** The contribution is real but narrower than the manuscript’s software-level framing sometimes suggests. M7.4 shows exact identity-limit nesting, information in native affine maps relative to permuted maps, and better planted-conflict localisation. It does not show overall benefit against the edge-local comparator. M7.5 shows ranking gains in two planted scenarios but fails the complete calibration gate. Frame the contribution as a falsifiable benchmark and a conditional representation result, not as validation of HydroSheaf as a whole.

**Response**

The title, abstract, Discussion and Conclusion now frame M7 as a falsifiable controlled-synthetic benchmark and conditional representation result. They explicitly reject validation of HydroSheaf as a whole or a general superiority claim.

**Verification evidence**

- `manuscript/sections/00-abstract/section.md`
- `manuscript/sections/04-discussion/section.md`
- `manuscript/sections/05-conclusion/section.md`

## M7-20260729-R30 — ADDRESSED

**Reviewer comment**

2. **[Major]** Development selected a local weight of 1.0 for both M7.5 hybrid arms, as reported at lines 856 to 860. The selected method therefore did not combine local and global residuals when both endpoints were observed. Its added capability was global fallback for missing endpoints plus map-sensitive conflict information. Use “local-first/global-fallback estimator” consistently and avoid presenting it as evidence that a general local/global blend is superior.

**Response**

The selected M7.5 method is called local-first/global-fallback throughout explanatory prose. The paper states that development selected local weight 1.0, so this test is not evidence for a general local/global blend.

**Verification evidence**

- `manuscript/sections/03-results/section.md`
- `manuscript/sections/04-discussion/section.md`
- `manuscript/sections/05-conclusion/section.md`

## M7-20260729-R31 — ADDRESSED

**Reviewer comment**

3. **[Major]** The study has no independent cross-generator replication of the M7.4/M7.5 representation result and no field truth for topology, age or reaction mechanisms. This does not erase the contribution, but it limits its likely impact to a method and capability paper. State that the next claim-bearing step is replication under a different generator family or a field dataset with independently measured connectivity, and do not infer temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance from M7.

**Response**

Limitations now identify independent replication under another generator family or field data with independently measured connectivity as the next claim-bearing step. They explicitly prohibit inference to temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance.

**Verification evidence**

- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R32 — ADDRESSED

**Reviewer comment**

1. **[Major]** The title, Introduction line 203 and Limitations lines 1134 to 1136 say or imply one generator and a fixed twelve-case locked test, while Methods lines 348 to 390 describe a separate scalar generator with 64 and 128 locked cases. Revise every scope statement to distinguish M7.3 from M7.4/M7.5.

**Response**

All scope statements now distinguish the six-development/twelve-test M7.3 generator from the 64-case M7.4 and 128-case M7.5 scalar generator tests.

**Verification evidence**

- `manuscript/sections/00-abstract/section.md`
- `manuscript/sections/01-introduction/section.md`
- `manuscript/sections/04-discussion/section.md`

## M7-20260729-R33 — ADDRESSED

**Reviewer comment**

2. **[Major]** Introduction line 207 says “seven linked audits,” whereas Conclusion line 1180 says “six linked audits.” There are seven if the four original experiments, public-pipeline audit, M7.4 and M7.5 are counted. Use seven throughout and enumerate them once.

**Response**

The paper now uses seven audits consistently and enumerates them in the Introduction and design table.

**Verification evidence**

- `manuscript/sections/01-introduction/section.md`
- `tables/publication/table1_benchmark_design.md`
- `manuscript/sections/05-conclusion/section.md`

## M7-20260729-R34 — ADDRESSED

**Reviewer comment**

3. **[Critical]** Results line 856 reports that the M7.5 provenance gate passed, but the current runner fails the confirmatory lock. Restore the exact source or qualify the statement until a new locked run is completed.

**Response**

The exact M7.5 source was restored and independently hash-checked against the confirmatory lock. The Results report the matching hash and no-rerun status.

**Verification evidence**

- `scripts/run_m7_robust_hybrid_sheaf.py`
- `manuscript/sections/03-results/section.md`

## M7-20260729-R35 — FACTUAL_CORRECTION_COMPLETE_RELEASE_PENDING

**Reviewer comment**

4. **[Critical]** Code availability cites commit `d336e87` for M7.4/M7.5, while both later manifests report `53beb460...` and neither commit contains the later files. Replace the single commit claim with accurate, experiment-specific release identifiers.

**Response**

The false single-commit claim was replaced with experiment-specific run identifiers, the distinct M7.3 freeze commit, the historical M7.4/M7.5 manifest revision and an explicit warning that the latter revision does not contain the files. The final commit/tag/DOI must be inserted only after release publication.

**Verification evidence**

- `manuscript/sections/06-availability-statements/section.md`

## M7-20260729-R36 — ADDRESSED

**Reviewer comment**

5. **[Minor]** The manuscript calls the M7.5 arm a hybrid while reporting a selected local weight of 1.0. The methods explain the missing-endpoint fallback, but figure and table labels do not. Rename the selected arm in explanatory prose or add “selected local weight = 1.0” to Figure 7 and Table 9 captions.

**Response**

Figure 7 and main Table 4 state that the selected nominal hybrid had local weight 1.0 and is therefore local-first/global-fallback. Machine-readable arm names are retained only where necessary to preserve provenance.

**Verification evidence**

- `manuscript/Manuscript-Final-Revised.md`
- `figures/publication/figure_source_manifest.csv`
