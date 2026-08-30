# Peer review of the M7 manuscript

**Manuscript:** *Conditional evidence integration and the incremental contribution of sheaf structure in a single-generator groundwater benchmark*  
**Files reviewed:** `manuscript/Manuscript-Final.md`, `manuscript/Manuscript-Final.docx`, `manuscript/supplementary/Supplementary-Information.md`, `manuscript/supplementary/Supplementary-Information.docx`, the locked M7.3, M7.4, M7.5 and public-pipeline result records, protocol locks, source files, tests, tables, figures, manifests, and availability statements  
**Review date:** 29 July 2026  
**Reviewer lens:** groundwater hydrology, computational hydrogeology, probabilistic evidence integration, graph and sheaf methods, and reproducible research software  
**Decision:** Return for major revisions  
**Score:** 59/100

## 1. Section-by-section critique

### 1.1 Title and abstract

1. **[Major]** The phrase “single-generator groundwater benchmark” is not accurate for the assembled paper. The M7.3 experiments used the MODFLOW 6, MODPATH 7, tracer and chemistry generator described at main-text lines 221 to 242, whereas M7.4 used a separate scalar graph-case generator that imported no HydroSheaf code, as stated at lines 348 to 358 and Supplementary lines 352 to 361. The title, Introduction line 203 and Limitations lines 1134 to 1136 therefore collapse two distinct generator systems and three different locked-test sizes into one. This matters because the claimed scope and independence of the evidence depend on which generator produced which result. Revise the title to refer to “controlled-synthetic benchmarks” or “two controlled-synthetic generator systems,” then identify M7.3 and M7.4/M7.5 separately in the abstract and limitations.

2. **[Minor]** The abstract reports the important positive, null and adverse findings, but it compresses seven audits and several estimands into one dense paragraph. At lines 8 to 26, the evidence-panel results, age results, reaction results, M7.4 representation test and M7.5 estimator test arrive without an explicit primary-versus-secondary hierarchy. A reader can miss that M7.4 rejected overall benefit against the claim-bearing comparator and that M7.5 also failed its complete three-outcome gate. Divide the abstract into purpose, design, primary representation result, follow-up result and scope conclusion, and use the words “failed the prespecified complete gate” for both M7.4 and M7.5.

### 1.2 Introduction

1. **[Major]** The originality claim at main-text lines 102 to 109 is not supported by a reproducible literature search. `manuscript/methods/literature_search.json` records general web queries but no bibliographic databases, search dates by database, returned counts, deduplication, inclusion criteria or screening record for the claim that no groundwater benchmark combined an external generator, a locked split and adverse controls. The claim may be correct, but the current evidence cannot establish its completeness. Either replace it with a narrower statement that the authors did not identify such a benchmark in a targeted, non-systematic search, or archive a reproducible search across at least Web of Science or Scopus, Crossref, Google Scholar and a preprint index, with queries, dates, counts and screened records.

2. **[Major]** The objective structure obscures which experiments were confirmatory. Main-text lines 159 to 190 list five objectives, but objective five contains M7.4 and a second fresh-seed diagnostic, while the manuscript later describes seven linked audits. M7.3, the public-pipeline run, M7.4 and M7.5 also used different generators, case counts, primary outcomes and claim gates. This makes it difficult to distinguish the original scientific questions from later estimator diagnosis. Add a one-page design table at the end of the Introduction that maps each audit to its generator, development and test counts, primary comparator, primary metrics, lock point, multiplicity family and permitted claim.

3. **[Minor]** The practical significance is framed more strongly than the decision analysis supports. Lines 198 to 203 acknowledge that translation to a field decision protocol remains future work, yet the surrounding text says the benchmark identifies which evidence was “worth collecting.” No value-of-information threshold, cost function or field decision loss was defined. Replace “worth collecting” with “improved the prespecified synthetic-benchmark outcomes,” and reserve collection recommendations for a future study that defines measurement cost and decision consequences.

### 1.3 Methodology

1. **[Critical]** The current M7.5 runner does not match the confirmatory lock. `m7_robust_hybrid_confirmatory.lock.json`, line 6, requires SHA-256 `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`; the current `scripts/run_m7_robust_hybrid_sheaf.py` hashes to `994e36954775c0577dd6f7a7655e9cf6d90d4cec2463a2cf4fc19891d8af7c12`. Direct execution stops at source line 150 with `Confirmatory-bound artifact changed: runner_sha256`. Stored output hashes remain intact, but the exact code that passed the one-time lock is not recoverable from the present checkout. Recover the exact locked runner from backup or an earlier task snapshot and archive it under its content hash. If it cannot be recovered, state that the stored outputs are an unverifiable locked-run record, create a new protocol and confirmatory lock, and run one new test set with previously unused seeds.

2. **[Critical]** The software availability statement is factually wrong for M7.4 and M7.5. Main-text lines 1249 to 1259 say the later protocols, runners, tests and result generators are held in the repository at protocol-freeze commit `d336e87`. Both M7.4 and M7.5 manifests instead record Git revision `53beb46034d5230c1a061341a5cf2175d9af858e`, and the M7.4/M7.5 protocol, runner and result files are absent even from that commit and from the current remote `main`. The recorded revision therefore cannot recreate either later experiment. Commit the complete M7.4/M7.5 protocol, code, tests, manifests and immutable outputs, tag a release, deposit that release in a persistent repository, and cite separate M7.3 and M7.4/M7.5 freeze identifiers in the manuscript.

3. **[Major]** The paper needs to separate the construct validity of its two generator systems. M7.3 used a flow, tracer and chemistry generator, whereas M7.4/M7.5 used scalar affine graph cases in four planted scenarios. Both are legitimate capability tests, but success in the scalar generator does not validate performance in the MODFLOW-based generator, and neither validates a field aquifer. Add a methods subsection that names the two generator systems, states which objectives each can test, and defines the non-transfer boundary between them.

4. **[Major]** No prospective sample-size or precision rationale is reported for six development cases and twelve locked M7.3 cases, 64 M7.4 test cases or 128 M7.5 test cases. Bootstrap intervals quantify the realized uncertainty but do not show why these sample sizes were adequate to detect a scientifically relevant difference. This is especially important because several findings rest near zero and the M7.3 exact tests are discrete at 12 blocks. Add a simulation-based precision analysis using development-only effect distributions, state the minimum difference of interest for each primary metric, and report the probability that each locked design would exclude the prespecified null or non-inferiority margin.

5. **[Major]** M7.3 received an explicit 24-test Benjamini-Hochberg analysis, but the M7.4 and M7.5 scenario claims rely on unadjusted percentile intervals drawn from much larger contrast matrices. Supplementary Table S7 contains 120 rows and Table S9 contains 560 rows. Prespecified directions reduce researcher discretion but do not remove family-wise selection risk when multiple scenarios and metrics are interpreted. Define the exact confirmatory family for M7.4 and M7.5, then provide simultaneous bootstrap intervals or adjusted permutation tests for the scenario statements that are retained in the abstract, results and conclusion.

6. **[Major]** The reaction-family experiment uses the same six-family vocabulary for generation and scoring. Supplementary lines 172 to 186 show that predicted reactions were mapped to the six planted process families, while the generator planted those same archetypes. Code independence prevents leakage, but a shared ontology can inflate recovery for well-separated families and does not test out-of-dictionary chemistry. Add a mechanism-mismatch sensitivity that perturbs stoichiometry, mineral assemblage and reaction combinations outside the scoring dictionary, or restrict the claim to discrimination among the six planted archetypes.

7. **[Minor]** The Supplementary Methods describe unweighted logistic regression as “the only choice consistent with reporting calibration honestly.” That statement is too absolute. Weighting changes the fitted target and can harm raw calibration, but weighted estimators can still be assessed and recalibrated on an untouched development fold. Rewrite this as a prespecified design choice, report the class prevalence in each generator, and explain why the unweighted target matched the intended probability estimand.

### 1.4 Results and discussion

1. **[Major]** The statement at main-text lines 856 to 865 that M7.5 passed its execution and provenance gate cannot remain unqualified while the current runner fails the confirmatory hash. The stored case tables support the numerical contrasts, but the present archive does not support the stronger provenance statement. Suspend that sentence until the exact runner is restored, or revise it to state that the stored manifest reports a pass but the current source archive does not verify the runner hash.

2. **[Major]** Statistical precision is not paired with practical magnitude. Complete topology reduced mean absolute age error by 0.062 years relative to a 2.764-year baseline and by 0.164 years relative to a 4.750-year baseline, approximately 2.2% and 3.5%. Interval-width reductions were also small relative to the baseline widths. Similarly, the M7.5 PR-AUC difference of 0.0200 is about 3.1% of the edge-local mean, while the complete calibration gate failed. Add relative effects, baseline dispersion and a prespecified minimum practically important difference to the Results, then state whether each result cleared that threshold.

3. **[Major]** The incompatible-cycle and noisy/missing findings at lines 892 to 900 are presented as supported scenario gains without multiplicity control across the full scenario-metric matrix. These are useful diagnostic observations, but their confirmatory status is weaker than the prose implies. Apply the family correction requested under Methodology, or relabel these findings as prespecified exploratory diagnostics and remove them from the abstract until corrected intervals are available.

4. **[Major]** The heading “Carbonate reactions remain non-identifiable regardless of panel richness” and related discussion extend beyond the tested panels. The study compared one core and one enhanced indicator set under one reaction dictionary and one noise model. The zero recovery is a valid finding for those conditions, but it does not establish non-identifiability regardless of isotopic, mineralogical or thermodynamic measurements. Rename the subsection “Carbonate reactions were not recovered under either tested indicator panel” and replace every broader formulation with the same boundary.

5. **[Major]** Figure 5 mixes the M7 Ghana audit with a companion M6 field-transfer experiment. The caption discloses this, but placing the M6 tier ablation inside an M7 result figure invites readers to treat all four panels as evidence from the same protocol. Move the M6 panel to a clearly labelled contextual figure in the Supplement, or present it in a boxed comparison with its own protocol, sample, outcome and citation. The M7 field claim must rest only on the truth-free Ghana scope audit and hold-forward analysis.

6. **[Minor]** The public-pipeline table reports the same selected F1 value, 0.4222, for all four arms, despite large changes in probability scores and log loss. Candidate recall is also 0.9815 rather than complete. These observations suggest that the selected threshold was insensitive in six cases and that the system result is conditional on candidate generation. Add thresholds and confusion counts for each arm, explain the identical F1 result, and state explicitly that the system audit tested scoring only on the candidates that were recovered.

7. **[Minor]** The conclusion says the sheaf layer “has therefore earned a bounded scientific-workflow role” at lines 1203 to 1205. “Earned” is evaluative language rather than a test result. Replace it with: “The experiments support use of the sheaf layer as a prespecified model of non-identity relations and as a global-compatibility diagnostic under the tested scalar scenarios, with global fallback when endpoint evidence is missing.”

### 1.5 Figures and tables

1. **[Major]** Figure 5 combines data and claims from two modules, as described above. Its provenance is disclosed in the caption and artifact registry, but the visual grouping still implies a common experiment. Separate the M6 panel or mark it inside the panel title and legend as external companion evidence, including the independent sample size and protocol identifier.

2. **[Minor]** The Word render is clean, but the main paper contains seven figures and nine tables in 30 pages. Table 7 spans three pages and interrupts the transition to the representation results. Retain the design table, the primary M7.3 decision table, Table 8 and Table 9 in the main paper, and move detailed metric tables to the Supplement, where the complete CSVs are already cited.

3. **[Minor]** Figures 6 and 7 are legible at full-page width but their four-panel labels and several axis labels become small in the 30-page Word layout. Increase label and tick sizes for final journal dimensions, and include the intended printed width in the figure-generation manifest so that legibility can be checked automatically.

### 1.6 References

1. **[Major]** Eighteen unique references are too few to support the paper’s broad positioning across groundwater joint inversion, age dating, reactive transport, probabilistic calibration, graph inference and sheaf methods. The present search record is strongest for software documentation and selected method citations but weak for competing groundwater graph or structured-residual approaches. Expand the literature review using the reproducible search requested above, and add the closest non-sheaf alternatives rather than only mathematical sheaf sources and one recent neural-network preprint.

2. **[Minor]** Davis and Goadrich are cited at main-text lines 416 to 420 and Supplementary lines 231 to 236 as showing that PR-AUC is more informative under imbalance. The paper’s primary statement concerns precision-recall curves under highly skewed data and the relation between ROC and PR spaces, not an unconditional ranking of scalar AUC summaries. Change the wording to “precision-recall curves can give a more informative view under class imbalance,” and retain ROC-AUC as the manuscript already does.

### 1.7 Editorial and language quality

1. **[Major]** The declarations are submission blockers. Main-text lines 1213 to 1222 contain placeholders for author contributions and competing interests; lines 1233 to 1240 and 1257 to 1260 contain data and software DOI placeholders. There is no funding statement. Replace every placeholder with final contributor roles, competing-interest and funding declarations, dataset DOI, software DOI, licence and versioned release before submission.

2. **[Minor]** Main-text line 318 reads, “The reaction solver was evaluated inference, not part of the synthetic generator.” The sentence is ungrammatical and obscures an important independence claim. Replace it with: “The reaction solver was the inference method under evaluation and was not part of the synthetic generator.”

3. **[Minor]** “Open Research” at line 1211 contains authorship, interests, data and code statements rather than only open-research material. Rename the section “Declarations and open research,” and order it as author contributions, funding, competing interests, data availability and code availability.

4. **[Major]** The current `Manuscript-Final.docx` fails LibreOffice headless conversion with `libpng error: Write Error`, although all seven embedded PNGs can be read and Microsoft Word exports a complete 30-page PDF. The supplement exports to 17 pages in Word. This is an interoperability defect in the submission package, not a scientific flaw. Rebuild both DOCX files from clean source, test them in Word and LibreOffice, and keep the package only when both applications produce complete PDFs without repair prompts or conversion errors.

### 1.8 Novelty, contribution and significance

1. **[Major]** The contribution is real but narrower than the manuscript’s software-level framing sometimes suggests. M7.4 shows exact identity-limit nesting, information in native affine maps relative to permuted maps, and better planted-conflict localisation. It does not show overall benefit against the edge-local comparator. M7.5 shows ranking gains in two planted scenarios but fails the complete calibration gate. Frame the contribution as a falsifiable benchmark and a conditional representation result, not as validation of HydroSheaf as a whole.

2. **[Major]** Development selected a local weight of 1.0 for both M7.5 hybrid arms, as reported at lines 856 to 860. The selected method therefore did not combine local and global residuals when both endpoints were observed. Its added capability was global fallback for missing endpoints plus map-sensitive conflict information. Use “local-first/global-fallback estimator” consistently and avoid presenting it as evidence that a general local/global blend is superior.

3. **[Major]** The study has no independent cross-generator replication of the M7.4/M7.5 representation result and no field truth for topology, age or reaction mechanisms. This does not erase the contribution, but it limits its likely impact to a method and capability paper. State that the next claim-bearing step is replication under a different generator family or a field dataset with independently measured connectivity, and do not infer temporal, three-dimensional, vadose-zone, vector-stalk or active-learning performance from M7.

## 2. Cross-section consistency audit

1. **[Major]** The title, Introduction line 203 and Limitations lines 1134 to 1136 say or imply one generator and a fixed twelve-case locked test, while Methods lines 348 to 390 describe a separate scalar generator with 64 and 128 locked cases. Revise every scope statement to distinguish M7.3 from M7.4/M7.5.

2. **[Major]** Introduction line 207 says “seven linked audits,” whereas Conclusion line 1180 says “six linked audits.” There are seven if the four original experiments, public-pipeline audit, M7.4 and M7.5 are counted. Use seven throughout and enumerate them once.

3. **[Critical]** Results line 856 reports that the M7.5 provenance gate passed, but the current runner fails the confirmatory lock. Restore the exact source or qualify the statement until a new locked run is completed.

4. **[Critical]** Code availability cites commit `d336e87` for M7.4/M7.5, while both later manifests report `53beb460...` and neither commit contains the later files. Replace the single commit claim with accurate, experiment-specific release identifiers.

5. **[Minor]** The manuscript calls the M7.5 arm a hybrid while reporting a selected local weight of 1.0. The methods explain the missing-endpoint fallback, but figure and table labels do not. Rename the selected arm in explanatory prose or add “selected local weight = 1.0” to Figure 7 and Table 9 captions.

## 3. Research integrity red-flag scan and adversarial verification ledger

No evidence of data fabrication, data falsification, citation fabrication or figure manipulation was found. This conclusion is based on direct checks rather than the manuscript’s own validation labels. The stored M7.4 and M7.5 output files matched every SHA-256 recorded in their run manifests. Case tables contained 64 cases by four models for M7.4 and 128 cases by seven models under each of two calibration regimes for M7.5, with no duplicate case-model keys. Bounded metrics remained within valid ranges. Recalculated model means reproduced Tables 8 and 9 exactly. The focused M7 test suite passed 12 of 12 tests. Microsoft Word rendered all 30 main-manuscript pages and all 17 supplementary pages without clipped figures or missing tables.

| Claim checked | Independent check | Verdict |
|---|---|---|
| M7.4 affine sheaf minus identity PR-AUC = +0.0854 | Recomputed from the 64-case metrics and matched the locked contrast record | Verified |
| M7.4 affine sheaf minus edge-local PR-AUC = +0.0097 with interval crossing zero; log loss = +0.0117 with adverse interval | Matched the locked contrast file and Table S7 | Verified |
| M7.5 LOO hybrid minus edge-local PR-AUC = +0.0200; Brier interval crosses zero; mean log loss is adverse | Matched the 128-case locked contrast file and Table S9 | Verified as a stored numerical result |
| M7.5 native maps outperform the frozen permuted control on the three primary metrics | Matched the locked claim decision and contrast file | Verified as a stored numerical result |
| M7.5 passed a currently reproducible confirmatory source lock | Current runner hash differs from the confirmatory lock; direct execution stops before test generation | Not verified |
| M7.4/M7.5 code is available at the cited repository commit | Files are absent from the cited commit, current HEAD and remote `main` | False in the present availability statement |

The runner mismatch is a serious provenance and transparency problem, but it is not evidence of misconduct by itself. The numerical outputs show no impossible values, unexplained duplication or selective removal. The proper response is to recover and archive the exact source, or to repeat the confirmation under a new lock and unused seeds.

## 4. Reporting standards compliance

No single CONSORT, STROBE or PRISMA-style standard governs this computational groundwater benchmark. Compliance is therefore assessed against reproducible simulation and environmental-model reporting practice. Model names, versions, seeds, case splits, estimands, metrics, blinding rules, software environment and adverse controls are reported. Official MODFLOW 6, MODPATH 7 and PHREEQC documentation is cited. Compliance is partial because exact M7.5 executable provenance is broken, M7.4/M7.5 are not in the cited repository revision, sample-size precision is not justified, later scenario families lack multiplicity correction, and persistent data and software identifiers are absent.

## 5. Reproducibility assessment

| Category | Verdict | Basis |
|---|---|---|
| Scientific design and methods | Present | Main and Supplement describe generators, evidence streams, controls, splits, estimators and metrics |
| Seeds, grids and decision gates | Present | Protocol locks and Supplement list seeds, folds, grids and claim rules |
| M7.3 predeclaration | Present | Commit `d336e87` is identified for the M7.3 protocol |
| M7.4 source lock | Present locally | All M7.4 bound file hashes match, but the files are not in the recorded Git revision or remote release |
| M7.5 exact confirmatory runner | Absent | Current source does not match the confirmatory runner hash |
| Stored output integrity | Present | Every artifact hash in the M7.4 and M7.5 manifests matches |
| Numerical regeneration from stored case tables | Present | Tables 8 and 9 means and headline contrasts reproduce |
| Automated tests | Present | Focused M7 suite: 12 passed |
| Persistent software environment | Partially present | Versions are listed, but no locked environment file and executable archive is cited for the full package |
| Document rendering | Partially present | Word renders complete output; LibreOffice conversion fails for the current main DOCX |

The overall reproducibility verdict is **partially reproducible, with one essential source artifact missing from the verifiable chain**.

## 6. Data availability and transparency assessment

| Category | Verdict | Required action |
|---|---|---|
| Data and stored result tables | Partially adequate | Deposit all claim-bearing tables, per-case files and provenance manifests with a dataset DOI |
| Analysis code | Inadequate | Recover the exact M7.5 runner, commit all M7.4/M7.5 files, tag a release and deposit a software archive |
| Supplementary materials | Adequate with formatting revisions | Retain complete tables, add corrected family-wise tests, and repair cross-application DOCX rendering |
| Author contributions | Absent | Add named CRediT roles for every author |
| Competing interests | Absent | Add the final declaration |
| Funding | Absent | Add funding sources and grant identifiers, or state that there was no specific funding |
| Data rights and licences | Partially adequate | State ownership, redistribution rights and licences for the Ghana workbook, code and generated outputs |

## 7. Citation integrity check

Four key sources were checked against their primary records. [Davis and Goadrich (2006)](https://doi.org/10.1145/1143844.1143874) support use of precision-recall curves for highly skewed classification but do not establish a universal scalar PR-AUC advantage, so the manuscript needs the wording correction identified above. [Hansen and Ghrist (2019)](https://arxiv.org/abs/1808.01513) support the description of cellular sheaf Laplacians as an extension of graph spectral constructions. [Caralt et al. (2026)](https://arxiv.org/abs/2603.05395) report competitive identity-sheaf baselines across five heterophilic graph-learning benchmarks, and the manuscript correctly labels this as a non-peer-reviewed, non-groundwater comparison. [Vehtari et al. (2021)](https://arxiv.org/abs/1903.08008) concern MCMC convergence and effective sample size, not single-set importance sampling; the manuscript accurately discloses that its ESS threshold is only an analogy.

No fabricated citation, citation cartel or clear salami-slicing pattern was detected. Self-citation could not be assessed because the manuscript does not yet name its authors. Figure 5 transparently identifies the M6 companion result, but the cross-module reuse should be separated visually as requested.

## 8. Scoring

| Section | Score | Maximum | Reason |
|---|---:|---:|---|
| Title and abstract | 3 | 5 | Accurate negative-result reporting, but the single-generator title is wrong and the abstract is too compressed |
| Introduction | 6 | 10 | Clear problem and aims, but originality search and audit hierarchy are insufficiently documented |
| Methodology | 11 | 25 | Strong controls and truth separation, offset by a failed confirmatory source lock, incorrect code provenance, missing precision rationale and later multiplicity gaps |
| Results and discussion | 20 | 30 | Stored numbers reproduce and adverse outcomes are reported, but practical magnitude, scenario multiplicity and scope require revision |
| Figures and tables | 3 | 5 | Visually clean in Word, but dense, cross-module Figure 5 and small multi-panel labels need correction |
| References | 3 | 5 | Key citations exist, but coverage and one PR-curve claim need improvement |
| Editorial and language quality | 3 | 5 | Prose is mostly clear, but declarations, count inconsistency, one grammar error and DOCX interoperability remain |
| Novelty, contribution and significance | 10 | 15 | A meaningful conditional benchmark contribution, limited to controlled scalar scenarios and not a general HydroSheaf validation |
| **Total** | **59** | **100** | **Major revision required** |

## 9. Critical flaw assessment

No fabrication override was triggered. No design flaw was found that permanently invalidates the stored numerical results, and the bounded conclusion generally follows from those results. The exact M7.5 confirmatory runner is missing from the verifiable chain, but this is addressable by restoring the content-hashed source or by repeating confirmation under a new protocol with unused seeds. The manuscript also makes a meaningful contribution by testing identity nesting, native-map information, conflict localisation, negative transfer and a strong edge-local comparator. The critical methodology, unsupported-conclusion and insufficient-contribution rejection overrides are therefore not triggered. Publication must still stop until the source-lock and repository claims are corrected.

## 10. Decision

**Return for major revisions.** The study is salvageable and contains a defensible conditional contribution, but the present package cannot be submitted while the M7.5 runner fails its lock, the cited repository commit omits M7.4/M7.5, the generator scope is internally inconsistent, and the required declarations and persistent archives are absent.

## 11. Structured recommendation answers

1. **Recommendation:** (d) Return for major revisions.

2. **Study design appropriateness:** (b) No, but addressable with revision. The controls are strong, but the authors must separate the two generator systems, justify precision, correct later multiplicity families and restrict reaction claims to the planted ontology.

3. **Methods reproducibility:** (b) No, but addressable. The exact M7.5 runner must be restored or a new locked confirmation must be run, and the complete package must be committed and archived.

4. **Statistics and uncertainty treatment:** (b) No. M7.3 multiplicity is handled, but M7.4/M7.5 scenario families and minimum important effects are not.

5. **Guidance on overstated claims:** (a) Yes, specific replacements and scope corrections are provided above.

6. **Presentation clarity:** (c) Needs language corrections before publishing.

7. **Editorial and language quality:** (c) Requires substantial editorial revision before publication because declarations, availability statements and internal scope statements are unfinished.

8. **Novelty and contribution:** (c) Marginal-to-adequate contribution that needs stronger, accurate framing and independent replication before broader claims.

9. **Reproducibility:** (b) Partially reproducible, with the exact M7.5 runner and persistent release missing.

10. **Research integrity:** (b) A provenance concern is flagged, but no evidence of fabrication or manipulation was found.

11. **Data availability and transparency:** (c) Inadequate until data, code and results are deposited under persistent identifiers and the exact source chain is restored.

12. **Citation integrity:** (b) Minor concerns. One methodological claim needs more precise wording, the originality search needs a reproducible record, and no manipulation pattern was found.

## 12. Summary letter to the authors

The manuscript does something scientifically useful: it does not convert every reduction in uncertainty into a success claim. The adverse controls, strong edge-local comparator, identity-limit equality test, stored negative results and explicit field-data boundary make the central conditional conclusion more credible than a simple demonstration paper. The stored M7.4 and M7.5 numerical results reproduce from the case tables, and the paper correctly rejects general superiority.

The current submission package nevertheless has one non-negotiable provenance defect. The M7.5 runner no longer matches its confirmatory lock, and the repository commit cited for M7.4/M7.5 does not contain those experiments. Recover the exact runner or repeat the confirmation under a new lock and unused seeds, then publish the full code and result package under a versioned DOI. Until that is done, the paper cannot claim a presently reproducible one-time confirmation.

The scientific framing also needs correction. The assembled manuscript contains two generator systems rather than one, the later scenario claims need family-wise statistical treatment, the practical size of the age and PR-AUC gains must be reported, and reaction conclusions must be limited to the two tested panels and six planted archetypes. Separate the M6 evidence from the M7 Ghana audit and replace the six-versus-seven audit inconsistency.

Finally, complete every declaration, deposit the data and software, correct the availability statement, and rebuild the DOCX package so that it renders in both Word and LibreOffice. With those revisions, the manuscript can support a publishable claim: affine maps encode relation information beyond identity connectivity and help in selected planted conflict or missing-endpoint settings, but they have not established general predictive superiority over a strong weighted graph.
