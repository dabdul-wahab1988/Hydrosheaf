# Peer review: M2.3, HydroSheaf claim-gated evidence-integration framework

**Target journal:** Computers & Geosciences
**Manuscript:** `M2.3/manuscript/Manuscript-M2.3.md` (7,491 prose words, 7 figures, 4 tables, 28 references)
**Review date:** 2026-08-02

---

## 1. Title and Abstract

1. **[Major] The abstract runs to 278 words and is dominated by numbers.** Computers & Geosciences expects roughly 250. More importantly, the abstract reports eleven separate quantities before it reaches the claim boundary, and a reader scanning it will struggle to extract the one idea that makes this paper worth reading: that the framework locates its own identifiability boundary. Revision instruction: cut to 250 words by removing the transport median absolute error and the 91.7% model-selection figure, both of which are secondary, and move the identifiability framing into the second sentence rather than the final one.

2. **[Major] "Claim-gated" appears in the title without appearing anywhere in the indexing vocabulary of the field.** The term is defined nowhere in the abstract, and a reader searching for groundwater inverse modelling, equifinality or uncertainty calibration will not find this paper by that phrase. Revision instruction: retitle along the lines of "Identifiability-aware evidence integration for groundwater connectivity, residence-time and reaction inference in data-limited aquifers", which keeps the contribution and uses terms the field actually searches on, and retain "claim-gated" as a defined term inside the Introduction.

3. **[Minor] The abstract states the software name and version nowhere, and gives no availability pointer.** For a software paper at this journal that is an omission a desk editor will notice. Revision instruction: add "HydroSheaf v0.5.1" at first mention and one clause pointing to the archived repository.

---

## 2. Introduction

1. **[Major] No hypotheses are stated, only objectives.** The paper poses a central question at line 27 ("Which groundwater inference targets are recoverable...") and then lists six objectives, but never commits to a predicted answer. This matters because the paper's most interesting result is a negative one, and a reader cannot tell whether the authors expected reaction-extent recovery to fail or discovered it and reframed the paper around it. That is the HARKing question, and the manuscript currently invites it. Revision instruction: state explicitly what was predicted before the benchmark was run, and if the reaction-extent failure was not predicted, say so plainly in the Results. The archived analysis plan can support that statement.

2. **[Major] The literature is anchored in 1994 to 2001 and thin on the last five years.** Beven, Kennedy and O'Hagan, and Oreskes are correctly used, but the only reference after 2021 is GraphFlow. Work on machine-learning surrogates for groundwater flow, graph neural networks for subsurface transport, and conformal prediction applied to environmental models has grown considerably since 2022, and a reviewer at this journal will read its absence as a gap. Revision instruction: add three to five references from 2022 onward covering data-driven groundwater surrogates and uncertainty quantification, and state in one sentence why this framework is not a competitor to them.

3. **[Minor] The framing of "sheaf-inspired" is deferred to the Discussion.** A reader encountering the software name in the title will expect to learn in the Introduction what the sheaf connection is and is not. Revision instruction: move the sheaf-inspired qualifier from Section 5.2 to the paragraph introducing the framework, and keep it to one sentence.

---

## 3. Methodology

1. **[Critical] The 100-realisation benchmark is not what the manuscript implies. It is 100 measurement-noise replicates of a single fixed scenario.** Inspection of the shipped records shows ten unique edges, the same ten in every realisation, and true reaction extents that are identical across all 100 realisations for all 70 edge-by-reaction combinations. True transport values take only six distinct values across the entire benchmark, fixed per edge and parameter. The manuscript nonetheless reports "2,100 active reaction terms" and "550 scored parameter rows" as though these were independent observations. They are not: the effective sample size for anything concerning the truth configuration is one scenario, and the 100 replicates inflate n by a factor of 100 without adding information about generalisation. This directly affects the headline negative result, because the coefficient of determination is computed over a truth vector containing roughly 21 distinct values repeated 100 times, and its denominator is entirely between-edge variance from a single aquifer. Revision instruction: state the design explicitly in Section 3.4 (ten edges, one locked truth specification, 100 noise replicates, seeded per realisation), report n as the number of distinct truth values alongside the number of scored rows, and either compute uncertainty by clustering on edge-by-reaction rather than on row, or drop numerical precision claims that depend on the inflated n. Do not report the benchmark as "100 realisations" without the qualifier.

2. **[Critical] The regularisation weights governing the reaction inversion are never given, and no evidence is presented that they were tuned.** Supplementary Methods S4 states the objective and then refers to "the declared regularisation weights" without providing values, a selection procedure, or a sensitivity analysis. Since the paper's pivotal claim is that this inversion fails, the reader cannot distinguish a structural non-identifiability from a poorly chosen penalty. An L1 weight set too high will drive extents toward zero and produce exactly the pattern reported: a negative coefficient of determination on active terms combined with heavy leakage. Revision instruction: report the numerical values of both penalties, state how they were selected, and add a sensitivity analysis across at least a decade of each penalty showing that the negative recovery persists. Without that analysis the central claim is not established.

3. **[Major] Measurement uncertainty for the field data was assigned rather than estimated, and the assigned values are not stated.** Section 2.5 notes that no replicate analyses or laboratory covariance metadata are available, and Section 3.1 refers to "declared analytical conventions", but no numbers appear. Every calibrated interval in the field application depends on these. Revision instruction: tabulate the assumed relative and absolute uncertainty per determinand in the Supplementary Methods, and cite the source of each convention.

4. **[Major] The competence-matched baselines are described in the abstract as comparators but never specified in the main text.** The reader is told the age baseline achieved a mean absolute error of 7.65 years and the reaction baseline a log loss of 2.852, but not what those baselines are. A comparator that is weak by construction makes any contrast meaningless. Revision instruction: describe both baselines in Section 3.2 in two or three sentences each, name the artefact registry entries that define them, and state what each baseline was and was not permitted to see.

5. **[Major] The six locked cases are treated as the programme's evidence base without any power justification.** Section 3.4 states that six cases were held out and Section 4.3 reports gate results over them, but the manuscript never addresses whether six cases can discriminate a passing from a failing method for any of the four gates. Revision instruction: either report the minimum detectable difference for each gate at six cases, or explicitly frame the programme as a smoke test for the claim architecture rather than as a performance evaluation.

---

## 4. Results and Discussion

1. **[Critical] The conclusion that point reaction extents are "not recoverable from this evidence" is not supported by the evidence presented.** What the benchmark shows is that one implementation, HydroSheaf's own sparse inverse solver, with unreported regularisation weights, failed to recover extents on one synthetic scenario. The manuscript then generalises this to a property of the inference target itself, and builds Section 5.1 on that generalisation. The contrast with the kinetic case is instructive and works against the authors here: for kinetics they present a structural argument (collinear sensitivity columns, numerical rank falling to one under k times A), and that argument genuinely establishes non-identifiability. No equivalent argument is offered for reaction extents. Revision instruction: either establish non-identifiability structurally by reporting the rank or condition number of the reaction dictionary matrix A under the observation set, or run a second, independent inverse solver (PHREEQC inverse modelling is the obvious choice and is already a declared dependency) and show that it fails similarly. Absent one of these, reword the claim throughout to "the implemented sparse inversion did not recover point extents under the tested generator" and revise Section 5.1 accordingly.

2. **[Critical] The 54.1% leakage figure has no baseline and is therefore uninterpretable.** The manuscript reports that 54.1% of truth-zero terms received a recovered magnitude above 0.05 mmol/L, and treats this as evidence of failure. But no comparison is given to what a null procedure would produce, and the threshold of 0.05 mmol/L is not justified relative to the simulated measurement noise (major-ion relative sigma of 0.035 in the benchmark configuration). If propagated analytical noise alone produces apparent extents near 0.05 mmol/L, then a large leakage fraction is expected and is not evidence about the method. Revision instruction: report the leakage fraction for a null model that assigns extents from noise alone, justify the 0.05 mmol/L threshold against the propagated measurement uncertainty, and report leakage as a function of threshold rather than at a single cut.

3. **[Major] The Talensi charge imbalance is diagnosed without reference to the source publication.** Section 2.3 attributes a median charge-balance error of −29.9% to unmeasured cations or an unmatched alkalinity determination, having tested two alternative unit interpretations. That analysis is sound as far as it goes, but the source study is cited in the manuscript and reports its own hydrochemistry; the authors did not check whether the original publication reports charge balance, whether the discrepancy arises in transcription into the project file, or whether the original reports species this dataset omits. A reviewer will ask, and if the imbalance is a transcription artefact the entire Talensi discussion changes. Revision instruction: compare the project file against the source publication's reported tables, state whether the imbalance is present in the original, and correct or annotate accordingly.

4. **[Major] The field candidate-edge results silently switch denominators between 208 and 161.** Section 4.5 reports 208 candidate edges and an overall median closure of 0.713, then reports per-site medians over 95 Talensi and 66 Lower Anayari edges, which total 161. The 47 excluded edges are not described, and inspection shows they are exactly the edges assigned an evaporation rather than a mixing transport model, which is not a random exclusion. Revision instruction: state why 47 edges carry no site-identifying end-member, report whether their closure distribution differs from the assigned edges, and use one denominator consistently or justify the switch at each use.

5. **[Major] The external age comparison discards nearly half its data without comment.** Section 4.4 reports 675 sites with finite paired values, and the Supplementary Methods note that this is 675 of 1,272 rows. Almost 47% of the reference set produced a non-finite estimate or reference value and was dropped. If those failures are not missing at random, and age-inference failures very plausibly concentrate in the hardest sites, the reported R² of 0.482 is optimistic. Revision instruction: characterise the 597 dropped rows against the retained ones on depth, aquifer group and reported age, and state whether the framework or the reference produced the non-finite value in each case.

6. **[Minor] The "mixed" residence-time class carries the worst relative error in both modes and is not discussed.** Table values show median absolute relative error of 0.574 and 0.344 for mixed waters, roughly triple the young-water figures, and the Results note this only in passing as "expected". Since binary mixing is the common real-world case, this deserves more than a clause. Revision instruction: add two sentences in Section 5.1 on what the mixed-class degradation implies for field deployment.

---

## 5. Tables and Figures

1. **[Major] Figure 2(a) is described as showing sampling locations but is a bare longitude-latitude scatter.** For a geoscience journal this is not a map. There is no basemap, no scale bar, no north arrow, no administrative or geological context, and no indication of the Birimian and Voltaian formations the Data section relies on. A reader cannot relate the three datasets to the hydrogeology. Revision instruction: replace with a proper map including a geological or administrative base layer, a scale bar, a north arrow and an inset locating northern Ghana within West Africa.

2. **[Minor] Figure 4(b) plots 2,800 points with heavy overplotting at alpha 0.30.** The two methods substantially occlude one another in the dense mid-range, which is exactly where the comparison matters. Revision instruction: add marginal error distributions or replace the point cloud with binned summaries, retaining points only in sparse regions.

3. **[Minor] Table 4 discloses disagreements with documents the reader cannot see.** Publishing a comparison against "M2 abstract" and "M2 summary", neither of which is a published artefact, will confuse a reviewer who has no access to them. The disclosure instinct is right, but the framing is wrong. Revision instruction: reframe as a reproducibility note stating that all values were recomputed from primary records, and move the specific comparisons to the supplementary material with a clear statement that the earlier documents were internal and unpublished.

---

## 6. References

1. **[Major] Twenty-eight references is thin for a paper of this scope.** The manuscript spans groundwater age inference, inverse geochemistry, graph methods, uncertainty calibration, selective prediction, experimental design and West African hydrogeology, and cites roughly four sources per topic. Revision instruction: expand to 45 to 60 references, with the additions concentrated in post-2022 work as noted under the Introduction.

2. **[Minor] West African hydrogeology is supported by two references.** For a paper whose field application is entirely Ghanaian, MacDonald and Banoeng-Yakubo are not enough. Revision instruction: add regional literature on Birimian and Voltaian aquifer properties, recharge estimation in the Sudano-Sahelian belt, and fluoride occurrence in northern Ghana.

---

## 7. Editorial and Language Quality

1. **[Minor] Several paragraphs carry too many clauses.** The paragraph beginning "The locked programme is a controlled-synthetic run" packs three generator families, case splits, five stress transformations, calibration restrictions and truth-release timing into one unit. Revision instruction: split into design, leakage control and scoring paragraphs.

2. **[Minor] The paper uses "closure", "in-sample closure" and "chemistry closure" for the same quantity.** Revision instruction: fix on "in-sample chemistry closure" at first use, then "closure" thereafter, and apply the same term in figure captions.

3. **[Minor] Section 5.1 opens with a sentence that states the paper's thesis and is easy to miss.** "The results do not divide into successes and failures of the software" is the most important sentence in the Discussion and it sits mid-paragraph in effect. Revision instruction: give it its own short opening paragraph.

---

## 8. Novelty, Contribution and Significance

1. **[Major] The individual components are established; the assembly is what is new, and the manuscript does not always defend that distinction.** Sparse inverse geochemistry, lumped-parameter age models, model discrepancy, discrete model averaging, conformal-style calibration and Bayesian experimental design are all prior art, and the manuscript cites them correctly. The contribution is the claim-gating architecture and the identifiability boundary it exposes. Section 4.1, however, reads as a capability inventory and risks presenting assembly as innovation. Revision instruction: cut Section 4.1 to a short paragraph, move the capability inventory entirely to Table 2, and use the recovered space for the sensitivity analysis requested above.

2. **[Major] The most transferable result is under-sold.** The finding that strontium and dissolved silica move the non-identifiable fraction from 51.9% to 0.6% while stable isotopes and fluoride move it not at all is a concrete, actionable survey-design result that any hydrogeologist working in comparable terrain could apply tomorrow. It currently sits in one paragraph of Section 4.5 and one clause of Section 5.6. Revision instruction: promote this to a named subsection in the Discussion and state the implication for survey budgeting explicitly.

3. **[Minor] The negative result is the paper's strongest claim to originality and is framed defensively.** Papers reporting that a widely used class of inversion does not recover its target are rare and valuable. Subject to the two Critical issues above being resolved, this should be foregrounded. Revision instruction: once the sensitivity analysis and structural argument are in place, state the negative finding in the abstract's opening rather than its middle.

---

## Cross-Section Consistency Audit

The audit was performed across the abstract, all main sections, the tables, the figure captions and the Supplementary Methods.

Headline values agree throughout. Age coverage 0.9655, specialist MAE 4.24 against baseline 7.65, reaction log loss 0.896 against 2.852, coverage 0.739, false commitment 0.038, kinetic RMSE 0.119, conditional identification 1.0, overall identification 0.167, abstention 0.833, integrated utility 1.478 against 0.246 and −0.004, held-out uncalibrated R² 0.482, calibrated R² 0.962, topology recall 0.845 and precision 0.487, and the measurement-value figures 51.9% and 0.6% are consistent everywhere they appear. All were checked against the exported records and reproduce.

Three inconsistencies were found. First, the candidate-edge denominator switches between 208 and 161 within Section 4.5 without explanation, as noted above. Second, the abstract states "264 samples from northern Ghana" while Table 1 lists 160, 41 and 63, which do sum to 264 but only after the Northern Ghana panel has been reduced to its primary panel; a reader comparing the abstract with the earlier internal figure of 320 will be confused without the Section 2.2 explanation, which comes later. Third, Section 4.2 reports the transport median absolute error as 0.058 while the exported summary carries 0.0577, and the conclusion rounds to 0.058; this is consistent but the manuscript should fix on a single precision convention.

The claim boundary is consistent: field validation, measured seasonality, universal superiority and formal sheaf theory are denied wherever they could be inferred.

## Research Integrity Red Flag Scan and adversarial verification

No evidence of fabrication, falsification, image manipulation or citation fabrication was found. The manuscript reports negative and limiting outcomes prominently rather than burying them, including a negative coefficient of determination, a 54.1% leakage rate, a 0.167 overall identification rate, a precision below 0.5, and a held-out R² roughly half the calibrated value.

Two features deserve positive note. The Northern Ghana seasonal reconstruction is disclosed in the Data section, the Statements and the Supplementary Methods, with the specific diagnostic evidence given, and the affected analyses are excluded rather than caveated. Table 4 discloses that earlier internal reporting of the same benchmarks is not reproducible. Both disclosures work against the authors' immediate interest, which is the correct signal.

The adversarial verification ledger is as follows. Age, reaction, kinetic and integrated gate values were traced to the locked performance-gate record and reproduce to the reported rounding: **verified**. The external age parity metrics were recomputed independently from the parity record and reproduce the archived summary exactly at n = 675: **verified**. Topology precision, recall and F1 were recomputed from the confusion counts and reproduce: **verified**. Field charge-balance tiering was recomputed from the raw ion concentrations and reproduces the previously reported 294/19/7 split on the full workbook: **verified**. The measurement-value ablation reproduces from the source record: **verified**. The claim that point reaction extents are not recoverable is **unsupported as stated**, for the reasons given under Results issue 1; the underlying numbers are correct but the inference from them is not established. The reported sample sizes of 2,100 and 550 are **contradicted** as measures of independent information, for the reasons given under Methodology issue 1.

No integrity concern is raised. The two adverse findings are analytical overreach and inadequate design disclosure, both correctable, not misconduct.

## Reporting Standards Compliance

No clinical, observational or laboratory reporting standard applies; this is a computational software-validation study with a secondary audit of existing hydrochemical records. STROBE is not appropriate because no new sampling or epidemiological inference was performed. The manuscript meets most computational expectations: software version, generator families, development and locked splits, truth-blind selection, observation stresses, metrics, costs, source hashes and an explicit claim boundary are all stated. It does not meet the expectation that hyperparameters governing a headline result be reported, as noted under Methodology issue 2.

## Reproducibility Assessment

Software identity and package boundary: **Present.** Version, entry points and the separation between package code, control scripts and generated evidence are stated.

Analysis provenance: **Present.** Every reported value is generated by named scripts from primary records, and the derivation and figure layers are separated with Python owning computation and R owning rendering.

Execution provenance: **Partially present.** The locked run's manifest records the generating revision, 31 source hashes and an uncommitted-tree flag. The manuscript discloses this rather than concealing it, but a clean-tree regeneration has not been performed.

Parameter reporting: **Absent for the reaction inversion.** The regularisation weights driving the paper's pivotal result are not given.

Statistical reporting: **Partially present.** Metrics, calibration, selective risk and bootstrap settings are recorded, and the six-case limitation is explicit. The pseudo-replication in the 100-realisation benchmark is not disclosed.

Data availability: **Partially present.** Two of three field datasets trace to publications; the third is disclosed as author-compiled with a reconstructed attribute.

## Data Availability and Transparency Assessment

Raw field data: **Partially adequate.** Files, inventory and crosswalk are present, and the Northern Ghana provenance limitation is stated. Redistribution terms for the two published datasets are not confirmed.

Analysis code: **Partially adequate.** The derivation, figure and assembly scripts are present and named. The repository DOI and URL are placeholders, which for this journal is a submission blocker rather than a stylistic gap: Computers & Geosciences requires a public, documented, non-archived repository at submission and treats private or incomplete repositories as grounds for rejection.

Supplementary materials: **Adequate,** with the exception of the missing hyperparameters.

Author contributions: **Inadequate.** The Statements section carries a placeholder rather than named CRediT roles.

## Citation Integrity Check

Three citations were verified against their sources. GraphFlow v1.0 was confirmed as Geoscientific Model Development volume 18, pages 7147 to 7163, 2025, with the DOI as given, and the manuscript's characterisation of it as a graph-based approximation of contaminant transport is accurate. Chegbeleh, Akurugu and Yidana 2020 was confirmed in The Scientific World Journal with the DOI as given; the source reports silicate and carbonate weathering as the main controls and mixed cation water types, which corroborates rather than contradicts the manuscript's facies finding of 65% mixed-bicarbonate at Talensi. Kennedy and O'Hagan 2001 was confirmed in Journal of the Royal Statistical Society Series B with the stated volume and pages, and is used correctly for model discrepancy.

Every citation key in the manuscript and supplement resolves to a bibliography entry, and the bibliography contains no uncited entries. No citation cartel, excessive self-citation, or salami-slicing signal was detected. Self-citation is effectively absent, which given the existence of earlier internal manuscripts on the same software is appropriate, since those were never published.

## Scoring

| Section | Score | Maximum | Basis |
|---|---:|---:|---|
| Title and Abstract | 3 | 5 | Over length, dense, non-discoverable title term. |
| Introduction | 8 | 10 | Clear problem and boundary; thin recent literature, no stated hypotheses. |
| Methodology | 15 | 25 | Two critical issues: undisclosed pseudo-replication and unreported hyperparameters. |
| Results and Discussion | 19 | 30 | Traceable and candid, but the pivotal negative claim overreaches its evidence and the leakage figure lacks a baseline. |
| Tables and Figures | 4 | 5 | Well constructed; the location panel is not a map. |
| References | 3 | 5 | Sound sources, too few, and weighted to pre-2002. |
| Editorial and Language Quality | 4 | 5 | Clear UK English with some dense paragraphs. |
| Novelty, Contribution and Significance | 12 | 15 | Real architecture-level contribution; the transferable survey-design result is under-sold. |
| **Total** | **68** | **100** | **Moderate revisions needed.** |

## Critical Flaw Assessment

Three conditions were evaluated explicitly.

On critical methodology flaw: the pseudo-replication in the 100-realisation benchmark is a serious reporting and analysis flaw, but it does not require new data collection. The recorded data support a correct analysis clustered at the level of distinct truth values, and the qualitative direction of the transport and age results is unlikely to change. The reaction-extent result is more exposed, because its strength genuinely depends on the inflated n and on unreported hyperparameters. This is correctable by reanalysis and additional sensitivity work using data already in hand, so it does not trigger the override, but it is Critical and the paper must not proceed without it.

On unreliable conclusions: the claim that point reaction extents are not recoverable overreaches the evidence, since only one solver on one scenario was tested. This is fixable by rewording plus either a rank analysis of the reaction dictionary or a comparison against an independent inverse solver. Because the fix is available and the underlying numbers are sound, the override is not triggered, but Section 5.1 currently rests on an unestablished premise.

On insufficient contribution: the contribution is real. An implemented architecture that reports its own identifiability boundary, together with a concrete survey-design result showing which two determinands buy identifiability, is a meaningful addition. The override is not triggered.

No critical flaw requiring rejection was identified.

## Decision

**Return for moderate revisions.** The core study is sound and the claim discipline is unusually good, but two Critical issues must be resolved before the paper can be considered: the benchmark design must be reported honestly and reanalysed at the correct level of independence, and the pivotal negative claim must be supported by a hyperparameter sensitivity analysis together with either a structural rank argument or an independent solver comparison.

## Structured recommendation answers

1. **Recommendation:** Return for moderate revisions.
2. **Study design appropriateness:** No, but addressable. The benchmark's pseudo-replication and the absence of a second inverse solver can both be handled with data in hand.
3. **Methods reproducibility:** No, but addressable. Reaction-inversion hyperparameters and field measurement-uncertainty conventions must be reported, and the clean-tree regeneration completed.
4. **Statistics and uncertainty treatment:** No. Uncertainty is computed at the row level over data that are replicated at the scenario level, which overstates precision.
5. **Guidance on overstated claims:** Yes, specific rewrites suggested for the reaction-extent claim and the benchmark sample sizes.
6. **Presentation clarity:** Yes, with light editing.
7. **Editorial and language quality:** Generally well written, needs a light editing pass.
8. **Novelty and contribution:** Adequate contribution with minor gaps.
9. **Reproducibility:** Partially reproducible; key parameters missing.
10. **Research integrity:** No concerns. Disclosure practice is better than typical.
11. **Data availability and transparency:** Partially transparent. The placeholder repository is a submission blocker at this journal.
12. **Citation integrity:** No concerns.

## Summary letter to the authors

The strongest thing about this manuscript is that it tells you where it stops working. The claim architecture is not decoration: reaction-family inference passes a calibrated gate while point extents fail, kinetics predicts while abstaining on confounded parameters, and the external age comparison reports its honest number rather than its flattering one. The disclosure that calibrating to a reference measures emulation rather than agreement, illustrated with an R² that moves from 0.482 to 0.962, is a point many benchmarking papers get wrong and should be read widely. The Northern Ghana provenance disclosure and Table 4 both work against your immediate interest, and that is to your credit.

Two things must be fixed before this can go out. First, the 100-realisation benchmark is not 100 independent scenarios. It is ten edges with fixed truth values and 100 measurement-noise replicates, and reporting 2,100 active terms and 550 parameter rows as though they carried 2,100 and 550 units of information overstates your precision by roughly an order of magnitude. Report the design plainly and cluster your uncertainty on distinct truth values.

Second, and more consequentially, your pivotal claim is that point reaction extents are not recoverable, and Section 5.1 is built on it. What you have shown is that your own sparse solver, with regularisation weights you do not report, failed on one scenario. That is not the same thing. You know this, because for the kinetic case you did it properly: you gave a structural argument about collinear sensitivity columns and rank, and that argument is convincing. Do the equivalent for reactions. Report the condition number or rank of the dictionary matrix under your observation set, or run PHREEQC inverse modelling as a second solver and show it fails too. Add a sensitivity analysis over the regularisation penalties. If the negative result survives, you have a genuinely valuable paper and you should lead with it. If it does not survive, you have learned something more important than a publication.

Beyond those, tighten the abstract to 250 words, replace the longitude-latitude scatter with a real map, expand the references and rebalance them toward the last five years, and reconcile the 208 versus 161 edge denominators. Promote the strontium and silica ablation: moving non-identifiability from 51.9% to 0.6% with two determinands, while isotopes and fluoride buy nothing, is the most directly useful result in the paper and it is currently buried. Finally, the placeholder repository DOI has to be replaced with a live, documented, public repository before submission, because this journal treats that as a condition of review rather than a formality.
