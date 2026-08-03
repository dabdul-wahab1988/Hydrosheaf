# M7 manuscript decisions (r2m Track B)

## 2026-07-31 — M7.6 auxiliary M3 mechanism claim decision

- Executed the separately locked M7.6 auxiliary controlled-synthetic diagnostic
  under `RUN-M7-6-M3-MECHANISM-20260731-01` after freezing the v2 generator,
  analysis engine, runner and pre-test amendment in commit `b06523e`.
- The truth-blind execution gate passed. The environmental-isotope age control
  was exactly zero: adding E to N for T2 changed MAE by 0.0000 years with a
  95% interval of [0.0000, 0.0000].
- Severe shared nuisance increased full-panel synthetic infeasibility by
  +0.2882 (95% CI +0.2118 to +0.3646), but the redox-stratified CFC-11
  contrast failed its predeclared CFC-12 specificity control (+0.7188 versus
  +0.7396). The selective CFC-11 mechanism is therefore **not supported under
  this generator**.
- The result is retained as an auxiliary controlled-synthetic M3 diagnostic.
  It is not field validation, does not resolve the USGS cause, and does not
  revise any M7.3–M7.5 result or general claim boundary. The manuscript was
  updated only after this decision; the new text reports the negative
  mechanism decision explicitly.

## 2026-07-30 — Replace milestone codes with scientific terminology

- Removed internal labels such as M7.2--M7.5 from reader-facing manuscript,
  supplement, table and figure language.
- Adopted the terminology contract in
  `manuscript/READER_FACING_TERMINOLOGY.md`.
- Retained exact run IDs, locked artifact filenames and paths solely for
  reproducibility; these identifiers are not used as scientific names.

## 2026-07-30 — Full technical response to independent manuscript review

- Restored the exact M7.5 runner required by the confirmatory lock (SHA-256
  `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`) and
  archived it by content hash. Locked-test results were not rerun.
- Treated all 120 published M7.4 contrasts and all 560 published M7.5 contrasts
  as complete families using 10,000 shared case-block bootstrap resamples and
  two-sided max-z simultaneous 95% intervals.
- Withdrew overall and scenario-ranking superiority claims against edge-local.
  Retained exact identity nesting, native-versus-permuted map information,
  family-wise supported conflict localisation and missing-endpoint fallback.
- Recorded post-review, non-field-validated planning margins and 20,000-
  simulation precision audits without presenting them as prospective power.
- Rebuilt Figure 5 from M7 evidence only; no M6 result appears in the main
  manuscript. Figures 6 and 7 use the simultaneous intervals and record a
  7.08-in intended width with minimum 8-pt labels.
- Reduced the main paper to four tables and moved detailed results to
  Supplementary Tables S1--S13.
- Restricted reaction claims to discrimination among the six planted
  archetypes under the two tested indicator panels and one noise model.
- Separated the M7.3 MODFLOW/MODPATH generator from the M7.4/M7.5 scalar affine
  generator throughout the title, abstract, methods and limitations.
- The M7.4/M7.5 versioned release, data/software DOI, author CRediT roles,
  funding statement and competing-interest declaration remain submission
  blockers requiring external author or repository action; none was invented.

## 2026-07-24 — Revision round 1: pre-submission peer review and response

A full pre-submission peer review (`manuscript/M7_Manuscript_Review_Report.docx`,
score 47/100, major revision, no critical-flaw override triggered) identified
several issues that were addressed in this revision:

- **Unverifiable predeclaration claim (Critical).** Added the protocol-freeze
  commit reference (d336e87) to Methods so "predeclared"/"locked" is now a
  checkable claim rather than an assertion.
- **Statistical robustness of the headline age/hydraulics PR-AUC contrasts
  (Critical).** Wrote `scripts/multiplicity_check.py`, which recomputes an
  exact sign-flip permutation test directly from the locked per-case metrics
  (`results/m7_3_locked/evidence_case_metrics.csv`, no locked value altered)
  and applies a Benjamini-Hochberg correction across the full family of 24
  reported contrasts. Result: the native age and native hydraulics PR-AUC
  contrasts do not survive correction (adjusted p = 0.070 both), while their
  entropy contrasts (and, for hydraulics, calibration contrasts) do. This is
  reported as new Supplementary Table S6 and folded into Results, Discussion,
  Abstract, and Conclusion, refining rather than reversing the paper's
  central claim.
- **Two misapplied citations.** `macdonald2012africa` no longer supports a
  claim about data/evidence heterogeneity it does not make (reworded to the
  resource-distribution claim it actually supports). `vehtari2021rhat` is no
  longer cited as direct endorsement of the importance-sampling ESS >= 400
  threshold; `kish1965survey` was added for the ESS formula itself, and the
  Vehtari citation is now explicitly scoped as an analogy to a related but
  distinct multi-chain MCMC diagnostic, with the unresolved gap recorded as
  Discussion Limitation 5.
- **Single-generator overgeneralisation.** Added "in this generator"/"under
  the tested conditions" qualifiers throughout Discussion and Conclusion;
  added a fifth Limitations point on the untested sensitivity of fixed
  constants (conflict threshold, order scale, noise level, ESS rule).
- **Missing transparency information.** Added a new `06-availability-statements`
  section (Author Contributions, Competing Interests, Data Availability, Code
  Availability with verified software versions: MODFLOW 6.7.0, MODPATH
  7.2.001, FloPy 3.10.0, scikit-learn 1.9.0, numpy 2.4.6, scipy 1.17.1,
  pandas 2.3.3, Python 3.14.6) with placeholders flagged for the submitting
  authors to complete (author names/affiliations and a persistent
  repository DOI do not exist yet).
- **"Summary decision table" contained no table.** Added Table 7
  (`tables/publication/table7_summary_decision_table.csv`), synthesising
  every contrast's classification and reflecting the multiplicity-corrected
  results.
- **Introduction scope and citation diversity.** Softened the "no widely
  used benchmark" absence claim to name the specific gap; added a Bayesian
  network/causal-discovery positioning paragraph; qualified the Linde et al.
  (2006) analogy (two continuous geophysical fields vs. this study's
  discrete evidence-panel fusion); noted the contested status of the
  Davis-Goadrich PR-AUC-over-ROC-AUC claim.
- **Editorial fixes.** Standardised PR-AUC decimal precision between prose
  and Table 3 (3 dp, matching the table); standardised isotope notation
  ("³H + ³⁹Ar"); reduced repeated "central guardrail finding" phrasing across
  Introduction/Results/Discussion; updated Figure 2/3 captions to state
  which panels carry 95% CI error bars (verified against the actual
  matplotlib `errorbar` calls in `scripts/make_m7_3_publication_assets.py`
  rather than assumed).

Not addressed in this round (flagged for a future revision or for the
submitting authors): replication across additional generator geometries;
an independent per-replicate PHREEQC-bound sensitivity check for the
reaction bootstrap; a full sensitivity analysis across the adverse-control
and decision-rule constants; depositing code/data at a persistent DOI;
completing author names, affiliations, and the CRediT statement.

## 2026-07-24 — Journal, title, and word budget

Selected **Water Resources Research** as the primary target journal, matching
M5's positioning and rationale (AGU publication-unit allowance comfortably
fits a 6,500-word main text plus six main figures and one main table; WRR is
the correct comparison literature for an identifiability/equifinality and
evidence-integration benchmark). Locked the 6,500-word main-text budget
(Introduction 1,050 / Methods 1,550 / Results 2,100 / Discussion 1,550 /
Conclusion 250) per explicit user instruction, mirroring M5's revised
budget. Full rationale recorded in `Outline.md`.

## 2026-07-24 — Scope: registering already-locked results, not new computation

M7.3 (`results/m7_3_locked/`) and its M7.2 supporting-validation audit trail
were already executed and locked (commits `d336e87`, `771388a`, `8ca2036`)
before this manuscript task began. `manuscript/analysis_plan.json`,
`artifact_validation.json`, and `artifact_registry.json` therefore document
already-executed, already-validated artifacts retrospectively rather than
driving new B1/B2 generation. `generation_script.py` is a thin wrapper over
the existing `scripts/run_m7_3_nonuniqueness.py` and
`scripts/make_m7_3_publication_assets.py` entrypoints, added only to satisfy
the B1 gate contract; it performs no new analysis.

## 2026-07-24 — Methods word count under budget

Main Methods came in at 928 words against a 1,550-word outline target (gate
requires `count <= limit`, not a minimum). All seven experiment-design
elements, the primary fusion equation, and required citations were retained;
the shortfall was not from omitted content but from Supplementary Methods
absorbing every derivation, secondary validation, and implementation detail
per `methods-writing.md` allocation rules. No content was cut to fit; the
section was written to completeness and happened to need fewer words.

## 2026-07-24 — Supplementary figures/tables kept out of main assembly numbering

`assemble_manuscript.py --require-all-artifacts-referenced` was **not**
used, by design. Registering Figure S1 and Tables S1-S5 in the same
`artifact_registry.json` as the main items (for B3 traceability) but
embedding them via the shared `[[FIG:...]]`/`[[TAB:...]]` token system would
have caused the assembler to assign them sequential main-text numbers
(e.g. "Figure 6") instead of their intended S-prefixed supplementary
numbering. Instead:
- Main assembly (`Manuscript-Final.md`/`.docx`) references only
  FIG-1..FIG-5 and TAB-1..TAB-6 via tokens, giving clean Figure 1-5/Table
  1-6 numbering.
- `manuscript/supplementary/Supplementary-Figures-and-Tables.md` was
  authored directly (not via the token system) with explicit "Figure S1"/
  "Table S1"-"S5" labels, converted to its own DOCX.
- Results-section prose refers to the supplementary items by plain-text
  cross-reference ("Supplementary Table S1", etc.), not by token.

## 2026-07-24 — Bibliography trimmed to actually-cited sources

Removed three bibliography entries reused from the M5/M6 template
(`sreekanth2017monitoring`, `anku2008ghana`, `banoengyakubo2011ghana`) that
were not cited anywhere in the M7 manuscript, to avoid an uncited-reference
warning. Kept `vehtari2021rhat`, which is cited only in Supplementary
Methods (outside `citation_audit.py`'s `manuscript/sections/**` scan scope),
because it is genuinely used there.

## 2026-07-24 — Figure/table token placement fix

Initial drafts embedded `[[FIG:...]]`/`[[TAB:...]]` tokens inline inside
parenthetical prose (e.g. "the benchmark architecture ([[FIG:FIG-1]])"),
which produced broken sentence flow once hydrated (the full image+caption
block was substituted mid-parenthesis). Corrected by using
`[[FIGREF:...]]`/`[[TABREF:...]]` for inline textual mentions (renders as
plain "Figure N"/"Table N") and placing bare `[[FIG:...]]`/`[[TAB:...]]`
tokens as their own paragraph immediately after the introducing sentence,
consistent with the tool's first-mention/reference-elsewhere design.

## Canonical field-data correction (2026-07-26)

Following the same correction already applied to M5 and M6, M7's Experiment
4 (Ghana data-scope audit) and the M7.2 supporting-validation field
component depended on a retired antecedent-study workbook
(`data/NorthenGhana/Aquifers_Dataset_Mendeley.xlsx`) rather than this
project's own canonical raw field data
(`data/FieldData/NorthenGhana/NorthernGhana.xlsx`, Dry/Wet sheets). Reading
`data/FieldData/NorthenGhana/SI.pdf` (already established during the M5
correction) confirmed the Mendeley workbook is the supplementary
information for a separate publication, "Graph-inverted Ghanaian aquifers
under aridification" — a different study's own derived reprocessing of the
same 160 boreholes (aquifer/geology/land-use classification, synthetic
recharge/aridity/risk-score columns, a processed graph-edge sheet, and a
fabricated per-record `Sampling_Date` field), not independently observed
field data and not part of this project.

**Scope-changing finding.** Unlike M5/M6, this was not a pure relabeling
fix. `field_prequential.py` (the M7.2 within-campaign hold-forward
experiment, whose output CSVs M6's Table 5/Figure 5/Supplementary S11 also
reuse) depends on a per-record sampling date to define "20 sequential August
issue batches" for a leakage-audited prequential test. The canonical
workbook records exactly one dry-season and one wet-season observation per
well with **no intra-season date field at all** (SI.pdf documents only
campaign-level date ranges — 5-24 March 2025 dry, 5-24 August 2025 wet —
not a per-well date), so the literal experiment as designed cannot be
reconstructed from canonical data. Presented to the user as a three-way
choice (redesign with a disclosed arbitrary order / drop the experiment /
collapse to a single non-sequential comparison); the user chose to redesign
with a disclosed arbitrary order.

- **`field_prequential.py` rewritten**: reads the canonical `Dry`/`Wet`
  sheets directly (no more `Wells_Nodes`/`Hydrochemistry_Seasonal` split);
  wells are revealed in a fixed, seeded pseudo-random permutation split into
  20 batches (`batch_order_seed = 2025`) instead of real chronological
  order, with the audit output (`batch_assignment`, `batch_order_seed`)
  stating this explicitly rather than implying a real sampling sequence.
  `Data_Class` (retired workbook only) was replaced with an independently
  computed charge-balance-error quantitative screen (|CBE| <= 5%, same
  formula convention as M6). `Aquifer_Type`/`Geology_Group`/`Land_Use` and
  the four synthetic Mendeley risk-score columns
  (`Recharge_Zone_Score_0_1`, `Aridity_Stress_0_1`,
  `Fluoride_Geogenic_Risk_0_1`, `Anthropogenic_Risk_0_1`) were removed from
  the design matrix and the upstream-graph "same unit" preference, which now
  uses `Region` instead of the removed `Aquifer_Type` (no independent
  aquifer-type classification exists for these boreholes). Result: 140
  complete quantitative pairs (vs 138 previously — an artefact of the
  independent CBE screen replacing the retired `Data_Class` field, not a
  data change); graph ridge vs persistence paired MAE difference -0.173
  (95% CI [-0.188, -0.158]); graph ridge vs expanding-mean-delta -0.005 (95%
  CI [-0.012, +0.002]) — closely reproducing the previous qualitative
  finding (graph ridge clearly beats persistence, is not distinguishable
  from the simple seasonal baseline) under the new, honestly-disclosed
  batch design.
- **`audit_ghana_workbook()` (`m7_3_analysis.py`) rewritten** for the
  canonical Dry/Wet schema: no `Graph_Edges` or `Coordinate_Masking_Note`
  sheets to read (`processed_graph_edges_available` is now correctly
  `False`, `coordinates_masked` kept `True` per SI.pdf's own masking
  statement, independently of any workbook sheet), no `Sampling_Date`
  field (`sampling_date_field_available: False`, with a
  `sampling_granularity` string replacing the old
  `first_sampling_date`/`last_sampling_date` fields), and a new
  `independent_aquifer_type_classification_available: False` flag. The six
  boolean flags this audit contributes to `manifest.json`'s `"ghana_scope"`
  block are unchanged in value between the old and new implementations
  (verified before deciding not to re-run the expensive M7.3 MODFLOW
  benchmark — see below).
- **Data paths fixed**: `GHANA_WORKBOOK` in `run_m7_3_nonuniqueness.py` and
  the `run_field_prequential(...)` call in `run_supporting_validation.py`
  both now point at `data/FieldData/NorthenGhana/NorthernGhana.xlsx`.
- **Compute scope, chosen deliberately.** `run_supporting_validation.py
  --confirmatory` (the M7.2 confirmatory suite: 12 locked test seeds + 6
  development seeds of MODFLOW 6/MODPATH 7 simulation, topology MCMC, and
  10,000-replicate bootstraps) **was re-run in full** because its only
  Ghana-data touchpoint (`field_prequential`) sits inside that same
  orchestrator and its manifest embeds `field_prequential_claim`/
  `field_n_pairs`; a full re-run kept provenance simple. By contrast,
  `run_m7_3_nonuniqueness.py`'s only Ghana-data touchpoint
  (`audit_ghana_workbook`) writes to its own standalone
  `ghana_data_scope_audit.json` and contributes only the six unchanged
  boolean flags to `manifest.json`; this was verified by code inspection
  before choosing **not** to re-run the far more expensive M7.3 benchmark
  (600 age particles, 64 reaction-bootstrap replicates, 10,000 paired
  bootstrap draws per case, deterministic and therefore guaranteed
  bit-identical for everything except this one already-regenerated file).
  `ghana_data_scope_audit.json` in `results/m7_3_locked/` was regenerated
  directly instead.
- **Pre-existing bug found and fixed in passing**: `make_m7_3_publication_
  assets.py`'s `verify_outputs()` hardcoded "expect 11 publication CSV
  tables", stale since Table 7 and Supplementary Table S6 were added in the
  2026-07-24 revision round by scripts outside this file's own
  `build_tables()` (7 main + 6 supplementary = 13 files now share the
  output directory, of which this script itself produces 11). Confirmed via
  `git diff HEAD` that this script already carried other uncommitted,
  unrelated path fixes from before this session (the `FIELD_RESULTS` path
  had already been updated to the renamed `m7_nonuniqueness_benchmark`
  layout but never successfully run), so the stale count was pre-existing
  breakage, not something this session introduced. Updated the check to 13
  rather than deleting the two independently-produced tables.
- **Manuscript and docs propagated**: main Results (`03-results/section.md`)
  and the Data Availability statement
  (`06-availability-statements/section.md`) no longer claim processed graph
  edges or an external Mendeley citation are needed; `Supplementary-
  Methods.md`'s Ghana audit-criteria section, `docs/m7_3_protocol.md`,
  `docs/m7_3_results.md`, and `docs/supporting_validation_protocol.md`/
  `supporting_validation_results.md` were rewritten to match. `table6_
  ghana_claim_boundary`'s hand-typed "Processed graph edges: Available" row
  was corrected to "Absent", with new rows for intra-season sampling dates
  and independent aquifer-type classification.
- **Tests fixed**: `tests/test_m7_nonuniqueness.py` and `tests/test_m7_
  supporting_validation.py` referenced the retired workbook path and
  Mendeley sheet names directly; the leakage-audit test (proving later
  batches cannot influence earlier predictions) was rewritten to alter
  chemistry for wells in later *batches* rather than dates after a cutoff,
  exploiting the fact that batch assignment depends only on well eligibility
  (scale-invariant under the charge-balance screen), not chemistry
  magnitude, so eligibility — and therefore batch membership — is provably
  unaffected by the alteration. All 14 M7 tests and all M7 gates (methods
  draft, citation audit, abstract, assemble) pass.

## 2026-07-28 — Strict public-pipeline system acceptance

- Added and hash-locked `RUN-M7-SYSTEM-20260728-01` before generating six
  fresh independent MODFLOW 6/MODPATH 7 cases (seeds 6301--6306).
- The public pipeline completed every requested nuclear-age,
  sheaf-refinement, and network-fit stage; mean candidate recall was 0.9815
  and the generator retained `imports_hydrosheaf=false`.
- Full global-sheaf versus local-age PR-AUC difference: +0.0586, 95% paired
  bootstrap interval [+0.0386, +0.0777]. Brier difference: -0.00849
  [-0.00983, -0.00723]. This is evidence that the global section added value
  relative to the local age gate in this controlled setup.
- The preregistered overall incremental claim failed: full-sheaf PR-AUC was
  0.3075 versus 0.3272 for hydraulic-only, and the full-versus-permuted-age
  PR-AUC interval crossed zero. Selected-edge F1 was unchanged across all
  conditions because the fixed out-degree rule selected the same count.
- Decision: retain the system-execution and conditional global-section
  finding; prohibit an overall topology-superiority or field-validation claim.

## 2026-07-29 — M7.4 competence-matched sheaf-versus-graph benchmark

- Added a separate prospectively locked comparator rather than altering the
  M7.3 or system-acceptance protocols. The M7.4 lock fixes 32 development and
  64 held-out cases, four scenario strata, 10,000 paired case-block bootstrap
  replicates, equal optimisation/tuning budgets, and an identity-limit
  equivalence gate.
- Compared an edge-local weighted graph, an ordinary identity-coupled graph
  Laplacian, the production affine directed-section solver, and a within-case
  permuted-map adverse control. The independent generator imports no
  HydroSheaf code; truth does not enter inference features.
- The execution/equivalence gate passed. Graph and sheaf raw residuals and
  predictions were exactly equal in all 16 identity-limit test cases.
- Relative to the identity Laplacian, the affine sheaf improved PR-AUC by
  +0.0854 (95% CI [+0.0666, +0.1050]) and Brier score by -0.0193
  [-0.0235, -0.0152]. Native maps also beat the permuted-map control by
  +0.0909 PR-AUC [+0.0705, +0.1117].
- Relative to the stronger edge-local weighted graph, overall PR-AUC difference
  was +0.0097 [-0.0054, +0.0248] and Brier difference +0.00053
  [-0.00330, +0.00443]; the full incremental-superiority rule failed.
- In the planted incompatible-cycle stratum, the sheaf improved conflict-
  localisation PR-AUC over the edge-local graph by +0.0689
  [+0.0467, +0.0915] and over the identity Laplacian by +0.1098
  [+0.0917, +0.1278].
- Decision: allow a controlled-synthetic conditional claim that affine
  restriction maps represent non-identity relationships and help localise
  global inconsistencies. Prohibit claims of overall weighted-graph
  superiority, field validation, or temporal/3-D/vadose performance.
- Rebuilt the printable supplement from the authoritative publication CSVs.
  The former 80-page appendix printed 10--20-column per-case and per-edge
  tables that were unreadable in portrait Word/PDF output. The revised
  supplement prints compact claim-bearing views and points to every complete
  machine-readable CSV, reducing the document to 14 legible pages without
  deleting a row or column from the archived data products.

## 2026-07-29 — M7.5 robust/hybrid diagnostic protocol freeze

- Froze `M7_ROBUST_HYBRID_SHEAF_PROTOCOL.md` and
  `m7_robust_hybrid_protocol.lock.json` before implementing the new runner or
  generating development/test cases.
- Fixed 64 development cases (8401--8464), 128 locked-test cases
  (8501--8628), four unchanged M7.4 scenarios, three LOO iterations, fixed
  hybrid/calibration grids, eight seed-group folds, and 10,000 paired
  case-block bootstrap replicates.
- Prespecified separate tests of candidate self-influence, loss of local
  endpoint evidence, calibration mismatch, and a permuted-map adverse control.
- The locked-test run will occur once after a second lock freezes the runner
  and development-only selections. M7.4 will remain reported irrespective of
  whether M7.5 succeeds or fails.

## 2026-07-29 — M7.5 locked claim decision

- The development stage selected a local weight of 1.0 for both hybrids;
  global compatibility was therefore retained only as a fallback when a local
  endpoint residual was unavailable. The runner, settings, development
  manifest and protocol were then frozen in
  `m7_robust_hybrid_confirmatory.lock.json`.
- The single authorised 128-case test run passed all execution and provenance
  gates. Robust hybrid minus edge-local PR-AUC was +0.0200 (95% CI
  +0.0073 to +0.0324), but Brier was -0.00151 (-0.00419 to +0.00105) and
  log loss was +0.00333 (-0.00341 to +0.01009). The prespecified overall
  superiority gate failed.
- Restoring direct local evidence improved the original global estimator on
  PR-AUC, Brier score and log loss. Stronger leave-one-edge-out robustification
  instead worsened Brier and log loss, both as a global estimator and relative
  to the original hybrid.
- Conditional PR-AUC gains over edge-local were supported in incompatible-
  cycle (+0.0437) and noisy/missing (+0.0335) strata. The heterogeneous-affine
  log-loss penalty remained (+0.0251, CI excluding zero).
- Native maps beat the permuted-map control on PR-AUC, Brier and log loss.
  Decision: retain a bounded local-first/global-fallback contribution and
  map-semantic claim; prohibit general superiority and prohibit describing the
  LOO robust estimator as excellent.
