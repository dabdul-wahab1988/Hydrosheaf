# M3.1 decisions and audit trail

## Status: draft complete, pending human review and journal-specific formatting

`manuscript/Manuscript-M3.1.docx` (~5,150 prose words) and
`manuscript/Supplementary-Information-M3.1.docx` are drafted, cite-checked
(`citation_audit.py` PASS), and assembled with all figures/tables regenerated
in R from the current codebase. Known open items before submission:

- Equation numbering (S1, S2, ...) may not render visibly next to each
  equation in the DOCX; pandoc's `\tag{}` support for Word math (OMML) is
  unconfirmed and could not be visually verified (no PDF-rendering tool
  available in this environment — see below). Recommend a manual visual pass
  in Word before submission, adding equation-number table cells if missing,
  matching the convention already used in `M3/M3_geochemistry/Manucript_3.docx`.
- Visual QA was not completed to the same standard as the M3
  `M3_FINAL_ACCURACY_AUDIT.md` (page-by-page PNG inspection): this environment
  had no `pdftoppm`/poppler or ImageMagick installed, and installing new
  system software was judged out of scope for an unprompted QA convenience.
  The DOCX conversion itself reported PASS and a markdown round-trip
  (`pandoc docx -> markdown`) confirmed equations, citations, and structure
  survived conversion; page-level visual inspection is still recommended
  before submission.
- Author affiliation for Ebenezer Aquisman Asare vs. Samuel Ganyaglo is a
  guess, not a confirmed fact (D5) — confirm with the authors before
  submission.
- SHA-256/manifest provenance lock (mirroring `M3_FINAL_ACCURACY_AUDIT.md`)
  has not been produced for M3.1; recommend generating one once the text is
  final and no further regeneration is expected.


M3.1 upgrades `M3` (accuracy-locked 2026-07-28, commit `47fbf7c`, submission
package in `M3/M3_geochemistry/`). M3 is retained unchanged; it is not deleted
or overwritten. M3.1 supersedes it as the version to submit once regeneration
and review are complete.

## D1. Reason for the upgrade

The user's premise was that codebase changes since the M3 lock might have
changed what the manuscript should report. Verified by reading every post-lock
document and diffing every changed file under `hydrosheaf/nuclear`,
`hydrosheaf/graph`, and `hydrosheaf/validation` against the lock commit
(`47fbf7ca385b3210547673519ac2b30f910a80d9`).

Finding: every post-lock document that discusses the design-matrix numbers
(`m3_correctness_rerun_full_qa.md`, the identified-TTD QA/audit docs, the CFC
reconciliation and redox-exclusion notes) states explicitly that it does not
change a currently printed manuscript number. One code-level change was found
that could, in principle, have mattered: `calculate_tracer_reliability_weights`
in `hydrosheaf/nuclear/multi_tracer.py` had its `reference_age` parameter
removed (commit `2e73d51`) so that reported/reference ages can no longer reach
tracer-weighting — a leakage fix. Whether the locked run ever exercised that
parameter with a non-NaN value (as opposed to always calling it with the
default) determines whether this changes any locked number. **This is settled
empirically, not by reading the diff**: the full manuscript-analysis pipeline
was rerun against current `HEAD` in an isolated copy
(`M3.1/m3_age_benchmark/`, run started from repo root, log at
`run_m3_1_full_rerun.log`) and every reportable statistic is diffed against the
M3 lock's manifests. Result recorded in D2 once the run completes.

Independently of whether any number moved, two genuinely new results exist
that M3 (locked 2026-07-28) could not have reported:

- an independent set-valued ("identified-TTD") compatibility audit of the
  graph-prior question, reaching the same negative conclusion as the main
  benchmark by a different method (`d287221`);
- a 27.85% local-infeasibility rate in CFC tracer reconciliation, with seven
  candidate explanations tested and rejected (`57b0b00`, `8a4a354`, `d287221`),
  and a synthetic control on an independent MODFLOW/MODPATH generator (M7.6)
  that did not support the leading hypothesis either.

Both are explicitly labelled `development`/`implementation_only` in
`configs/identified_ttd_protocol.yaml` and in their own QA docs, i.e. exploratory
diagnostics, not confirmatory validation. See `Outline.md` for how they are
scoped into M3.1 (compact main-text mention, full detail in Supplementary,
never described as validated or as explaining the infeasibility cause).

## D2. Numerical reproduction — result

The full locking command (`run_m3_manuscript_analysis.py --full --age-steps 90`)
was rerun against current `HEAD` in `M3.1/m3_age_benchmark/` and every
regenerated result file was diffed line-by-line against the M3 lock
(`M3/m3_age_benchmark/results/`). Verdict: **the manuscript's reportable
numbers are unchanged.** Specifically:

- `m3_design_matrix_summary.csv` (Table 2 / abstract numbers): identical for
  every headline scenario (`tracerlpm_strict_parity`, N=329;
  `tracerlpm_parity_agefractions`, N=289; `hydrosheaf_selection_corrected`,
  N=309; and all others). Three non-headline scenarios differ only in the
  15th significant digit (e.g. `-3.1958677627964587` vs.
  `-3.195867762796456`) — ordinary floating-point noise, not a result change.
- All five `m3_cv_benchmark_{3H,SF6,14C,CFC11,CFC12}.csv` files (the
  leakage-guarded tracer-withholding test, Figure 6/Table in the manuscript):
  **byte-for-byte identical**, zero diff lines.
- `m3_real_usgs_graph_benchmark.csv` (Table 4, the graph-regularisation
  test): the one row that matters for the paper's claim —
  `hydraulic_proxy_constrained`, the only family with
  `improved_vs_single = True` — is byte-for-byte identical. Other families
  (`coordinate_global`, `depth_constrained`, `parameter_smooth_aicc`,
  `wrong_direction_negative_control`, `randomized_negative_control`) show
  small numeric drift (single-digit-percent in `rmse_graph_log10`,
  `n_violations_before/after`) but **no row changes its
  `improved_vs_single` classification**. The likely cause is that the design
  matrix runs 16 parallel worker processes (`m3_design_matrix_results.csv`
  row order is not guaranteed identical run-to-run), and the eight-iteration
  graph-regularisation update (manuscript Eq. 5) compounds tiny
  floating-point ordering differences across a connected graph; it does not
  reflect a scientific change and is present as an inherent property of the
  existing pipeline, not something newly introduced. This is disclosed for
  transparency; it does not alter any number cited in manuscript prose,
  which uses only the unchanged `hydraulic_proxy_constrained` result and the
  unchanged design-matrix/cross-validation numbers.

**Conclusion for the upgrade rationale:** contrary to the premise that motivated
this upgrade, none of the post-lock codebase changes altered a number M3
reports. The leakage fix to `calculate_tracer_reliability_weights`
(`reference_age` removed) turned out not to change the strict-parity,
age-fraction, or model-selection scenarios at all. M3.1 is therefore not a
correction of stale numbers; it is (a) an empirically verified confirmation
that the locked numbers still hold, (b) the two new exploratory findings
(D1), (c) a reader-friendly rewrite at a 6,000-word budget, and (d) new R-made
figures and a study-area map, none of which existed in M3.

## D2b. Identified-TTD diagnostic and network-dating demo: not rerun, and why

`M3.1/m3_age_benchmark/results/m3_infeasibility_diagnostics.json` (feeds
FIG-7) and `m3_network_dating_demo.csv`/`_summary.csv` (feeds Supplementary
Figure S4) were carried forward from the copied M3 tree rather than rerun.
Checked before relying on them: no commit after each file's generating
commit (`d287221` for the infeasibility diagnostic; predates the M3 lock
entirely for the network-dating demo) touches any file the respective script
imports (`hydrosheaf/nuclear/ttd_identified.py`, `ttd_design.py`,
`ttd_graph.py`, `run_m3_infeasibility_diagnostics.py`; `hydrosheaf/nuclear/lpm.py`,
`input_history.py`, `run_m3_network_dating_demo.py`), and the working tree is
clean. A rerun is therefore expected to reproduce byte-identically (as most
of D2's rerun did) and was not spent on given the time cost of the LP-based
infeasibility sweep; this is disclosed rather than silently assumed.

## D3. Out-of-scope material confirmed by audit

The following post-lock work touches no file M3 imports and reports no M3
number; it is excluded from M3.1:

- M7.6 (`e247318`, `b06523e`): self-contained M7 nonuniqueness benchmark on an
  independent synthetic generator. Reuses an M3-style diagnostic as a method
  only; explicitly forbids revising any M3 number.
- `hydrosheaf/validation/reaction_rapm.py`, `synthetic_claims.py`, and the
  wider Aug 1-2 "RAPM reaction layer" / "synthetic validation gates" additions
  (`8718d66`, `2d4b8af`): reaction/kinetics chemistry evidence models, tested
  only on independent synthetic generators. Never touch the M3 USGS age-dating
  dataset.
- M8 Bayesian active learning (`53beb46`): measurement-design scope, no M3
  file touched.

## D4. Computation and figure authority

Following the M2.3 precedent in this repository: the existing Python benchmark
pipeline under `m3_age_benchmark/` remains the computational authority (it is
proven, scenario-matrix-based, and already produces SHA-256-manifested
results). It is copied to `M3.1/m3_age_benchmark/` and rerun against current
`HEAD` rather than rewritten. R is added as the figure/map authority
(`M3.1/analysis/r/`), consuming CSV exports only and recomputing no reported
statistic, matching `M2.3/analysis/r/` conventions (`_theme.R`, `_map.R`,
one script per figure).

## D5. Author information source and an unresolved affiliation discrepancy

The user pointed to `AuthorsInformation.txt` (repository root) as the source
for author metadata. It confirms the same five authors and ORCID as the
current M3 docx, but reassigns affiliations: it lists **both** "Nuclear
Chemistry and Environmental Research Centre, NNRI, GAEC" and "Water Resources
Research Centre, NNRI, GAEC" under **Samuel Ganyaglo**, and gives no
institutional affiliation at all for **Ebenezer Aquisman Asare** (email only).
The current locked M3 docx instead splits these: Asare carries the Nuclear
Chemistry and Environmental Research Centre affiliation, and Ganyaglo carries
University of Ghana plus the Water Resources Research Centre affiliation.

Resolution adopted for drafting (not verified against the authors
themselves): keep `AuthorsInformation.txt` verbatim for contact
details/ORCID/Adomako's dual email, and — since removing Asare's institutional
affiliation entirely would look like an error rather than a deliberate change —
retain his existing Nuclear Chemistry and Environmental Research Centre
affiliation from the current M3 docx alongside Ganyaglo's now-dual
affiliation, i.e. treat the two as sharing that affiliation rather than it
transferring. **This is a guess, not a confirmed fact, and is flagged to the
user for correction before submission.**

## D6. Map data source

M3's benchmark uses the USGS national public-supply aquifer release (sites
across the continental United States), not the Ghanaian field datasets used
elsewhere in this research programme. The new site map (FIG-1) therefore uses
US state boundaries via `geobounds::gb_get_adm1("USA")` (same package/cache
convention as `M2.3/analysis/r/fetch_boundaries.R`), not the Ghana boundary
GeoPackage already cached under `M2.3/analysis/r/data/`.
