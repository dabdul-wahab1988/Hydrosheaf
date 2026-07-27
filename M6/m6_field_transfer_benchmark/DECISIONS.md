# M6 manuscript — material decisions

## Revision round (manuscript-reviewer report, 2026-07-24)

A full peer-review pass was run using the manuscript-reviewer skill
(`manuscript/review/M6_manuscript_review.md` / `.docx`), scoring 59/100
(major revision required, no critical-flaw override triggered). The
following fixes were applied directly to the manuscript in response:

- Abstract rewritten to disclose Talensi's charge-balance screening failure
  and to attribute the tier-3-to-2 identifiability collapse to the
  evidence-gating prior rather than implying "fit quality" alone reveals it.
- Introduction: hedged the "no existing study" gap claim with "to the
  authors' knowledge" and anchored it to the cited regional studies;
  disclosed that the Predictions in Ungauged Basins citations are drawn from
  streamflow literature applied by analogy; fixed the five-objectives-vs-six-
  experiments inconsistency by adding a sixth objective; added a
  locked-analysis-plan disclosure sentence.
- Methods: added a self-contained summary of how the frozen M5 calibration
  was originally fitted and validated; disclosed that the M5/M7.2 components
  this study depends on are not yet independently citable and that a data
  and code availability statement will accompany submission; disclosed that
  bootstrap resample counts vary by experiment; disclosed that evidence-gate
  thresholds were not sensitivity-tested (see new Discussion limitation).
- Results: added an explicit note that Table 3 (seasonal-pair residuals) and
  Table 4 (edge residuals) report non-comparable "Tier 1" Northern Ghana
  constructions, mirrored in both table captions; added a chance-baseline
  caveat to the prior-label concordance claim; disentangled the Talensi
  charge-balance/tracer-scarcity conflation and noted Talensi's rate is
  *lower* than Northern Ghana's own seasonal-pair Tier 1 rate; reordered the
  hold-forward paragraph to lead with the informative (inconclusive)
  graph-vs-mean-delta comparison rather than the graph-vs-persistence
  comparison; softened the transport-correction "not sensitive to choice"
  claim (63.8% agreement means over a third of wells disagree).
- Discussion: added a new subsection disclosing the overlap between this
  study's Northern Ghana seasonal-transfer experiment and the antecedent
  project's own prior Northern Ghana field application; reframed the
  strontium/silica "next-best-measurement" finding as a quantified
  demonstration of established geochemical domain knowledge rather than a
  new discovery; removed the Sreekanth et al. (2017) citation (verified via
  web search to concern spatial monitoring-well placement, not tracer-panel
  selection, and therefore misapplied) and replaced the argument with the
  authors' own reasoning; expanded "Three limits" to five, adding the
  gate-threshold-sensitivity gap and the sample-size-dependent kNN-graph
  confound between external and reference datasets.
- Study design: corrected the Chegbeleh et al. (2020) and Zakaria et al.
  (2021) citations, verified via web search to be overstated in the original
  drafting (Chegbeleh's own emphasis is weathering plus agrochemical/
  wastewater nitrate and lead, not mining pressure specifically; Zakaria
  reports 94% good-quality water with fluoride as one of several minor
  contributing factors, not a defining geogenic characteristic).
- Editorial: consolidated the honesty-boundary statement (previously
  repeated near-verbatim in four sections) to a single full statement in
  Study design with brief cross-references elsewhere; standardised
  "Mechanism Resolution Score" to "MRS" after first definition in Results.
- Bibliography: removed the now-unused `sreekanth2017monitoring` entry;
  updated `claim_source_map.csv` and `methods_evidence.csv` (METH-24 moved
  to `supplementary` destination) to match; reassembled with
  `--max-table-rows 20` to fix a silent truncation bug that had dropped 4 of
  Table 2's 16 variable rows (Sr, SiO2, Calcite_SI, Aquifer_Type) from the
  first-drafted DOCX.

All r2m gates (methods plan/draft/docx, citation audit, abstract gate) pass
on the revised manuscript. Main-body word count after revision: 6,126 words
(within the 6,500-word ceiling). Not fixed in this pass, and flagged
explicitly in the manuscript itself as outstanding: an actual data/code
availability statement with real repository locations, a citable reference
for M5/M7.2, and a numeric threshold-sensitivity analysis for the
evidence-lift gate all require actions (repository creation, rerunning the
pipeline) outside the scope of a text-only revision pass.


Decisions taken while drafting the Q1 manuscript package under `manuscript/`
(Track B, r2m skill). Dated 2026-07-24.

## Word budget

- User set an explicit 6500-word ceiling for the main manuscript (Introduction
  through Conclusion, prose only), superseding the Outline.md's original
  ~5,000-word Communications Earth & Environment target. The section budget
  was rescaled proportionally from the outline's per-section split.
- Delivered main-body prose totals ~5,000 words (Introduction 807, Study
  design 642, Methods 742, Results 1,676, Discussion 888, Conclusion 249),
  leaving headroom under the 6,500-word ceiling rather than padding to it.
- Main Methods is capped separately at 1,300 words per
  `manuscript/methods/methods_manifest.json` (delivered: 742 words); all
  derivations and secondary validation moved to Supplementary Methods
  (~4,800-word target, 15 subsections).

## Scope of the Methods evidence ledger

- `manuscript/methods/methods_evidence.csv` was restricted to statements
  describing the Methods workflow itself (destinations `main`/`supplementary`
  map to main-Methods/Supplementary-Methods only, per the methods-writing.md
  contract). Citations that support Introduction or Discussion prose only
  (e.g. `macdonald2012africa`, `plummer1994netpath`, `plummer1980massbalance`,
  `podgorski2022fluoride`, `onipe2020fluoride`, `tredoux2006nitrate`,
  `gibbs1970mechanisms`, `piper1944graphic`, `sreekanth2017monitoring`) were
  removed from that ledger and tracked instead via
  `manuscript/citations/claim_source_map.csv`, which covers the whole
  manuscript body. This was necessary because the Methods draft gate checks
  that every ledger citation actually occurs inside the declared
  main-Methods/Supplementary-Methods text, and these citations legitimately
  belong to other sections.

## Cross-milestone provenance anchoring

- M6 genuinely reuses frozen M5 artefacts (inverse solver, MRS calibration)
  and an unmodified M7.2 supporting-validation benchmark (field prequential
  hold-forward). The r2m "keep paths local" rule requires evidence-ledger
  `project_sources` to stay contained under the M6 project root
  (`M6/m6_field_transfer_benchmark/`), so cross-milestone files (M5/M7) are
  cited in evidence-ledger *statements* (prose) but the ledger's *source
  anchors* point to the local M6 script that performs the import/reuse
  (e.g. `scripts/m6_common.py`'s import of the frozen M5 calibration;
  `scripts/make_objective6_prequential_figure.py`'s read of the M7.2 CSVs),
  not to the upstream files themselves.

## Figure/table numbering

- The assembly tool (`assemble_manuscript.py`) numbers figures and tables by
  order of first mention in the assembled document, not by artifact ID. To
  preserve the package's existing figure numbering (matching the R script
  output filenames `figure1_..figure6_...` and all prior docs/*.md
  descriptions), a forward reference to Figure 1 (`[[FIGREF:FIG-1]]`) was
  added in "Study design and datasets" so it is mentioned before Figure 2
  (the Methods workflow schematic), even though Figure 1 is embedded later,
  in Results.
- Extended Data Figures 1-3 are embedded as plain Markdown images with
  hand-written "Extended Data Fig. N" captions rather than through the
  `[[FIG:...]]` token system, because that system uses one shared sequential
  counter per artifact kind and would otherwise renumber them as Figure 7-9
  instead of preserving the Nature-portfolio Extended Data numbering
  convention.
- Methods no longer embeds Figure 1 (dataset/tier design), Figure 5 (field
  prequential) or Table 2 (variable availability): those are Results-owned
  artifacts per Outline.md's own display-item plan (§4.1), and embedding them
  in Methods both violated that plan and produced out-of-order numbering.

## Title and proposal metadata

- No `proposal.normalized.json` existed for M6. `assemble_manuscript.py`
  falls back to the first `#`-heading in `Outline.md` ("M6 — Q1
  Nature-Portfolio Manuscript Outline") for the document title if this file
  is absent, which is not a manuscript title. A minimal
  `proposal.normalized.json` was created at the M6 project root using the
  Outline's own stated working title and RQ1-RQ6 research questions as
  objectives (no invented content) so the assembled manuscript carries the
  correct title.

## Human review attestation

- `manuscript/methods/human_review.json` records a self-review performed by
  the drafting agent against all eight required checks, since no separate
  independent human reviewer was available in this session. This is recorded
  honestly as an agent self-review, not represented as independent human
  sign-off, and the equation-accuracy check flags that visual inspection of
  the rendered OMML equations in the DOCX remains an outstanding manual step
  for the user.

## Supplementary Methods DOCX rendering

- `Supplementary-Methods.md` legitimately keeps two `[[TAB:TAB-S7]]`/
  `[[TABREF:TAB-S8]]`-style tokens (required by the Methods draft gate's
  artifact-reference check), but that file is never hydrated by
  `assemble_manuscript.py` (which only processes `manuscript/sections/`), so
  converting it to DOCX directly would leave literal bracket text in the
  reader-facing document. A token-free copy was generated on the fly
  (tokens replaced with plain "Table S7"/"Table S8" text), converted to
  `Supplementary-Methods.docx` with Pandoc, and the temporary copy deleted.
  The canonical `Supplementary-Methods.md` (with tokens intact) remains the
  gate-checked source of truth.

## Supplementary Information docx

- `manuscript/supplementary/Supplementary-Figures-Tables.docx` was built
  outside the r2m token/registry pipeline (which only assembles
  `manuscript/sections/`) by programmatically concatenating each
  supplementary figure's registry caption/technical description/
  interpretation with its image, and each supplementary table's locked,
  already-generated Markdown (`tables/tableSx_*.md`) verbatim, then
  converting with Pandoc. Table content was read directly from the locked
  table files rather than retyped, to avoid transcription risk.

## Canonical field-data correction (2026-07-26)

The user reported that M6 (and M5/M7) had been analysing field data from a
retired antecedent-study workbook (`Aquifers_Dataset_Mendeley.xlsx`, an SI
export of a *different* study's own derived/inferred products — aquifer-type
labels, dominant-process priors, charge-balance-error values and a supplied
graph-edge sheet) rather than the project's own canonical raw field data at
`data/FieldData/` (`LowerAnayari/manu.csv`, `Talensi_MiningArea/talensi.csv`,
`NorthenGhana/NorthernGhana.xlsx`). This was corrected end-to-end:

- **Data source.** `m6_common.py` now loads exclusively from
  `data/FieldData/`. `load_northern_ghana()` was rebuilt from the raw
  `Dry`/`Wet` sheets only (no `Graph_Edges` sheet, no supplied `Aquifer_Type`,
  `Dominant_Process_prior`, `Data_Class` or `CBE_provided` fields). Talensi
  and Lower Anayari loaders had the same borrowed fields removed.
- **Saturation indices.** Rather than leaving Tier 4 unconstrained, real
  PHREEQC saturation indices (Calcite, Dolomite, Gypsum, Halite, Fluorite)
  are now computed for every Northern Ghana dry-season sample using the
  project's own `hydrosheaf.phreeqc.runner.run_phreeqc()` +
  `hydrosheaf.config.Config` interface (`phreeqpython`, in-process — no
  subprocess/executable dependency), per the user's request to reuse
  Hydrosheaf's existing, more efficient PHREEQC integration instead of a
  hand-rolled one.
- **Region replaces aquifer.** No independent aquifer-type classification
  exists in the canonical raw data for any of the three datasets. Every
  aquifer-stratified code path (experiments, figures, tables) was
  re-stratified by administrative region/district instead; hardcoded
  Aquifer_Type string labels for Talensi/Lower Anayari were kept only where
  they describe the dataset as a whole (mining area / regolith-alluvial),
  never as a per-sample classification.
- **Provided-graph edge set removed.** The RQ4 edge-uncertainty experiment
  compared a fourth, imported "provided graph" edge set against three
  self-generated ones (chemistry-kNN, geographic-nearest, random/perturbed).
  Per the user's explicit choice, the provided graph was dropped entirely
  rather than substituted; the experiment now compares only the three
  Hydrosheaf-generated sets, with chemistry-kNN as the reference for
  divergence/stability statistics.
- **Prior-label concordance removed.** All comparisons against the
  antecedent study's own process/aquifer labels (`prior_concordant`,
  `frac_concordant`, `_family_matches_prior`) were removed as a category, since
  those labels are that other study's inferences, not independent field
  truth, and treating them as a benchmark target would have been circular.
- **Numbers changed.** All six M6 experiment CSVs were regenerated from the
  canonical data and real SI values; results, figures and tables were
  rebuilt from these new numbers (e.g. the tier-3-to-2 identifiability
  collapse, region-stratified MRS/stability, and external-transfer rates all
  shifted from their previous Mendeley-sourced values). A pre-existing,
  unrelated bug in `experiment3_tier_ablation` (`_seasonal_wells` was called
  but had been deleted from the working tree, leaving a dangling name)
  surfaced during the rerun and was fixed by calling the still-present
  `m6.seasonal_well_pairs(ng)` instead; confirmed via `git diff HEAD` that
  this was pre-existing breakage, not something introduced this session.
- **Manuscript propagated.** Every main section (00-Introduction through
  06-Conclusion), `Supplementary-Methods.md`, `Supplementary-Figures-Tables.md`,
  `artifact_registry.json`, `methods_evidence.csv`, `methods_manifest.json`
  (S14 heading rename to "Removal of borrowed antecedent-study fields") and
  `claim_source_map.csv` (CLAIM-08 claim_text) were updated to match the new
  numbers and the region-based (not aquifer-based) framing.
- **Gates.** `methods_gate_check.py --stage draft`, `citation_audit.py
  --strict-bib --write-references-bib` and `abstract_gate_check.py` all PASS.
  `assemble_manuscript.py` was run **without**
  `--require-all-artifacts-referenced`: that flag cannot pass for this
  manuscript by design (see "Supplementary Methods DOCX rendering" above) —
  Extended Data and Supplementary figures/tables are deliberately numbered
  and cited outside the `manuscript/sections/` token pipeline, so 20 of the
  registry's artifacts (`FIG-ED1-3`, `FIG-S1-11`, `TAB-S2-9`) can never
  register as "used" by that scan; this is unrelated to, and predates, this
  session's corrections. The plain `assemble_manuscript.py --project-root .`
  call passed and regenerated `Manuscript-Final.md`.
- **Known remaining gap (closed 2026-07-26).** M6's Table 5/Figure 5/
  Supplementary S11 (the within-campaign hold-forward benchmark) is reused
  unmodified from M7.2. M7 was subsequently fixed the same day: the field
  data source was corrected and, because the canonical data carry no
  intra-season sampling date, the "20 sequential August issue batches"
  design was redesigned as a fixed, disclosed arbitrary batch order (user
  decision, M7 DECISIONS.md). `make_objective6_prequential_figure.py` and
  `make_m6_tables.py` were re-run to pick up the refreshed
  `results/supporting_validation/field_prequential_*` outputs (138 -> 140
  complete quantitative pairs; MAE persistence 0.342, mean-delta 0.174,
  graph-ridge 0.169; paired differences -0.173 [-0.188,-0.158] vs
  persistence and -0.005 [-0.012,+0.002] vs mean-delta; coverage90 0.631 /
  0.903 / 0.923) — closely reproducing the previous qualitative finding.
  `04-results/section.md`, `artifact_registry.json` (TAB-5), 
  `methods_evidence.csv` (METH-24), `docs/objective6_data_limited_synthesis.md`,
  and `Supplementary-Methods.md` S11 (which also had a pre-existing,
  unrelated March/August wet/dry season mislabelling, fixed in the same
  pass) were updated with the new numbers and the disclosed-arbitrary-order
  framing. M2 has the identical broken field-data path but remains out of
  scope per the user's module-by-module sequencing.
- **Evidence-ledger line-number audit.** `manuscript/methods/methods_evidence.csv`
  and `manuscript/methods/human_review.json` cite specific line ranges in
  `scripts/m6_common.py` and `scripts/run_m6_field_transfer.py` as source
  anchors; the canonical-data rewrite shifted most of these (e.g. the new
  `Config`/`run_phreeqc` import and `_compute_northern_ghana_si` function
  pushed everything below line ~90 down by 30-60 lines). Every anchor was
  re-verified against the current file and corrected. One CSV-quoting bug
  was introduced and fixed in the same pass: an edited `source_anchors` field
  (METH-04) gained an internal comma without surrounding quotes, which split
  the row into 9 columns and made `citation_placement()` misread a line-range
  fragment as a citation key, failing `methods_gate_check.py --stage draft`;
  wrapping the field in quotes restored 8-column parsing and the gate passed
  again. All four M6 gates (methods draft, citation audit, abstract,
  assemble) and `pytest tests/test_m6_q1_package.py` (3/3) were re-verified
  green after this fix.
## Saturation indices computed for Talensi and Lower Anayari (2026-07-27)

The user flagged a suspected overclaim: Supplementary Methods S3 stated that
Talensi and Lower Anayari's native tier ceilings were "fixed by data
availability rather than by choice." Checked directly: both datasets have
complete pH, temperature and major-ion records (63/63 and 41/41
respectively), the only inputs PHREEQC saturation-index computation needs —
Sr/SiO2 are not required. The tier ladder's own nested construction (Tier 4
= Tier 3 + saturation indices, and Tier 3 requires Sr/SiO2) was gating SI
out for both external datasets by design, not by a genuine data ceiling;
the manuscript's "data availability" framing was therefore inaccurate for
the SI component specifically (it remains accurate for Sr/SiO2 and season).

Given the user's explicit choice to fix this substantively rather than by
rewording alone:

- **`m6_common.py`**: `_compute_northern_ghana_si` generalised to `_compute_si`
  and reused in `load_talensi()`/`load_manu()`. Found and fixed a real bug in
  the same pass: the function passed literal `NaN` (not `None`) for Talensi's
  unmeasured fluoride column straight into the PHREEQC composition builder,
  which silently poisoned every saturation index for every Talensi sample
  (`phreeqc_ok=True` but all SI values NaN). `_add_value`'s `_parse_float`
  accepts `float('nan')` as a valid number, so nothing downstream caught
  this. Fixed by converting non-finite values to `None` before building the
  sample dict. After the fix: Talensi 61/63 samples got real SI (fluorite
  correctly NaN for all, since F is genuinely unmeasured there); Lower
  Anayari improved from 40/41 to 41/41.
- **`run_m6_field_transfer.py`**: `experiment5_external()` now passes
  `si=m6._upstream_si(tgt)` (previously hardcoded `None`) and includes
  `"si"` in both datasets' lifters list; the matched-tier Northern Ghana
  reference loop was given the same `si`/`"si"` treatment for a fair
  comparison. This is not only an evidence-lift-gate change: `fit_unit`
  selects thermodynamic-bounded fitting (`method="thermo"`) whenever `si`
  is non-empty, so the underlying inverse fit itself, not just the
  reported class, now differs for these edges.
- **Numeric effect, direction included.** Talensi (Tier 1): 36.0% -> 37.2%
  non-identifiable (small). Lower Anayari (Tier 2): 96.5% -> 95.3% (small).
  Matched-tier Northern Ghana reference: Tier 1 35% -> 54.2%, Tier 2 35% ->
  53.3% (large) — the NG reference edges are more often carbonate-dominant
  than Talensi/Manu's, so gaining SI corroboration mattered far more for
  the reference than for the external sets themselves. Net effect: Talensi
  now reads *below* its matched NG reference (37.2% vs 54.2%), reversing
  the previous "essentially the same rate" reading — read as evidence that
  per-dataset process composition, not evidence tier alone, drives the
  non-identifiable rate. Lower Anayari's own number barely moved because
  its edges are overwhelmingly silicate-dominant, a family SI does not
  corroborate (only Sr/SiO2 does); the "Sr/SiO2 is the next-best
  measurement for Lower Anayari" conclusion is unchanged.
- **E1 readiness, E6 limitation map** picked up the change automatically:
  `Calcite_SI` now shows "present" for Talensi/Lower Anayari in the
  variable-availability table (was "absent"), and carbonate weathering is
  now partially identifiable in all three datasets in the limitation map
  (was "non-identifiable in both external datasets"); silicate remains
  non-identifiable in both, unchanged, since SI does not corroborate it.
- **Propagated** to Study design, Methods, Supplementary Methods S3
  (rewrote the "fixed by data availability" sentence to state precisely
  that SI availability does not follow the ion-panel tier ceiling), Results
  E5 and E6 paragraphs, `artifact_registry.json` (TAB-4, FIG-ED1),
  `docs/m6_results_summary.md`, `docs/objective6_data_limited_synthesis.md`.
  All figures and tables regenerated from the fresh CSVs. All four gates
  and `pytest tests/test_m6_q1_package.py` (3/3) re-verified green.
- **M5 has the identical issue** (external field-transfer ELRI audit passes
  `upstream_si={}` unconditionally for Talensi/Lower Anayari, and does not
  even request pH/Temperature from their CSVs) and is being fixed in the
  same session; see `M5/m5_inverse_reaction_benchmark/DECISIONS.md`.

- **Internal docs (`docs/`, `Outline.md`, `README.md`,
  `proposal.normalized.json`) refreshed.** `docs/m6_data_readiness_audit.md`,
  `docs/02_figures.md`, `docs/03_tables.md`, `docs/objective6_data_limited_synthesis.md`
  and `docs/m6_robustness_diagnostics.md` had their aquifer/Mendeley/provided-graph
  language and (where affected) numbers corrected directly — the transport-correction
  agreement figure in `m6_robustness_diagnostics.md` moved from a stale 63.8%
  (66% vs 73% silicate share) to the current 64.4% (71% vs 70%), since regrouping
  the recharge endmember from Aquifer_Type to Region changed that specific
  computation; the gate-on/off, synthetic-validation and Talensi-CBE numbers in
  that same document were unaffected by the fix and left as-is after verifying them
  against the current result CSVs. `docs/m6_locked_analysis_plan.md` and
  `Outline.md` are historical Stage-B1 planning documents pre-dating the
  correction; rather than silently rewriting locked prose, each was given a
  dated amendment note at the top stating what changed and why, with the
  original as-drafted wording left intact below for audit purposes (`Outline.md`'s
  "Data & Analysis Status" results section, which states current facts rather
  than a plan, was updated directly with current numbers). `README.md` and
  `proposal.normalized.json` are user/tooling-facing rather than a locked
  historical record and were updated directly (dataset table, honesty-boundary
  statement, RQ2/RQ4 wording, `DATA-1` path and provenance, `RISK-1`).
