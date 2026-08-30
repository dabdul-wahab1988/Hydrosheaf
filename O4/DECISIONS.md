# O4 decisions and claim ledger

O4 is a new manuscript package for Objective 4 of Chapter 1 (HydroSheaf PhD):
"To assess the uncertainty and robustness of the framework under different
assumptions and data limitations" (`draftbychapter/PhDChapter1_refined.docx`,
Section 1.4.1). It harmonizes, but does not restate the software architecture
of, `M6` (field-transfer/robustness under data scarcity), `M7` (conditional
integration and non-identifiability, M7.3-M7.6) and `M8` (calibration is not
identification: parameter recovery, uncertainty calibration, optimal
measurement design). It is written after, and positioned as a sibling to,
`O3` (Objective 3, M3+M4+M5 benchmark harmonization) and as a companion to
`M2.3` (the HydroSheaf framework paper), all three targeting Computers &
Geosciences.

## Locked decisions

### D1. Package status and journal target

O4 is a new, standalone manuscript targeting Computers & Geosciences,
distinct from `M2.3` and `O3`. It replaces no prior submission plan for `M6`,
`M7` or `M8`; those three packages and their own journal-specific Word/
markdown deliverables (`M6/.../manuscript/Manuscript-Final-Revised.md`,
`M7/.../manuscript/Manuscript-Final-Revised-M7.6-20260731.docx`,
`M8/.../manuscript/Manuscript-Final.md`) are retained unchanged as historical
records of the underlying benchmarks and are not superseded as evidence
sources. This mirrors `O3/DECISIONS.md` D1.

### D2. No re-simulation

O4 performs no PHREEQC, MODFLOW/MODPATH, PEST-GLM, or Bayesian
active-learning re-run. Every number reported in O4 is read from an
already-locked result file under `M6/m6_field_transfer_benchmark/results/`,
`M7/m7_nonuniqueness_benchmark/results/` (specifically `m7_3_locked/`,
`RUN-M7-SHEAF-VS-GRAPH-20260729-01/`, `RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01/`,
`RUN-M7-6-M3-MECHANISM-20260731-01/`), or `M8/m8_calibration_benchmark/results/`
and `M8/m8_calibration_benchmark/manuscript/Manuscript-Final.md`'s own
registered tables, or is a harmonization-layer statistic (a ratio, a paired
contrast already computed by the source component, a re-grouping) computed by
`O4/analysis/python/` directly from those files with no change to the
underlying inference. This mirrors `O3` D2, itself mirroring `M2.3` D5.

### D3. Staleness verification against 2026-08-01/02 changes (evidence-based, not inferred)

Three commits made after `M6`, `M7` and `M8`'s results were locked
(`2bd5db0` "consolidate HydroSheaf validation and evidence", 2026-08-01;
`8718d66` "add synthetic validation and RAPM reaction layer", 2026-08-01/02;
`2d4b8af` "complete calibrated synthetic validation gates", 2026-08-02)
modified `hydrosheaf/inference/network_fit.py`,
`hydrosheaf/inference/topology_posterior.py`, `hydrosheaf/sheaf/topology_refine.py`,
`hydrosheaf/models/reactions.py`, `hydrosheaf/models/redox.py`,
`hydrosheaf/config.py`, `hydrosheaf/api.py`, `hydrosheaf/calibration/benchmark.py`,
`hydrosheaf/calibration/validation_workflow.py`, and the entire
`hydrosheaf/validation/` package -- the same integrated decision-utility/
validation-programme layer `O3` D3 already checked for `M3`/`M4`/`M5`.

Unlike `O3` D3, which reasoned only from module-dependency lists, this check
traced the actual per-script import graph of every `M6`/`M7`/`M8` script that
produces a number O4 cites, and inspected the literal diff of every changed
shared module:

- `network_fit.infer_edges` and `topology_refine.refine_edges_with_sheaf`
  gained new branches (`method="null_aware_sheaf"`, a Bayesian
  topology-posterior selection path), but every new branch is gated behind
  either a new opt-in `method` string or `config.topology_posterior_enabled`,
  which defaults to `False` in `hydrosheaf/config.py` and is not set to `True`
  anywhere in `M6`, `M7` or `M8`'s scripts (checked by direct grep). The body
  of `network_fit.fit_network`, which `M7.3`'s `run_m7_3_nonuniqueness.py`
  imports, is textually unchanged in the diff. `M7.4`/`M7.5` use
  `hydrosheaf.sheaf.directed_section`, which was not touched by any of the
  three commits (last changed in an unrelated, much earlier commit,
  `4b44459`).
- `models/reactions.py`'s `build_reaction_dictionary` changed how a missing
  (`None`) SO4/Cl/F/NO3 measurement is handled in its concentration-based
  phase-pruning gate (previously `sample.get("SO4") or 0.0`, silently
  treating a missing measurement as zero and triggering pruning; now
  `_sample_float` returns `None` and pruning is skipped when the measurement
  is absent). This function is imported only by `M5`'s
  `run_m5_inverse_reaction_benchmark.py` (whose 240 scenarios are
  fully-specified synthetic PHREEQC output with no missing ions, so the fix
  is inert there) and by `M6`'s `run_m6_chemistry_robustness.py` and
  `run_m6_chemistry_stress_tests.py`, which `M6/README.md` states are
  "deliberately excluded from `run_m6_q1.py` and from the M6 Q1 submission
  assets." `M6`'s actual Q1 pipeline (`run_m6_field_transfer.py`,
  `run_m6_robustness_diagnostics.py`, `run_m6_null_sensitivity.py`,
  `m6_common.py`) imports only `hydrosheaf.config.Config` and
  `hydrosheaf.phreeqc.runner.run_phreeqc`, plus `M5`'s own `m5_common.py`,
  which defines its own local `REACTIONS` dictionary independent of
  `hydrosheaf/models/reactions.py` (confirmed: `m5_common.py` contains no
  import of `hydrosheaf.models`). O4 therefore treats `M6`'s tier-ablation
  and external-transfer numbers as unaffected by this diff.
- `hydrosheaf/models/redox.py`, `hydrosheaf/calibration/benchmark.py`,
  `hydrosheaf/calibration/validation_workflow.py`, and every module under
  `hydrosheaf/validation/` are imported by none of `M6`, `M7` or `M8`'s
  scripts (checked by repository-wide grep restricted to
  `M6/.../scripts/*.py`, `M7/.../scripts/*.py`, `M8/.../scripts/*.py`).
  `M8`'s confirmatory scripts import only `hydrosheaf.calibration.{adapters,
  definitions, glm, active_learning, bayesian_active_learning}`, none of
  which appear in the changed-file list of any of the three commits.
- One genuine hit was found: `M7/m7_nonuniqueness_benchmark/scripts/
  run_m7_system_acceptance.py` imports `hydrosheaf.api.fit_network_pipeline`,
  and `hydrosheaf/api.py` was modified on 2026-08-01, three days after that
  run was locked (`RUN-M7-SYSTEM-20260728-01`, 2026-07-28). This is a
  secondary, non-load-bearing supporting check in `M7`'s own package (public-
  pipeline system acceptance), not part of `M7.3`'s core locked decision
  table. **O4 excludes this specific run's numbers from every claim**; it is
  not cited anywhere in `Manuscript-O4.md` or `Supplementary-O4.md`.
- As an independent empirical cross-check (not only code reading), three
  headline numbers were re-derived directly from raw locked CSVs rather than
  taken from prose summaries: `M6`'s tier-ablation table (recomputed from
  `results/m6_tier_ablation.csv`, 800 well-tier rows) reproduced
  `docs/m6_results_summary.md`'s E3 table exactly (T4 0.0%/71.00 MRS through
  T0 60.0%/68.43 MRS); `M7.3`'s native and adverse-control evidence-value
  contrasts (recomputed from `results/m7_3_locked/
  evidence_case_bootstrap_contrasts.csv`) reproduced `docs/m7_3_results.md`
  exactly (HAC-HC PR-AUC -0.0060 [-0.0122, -0.0011]; permuted-age HAC-HC
  PR-AUC -0.0754 [-0.1353, -0.0148]); and `M7.4`'s sheaf-vs-graph contrasts
  (recomputed from `results/RUN-M7-SHEAF-VS-GRAPH-20260729-01/
  paired_bootstrap_contrasts.csv`) reproduced `M7/README.md`'s summary
  exactly (affine sheaf vs. graph Laplacian PR-AUC +0.0854 [0.0666, 0.1050]).
  No discrepancy was found between any prose summary and its underlying raw
  result file.

**Conclusion.** Every headline number O4 draws from `M6`, `M7` and `M8` is
produced by code paths that were either untouched by the three post-lock
commits or touched only behind an opt-in flag left at its pre-existing
default; the one number chain found to be genuinely stale (`M7` public-
pipeline system acceptance) is excluded from every claim. This is disclosed
as a verification performed, not assumed; residual risk from a dependency
neither this nor any check can fully rule out is still carried in
Limitations (`Outline.md`; `manuscript/sections`), consistent with `O3` D3's
own residual-risk disclosure.

### D4. Computation and figure authority

Python owns computation and writes read-only tidy CSV exports under
`O4/manuscript/artifacts/data/`. R consumes those exports and owns all
figures under `O4/manuscript/artifacts/figures/`. No figure recomputes a
reported statistic; every plotted value traces to a Python export column,
which in turn traces to a named `M6`/`M7`/`M8` result file (D2/D3).

### D5. Non-duplication with M2.3 and O3

`M2.3` reports the integrated claim-gate architecture and an end-to-end
Ghana demonstration; it does not benchmark `M6`, `M7` or `M8` as its own
primary evidence. `O3` harmonizes `M3`+`M4`+`M5` (capture-vs-correctness
under independent evaluation, Objective 3); O4 harmonizes `M6`+`M7`+`M8`
(robustness, identifiability and calibration under stress, Objective 4). The
two O-series papers share a house style and a companion relationship to
`M2.3` but ask different questions of disjoint evidence: O3's central axis is
whether detection exceeds correctness; O4's central axis is whether an
internally-generated confidence signal (fit quality, entropy reduction,
nominal interval coverage, convergence) tracks an externally-verifiable truth
signal (true identifiability class, independent-generator predictive skill,
realised parameter-recovery error) once evidence is limited, misspecified, or
the calibrating model's form does not match the truth-generating process. No
number is claimed as a fourth demonstration of `M2.3`'s integrated claim, and
no figure or table from `O3` is reproduced.

### D6. Terminology carried over from the wider project

"Sheaf-inspired" or "sheaf-style" is used, never a claim of a formal sheaf
theorem (consistent with `M2.3` D7 and `O3` D6). "Robust" and "identifiable"
are used only in the specific, bounded sense each source component's own
`DECISIONS.md` defines them; O4 does not upgrade a conditional claim (for
example `M7`'s "conditional-integration, not universal-integration" finding)
into an unconditional one. "Noninferior" (`M8` C7) is never restated as
"superior."

### D7. Claims excluded by construction

Field validation of groundwater robustness, connectivity non-uniqueness, or
calibration accuracy beyond each component's own claim ledger; a claim that
the cross-component divergence between internal and external signals was
predicted in advance of seeing `M6`, `M7` and `M8`'s independently locked
results; a second software architecture description; a claim that `M8`'s
controlled-synthetic calibration findings transfer to any field dataset
(`M8` has no field-transfer component at all, unlike `M6`; this asymmetry is
reported, not smoothed over); description of `M7`'s conditional sheaf-vs-
graph finding as general sheaf superiority; description of `M6`'s Northern
Ghana/Talensi/Lower Anayari diagnostics as validated field robustness rather
than a data-limited transfer audit; any number from `M7`'s public-pipeline
system-acceptance run (D3).

### D8. Adversarial review and interval-aware revision (2026-08-03)

An internal adversarial review (`manuscript/review/O4_manuscript_review.md`)
found two issues before the manuscript was finalized. First, `M6`'s tier and
external-transfer percentages (proportions of a known well count) were
reported as bare point estimates in the first draft, despite a known
denominator making a Wilson score interval straightforward to compute; this
mirrors `O3` D7's identical finding for its own capture/correctness
proportions. 95% Wilson intervals were added to every `M6` proportion cited
in the main text (Table 2, Section 4.2) and to `derive_robustness_gradient.py`
(now emits `external_ci_low`/`external_ci_high` for both the tier-ablation
and external-transfer exports); the T3→T2 cliff and the overall T4→T0 range
both remain non-overlapping at 95% confidence, so the finding is strengthened
by the interval, not narrowed by it, unlike `O3`'s reaction-layer finding.
Second, the review found that a hand-written summary row in Table 2 (the
first draft's "M8 ... matched → independent truth" coverage/error row) had
silently substituted the independent-truth-arm's own no-new-measurement
baseline coverage (0.788/0.640) for the matched analytical model's actual
coverage (0.832/0.836) -- the two are numerically close but are not the same
quantity, and the row's label did not match what was actually being
contrasted. This was caught by re-deriving the three-way comparison directly
from `calibration_model_form_shift.csv` (which was already correct; only the
hand-written manuscript prose and table row were wrong) and is corrected in
Table 2, Section 4.2 and the Abstract; the underlying figure (FIG-5) was
never affected, since its R script read the same correct three-condition
export throughout. Two figure-script bugs were also found and fixed while
rendering figures, before this review examined them: `fig03_robustness_
gradient.R`/`fig02_central_divergence.R` initially plotted `robustness_
gradient.csv` without filtering on its `condition` column, so the tier-
ablation series and the evidence-gate-off negative-control series (which
share identical tier labels) were drawn as one zigzagging line; and
`fig02_central_divergence.R`'s three `scale_colour_manual`/`scale_fill_
manual` calls passed an unnamed `labels` vector, which ggplot2 matched to
alphabetically sorted factor levels rather than to the order the vector was
written in, silently swapping the internal/external legend labels in all
three panels. Both were caught by visually inspecting the rendered PNG
output against the intended pattern, not by code review alone, and both are
fixed in the committed R scripts. No other numerical or presentation
discrepancy was found in this pass.

## Claim ledger

| ID | Claim | Evidence | Status |
|---|---|---|---|
| C1 | A common taxonomy classifies retained `M6`, `M7` and `M8` experiments by stress axis (data limitation, evidence misspecification, model-form shift) and by whether an internal confidence signal and an external truth signal are both reported, without contradicting any component's own claim ledger. | Cross-component taxonomy table (TAB-1) | PASS |
| C2 | Under data limitation, evidence misspecification, or a model-form shift away from the calibrating model, an internally-generated confidence signal (fit quality/MRS, posterior-entropy reduction, nominal interval coverage, optimiser convergence) can remain stable or favourable while an externally-verifiable truth signal (true identifiability class, independent-generator predictive skill, realised parameter-recovery error, structural information rank) degrades or diverges, in all three layers. | `M6` tier-ablation (flat MRS, collapsing identifiability), `M7` adverse-control entropy/PR-AUC divergence, `M8` coverage collapse under independent numerical truth and kinetic rank-one confound (TAB-2, FIG-2) | PASS |
| C3 | The size and shape of this divergence differs by layer and by stress type; it is not one universal curve. | Within-component detail (FIG-4); `M6`'s divergence is gradual across an ordered tier ladder, `M7`'s appears only under explicit adverse controls (native evidence combination shows a small, calibrated effect), `M8`'s is parameter-specific and sign-dependent (dispersivity improves, decay worsens, under the same design and the same model-form shift) | PASS |
| C4 | Negative and structural controls bound each layer's finding: `M6`'s evidence-gate-off ablation shows 0% non-identifiable at every tier, proving the tier "collapse" is the conservative prior, not an artefact of the classifier; `M7`'s permuted-map/permuted-stream controls are beaten by every native condition on every primary metric; `M8`'s grid-convergence gate and doubled-product off-ridge check confirm the kinetic rank deficiency is structural, not a dead experiment. | `M6` Extended Data Figure ED3a; `M7` native-vs-permuted contrasts; `M8` Table 3, Section 3.3 | PASS |
| C5 | None of the three layers licenses a field-validated robustness, identifiability, or calibration claim; `M8` in particular has no field-transfer component. | `M6` DECISIONS.md honesty boundary; `M7` Ghana scope decision; `M8` scope statement (controlled synthetic only) | ABSTAIN (explicitly out of scope) |
| C6 | This is a second HydroSheaf framework or architecture paper, or a fourth validation of any `M6`/`M7`/`M8` finding beyond what each already claims. | Scope is a benchmark synthesis; see D1, D5 | ABSTAIN (explicitly out of scope) |
