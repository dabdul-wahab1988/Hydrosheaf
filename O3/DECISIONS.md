# O3 decisions and claim ledger

O3 is a new manuscript package for Objective 3 of Chapter 1 (HydroSheaf PhD).
It benchmarks, but does not restate the software architecture of, `M3`
(nuclear-tracer/age), `M4` (MODPATH-referenced topology) and `M5`
(PHREEQC-screened reaction). It is written after, and positioned as the
depth companion to, `M2.3` (the HydroSheaf framework paper, same target
journal).

## Locked decisions

### D1. Package status and journal target

O3 is a new, standalone manuscript targeting Computers & Geosciences,
distinct from `M2.3`. It replaces the earlier plan of submitting `M3`, `M4`
and `M5` as three fully separate journal articles (Applied Geochemistry,
Advances in Water Resources, Water Resources Research per each package's own
prior outline); those three packages and their journal-specific Word/markdown
deliverables are retained unchanged as historical records of the underlying
benchmarks and are not superseded as evidence sources.

### D2. No re-simulation

O3 performs no PHREEQC, MODFLOW/MODPATH, or TracerLPM re-runs. Every number
reported in O3 is read from an already-locked result file under
`M3/m3_age_benchmark/results/`, `M4/m4_topology_benchmark/results/` or
`tables/Manuscript_Ready/`, or `M5/m5_inverse_reaction_benchmark/results/` or
`tables/`, or is a harmonization-layer statistic (a ratio, a rank, a common
grouping label) computed by `O3/analysis/python/` directly from those files
with no change to the underlying inference. This mirrors how `M2.3` itself
reused `M3`'s pre-2026-08-01 results without rerunning them (see
`M2.3/DECISIONS.md` D5 and `M2.3/analysis/python/derive_model_evidence.py`).

### D3. Staleness check against the 2026-08-01/02 validation-layer changes

Commits `2bd5db0`, `8718d66` and `2d4b8af` (2026-08-01/02) added or grew
`hydrosheaf/validation/{reaction_rapm.py, age_competent_baseline.py,
kinetics_specialist_benchmark.py, decision_utility.py, performance_contract.py,
specialist_calibration.py}` -- the integrated programme/decision-utility layer
`M2.3` benchmarks as its own component and integrated gates (its claims C2-C5).
Checked directly: `M3`'s benchmark exercises `hydrosheaf/nuclear/*`
(frozen 2026-07-31, commit `2e73d51`, before the validation-layer changes);
`M4`'s independent-mode topology benchmark exercises head-gradient/
depth/hydrostratigraphic graph construction against the MODPATH archives, not
the new validation-layer modules; `M5`'s benchmark exercises
`hydrosheaf/models/reactions.py` and the sparse-inversion solver, last
changed 2026-08-01 in the same consolidation commit that also touched
`inference/edge_fit.py` and `inference/network_fit.py` used by the transport
correction step. The reaction dictionary itself (`models/reactions.py`) was
not rewritten, only extended; M5's own locked stoichiometry table
(`tableS1_reaction_stoichiometry.csv`) was not invalidated by that diff.
Residual risk that a downstream module changed in a way not caught by this
check is disclosed in Limitations rather than assumed away.

### D4. Computation and figure authority

Python owns computation and writes read-only tidy CSV exports under
`O3/manuscript/artifacts/data/`. R consumes those exports and owns all
figures under `O3/manuscript/artifacts/figures/`. No figure recomputes a
reported statistic; every plotted value traces to a Python export column.

### D5. Non-duplication with M2.3

O3's Discussion states its relationship to M2.3 explicitly (see
`Outline.md`, Positioning). O3 does not reproduce M2.3's architecture
figures, does not restate M2.3's integrated-decision-utility claim (C5) as
its own, and does not claim field validation beyond what each of M3, M4 and
M5 already permits in their own claim ledgers. Where a number appears in
both papers (for example, topology F1 0.618 or reaction phase F1 0.563), it
is cited to the same primary result file in both, and O3 states that M2.3
reports the same primary result at lower resolution.

### D6. Terminology carried over from the wider project

"Sheaf-inspired" or "sheaf-style" is used, never a claim of a formal sheaf
theorem (consistent with `M2.3` D7). "Independent" is reserved for a
benchmark row evaluated with no access to the reference used to score it, and
is never used for a calibrated, prior-informed, or emulation row without the
qualifying word attached in the same sentence.

### D7. Adversarial review and interval-aware revision (2026-08-03)

An adversarial review (`manuscript/review/O3_manuscript_review.md`, score
60/100, return for major revisions, no fabrication or critical-flaw override
triggered) found that the central capture/correctness comparison was
reported as bare point estimates despite source data supporting confidence
intervals in at least the reaction layer. Wilson score 95% intervals were
computed for every proportion in Table 2 from the underlying counts already
in each component's own result files (no new simulation) and revealed that
the reaction layer's capture and correctness intervals overlap while
topology's and age's do not. Every section stating the central claim was
revised to report this qualification (`manuscript/review/
O3_reviewer_issue_resolution.md` records the full resolution ledger). This
is treated as a structural finding, not an error correction: it narrows
claim C2 below to what two of the three layers, not all three, actually
support at conventional confidence.

### D8. Claims excluded by construction

Field validation of groundwater age, directed connectivity or reaction
mechanism beyond each component's own claim ledger; a claim that O3's
cross-component pattern was predicted in advance of seeing any of the three
benchmarks; validation of M3's identified-TTD development-stage path; a
second software architecture description; universal superiority over
MODFLOW, MODPATH, TracerLPM or PHREEQC.

## Claim ledger

| ID | Claim | Evidence | Status |
|---|---|---|---|
| C1 | A common evidence taxonomy classifies every retained row of the M3, M4 and M5 benchmarks without contradicting any component's own claim ledger. | Cross-component taxonomy table (TAB-1) | PASS |
| C2 | Under independent, prior-free, uncalibrated evaluation, the point estimate of recall/sensitivity exceeds precision/specificity in all three inference layers; with 95% Wilson intervals attached, the gap is statistically resolved (non-overlapping) for M4 and M3 but not for M5, whose intervals overlap. | M3 held-out parity, M4 independent topology metrics, M5 primary reaction metrics, all with 95% CIs (TAB-2) | PASS, qualified per interval check added after adversarial review (see `manuscript/review/`) |
| C3 | Negative and informed-structural controls bound each layer's independent result from below and above respectively. | M3 randomised-graph delta, M4 negative-control rows, M5 conventional-PHREEQC-baseline multiplicity | PASS |
| C4 | Calibrating age and reaction diagnostics against the same reference used for scoring narrows the apparent performance gap, and is reported as emulation rather than independent agreement. | M3 calibrated parity, M5 MRS calibration on held-out archetype (TAB-3) | PASS, non-independent by construction |
| C5 | The three benchmarks differ by roughly two orders of magnitude in scenario/fit count, and only M5 records per-fit runtime. | Benchmark-scale table (TAB-4) | PASS |
| C6 | Field-transfer scope is uneven across the three layers (chemistry: three field sites; topology: three public MODFLOW archives; age: public national release only, no independent field transfer). | FIG-6 / field-transfer inventory | PASS |
| C7 | Any component's inference is field-validated by virtue of appearing in this comparison. | No independent field truth for any of the three | ABSTAIN |
| C8 | This is a second HydroSheaf framework or architecture paper. | Scope is a benchmark synthesis; see D1, D5 | ABSTAIN (explicitly out of scope) |
