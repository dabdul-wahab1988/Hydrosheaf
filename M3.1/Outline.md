# M3.1 outline and word contract

## Target

Applied Geochemistry (Elsevier), Research Article — same venue as the locked
M3 submission package in `M3/M3_geochemistry/`. M3.1 is a version upgrade, not
a new submission track: it supersedes M3 once regenerated evidence and review
are complete. M3 is retained unchanged as the prior locked record.

## Why an upgrade rather than a resubmission from scratch

M3 was accuracy-locked on 2026-07-28 at commit `47fbf7c`. Since then, 11 commits
on `codex/m3-correctness` added a set-valued ("identified-TTD") age-inference
diagnostic, an independent graph-compatibility audit, a tracer-reconciliation
negative result, and a leakage fix to `calculate_tracer_reliability_weights`
(the reference age can no longer reach tracer-weighting, closing a path that
could otherwise leak the evaluation target into the fit). Every post-lock
document that discusses the design-matrix numbers states explicitly that they
are unchanged; M3.1 verifies that by rerunning the full locking command
(`run_m3_manuscript_analysis.py --full --age-steps 90`) against current `HEAD`
in an isolated copy (`M3.1/m3_age_benchmark/`) and diffing every reportable
statistic against the M3 lock, rather than trusting the prior claim.
Independently of whether any number moves, two genuinely new, honestly-scoped
findings exist that M3 could not have reported and that belong in an upgrade.

## Positioning

Unchanged core thesis: graph topology is evaluated as a falsifiable prior on
multi-tracer groundwater-age inference, not assumed to help. M3.1 adds a
second, independent falsification test (a set-valued compatibility audit) that
reaches the same negative conclusion by a different route, and reports a new,
unresolved tracer-identifiability limitation discovered while building it.

## Central argument

> Graph topology should be tested, not assumed. Two independent tests — a
> reference-agreement/negative-control benchmark and a set-valued compatibility
> audit — agree that no candidate graph earns its keep on this dataset. The same
> audit process also surfaces a tracer-level problem: for roughly a quarter of
> evaluable cases, no transit-time distribution reconciles the measured tracers
> within their stated uncertainty, and seven tested explanations all fail to
> account for it. Reporting that honestly is itself part of the framework's
> contribution.

Evidence, in order:

1. Reference-agreement emulation reproduces published USGS ages closely only on
   an identifiability-screened subset (329 of 1,272); wider claims are not
   supported.
2. No candidate graph prior meets the pre-declared robust-improvement rule on
   that subset; implausible controls (wrong-direction, randomised) are more
   harmful, showing the model responds to topology without confirming any
   candidate topology is correct.
3. Leakage-guarded withholding of five tracers shows no meaningful predictive
   graph gain, and CFC comparisons remain non-estimable for lack of eligible
   edges.
4. [NEW] An independent set-valued (interval-bound) compatibility audit,
   run as an exploratory diagnostic rather than a confirmatory analysis, likewise
   finds no support for graph-based tightening — a second method reaching the
   same conclusion.
5. [NEW] The same audit found that 27.85% of evaluable cases admit no feasible
   transit-time distribution at all under standard tracer physics and declared
   uncertainty. Seven candidate explanations (shared correction-parameter error,
   a common-mode scale error, a single culprit tracer, and others) were tested
   and none accounts for it; a follow-up synthetic control on an independent
   MODFLOW/MODPATH generator also did not support the leading hypothesis. The
   cause is reported as unresolved, not attributed.

## Word budget: 6,000 main-text prose words

Prose only; excludes title, headings, references, table content, figure
captions and equations. Reader-facing language throughout: define terms on
first use, minimise nested caveats per sentence, move derivations and
implementation detail to Supplementary Methods.

| Section | Target | Purpose |
|---|---:|---|
| Abstract | 220 | Problem, what was tested, headline result, honest boundary |
| 1. Introduction | 700 | Why groundwater age is ambiguous, what a graph prior could add, the test design |
| 2. Methods (main) | 1,300 | Datasets, tracer suite, model idea, graph idea, scenario matrix, metrics — in plain language |
| 3. Results | 1,750 | Coverage, baseline agreement, graph findings, tracer withholding, new identifiability/infeasibility finding |
| 4. Discussion | 1,650 | What it means, comparison to LPM practice, data-limited aquifer implications, limitations, recommendations |
| 5. Conclusions | 230 | Exact claim and what remains open |
| **Total** | **5,850** | leaves ~150-word margin |

Supplementary Information: technical methods (full mathematics), all
derivations, the identified-TTD formal definitions, the seven-explanation
exclusion protocol, and supplementary figures/tables. No prose cap; write it
completely rather than to a length target.

## Display items

| ID | Type | Content | Status |
|---|---|---|---|
| FIG-1 | Figure | Study-area map: USGS public-supply sites used, by Principal Aquifer / study unit (US states basemap) | **NEW** — addresses the current manuscript's lack of any map |
| FIG-2 | Figure | Atmospheric input histories (³H, SF₆, CFC-12) | carried over, replotted in R |
| FIG-3 | Figure | Active design-matrix reference agreement and reportability | carried over, replotted in R |
| FIG-4 | Figure | Strict reported-configuration emulation vs. USGS reference ages (N=329) | carried over, replotted in R |
| FIG-5 | Figure | Graph-age effects on the 329 strict reportable nodes | carried over, replotted in R |
| FIG-6 | Figure | Leakage-guarded target-tracer withholding | carried over, replotted in R |
| FIG-7 | Figure | Set-valued compatibility audit and infeasibility summary | **NEW**, exploratory, clearly labelled |
| TAB-1 | Table | Nuclear/dissolved-gas tracer properties | carried over |
| TAB-2 | Table | Design-matrix performance by scenario | carried over, regenerated |
| TAB-3 | Table | Graph-age benchmark (329 nodes) | carried over, regenerated |
| TAB-4 | Table | Modelling-mode comparison (full vs. common support) | carried over, regenerated |
| TAB-5 | Table | Statistical significance (paired tests) | carried over, regenerated |

Age-class diagnostics, residual diagnostics, edge-geometry, CFC withholding
detail, and the reproducibility inventory move to Supplementary Figures/Tables
to protect the main-text word budget. Exact supplementary numbering is fixed
at artifact-registry time (B3), after regenerated data is in hand.

## Claims that must not appear

- Field validation of true groundwater age, confirmed flow paths, or general
  graph benefit.
- The set-valued/identified-TTD extension described as confirmatory,
  validated, or as replacing the scalar M3 benchmark.
- An attributed cause for the 27.85% infeasibility rate (all seven tested
  explanations failed; the M7.6 synthetic control did not support the leading
  one either).
- Ghana-specific or semi-arid African field application (this remains a
  US public-data benchmark; transfer is a stated future hypothesis only).
- MRVA used as cross-aquifer graph-age replication.
- Any manuscript number not reproduced by the regenerated M3.1 pipeline run
  against current `HEAD`.

## Open decisions for the user

1. Whether to include claim 4-5 above (the new set-valued compatibility audit
   and the unresolved infeasibility finding) as a compact new Results/Discussion
   subsection with full detail in Supplementary — recommended, since it is
   genuinely new, honestly bounded, and reinforces rather than complicates the
   paper's thesis — or to leave it out of M3.1 entirely and hold it for a
   separate future paper once the cause is understood.
2. Title: recommend keeping the current title unchanged (the core claim it
   describes is unchanged); open to a light adjustment if the new content
   should be foregrounded.
