# O3 outline and word contract

## Target

Computers & Geosciences, research article. Companion paper to `M2.3` (the
HydroSheaf framework article, same target journal). O3 does not introduce new
software architecture and must not be read as a second framework paper.

## Positioning and non-duplication statement

Three components of HydroSheaf were independently benchmarked against three
physically distinct external references: `M3` compares nuclear-tracer/lumped-
parameter age inference against the USGS national groundwater-age release;
`M4` compares reduced-order graph topology against MODPATH particle-tracking
connectivity from three public MODFLOW archives; `M5` compares sparse linear
reaction inversion against a 240-scenario live-PHREEQC factorial synthetic
benchmark. Each has its own locked `DECISIONS.md`, claim ledger, and
manuscript-ready display items, and each was built, reviewed and locked
independently, on a different timetable, by a different evaluation design.

`M2.3` already reports a compressed, three-line version of this pattern
(coverage/error, log-loss, precision/recall) inside a 7,000-word
software-and-methods article whose primary contribution is the claim-gated
architecture itself. O3's contribution is different in kind, not degree: it is
the first place any of the three benchmark designs are placed side by side,
the first place their differing notions of "independent" and "calibrated" are
reconciled into one taxonomy, and the only place their computational scale and
full ablation/sensitivity detail is reported together. A reader of O3 does not
need architecture detail; a reader of M2.3 does not get benchmark depth. Each
paper cites the other as the companion piece, and neither claims the other's
result as its own primary contribution.

## Central argument

> Three physically unrelated groundwater inference problems -- where water
> moves, how long it has resided, and what it has reacted with -- were each
> benchmarked against an independent external reference under a common
> design: an uninformed or prior-free evaluation, one or more negative
> controls, and an explicit reportability or abstention rule. In all three,
> the independent, uncalibrated evaluation shows recall or sensitivity
> substantially exceeding precision or specificity: fractions detected are
> higher than fractions correct. In the two layers where a calibrated
> alternative exists (age, reaction), calibrating against the same reference
> the benchmark is judged on collapses most of the apparent gap -- and stops
# being an independent test the moment it does so.

The evidence, in order:

1. A common taxonomy classifies the retained rows of all three benchmarks by
   independence, control type and claim scope (Table 1); none of the three
   components' own claim rules is altered by the classification.
2. Under independent, prior-free, uncalibrated evaluation: topology recall
   0.845 against precision 0.487 (F1 0.618, MODPATH reference); age held-out
   agreement log10 R2 0.482 (56.4% within a factor of 2) with only 25.9-53.1%
   of fitted rows passing the reportability guard at all; reaction phase F1
   0.563 with false-discovery rate 0.361, and 55.0% of the lowest-residual
   quartile still below phase F1 0.80.
3. Negative and structurally informed controls bound each result: topology
   collapses to F1 <= 0.14 under every negative control and reaches F1 0.552
   from receptor-set knowledge alone with zero hydraulic information; age
   graph-regularisation fails a joint robust-improvement rule against
   single-node inference and randomised graphs increase error by 0.65 log10
   units; reaction's own conventional-PHREEQC baseline is feasible for only
   45.8% of scenarios at strict tolerance and returns 185.8 feasible models
   per scenario on average.
4. Calibrating age to the same held-out reference raises log10 R2 from 0.482
   to 0.962; calibrating reaction's Mechanism Resolution Score against a
   held-out archetype reaches only 48.9% four-class accuracy -- a smaller,
   more honest gain, reported beside the age result rather than in isolation,
   because the two calibration exercises are not equally informative and the
   paper does not let the more flattering one stand for both.
5. Benchmark scale differs by roughly two orders of magnitude across the
   three components (a few hundred graph edges; low thousands of TracerLPM
   fits; 21,600 factorial reaction fits), and only the reaction benchmark
   records per-fit runtime -- a comparison the paper reports rather than
   silently normalises away.
6. Field-transfer scope (Northern Ghana, Talensi, Lower Anayari for
   chemistry; Savage, Great Miami, Long Island archives for topology; no
   field transfer for age beyond the public USGS release) is uneven across
   the three layers and is reported as such.

## Word budget: 7,800 main-text prose words

Prose only; excludes title, headings, references, table content, figure
captions and equations. Follows the project's `methods_manifest.json`
word-count policy.

| Section | Target | Purpose |
|---|---:|---|
| Abstract | 250 | Problem, contribution, headline pattern, scope (excluded from the 7,800 total) |
| 1. Introduction | 1,000 | Three inference problems, why they are usually benchmarked separately, the harmonization gap, contribution, objectives |
| 2. Data | 1,300 | Full reader-facing description of every dataset the three benchmarks consume |
| 3. Methods | 1,600 | Harmonization taxonomy; each component's benchmark design summarised; what is and is not recomputed |
| 4. Results | 2,400 | Taxonomy application; cross-component recall/precision pattern; screening-to-calibration gap; per-component detail; benchmark scale; field-transfer scope |
| 5. Discussion | 1,250 | What the recurring pattern means; relation to M2.3 and to each standalone component; comparison with published practice; limitations |
| 6. Conclusion | 250 | Exact claim and next stage |
| **Total** | **7,800** | |

Supplementary Methods: approximately 3,000 words, sized to carry the full
per-component benchmark design detail (reference construction, reportability
guards, control definitions, hyperparameter selection) that Section 3 can
only summarise.

## Display items

| ID | Type | Content |
|---|---|---|
| FIG-1 | Figure | Three benchmark pipelines and the common evidence taxonomy |
| FIG-2 | Figure | Dataset and benchmark-scale overview (scenario/fit counts, external references, field sites) |
| FIG-3 | Figure | Cross-component recall/sensitivity versus precision/specificity, with negative and informed controls |
| FIG-4 | Figure | Screening-to-calibration/emulation gap for age and reaction, reported side by side |
| FIG-5 | Figure | Within-component identifiability detail: age graph falsification deltas, topology precision-recall by evidence level, reaction phase versus equivalence-class recovery by method |
| FIG-6 | Figure | Field-transfer scope map and scores across all six field/archive sites |
| TAB-1 | Table | Common evidence taxonomy applied to every retained benchmark row |
| TAB-2 | Table | Headline independent/uncalibrated recall-type and precision-type metrics, all three components |
| TAB-3 | Table | Screening/independent versus calibrated/prior-informed/emulated metric, with non-independence stated per row |
| TAB-4 | Table | Benchmark scale: scenario counts, fit counts, external reference type, recorded runtime |

## Claims that must not appear

- A fourth, independent validation of any component beyond what its own
  `DECISIONS.md` and claim ledger already permit.
- A new architecture, software contribution, or sheaf-theoretic result (that
  is M2.3's and M1's domain, not O3's).
- Field validation of groundwater age, directed connectivity, or reaction
  mechanism.
- Description of the calibrated/emulated numbers as independent agreement.
- A claim that the cross-component pattern was predicted in advance; it is
  reported as a retrospective synthesis of three independently locked
  benchmarks, not a preregistered joint test.
- Description of M3's identified-TTD development-stage path (frozen WP0,
  27.85% local infeasibility) as validated; it is out of scope for O3 and
  referenced only as a forward pointer to future work, consistent with
  `M3/m3_age_benchmark/README.md`.
