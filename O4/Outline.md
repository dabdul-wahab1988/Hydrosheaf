# O4 outline and word contract

## Target

Computers & Geosciences, research article. Sibling paper to `O3` (Objective
3) and companion to `M2.3` (the HydroSheaf framework article), same target
journal. O4 does not introduce new software architecture and must not be
read as a second framework paper.

## Positioning and non-duplication statement

Three components of HydroSheaf were independently benchmarked for
robustness, identifiability and calibration validity: `M6` field-transfers
the reaction-inference workflow to three real Ghanaian datasets under a
five-level evidence ladder and asks whether it degrades gracefully; `M7`
tests conditional evidence integration against fresh independent MODFLOW 6/
MODPATH 7 truth and asks when combining hydraulic, tracer and chemical
evidence helps, does nothing, or actively misleads; `M8` runs three
prospectively locked controlled-synthetic calibration protocols against
known parameters and asks whether optimiser convergence and nominal
uncertainty intervals actually track parameter recovery. Each has its own
locked `DECISIONS.md`, claim ledger and manuscript-ready display items, each
was built and locked independently, on its own timetable, and each already
has its own standalone manuscript (`M6/.../Manuscript-Final-Revised.md`,
`M7/.../Manuscript-Final-Revised-M7.6-20260731.docx`,
`M8/.../Manuscript-Final.md`).

`M2.3` reports the integrated claim-gated architecture and an end-to-end
Ghana demonstration; it does not attempt the cross-component comparison this
paper makes. `O3` harmonizes `M3`+`M4`+`M5` on a different axis (does
detection exceed correctness under independent evaluation). O4's
contribution is the first place `M6`, `M7` and `M8` are read side by side,
and the only place a single, disclosed question is asked of all three: does
an internally-generated confidence signal, computed without reference to
ground truth, track an externally-verifiable truth signal once the evidence
is limited, an evidence stream is misspecified, or the model used to
calibrate differs in form from the process that generated the data? A reader
of O4 does not need `M6`/`M7`/`M8`'s full ablation depth; a reader of any one
of those three does not get the cross-component pattern. Each paper cites
the others as companion pieces, and none claims another's result as its own
primary contribution.

## Central argument

> Three physically unrelated stress tests of the same evidence-integration
> and inverse-inference framework -- stripping measurement tiers from real
> field chemistry, combining or deliberately misspecifying hydraulic/tracer/
> chemical evidence streams against independent synthetic truth, and shifting
> from a matched calibration model to an independently generated numerical
> truth -- were each designed, executed and locked independently. In all
> three, an internally-generated signal of confidence or fit quality --
> mean mechanism-resolution score, posterior-entropy reduction, nominal 95%
> interval coverage, optimiser convergence -- can remain flat, improve, or
> stay nominally healthy while an externally-verifiable truth signal --
> true identifiability class against synthetic ground truth, predictive skill
> against an independent MODFLOW 6/MODPATH 7 generator, realised parameter-
> recovery error, structural Fisher-information rank -- degrades sharply or
> diverges in sign. Numerical health is not evidence of identifiability,
> integration value, or calibration validity.

The evidence, in order:

1. A common taxonomy classifies every retained `M6`, `M7` and `M8` experiment
   by stress axis (data limitation, evidence misspecification, model-form
   shift) and by whether it reports a paired internal/external signal (Table
   1); no component's own claim rule is altered by the classification.
2. Under a five-level evidence-tier ablation, `M6`'s mean mechanism-
   resolution score stays within 68.4-71.0 across all five tiers while the
   fraction of Northern Ghana wells correctly flagged non-identifiable rises
   from 0.0% (full tier) to 60.0% (majors only), with a sharp 51.3-point
   cliff at the single step that removes Sr/SiO2; a synthetic validation with
   known truth confirms the direction (structural rank rises 8 to 11, exact-
   mineral F1 stays 0.49-0.74 at every tier) rather than merely asserting it.
3. Under `M7`'s adverse controls, every one of permuted age, permuted
   hydraulics and jointly permuted evidence reduces posterior-edge entropy by
   more than the native full-panel condition (joint: -0.0706 nats) while
   simultaneously reducing PR-AUC (joint: -0.139) and, for two of the three
   controls, increasing Brier score or log-loss; native age addition alone
   already shows the same sign pattern at smaller magnitude (entropy -0.0006,
   PR-AUC -0.0060, both 95% CIs excluding zero).
4. Under `M8`'s shift from a matched analytical calibration model to an
   independently generated 240-cell numerical truth, the same 50 d
   observation that improves dispersivity recovery (median absolute log10
   error 0.826 to 0.167) collapses that parameter's own linearised 95%
   coverage from 0.788 to 0.02, and simultaneously degrades decay recovery
   (0.137 to 0.154) while its coverage falls from 0.64 to 0.004; separately,
   every calibration in the kinetic rate-law experiment converges and
   reports a well-defined objective value while the underlying Fisher
   information matrix is numerically rank one (infinite condition number),
   which only an independent surface-area observation resolves.
5. Structural and negative controls bound each layer: `M6`'s evidence-gate-
   off ablation returns 0% non-identifiable at every tier, proving the
   apparent "collapse" is the framework's conservative prior, not a
   classifier artefact; `M7`'s native evidence conditions beat every
   permuted-map or permuted-stream adverse control on every primary metric;
   `M8`'s grid-convergence gate and off-ridge product-doubling check confirm
   the kinetic rank deficiency is a structural property of the rate law, not
   an under-powered experiment.
6. The three layers' field- and archive-transfer scope is uneven and is
   reported as such: `M6` transfers to three real Ghanaian datasets under
   disclosed data limitations; `M7` uses fresh independent MODFLOW 6/MODPATH
   7-generated truth throughout, with the Northern Ghana workbook audited
   only for component-diagnostic readiness, not topology or age validation;
   `M8` has no field-transfer component at all, a limitation reported
   directly rather than implied away.

## Word budget: approximately 6,000 main-text prose words

Prose only; excludes title, headings, references, table content, figure
captions and equations.

| Section | Target | Purpose |
|---|---:|---|
| Abstract | 220 | Problem, contribution, headline pattern, scope (excluded from the 6,000 total) |
| 1. Introduction | 900 | Three stress tests, why fit quality is not identifiability, the harmonization gap, contribution, objectives |
| 2. Data | 850 | Full reader-facing description of every dataset the three benchmarks consume |
| 3. Methods | 1,300 | Harmonization taxonomy; each component's own stress design summarised; what is and is not recomputed |
| 4. Results | 1,900 | Taxonomy application; cross-component internal-vs-external divergence pattern; within-component detail; scale and transfer scope |
| 5. Discussion | 700 | What the recurring pattern means; relation to M2.3 and O3; comparison with published identifiability/calibration practice; limitations |
| 6. Conclusion | 130 | Exact claim and next stage |
| **Total** | **~6,000** | |

Supplementary Information: approximately 3,000 words, carrying full
per-component stress-design detail (tier construction, adverse-control
definitions, the independent-truth generator, hyperparameter and bootstrap
detail) that Section 3 can only summarise.

## Display items

| ID | Type | Content |
|---|---|---|
| FIG-1 | Figure | Three stress-test pipelines and the common internal-vs-external signal taxonomy |
| FIG-2 | Figure | Central result: paired internal confidence signal vs external truth signal, one panel per component, under increasing stress |
| FIG-3 | Figure | M6 tier-ablation cliff: mean MRS (flat) vs percent non-identifiable (collapsing), five tiers, with the evidence-gate-off negative control |
| FIG-4 | Figure | M7 entropy reduction vs predictive-skill change, native and three adverse controls, with 95% bootstrap intervals |
| FIG-5 | Figure | M8 coverage-vs-recovery divergence under matched vs independent numerical truth, and the kinetic rank-one structural result |
| FIG-6 | Figure | Benchmark scale and field-/archive-transfer scope across all three components |
| TAB-1 | Table | Common stress/signal taxonomy applied to every retained experiment |
| TAB-2 | Table | Headline internal-signal vs external-signal pairs, all three components, with 95% intervals where the source component reports them |
| TAB-3 | Table | Negative- and structural-control results bounding each layer's central finding |
| TAB-4 | Table | Benchmark scale: replicate/case/bootstrap counts, external reference type, field-transfer scope |

## Claims that must not appear

- A fourth, independent validation of any component beyond what its own
  `DECISIONS.md` and claim ledger already permit.
- A new architecture, software contribution, or sheaf-theoretic result (that
  is `M1`'s and `M2.3`'s domain, not O4's).
- A claim that `M7`'s conditional sheaf-versus-graph or hybrid finding is a
  general superiority result.
- A claim that `M8`'s controlled-synthetic calibration findings transfer to
  any field dataset; `M8` has no field-transfer component, and this asymmetry
  with `M6` is reported directly.
- A claim that the cross-component divergence pattern was predicted in
  advance; it is reported as a retrospective synthesis of three independently
  locked stress tests, not a preregistered joint test.
- Description of `M6`'s Northern Ghana/Talensi/Lower Anayari diagnostics, or
  `M7`'s Ghana workbook audit, as validated field robustness or connectivity.
- Any number from `M7`'s public-pipeline system-acceptance run
  (`RUN-M7-SYSTEM-20260728-01`), which imports a `hydrosheaf` module changed
  after that run was locked (`DECISIONS.md` D3).
