# M8 — Manuscript Outline

**Working title**

> **Calibration is not identification: parameter recovery, uncertainty calibration and optimal measurement design in groundwater inverse modelling**

**Target journals (in order of fit)**

1. *Water Resources Research* (AGU, Q1) — the natural home; WRR publishes inverse-method
   evaluation and has a long PEST/identifiability tradition.
2. *Journal of Hydrology* (Q1) — broader, faster, accepts benchmark studies.
3. *Environmental Modelling & Software* (Q1) — if the framing leans toward the released
   software and reproducibility contract.

WRR is the recommendation: the result is a methods-evaluation result, not a software
announcement, and WRR readers already accept that non-uniqueness is the central problem.

---

## 1. Why this is M8 and not a section of M2–M7

M3–M7 established, for **mechanisms and topology**, that a low residual does not imply a
correct answer. Every one of those benchmarks was ultimately limited by the same thing:
the comparison target was another model's output (USGS/TracerLPM ages, MODPATH edges), so
the ceiling was *emulation*, not accuracy.

Calibration removes that ceiling. Synthetic parameters are known exactly, so recovery is
measurable directly. **M8 is therefore the one place in the programme where a legitimate
positive result is available** — and the pilot confirms one exists.

It also closes an audit gap: `hydrosheaf/calibration` is the single largest subsystem
(6,249 lines, PEST-style GLM engine, seven domain adapters, active learning) and is
exercised by **zero** manuscript benchmarks to date, while passing 77/77 of its unit tests.

---

## 2. Contribution claims (ranked, and each already evidenced or directly runnable)

**C1 — The core dissociation.** Optimiser success, objective value, and parameter
identifiability are three different things, and in groundwater inverse problems they come
apart routinely. *Pilot evidence: φ = 7.2×10⁻¹³ with 39.9% median parameter error and 7.5%
interval coverage.*

**C2 — Linearised uncertainty under-covers even when the problem is well-posed.** In the
fully identifiable design, nominal 95% intervals achieved 76–80% coverage at realistic
noise. This is a quantified, reportable miscalibration of standard practice.

**C3 — The covariance fallback detects *exact* rank deficiency only, and is therefore the
wrong flag; the condition number is the right one.** When `J^T J` is singular the engine
falls back from matrix inverse to SVD pseudoinverse. Pilot: 0% fallback in the identifiable
design, 22–43% in the exactly-singular designs — **but 0% across the entire conditioning
sweep, even at cond(FIM) = 6×10¹²**, because a near-singular matrix is still numerically
invertible. *This is a correction to the naive reading and a concrete recommendation:
report cond(FIM), not the covariance path.*

**C3b — The identifiability threshold is locatable and lies at cond(FIM) ≈ 10²–10³.**
Median relative parameter error jumps from 0.139 at cond = 116 to **0.546 at cond = 735** —
a fourfold degradation across less than one order of magnitude in conditioning. That is an
actionable number a practitioner can compute before trusting a calibration.

**C3c — Coverage rises to 1.00 as identifiability fails, and this is a trap.** In the sweep,
empirical coverage of nominal 95% intervals *increases* monotonically to 1.00 as the problem
becomes unidentifiable, because the intervals inflate until they trivially contain the truth
while the point estimate is 70% wrong. **Coverage alone is not a sufficient diagnostic; it
must be reported jointly with interval width.** This is a genuinely counter-intuitive result
and a strong candidate for the paper's second headline.

**C4 — Planned, not yet evidenced: a recovery envelope across seven hydrogeological inverse
problems.** Kinetics and transport now have confirmatory controlled-synthetic evidence.
Vadose, nitrate source, temporal age, topology and isotope mixing remain future work and
cannot be included in the present contribution claim.

**C5 — A front observation repairs a late-only transport design, but the anticipated
parameter-target split was not confirmed.** In the locked 250-replicate production-adapter
benchmark, dispersivity-targeted, decay-targeted, D-, A-, E- and balanced local criteria all
selected 50 d. Median absolute log10 error fell from 0.2154 to 0.0173 for dispersivity and
from 0.0182 to 0.0146 for decay. The analytic 60-fold result remains a positive control and
must not be described as the real-adapter effect.

**C5b — Kinetic sampling has a structural limit.** The production PHREEQC adapter produced
rank-one information for calcite rate constant and surface area with one or four residence
times because the rate law identifies `k*A`. An independent surface-area observation
restored rank two. This is evidence about the information required to identify the two
parameters, not field validation of a surface-area measurement programme.

**C6 — The existing active-learning module is heuristic, not information-theoretic.**
`_compute_expected_benefit` scores candidate measurements from categorical evidence classes
and never touches the Jacobian. Benchmarking it against E-optimal design, with random and
oracle controls, is a falsifiable test of an existing 821-line claim. All three outcomes are
publishable: it works (validates the module), it does not (honest negative, replace it), or
it works only in specific regimes (conditional result).

---

## 3. Structure

**1. Introduction** — inverse modelling in data-limited aquifers; PEST/GLM as standard
practice; the difference between fitting, estimating and identifying; why the M3–M7
emulation ceiling motivates a known-truth design.

**2. Methods**
- 2.1 The GLM engine (`PESTGLM`): LM via trust-region reflective, log transforms, bounds,
  parallel Jacobian, covariance by inverse with SVD-pseudoinverse fallback.
- 2.2 Diagnostic suite: φ, optimiser status, true parameter error, 95% interval half-width,
  empirical coverage, Jacobian condition number, covariance path, multistart spread.
- 2.3 **Three canonical designs** (identifiable / correlated / over-parameterised) — the
  controlled core.
- 2.4 **Seven domain adapters** — the applied breadth.
- 2.5 Noise ladder 0 / 1 / 3 / 10%; ≥200 replicates per cell; seed 20260727.
- 2.6 Pre-registration: protocol and code committed before test-case generation, as in M7.

**3. Results**
- 3.1 Recovery in the identifiable case (**the positive result**).
- 3.2 The dissociation: perfect fit, wrong parameters (**the headline**).
- 3.3 Uncertainty calibration: coverage vs nominal, by design and noise.
- 3.4 Rank-deficiency detection: sensitivity/specificity of the covariance flag.
- 3.5 Multistart and equifinality: distinct optima with indistinguishable φ.
- 3.6 Domain envelope across the seven adapters.
- 3.7 Active learning: does the recommended measurement actually help?

**4. Discussion** — what a modeller should report alongside a calibration; why "converged"
and "good fit" are not evidence; identifiability as a first-class output; relation to the
equivalence-class language of M5/M6 and the falsification framing of M3/M4.

**5. Conclusions** — with the negative results stated as prominently as the positive.

---

## 4. Figures

**Figure 1 — Conceptual: three ways an inverse problem can fail.**
Three-panel schematic of objective-function surfaces in 2-parameter space: (a) a single
sharp minimum (identifiable); (b) a flat ridge (correlated — every point on the ridge fits
equally well); (c) a plateau (over-parameterised). Overlay the true parameter point and
three optimiser trajectories from different starts. *Sets the whole argument up visually in
one glance; no data needed.*

**Figure 2 — The dissociation (HEADLINE FIGURE).**
Scatter of **final objective φ (x, log) vs true relative parameter error (y, log)**, one
point per run, coloured by design, marker by noise level. The identifiable design forms a
diagonal band (better fit → better parameters). The correlated and over-parameterised
designs form a **vertical stripe at φ → 0** spanning the full error range. One annotated
exemplar: "φ = 7.2×10⁻¹³, parameter error 40%".
*This single panel is the paper. It should be readable from across a room.*

**Figure 3 — Uncertainty calibration.**
(a) Coverage of nominal 95% intervals vs noise level, one line per design, with a dashed
0.95 reference and binomial CIs. (b) Calibration curve: nominal vs empirical coverage
across nominal levels 50–99%. (c) Interval half-width distributions (log) by design —
showing that the correlated design's intervals blow up yet still fail to cover.

**Figure 4 — Rank-deficiency detection.**
(a) Jacobian condition number distributions by design (log, violin). (b) Confusion matrix
of the covariance-path flag (`inverse` vs `svd_pseudoinverse`) against ground-truth
identifiability. (c) ROC of condition number as a continuous identifiability predictor,
with the covariance flag marked as a single operating point. *Reports the flag honestly as
a partial detector.*

**Figure 5 — Equifinality under multistart.**
For the correlated design: parameter-space scatter of optima from N random starts, with the
true value marked and the analytic ridge (a·b = const) drawn. Colour by φ to show that
every point on the ridge is an equally good fit. *Makes non-identifiability tangible rather
than abstract.*

**Figure 6 — Domain envelope across seven adapters.**
Heatmap: adapter (rows) × noise level (columns), cell colour = median relative parameter
error, cell annotation = 95% coverage. Adapters ordered by identifiability. Immediately
answers "which of my inverse problems can I trust?"

**Figure 7 — The identifiability threshold (SECOND HEADLINE).**
Twin y-axes against cond(FIM) on a log x-axis: (left) median relative parameter error,
(right) empirical 95% coverage. Error rises sharply through cond ≈ 10²–10³ while coverage
*rises to 1.00*, with interval half-width plotted as a shaded band to show why. Annotate the
threshold. *This panel makes the counter-intuitive coverage trap immediate.*

**Figure 8 — Optimal experimental design works (PRESCRIPTIVE FIGURE).**
(a) Candidate observations scored by D-, A- and E-optimality, ordered by E. (b) Verification:
median parameter error after actually adding each candidate — optimal, random, none, worst —
with the 60-fold gain and the "redundant is worse than nothing" bar annotated. (c) Coverage
per strategy against the 0.95 reference, showing that good design restores calibration.

**Figure 9 — Heuristic versus information-theoretic design.**
Head-to-head: the existing `active_learning` recommendation against E-optimal, random and
oracle. Report whatever it shows.

**Supplementary:** S1 convergence traces; S2 sensitivity to LM settings (ftol/xtol/max_nfev);
S3 robust-loss comparison (linear/huber/soft_l1/cauchy); S4 parallel-Jacobian
reproducibility; S5 full per-adapter parameter tables.

---

## 5. Tables

**Table 1 — Benchmark design contract.** Design | parameters | observations | expected rank |
identifiability by construction | what it tests | allowed claim | required guardrail.
*(Mirrors M4's evidence-ladder table, which reviewers responded well to.)*

**Table 2 — Recovery and calibration summary (MAIN RESULT).**
Design × noise → n, median relative error, mean φ, optimiser success rate, 95% coverage,
median interval half-width, pseudoinverse fraction. **This table already exists from the
pilot** (`results/m8_pilot_summary.csv`).

**Table 3 — Identifiability diagnostics.** Design | median condition number | covariance-path
distribution | detector sensitivity | specificity | AUC of condition number.

**Table 4 — Domain adapter envelope.** Adapter | n parameters | n observations | median
relative error at 3% noise | coverage | fraction rank-deficient | verdict
(identifiable / conditionally / not).

**Table 5 — Multistart equifinality.** Design | n starts | distinct optima | φ spread across
optima | parameter spread | ridge dimension.

**Table 6 — Identifiability threshold.** eps | cond(FIM) | lambda_min | median relative error |
95% coverage | median interval half-width | pseudoinverse fraction.
**Already produced** (`results/m8_conditioning_summary.csv`).

**Table 7 — Optimal design: selection and verification.** Candidate | D | A | E scores |
selected by | median relative error after adding | coverage | fold-change vs no measurement.
**Already produced** (`results/m8_oed_candidate_scores.csv`, `m8_oed_summary.csv`).

**Table 8 — Heuristic vs information-theoretic design.** Strategy | measurement chosen |
Δ error | Δ coverage | n improved | n degraded.

**Table 9 — Reproducibility.** Seeds, engine settings, versions, commit hash, runtime,
artifact checksums. *(Carry over the M7 pre-registration statement.)*

---

## 6. Abstract (draft — updated after locked real-adapter confirmation)

> Groundwater inverse modelling can produce a small residual without recovering the
> parameters that generated the observations. We evaluated this dissociation using a
> PEST-style Gauss–Levenberg–Marquardt engine, controlled canonical problems and production
> HydroSheaf transport and PHREEQC kinetic adapters. In 250 paired transport replicates, a
> noise-whitened log-parameter Fisher analysis selected the same 50 d front observation under
> parameter-specific, D-, A-, E- and balanced criteria. Relative to a three-point late-time
> design, the added observation reduced median dispersivity absolute log10 error from 0.2154
> to 0.0173 and decay error from 0.0182 to 0.0146. The anticipated split between
> dispersivity- and decay-targeted next observations was therefore not confirmed. In the
> kinetic adapter, one and four residence-time experiments both yielded rank-one information
> for calcite rate constant and reactive surface area because the rate law identified their
> product; an independent surface-area observation restored rank two. These controlled
> results show both the value and the limit of additional sampling: measurement placement can
> repair a poorly timed transport design, whereas structural product confounding requires
> independent parameter information.

---

## 7. What is done, and what must still be run

| item | status |
|---|---|
| GLM engine | exists, 77/77 unit tests pass |
| Three canonical designs | **pilot complete**, 40 replicates × 4 noise levels |
| Table 2 | **already produced** from the pilot |
| Figure 2 (headline) | plottable from existing pilot data **today** |
| Conditioning sweep + Table 6 | **complete** — threshold located at cond ≈ 10²–10³ |
| Figure 7 (threshold) | plottable from existing data **today** |
| OED selection + verification + Table 7 | **confirmatory complete** — common 50 d optimum; target split not supported |
| Figure 8 (OED) | **reviewed 2x2 confirmatory figure produced** |
| Figures 3, 5 | need nominal-level sweep and multistart runs (hours) |
| Figure 4 | condition numbers now logged; needs the detector ROC assembled |
| Figures 6, Table 4 | transport and kinetics complete; five-adapter envelope remains future work |
| Figure 9, Table 8 | needs the heuristic active-learning module wired to the same harness |
| Prospective confirmation | protocol/source hashes locked before fresh run; exploratory results explicitly excluded |

**Critical sequencing note.** Register the protocol *before* generating the final test
cases, exactly as M7 did at `d336e87`. That single practice is what makes M7 the most
defensible chapter in the thesis, and it costs nothing.

---

## 8. Risks and how the design already answers them

- *"You engineered the failure cases."* The correlated design is the textbook
  product-of-parameters case that occurs constantly in real reaction and transport models;
  the over-parameterised case is the normal condition in data-limited aquifers. Both are
  disclosed by construction in Table 1, and the identifiable control shows the engine works
  when the problem is well-posed.
- *"This is just a scipy wrapper."* The contribution is the diagnostic suite and the
  measured dissociation, not the optimiser. Say so in the introduction.
- *"Linearised uncertainty is known to be approximate."* Known qualitatively; the
  contribution is quantifying the miscalibration and showing the engine can flag it.
- *"Synthetic only."* Follow with the Ghana application under M6's claim discipline, or
  scope explicitly to controlled evaluation as M7 does.
