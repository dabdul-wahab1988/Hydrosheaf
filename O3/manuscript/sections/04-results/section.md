## 4. Results

### 4.1 Applying the common taxonomy

Figure 1 summarizes the three pipelines this comparison harmonizes.
Applying the taxonomy of Section 3.3 to every retained row of the age,
topology and reaction benchmarks (Table 1) yields 25 classified rows: 16
independent, one non-independent, one informed control, four negative
controls, one sensitivity diagnostic, three prior-informed, and one
conventional baseline. No row's classification contradicts the claim
already recorded for it in its own source benchmark's decisions or
manuscript-ready tables; the taxonomy adds a shared vocabulary across
benchmarks, not a re-interpretation of any single result.

![Figure 1. Three independently designed HydroSheaf benchmarks and the common evidence taxonomy this paper adds. The age benchmark compares nuclear-tracer/lumped-parameter age inference against the USGS national groundwater-age release; the topology benchmark compares reduced-order graph topology against MODPATH particle-tracking connectivity from public MODFLOW archives; the reaction benchmark compares sparse linear reaction inversion against a 240-scenario live-PHREEQC factorial synthetic benchmark. Each pipeline was locked independently, on its own timetable, before this comparison was constructed.](artifacts/figures/FIG-1.png)

**Table 1.** Common evidence taxonomy applied to every retained benchmark row of the age, topology and reaction benchmarks.

| Benchmark | Scenario | Row type | Independent | External reference | Headline metric | Value |
|---|---|---|:---:|---|---|---:|
| Topology | Spatial only | independent | yes | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Sink-aware baseline | informed control | no | MODPATH (Savage archive) | F1 | 0.552 |
| Topology | Head gradient | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Hodge pruned (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Proj. gradient (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Head depth | independent | yes | MODPATH (Savage archive) | F1 | 0.616 |
| Topology | Hydrostratigraphic | independent | yes | MODPATH (Savage archive) | F1 | 0.614 |
| Topology | Negative: random | negative control | no | MODPATH (Savage archive) | F1 | 0.006 |
| Topology | Negative: wrong direction | negative control | no | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Negative: misrouted sink | negative control | no | MODPATH (Savage archive) | F1 | 0.138 |
| Topology | Negative: shortcut (inapplicable) | negative control | no | MODPATH (Savage archive) | F1 | n/a |
| Topology | Sparse node | sensitivity diagnostic | no | MODPATH (Savage archive) | mean F1, 20 trials | 0.445 |
| Topology | Prior: override / merge / only | prior-informed | no | MODPATH (Savage archive) | edges in/out | 302 -> 302/329/174 |
| Age / residence time | TracerLPM strict parity | independent | yes | USGS TracerLPM Table 4 | within-factor-2 | 0.875 (n=329/1272) |
| Age / residence time | TracerLPM age fractions | non-independent | no | USGS reported age fractions | within-factor-2 | 0.917 (n=289/1272) |
| Age / residence time | Hydrosheaf model selection | independent | yes | USGS national age release | within-factor-2 | 0.673 (n=309/1272) |
| Age / residence time | Old-water He-4 branch | independent | yes | USGS national age release | within-factor-2 | 0.697 (n=66/1272) |
| Age / residence time | Graph regularisation, randomised control | negative control | yes | Randomised candidate graph | delta log10 RMSE vs single-node | +0.655 |
| Reaction | Bounded LS | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.435 |
| Reaction | Lasso | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.554 |
| Reaction | Elastic Net | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Thermo Elastic Net (legacy primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Hydrosheaf Guarded (primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.563 |
| Reaction | Hydrosheaf-Core (evidence-gated) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.561 |
| Reaction | Conventional PHREEQC inverse (strict+relaxed fallback) | conventional baseline | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.510 |

Note. The Hodge-pruned and projected-gradient topology rows are diagnostic
post-processing variants of the Head-gradient scenario immediately above
them, not independent replications; all three converge on the same
underlying inference and are numerically identical (F1 0.618) by
construction, not by coincidence.

### 4.2 Capture exceeds correctness under independent evaluation, in all three layers

Figure 2 and Table 2 report the paper's central result. For
topology, the independent Head-gradient scenario recovers 84.5% of the true
MODPATH reference edges (recall 0.845) while only 48.7% of the edges it
infers are true (precision 0.487); the three other independent scenarios
(head-depth, hydrostratigraphic, and two head-gradient diagnostic variants)
recover the same pattern within 0.004 F1 of one another (F1 0.614-0.618).
For age, 25.9% of the 1,272 candidate TracerLPM fits clear the reference-free
reportability guard and are reported at all, while the fits that are
reported agree with an independent held-out cross-validation reference
within a factor of 2 in 56.4% of cases (log10 R2 0.482). For reaction, the
Hydrosheaf Guarded model recovers 59.9% of the true active reaction phases
across the 240-scenario benchmark, pooled across 16 reactions (recall
0.599), while 63.9% of the phases it selects as active are true (macro-mean
precision, 1 minus a false-discovery rate of 0.361). In every one of the
three physically distinct inference problems, the point estimate of the
capture-type axis exceeds the point estimate of the matched correctness-type
axis.

That point-estimate ordering does not hold with equal statistical force in
all three layers. Adding 95% Wilson score intervals to each proportion
(Table 2) shows that topology's recall and precision intervals do not
overlap (0.784-0.891 versus 0.431-0.543, n=174 and n=302 respectively) and
neither do age's reportability and within-factor-2 intervals (0.235-0.283
versus 0.527-0.601, n=1,272 and n=675): in both layers the capture-exceeds-
correctness ordering is resolved at conventional confidence. Reaction is
different: its pooled recall interval (0.562-0.635, n=681) and macro
precision interval (0.598-0.681) overlap, so the small numerical gap between
0.599 and 0.639 for this layer is not distinguishable from sampling
variation at the 95% level, even though the point estimates run in the same
direction as the other two layers. The pattern that survives interval-aware
scrutiny is therefore narrower than the point estimates alone suggest:
recall/reportability exceeds precision/agreement with resolved confidence in
the topology and age layers, and in the same direction but not resolved in
the reaction layer.

The magnitude of the gap also differs across the two layers where it is
resolved. Topology shows the larger separation (0.358 between recall and
precision); age shows a smaller separation (0.305) on a pair of axes that,
as Section 3.4 notes, are not the same kind of rate as topology's or
reaction's. Among the lowest-residual quartile of reaction
fits -- the cases a residual-only read would judge most trustworthy -- 55.0%
still show phase F1 below 0.80, and the conventional PHREEQC inverse
baseline, run in parallel at strict 5% uncertainty, is feasible for only
45.8% of the 240 scenarios (99.6% once a relaxed 20% fallback is allowed),
returning a mean of 185.8 thermodynamically feasible models per scenario
where it is feasible at all.

Across the seven compared reaction-inversion methods, correctness and
convergence move together but not in lock-step with capture. Unregularized
bounded least squares converges least often (35.4% of fits) and shows the
highest extent error (0.091 mmol/L median) and the highest false-discovery
rate (0.697) despite fitting concentrations closely; the four regularized
sparse methods (LASSO, elastic net, thermodynamic elastic net, Hydrosheaf
Guarded) converge in 76.7-82.5% of fits and cut extent error roughly in half
(0.050-0.052 mmol/L); Hydrosheaf Guarded and Hydrosheaf-Core are the only two
methods with zero recorded thermodynamic-bound violations. The leakage-guarded
tracer-withholding design applied to age produces the same qualitative
result as the graph-regularization comparison: refitting every graph-node
age with its target tracer withheld shows no candidate graph improving on
baseline prediction of the withheld tracer for any of the three tracers with
enough retained fits to compare (tritium: baseline RMSE 20.818 TU against
20.841 for a depth-based graph and 21.567 for a randomized graph, 121 of 794
eligible rows reportable; sulphur hexafluoride: 2.844 pptv baseline against
2.929 randomized, 75 of 262 reportable; carbon-14: 26.830 pmC baseline
against 27.819 randomized, 169 of 1,103 reportable).

![Figure 2. Capture-type versus correctness-type metrics under independent, prior-free, uncalibrated evaluation, one panel per benchmark on its own native scale, with 95% Wilson score error bars. Topology: recall and precision of edge classification against the MODPATH reference (Head-gradient scenario). Age: reportability rate (fraction of fitted rows admitted by the reference-free guard) and within-factor-2 agreement on a held-out cross-validated comparison. Reaction: recall and precision of active reaction-phase detection, Hydrosheaf Guarded method, against the PHREEQC benchmark; note the overlapping intervals for this layer.](artifacts/figures/FIG-2.png)

**Table 2.** Headline independent/uncalibrated capture-type and correctness-type metrics, all three benchmarks.

| Benchmark | Capture-type metric | Value [95% CI] | Correctness-type metric | Value [95% CI] | CIs overlap? |
|---|---|---|---|---|:---:|
| Topology | Recall (true MODPATH edges recovered), n=174 | 0.845 [0.784, 0.891] | Precision (inferred edges that are true), n=302 | 0.487 [0.431, 0.543] | No |
| Age / residence time | Reportability rate, n=1272 (fitted rows admitted by the reference-free guard) | 0.259 [0.235, 0.283] | Within-factor-2 agreement, n=675, held-out cross-validated, uncalibrated | 0.564 [0.527, 0.601] | No |
| Reaction | Recall, pooled across 16 reactions, n=681 (true active phases recovered) | 0.599 [0.562, 0.635] | Precision, macro mean across 240 scenarios (1 minus false-discovery rate) | 0.639 [0.598, 0.681] | Yes |

Notes. Intervals are 95% Wilson score intervals computed from the underlying counts (`O3/analysis/python/derive_headline_metrics.py`); the reaction-layer precision interval is inverted from the reaction benchmark's own published false-discovery-rate 95% CI rather than recomputed from scratch. Topology: independent Head-gradient scenario against the USGS Savage MODPATH reference; the informed structural floor (receptor-set only, zero hydraulic information) reaches F1 0.552 and is not independent evidence. Age: capture and correctness are two different kinds of rate, not a literal recall/precision pair -- reportability is the fraction of attempted fits that produce an actionable claim at all (an abstention design, not a defect: see Discussion); within-factor-2 is the accuracy of the claims that are produced, scored on an independent held-out cross-validation split. Reaction: recall is pooled (micro-averaged) from `tableS6_reaction_confusion_matrices.csv` because no macro-averaged recall is published; the macro-averaged precision (0.639, from `table1_comparative_inverse_performance.csv`) is the reaction benchmark's own published headline correctness value. A pooled-precision cross-check (0.586, 95% CI [0.549, 0.622]) is reported in `manuscript/artifacts/data/evidence_discrepancies.csv` and differs from the macro value only because of the averaging convention, not a computational disagreement. **The reaction layer's capture and correctness intervals overlap; only the topology and age layers show a gap that is clearly resolved at the 95% level.**

### 4.3 The gap narrows only when independence is given up

Figure 3 and Table 3 compare each independent or screening
condition against its calibrated, prior-informed, or emulated counterpart.
Calibrating age estimates with a ridge model fit on the same held-out folds
used to score them raises log10 R2 from 0.482 to 0.962 and within-factor-2
agreement from 0.564 to between 0.769 and 0.993 depending on the reported
statistic; because the calibration is fit and scored on the same folds, this
is reported as emulation of the reference rather than independent agreement
(Section 2.1). The reaction layer's calibrated alternative is built
differently: the Mechanism Resolution Score is trained on the carbonate,
crystalline and evaporitic archetypes and evaluated, untouched, on the mixed
archetype, reaching 48.9% four-class accuracy against a 25.0% chance level
for four resolution classes -- a smaller absolute gain than age's calibrated
figure, but a genuine transfer test rather than a same-fold emulation.
Topology's counterpart is neither a calibration nor an emulation: entering
the 174 MODPATH reference edges as priors changes the output graph's edge
count (302 independently inferred edges become 302 under override, 329 under
merge, or 174 under only-prior mode) rather than raising an accuracy metric
comparable to the independent F1 of 0.618, and this mode is reported on its
own edge-count axis rather than forced onto the same scale.

![Figure 3. Independent/screening versus calibrated, emulated, or prior-informed evaluation. Left: age benchmark log10 R2, held-out uncalibrated versus calibrated emulation (fit on the same folds it is scored against). Centre: reaction benchmark Mechanism Resolution Score four-class accuracy, chance level versus a genuine held-out-archetype transfer test. Right: topology benchmark, independent F1 versus prior-informed output-graph edge counts (a different metric, plotted on a separate axis).](artifacts/figures/FIG-3.png)

**Table 3.** Screening/independent versus calibrated, emulated, or prior-informed metric, with non-independence stated per row.

| Benchmark | Independent / screening condition | Value | Calibrated / prior-informed / emulated condition | Value | Independence lost |
|---|---|---:|---|---:|---|
| Age / residence time | Held-out uncalibrated log10 R2 | 0.482 | Calibrated emulation log10 R2 | 0.962 | Ridge calibration fit on the same held-out folds it is scored against |
| Age / residence time | Held-out uncalibrated within-factor-2 | 0.564 | Calibrated emulation within-factor-2 | 0.769-0.993 | Same non-independence as above |
| Reaction | Chance level, 4-class MRS | 0.250 | Held-out archetype transfer, 4-class MRS accuracy | 0.489 | Trained on 3 archetypes, tested on 1 untouched archetype (a genuine transfer test, not an emulation) |
| Topology | Independent F1 (no MODPATH access) | 0.618 | Prior-informed output-graph edges (override / merge / only) | 302 / 329 / 174 | MODPATH edges entered the graph as prior information; not independent validation |

Notes. The two calibration exercises are not equally rigorous and are not presented as such: the age benchmark's calibrated value is an emulation fit and scored on the same data, which is why its improvement (0.482 to 0.962) is large; the reaction benchmark's Mechanism Resolution Score is tested on a fourth, untouched archetype it never trained on, a harder and more independent test, which is why its absolute accuracy (0.489) is far lower and still only doubles chance. The topology benchmark's prior-informed row is not a calibration exercise at all -- it is a different mode of use (MODPATH edges as priors, not as an accuracy target) and is reported on its own edge-count scale rather than forced onto the same axis.

### 4.4 Within-component detail

Figure 4 plots every retained row of each benchmark on its own
headline metric, coloured by taxonomy row type. In topology, all four
negative controls fall at or below F1 0.138, the informed structural
baseline reaches F1 0.552 using the reference receptor set alone with zero
hydraulic information, and the sparse-node sensitivity diagnostic falls to a
mean F1 of 0.445 across 10-100% node-density trials. In age, the
non-independent age-fraction scenario (which shares provenance with its own
reference) shows the highest within-factor-2 agreement of the four
design-matrix scenarios plotted (0.917), above every independent scenario
including strict parity (0.875); the randomized-graph negative control
increases log10 RMSE by 0.655 relative to single-node inference, and the
smallest candidate improvement from any tested graph-regularization
condition is a decrease of 0.001 log10 units, accompanied by a decline in
within-factor-2 agreement, which does not meet the benchmark's own joint
robust-improvement criterion. In reaction, the seven compared methods span
phase F1 from 0.435 (bounded least squares, unregularized) to 0.563
(Hydrosheaf Guarded), with the conventional PHREEQC inverse baseline (0.510)
falling within this range on accuracy alone despite its far lower feasibility
and higher model multiplicity (Section 4.2). A data-tier ablation that adds
optional diagnostic ions on top of the core panel shows the same evidence
class present across all three tiers: class F1 rises only from 0.606 (core,
major-ion panel only) to 0.610 (plus optional SiO2, Sr and water-isotope
diagnostics) to 0.614 (plus Br, dissolved oxygen, dissolved organic carbon,
and sulphate/nitrate isotopes), while the false-discovery rate falls only
from 0.383 to 0.364 to 0.361 across the same tiers -- additional analytical
information narrows, but does not close, the gap reported in Section 4.2.

![Figure 4. Within-component detail: every retained benchmark row for the age, topology and reaction benchmarks, coloured by evidence-taxonomy row type (Table 1), plotted on each benchmark's own headline metric. The three panels use three different metrics on visually similar axes and are not comparable to one another; only within-panel bar heights may be compared directly.](artifacts/figures/FIG-4.png)

### 4.5 Benchmark scale and computational cost

Figure 5 and Table 4 report the scale of the three benchmarks
on their own terms. The reaction benchmark is the largest by two measures:
240 PHREEQC scenarios, expanding to 21,600 factorial inverse fits once every
comparator method and noise/panel/archetype combination is counted. The age
benchmark spans 13 design-matrix scenarios totalling 16,536 candidate rows,
plus 375 additional rows from a leakage-guarded tracer-withholding
cross-validation exercise run separately for five tracers. The topology
benchmark is the smallest by scenario and fit count: 12 independent
graph-inference rows compared against 242 endpoint-derived reference edges
pooled across the two active MODFLOW archives. Only the reaction benchmark
records per-fit runtime across every comparator method, with a median of
20.8 ms per fit; neither the age nor the topology benchmark's result files
carry a comparable runtime field.

![Figure 5. Benchmark scale and computational cost. Top: scenario counts and total fit counts for the age, topology and reaction benchmarks, each on its own axis (values differ by roughly two orders of magnitude). Bottom: median per-fit runtime, recorded only for the reaction benchmark.](artifacts/figures/FIG-5.png)

**Table 4.** Benchmark scale: scenario counts, fit counts, external reference type, and recorded runtime.

| Benchmark | External reference | Scenarios | Fits | Cross-validation rows | Median runtime per fit | Runtime recorded |
|---|---|---:|---:|---:|---:|:---:|
| Age / residence time | USGS national groundwater-age release (TracerLPM Table 4) | 13 design-matrix scenarios | 16,536 (10 scenarios x 1,272 rows, plus 3 ablations) | 375 (leakage-guarded tracer-withholding CV, 5 tracers) | not recorded | no |
| Topology | MODPATH particle-tracking connectivity (3 public MODFLOW/MODPATH archives) | 12 independent graph-inference rows | 242 endpoint-derived reference edges (Savage + Great Miami; Long Island is a documented fallback stub) | not applicable | not recorded | no |
| Reaction | 240-scenario live-PHREEQC factorial synthetic benchmark | 240 PHREEQC scenarios | 21,600 factorial inverse fits (240 x 5 comparator methods x noise/panel/archetype grid) | not applicable | 20.8 ms (median across methods) | yes |

Notes. Counts are read directly from each benchmark's own manifests (`m3_design_matrix_summary.csv`, `independent_graph_vs_modpath.csv` plus the three `tier_*_archive_summary.csv` files, `analysis_summary.json`); no scenario was re-run to produce this table. The reaction benchmark is the only one to record per-fit runtime across every comparator method, which this paper reports as a genuine difference in benchmark design rather than normalising it away.

### 4.6 Field- and archive-transfer scope is uneven

Figure 6 and Table 5 report where each layer has been
tested against something beyond its primary synthetic or archive reference.
The reaction layer has the broadest field-transfer footprint: 160 candidate
edges from the Northern Ghana workbook (median evidence-lifted resolution
index, ELRI, 0.072), 85 from the Talensi District (median ELRI 0.147), and
49 from the Lower Anayari catchment (median ELRI 0.072). The topology layer's
transfer footprint is archive-to-archive rather than archive-to-field: the
Savage and Great Miami MODFLOW archives both return endpoint/pathline edge
F1 of 1.000 as a self-consistency check, and the Long Island archive remains
an inactive fallback stub with no rows in this comparison. The age layer has
no field- or archive-transfer benchmark beyond its primary comparison against
the public USGS national release; no dataset in Table 5 extends the
age layer's evaluation beyond that single reference.

Within the reaction layer's field-transfer component, a retrospective
next-best-measurement analysis ranks candidate additional ions and tracers by
their expected reduction in reaction-class ambiguity; fluoride returns the
highest measurement-value score of any candidate in the Northern Ghana
transfer analysis, and all 160 Northern Ghana wet-to-dry pairs are classified
as partially identifiable, rather than fully identifiable, by the calibration
transferred from the synthetic benchmark. Median Hydrosheaf-Core evidence and
total-dissolved-solids consistency scores for the Northern Ghana panel are
0.690 and 0.942 respectively.

![Figure 6. Field- and archive-transfer scope. Top: MODPATH archive endpoint/pathline edge F1 (Savage, Great Miami; Long Island is a documented fallback stub with no active rows). Bottom: median evidence-lifted resolution index (ELRI) for the three field chemistry datasets (Northern Ghana, Talensi, Lower Anayari). No field-transfer benchmark exists for the age layer.](artifacts/figures/FIG-6.png)

**Table 5.** Full dataset inventory: every synthetic and field dataset any of the three benchmarks consumes, with source, size, and role.

| Dataset | Type | N | Source | Role in this comparison | Access |
|---|---|---:|---|---|---|
| USGS national groundwater-age release (TracerLPM Table 4 reported configuration) | Published field tracer data | 1,272 candidate site-fits; tracers 3H, SF6, CFC-11, CFC-12, CFC-113, 14C | [@jurgens2012tracerlpm] | External reference for the age/residence-time benchmark | Public (USGS) |
| USGS Savage Municipal Water-Supply Well archive | Published synthetic MODFLOW-2005/MODPATH 5 archive | 3,000 particles; 174 endpoint reference edges | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/F7J102FK | Primary external reference for the topology benchmark | Public (USGS) |
| Great Miami River Basin archive | Published synthetic MODFLOW 6/MODPATH 7 archive | 154 particles; 68 endpoint reference edges | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/P9X4C9R6 | Secondary external reference for the topology benchmark (self-consistency check) | Public (USGS) |
| Long Island Nitrogen archive | Published synthetic MODFLOW 6/MODPATH 7 archive | Fallback stub; no active rows in this comparison | U.S. Geological Survey ScienceBase, https://doi.org/10.5066/P97VFXZ4 | Tertiary archive for the topology benchmark, documented as inactive | Public (USGS) |
| 240-scenario live-PHREEQC factorial benchmark | Synthetic, HydroSheaf-generated | 240 scenarios; 21,600 factorial inverse fits; 4 hydrochemical archetypes (carbonate, crystalline, evaporitic, mixed) | Generated with [@parkhurst2013phreeqc] (PHREEQC 3.7.3/3.8.6) | Primary external reference for the reaction benchmark | Reproducible from the reaction benchmark's own repository folder |
| Northern Ghana workbook | Published field chemistry, reconstructed seasonal split | 160 boreholes, 4 regions; Dry panel is the primary measured unit | Author-compiled regional dataset; see Methods for the reconstruction disclosure | Field-transfer chemistry audit for the reaction benchmark | Project-held; see data-availability statement |
| Talensi District | Published field chemistry | 63 samples | [@chegbeleh2020talensi] | Field-transfer chemistry audit for the reaction benchmark (screening-level) | Published source |
| Lower Anayari catchment | Published field chemistry | 41 samples | [@zakaria2021anayari] | Field-transfer chemistry audit for the reaction benchmark (screening-level) | Published source |

Regional hydrogeological context for the three Ghanaian sites is drawn from [@banoengyakubo2011ghana] and [@macdonald2012africa]; neither is a source for the Northern Ghana workbook itself.
