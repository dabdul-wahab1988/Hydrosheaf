# M3 & M4 Manuscript Accuracy, Figure and Reproducibility Audit (2026-07-27)

> **STATUS: corrections applied.** See **PART Z — WHAT WAS CHANGED** at the end for the
> full list of code, figure, table and manuscript edits, and for two further defects
> (Z-1, Z-2) found while applying them. M4 is complete. M3 required a full benchmark
> re-run because the tritium input function was found to be broken (B-6).


Scope: `M3/M3_geochemistry/Manucript_3.docx` (+ Supplementary_Information_M3.docx) and
`M4/Re_ Manuscript 4/Manuscript-Advances of Water Resources.docx` (+ Supplementary Informationsv-AWR.docx),
checked line-by-line against the committed result CSVs, the generating scripts, and the
manuscript-ready figures.

**Headline:** every *reported* number in both manuscripts reproduces exactly from the
committed artifacts. The problems are not arithmetic — they are **four unreported results
that contradict or materially weaken the two papers' headline claims**, plus a set of
figure defects. Both headline claims survive, but only in weakened form.

---

## PART A — VERIFIED CORRECT

All recomputed from the CSVs; no discrepancies found.

**M3** — parity full-available (n=1,249; MdAE 0.167; RMSE 0.740; R² 0.681; wf2 0.614;
wf10 0.859); strict parity (n=1,254; 0.221/0.841/0.590/0.573); Hydrosheaf selection
(n=493; 0.487/0.980); high-N common support N=655 (0.157/0.761/0.709/0.626 vs
0.217/0.865/0.624/0.582); design common support N=272; 1,272 harmonised rows; graph
Table 4 (weak parameter smoothing ΔRMSE +0.003, violations 62→46; randomised strong
+0.376); withheld-tracer RMSEs (³H 12.42→7.97 = −35.87%; SF₆ 4.435→4.412;
¹⁴C 25.476/26.464/30.095/39.639); age-class S2 (modern MdAE 0.47, wf2 0.38–0.40 vs
0.67 old / 0.82 very old); old-water S3 (71/67/88/768/278).

**M4** — head-gradient 147 TP / 155 FP / 27 FN → P 0.4868, R 0.8448, F1 0.6176; spatial-only
306 FP / 174 FN / F1 0.000; random F1 0.00575; head-depth 0.6164 (156 FP);
hydrostratigraphic 0.6138 (158 FP); sparsity 0.149 / 0.251 / 0.539 / 0.617 / 0.618 with
4/20 and 9/20 successful trials; candidate universe 153×152 = 23,256 (TP+FP+FN+TN closes
exactly); Tier 1 Savage 3,000 / 410,769 / 3,000 / 174 / 174; Tier 2 Great Miami
154 / 1,494 / 154 / 68 / 68.

**The provenance claim in M3 §2.3 is real.** All 1,272 rows carry
`input_history_source = wiser_north_america_nearest`. The model forcing genuinely used
WISER. (But see A-2 — the *figure* does not.)

---

## PART B — CRITICAL: UNREPORTED RESULTS THAT UNDERCUT THE HEADLINE CLAIMS

### B-1 (M3, CRITICAL). The randomised negative control also improves withheld-³H RMSE — and the paper's own acceptance rule says that disqualifies the headline result.

M3's single "genuinely independent" positive result is the depth-constrained graph's
35.9% withheld-tritium RMSE reduction (abstract, §3.4, §3.7, §4.3, Conclusions).

Recomputed from `results/m3_cv_benchmark_3H.csv` (n=794):

| scenario | RMSE (TU) | vs single |
|---|---|---|
| single-well baseline | 12.425 | — |
| hydraulic-proxy graph | 12.390 | −0.3% |
| **depth-constrained graph** | **7.967** | **−35.9%** ✅ reported |
| **randomised negative control** | **9.979** | **−19.7%** ❌ **never reported** |

§2.9 pre-specifies: *"A graph-regularised result was interpreted as meaningfully improved
only where RMSE was lower than the single-node baseline, within-factor agreement did not
deteriorate, and **randomised negative-control graphs did not show comparable
improvement**."* A 19.7% improvement from a random graph is comparable. By the
manuscript's own stated rule the ³H result does not pass.

The internal QA (`docs/m3_results_summary.md`: *"Randomized negative-control rows with
apparent improvement: 0"*) is scoped **only to the Table-4 log-age RMSE benchmark** and
never audited the tracer-withholding CV. That is the gap.

**Figure 6a already plots this** — the randomised ³H bar sits visibly below the
single-well baseline — so a reviewer sees the contradiction in the figure while the text
and caption are silent.

### B-2 (M3, CRITICAL). The 35.9% gain is an artifact of ~5 pathological rows, and median error actually gets *worse*.

Same file, decomposed:

| metric | single | depth-constrained | randomised |
|---|---|---|---|
| RMSE | 12.425 | **7.967** | 9.979 |
| MAE | 2.531 | 1.999 | 2.467 |
| **Median abs error** | **0.395** | **0.513 (worse)** | 0.747 (worse) |
| P90 abs error | 3.807 | **4.125 (worse)** | 5.043 (worse) |
| RMSE after dropping 5 worst baseline rows | 4.959 | 4.728 (−4.7%) | **3.370 (best)** |

Only 179 of 794 rows change at all. The reduction comes almost entirely from a handful of
rows where the single-well LPM returned τ = 0 yr → predicted ³H = 154 TU against observed
0.20–1.30 TU (e.g. `METXPAS1-37`, `-18`, `-13`, `-38`). The depth graph ages three of
those, collapsing the prediction to ~0. Remove five rows and the effect is 4.7%, and the
**randomised control becomes the best method**.

Mechanistically, every graph "improvement" here is just *aging the water*: median τ goes
207 yr (single) → 621 (depth) → 10,112 (random). Aging pulls predicted ³H down toward the
mostly-low observations, so RMSE falls while typical-case accuracy degrades. The
randomised graph ages most aggressively — which is exactly why it "wins" after outlier
removal.

**Two things needed:**
1. Report the randomised-control ³H result and the median/MAE alongside RMSE.
2. Investigate the τ = 0 → 154 TU baseline fits. A piston-flow fit of zero age at a
   ~2010 sampling year should predict modern (~5–10 TU), not the bomb-peak value. This
   looks like an input-history lookup returning a clamped/extreme value, and it is
   driving the paper's headline number.

**Suggested reframing (defensible, and arguably a better paper):** the honest finding is
that *withheld-³H RMSE is dominated by a small number of degenerate young-water LPM fits,
and any age-inflating prior — physical or random — suppresses them.* That is a real
diagnostic result about non-identifiability and fits M3's falsification framing better
than a fragile 35.9%.

### B-3 (M4, CRITICAL). The MODPATH reference graph has only **three** distinct target cells, and 301 of 302 inferred edges already point at them.

`results/edge_classification.csv`: the 174 reference edges run from **150 unique source
cells into 3 unique target cells** (`cell_15098`, `cell_16745`, `cell_21316` — the supply
wells). It is a fan-in/star graph, not a general connectivity network. Fig. 4b already
draws exactly three capture envelopes; the manuscript never states the cardinality.

Consequences, all computable from the committed file:

- **301 of the 302 head-gradient inferred edges terminate at one of those 3 sinks.** The
  elevation-as-head downhill rule is, in this domain, largely rediscovering "everything
  drains to the three lowest cells."
- **A trivial sink-aware baseline with no hydraulics at all** — connect every node to each
  of the 3 sinks — scores **TP 174, FP 282, FN 0, P 0.382, R 1.000, F1 0.552**, versus the
  reported F1 0.618. The entire skill attributable to hydraulic directionality over a
  structurally-informed baseline is **ΔF1 ≈ 0.066**, not the 0.618-vs-0.000 gap the paper
  presents.

Every control in the paper (spatial-only 0.000, random 0.006, wrong-direction 0.000) is
*uninformed* — none knows the sink structure — so all of them are guaranteed to fail and
none of them tests the claim that actually matters.

**Required:** disclose the 150→3 reference structure in §2.2 and Fig. 1, add the
sink-aware baseline as the primary control, and restate the abstract's claim
("gave the strongest independent topology recovery") relative to *that* baseline. Also
note that the reference edges are long-range source→sink links (median 9,598 m) while
spatial-only edges are short (median 2,235 m) — "proximity fails" is a statement about
length scale, not about proximity as a hypothesis.

### B-4 (M4, HIGH). The shortcut control generated zero edges; reporting it as a control that "failed" is not supportable.

`independent_graph_vs_modpath.csv` → `negative_shortcut`: `n_inferred_edges = 0.0`.

`phase2b_independent_validation.py:676-687` builds shortcuts as u→v→w two-hop paths. Since
no node is both a source and a target (verified: |sources ∩ targets| = 0), **no two-hop
path exists** and the shortcut list is empty by construction.

Affected text:
- **Abstract**: "Spatial-only, wrong-direction and shortcut controls failed, confirming
  that topology recovery was not produced by proximity or chance."
- **§3.2**: "wrong-direction and shortcut controls produced F1 = 0.000 … visually
  plausible shortcuts were insufficient…"
- **Conclusion**: same claim repeated.
- **Table 2** row "Negative: shortcut … Shortcut-only topology failed to recover the
  reference graph" — reports 0/0/0 without disclosing n_inferred = 0.
- **Supp. Fig. S2 caption (j)**: "shortcut controls failed as expected, producing either
  unsupported edges, missed reference edges or both" — it produced neither.

A control that emitted no edges is not a failed control; it is an inapplicable one. Either
redefine the shortcut control so it is well-posed on a fan-in reference (e.g. skip the
nearest intermediate node in the *inferred* graph), or state plainly that the two-hop
shortcut control is undefined for this reference topology and drop it from the claim list.

### B-5 (M4, HIGH). The Hodge-pruned and projected-gradient scenarios prune nothing — `probability_threshold = 0.0`.

`head_gradient`, `head_gradient_bayesian_hodge` and `real_head_projected_gradient` report
**bit-identical** 147/155/27 and 302 inferred edges. §3.2 presents this as a result ("The
Hodge-pruned and projected-gradient cases reproduced the same aggregate metrics").

It is not a result — it is a configuration consequence. Both scenarios call
`infer_topology_map_edges(..., probability_threshold=0.0, ...)`
(`phase2b_independent_validation.py:196` and `:310`), so every candidate edge is retained
regardless of posterior probability. The MCMC runs and is discarded.

The results file's own `allowed_claim` for that row — *"Bayesian Hodge pruning suppresses
topologically inconsistent downhill shortcuts"* — is therefore unsupported by any reported
number.

**This is also the paper's best unexploited opportunity.** The posteriors are already
computed and stored: `posterior_n_edges_mean` = **171.9** (Hodge) and **63.8**
(projected-gradient), against 302 retained. The projected-gradient posterior wants ~64
edges — i.e. it is actively trying to solve the paper's central weakness (FP 155 > TP 147).
Re-running at a non-zero threshold (say 0.5) and reporting the precision/recall/F1
operating curve would turn a null diagnostic into a genuine contribution, and would
directly answer §4.3's own "higher precision would be needed" caveat.

---

## PART C — FIGURE DEFECTS

### C-1 (M3, HIGH) Figure 1a is the wrong curve, and the curve itself is buggy.

`make_publication_figures.py:229` builds panel (a) from
`build_default_tritium_input()` — the hard-coded stylised curve in
`hydrosheaf/nuclear/input_history.py:75-123` — **not** the WISER histories the benchmark
actually used. The caption states the panel "uses GNIP/WISER precipitation isotope records
as the source for ³H forcing (IAEA/WMO, 2026)." That is inaccurate on two counts: the
panel is synthetic, and it is not the forcing used in the analysis.

The generator is also internally inconsistent (`input_history.py:103-108`): the 1952–1963
rise is `5 + 100·exp((y−1952)/2)`, which reaches **≈19,000 TU at 1962.5**, then
*discontinuously drops* to the nominal 3,000 TU "peak" at 1963. The plotted maximum
(~2×10⁴ TU) is therefore ~2–3× any real Northern Hemisphere GNIP record and sits a year
*before* the arrow labelled "Peak nuclear testing ~1963."

**Fix:** regenerate panel (a) from an actual WISER station history (the loader is already
in the pipeline), or from a published GNIP composite; and repair or delete the synthetic
generator so it cannot be reused.

### C-2 (M3, MEDIUM) Figure 6.

- Panel (a) puts ³H (TU) and SF₆ (pptv) bars on two y-axes with different baselines — the
  in-figure note concedes they "should be compared only within each tracer," which is an
  admission the encoding is wrong. Split into two panels, or plot % change from baseline.
- The −35.9% annotation appears on the depth-constrained ³H bar with no counterpart on the
  randomised bar (−19.7%) — see B-1. As drawn, the figure emphasises the supporting result
  and leaves the contradicting one unlabelled.
- Panel (b) annotates the ¹⁴C randomised bar "+55.6% (random noise)"; verified (39.639 vs
  25.476). Correct — but it makes the asymmetry with panel (a) more conspicuous.
- Empty bottom-right quadrant in a 2×2 grid.

### C-3 (M4, HIGH) Figure 4a plots FP = 0 and FN = 0 as bars of height 1 on a log axis.

Tier 1 Savage projection is TP 174, FP 0, FN 0. On a log scale zero is unplottable, and
the renderer has drawn unit-height bars — the figure reads "FP = 1, FN = 1" and
understates a perfect projection result. Use a linear axis, a broken axis, or explicit
"0" labels.

Also in Fig. 4: panel (a) mixes particles, pathline points and edge counts on one axis
(not a meaningful comparison); panel (c) is labelled "MODPATH-derived Hydrosheaf edge
weight" vs "Endpoint particle time median" but reports rho = 0.963 / tau = 0.923, which
in `external_modpath_archive_summary.csv` are the *harmonised travel-time* rank
correlations — check the axis labels match the statistic. The caption should also state
that this correlation is MODPATH-vs-MODPATH and therefore circular by construction.

### C-4 (M4, MEDIUM) Figure 2.

- Panel (a) is fully redundant with panel (b) (F1 appears in both) and shows five visually
  identical bars (0.618/0.618/0.618/0.616/0.614) — three of which are identical only
  because of B-5.
- Panel (a) legend carries an "L0: Controls (incl. spatial-only)" swatch for bars that are
  all zero and therefore invisible; the dashed vertical line at F1 = 0.5 is unexplained in
  the caption.
- Panel (c): the x-axis title "Node fraction retained" **collides with** the rotated
  "n = successful/planned trials" annotations — unreadable as submitted.
- Panel (c): 95% CI whiskers at the 10% and 25% fractions extend above 1.0 for recall and
  below 0 for precision, i.e. outside the metrics' admissible range. Clip, or switch to a
  bootstrap/percentile interval.
- Four of nine x-categories in panel (b) and (d) carry no bars at all (Shortcut has none
  anywhere — see B-4).

### C-5 (M4, MEDIUM) Supp. Table S5 publishes an empty Tier 3 row.

`Supp_TableS5_External_MODPATH_Archive_Validation.csv` and the supplement (line 369/377)
carry a **Tier 3 Long Island** row that is entirely zeros/dashes with the note
*"Fallback scalability stub only; no active validation result."* — while still citing a
live DOI (10.5066/P97VFXZ4) and a full reference (Jahn & Walter, 2025). §3.4 of the
manuscript discusses Tier 1 and Tier 2 only and never mentions Tier 3.

Publishing a cited-but-unprocessed archive as a table row invites "did you actually run
this?" Either process it or remove the row and the reference.

---

## PART D — REPRODUCIBILITY DEFECTS (regeneration does not reproduce the submitted artifacts)

### D-1 (M3) `make_publication_tables.py` emits the wrong subset label.

`scripts/make_publication_tables.py:226, 271, 299, 304, 308` hard-code
`"Design Common Support (N=43)"`. The generated
`Manuscript_Table5_Mode_Comparison.csv` carries `N=43` on every row while the N column
reads 272. The docx was hand-corrected to `(N=272)`, so **re-running the pipeline
regenerates a table that contradicts the submitted manuscript.** Fix the string in the
script.

### D-2 (M3) Supplementary table numbering is off by one between repo and manuscript.

| repo file | manuscript / supplement |
|---|---|
| `Manuscript_Supp_TableS1_Age_Class_Performance.csv` | Supplementary **Table S2** |
| `Manuscript_Supp_TableS2_Old_Groundwater_Diagnostics.csv` | Supplementary **Table S3** |
| `Manuscript_Supp_TableS3_Gas_Correction_Audit.csv` | Supplementary **Table S4** |

(The supplement's S1 is the reproducibility inventory, which has no generated CSV.) A
reader following the repo lands on the wrong table. Renumber the outputs or add an
explicit mapping to `MANIFEST.md`.

### D-3 (M3) Undisclosed sample sizes in the CFC withholding results.

Supplementary Figure S3 reports CFC-11 and CFC-12 withholding RMSE without stating n.
Actual: **CFC-11 n = 28, CFC-12 n = 16**. At n = 16, the reported CFC-12 degradation
(140.1 → 184.6) is not interpretable without a confidence interval. Add n to the caption
and to §3.6, or demote to a qualitative statement.

### D-4 (M3) §2.11 promises archival material that is not in the repo.

*"A release tag, commit hash, dependency lockfile or environment file, and clean-checkout
execution order should be archived with the journal submission."* This is written as a
recommendation, not as something done. For a paper whose stated contribution is
reproducibility, pin it: tag the release, record the commit hash, and commit the lockfile
(`uv.lock` exists at repo root — reference it).

---

## PART E — TEXT DEFECTS CONFIRMED IN THE DOCX (not conversion artifacts)

Verified by reconstructing paragraph text directly from `word/document.xml` including
`<m:t>` math runs — these characters are genuinely absent from the file.

**M3:**
1. §2.7: *"no uphill movement under the elevation proxy ."* — the inline expression is
   missing and its math object is orphaned at the **end of the paragraph**
   (`elevation = −depth ≤ 25`). Reads as a dangling formula after the last sentence.
2. §2.7: *"the depth-constrained subset with edge distance  km"* — **the numeric
   threshold is missing.** A construction rule the reader cannot reproduce.
3. §2.10: *"with  and ."* — the factor-of-2 and factor-of-10 values are missing.
4. §2.10: *"Log-scale , age-class performance…"* — the R² symbol is missing.
5. §2.6 / eq. (1): *"where is sampling year, is the atmospheric or initial input function,
   and for radioactive tracers. For stable atmospheric tracers, ."* — every symbol in the
   where-clause has dropped out.
6. §4.3: **"₅⁵Kr"** should be **⁸⁵Kr**.
7. Stray "." on its own line between §5 Conclusions and Author contributions.
8. §2.3 cites "(Jurgens et al., 2022)" without the a/b suffix used elsewhere — ambiguous
   between the two 2022 references.

**M4:** no equivalent dropouts found. But note the mismatch between
*"Data and Code Availability"* (§2.6, public GitHub + USGS DOIs — correct) and
*"Availability of data and material"* under DECLARATIONS (*"available from the
corresponding author upon reasonable request"*). The second contradicts the first and
undercuts the paper's reproducibility claim; delete it or point it at the repository.

---

## PRIORITISED ACTION LIST

| # | Paper | Sev | Action |
|---|---|---|---|
| 1 | M3 | **Critical** | Report the randomised-control withheld-³H result (−19.7%) and apply §2.9's own acceptance rule to the headline claim (B-1) |
| 2 | M3 | **Critical** | Disclose that the 35.9% is tail-driven (median error worsens; 4.7% after removing 5 rows); investigate the τ=0 → 154 TU baseline fits (B-2) |
| 3 | M4 | **Critical** | Disclose the 150-source → 3-sink reference structure; add the sink-aware baseline (F1 0.552) as the primary control; restate the headline relative to it (B-3) |
| 4 | M4 | High | Fix or withdraw the shortcut-control claim in abstract, §3.2, Table 2, Conclusion, Supp. Fig. S2 caption (B-4) |
| 5 | M4 | High | Set a non-zero posterior threshold; report the projected-gradient operating curve (posterior wants ~64 edges vs 302 retained) — turns a null into a contribution (B-5) |
| 6 | M3 | High | Regenerate Fig. 1a from real WISER/GNIP data; fix the 19,000-TU spike bug in `build_default_tritium_input` (C-1) |
| 7 | M4 | High | Fix Fig. 4a — FP=0/FN=0 currently render as height-1 bars on a log axis (C-3) |
| 8 | M3 | Med | Fix `make_publication_tables.py` "N=43" label; renumber supplementary tables (D-1, D-2) |
| 9 | M3 | Med | Rebuild Fig. 6 (split dual axis; annotate the randomised bar); disclose CFC n=28/16 (C-2, D-3) |
| 10 | M4 | Med | Rebuild Fig. 2 (drop redundant panel a; fix axis-label collision and out-of-range CIs); resolve the Tier 3 stub row (C-4, C-5) |
| 11 | M3 | Med | Restore the 8 missing symbols/values in the docx; fix ₅⁵Kr → ⁸⁵Kr (E) |
| 12 | M4 | Low | Remove the contradictory "available … upon reasonable request" declaration (E) |
| 13 | both | Low | Pin release tag + commit hash + lockfile before submission (D-4) |

---

*Recomputed against: `M3/m3_age_benchmark/results/*.csv`,
`M4/m4_topology_benchmark/results/*.csv`, the generating scripts, and the
`figures/Manuscript_Ready/` PNGs.*

---

# PART Z — WHAT WAS CHANGED

## Z-0. Two further defects found while applying the fixes

### B-6 (M3, CRITICAL — root cause of B-1 and B-2). The tritium input function was fabricated for every recharge year after each station's record ends.

`InputHistory.interpolate` (`hydrosheaf/nuclear/input_history.py:53`) delegates to
`numpy.interp`, which **clamps out-of-range targets to the nearest endpoint value**.
The M3 benchmark forces ³H with the nearest WISER North America station history. Those
station records end between **1965 and 1999** (quartiles 1969 / 1978 / 1999), while the
benchmark samples were collected 2004–2020.

Measured across the benchmark: **1,269 of 1,272 samples (99.8%)** have a sampling year
after their nearest station's record end, by a **median of 26 years** (p90 46 yr, max
53 yr). For all of those, the tritium input over the most recent decades of recharge was
a flat clamp at whatever value that station last reported — e.g. **245 TU** for the
Mesilla/Rio Grande sites, whose nearest station stops in 1965.

This is the mechanism behind the reported headline:

- τ = 0 fits predicted **154 TU** against observed ~1 TU (the `METXPAS1` rows in B-2).
- Modern water (≤50 yr) had the worst error in the benchmark (§3.5 attributed this to
  "young-tracer discordance"; the actual cause is the input function).
- **Any** age-inflating prior — depth-constrained *or* randomised — improved withheld-³H
  RMSE by moving predictions off the clamped bomb-era plateau. That is why the negative
  control "worked" (B-1) and why the gain was tail-driven (B-2).

**Fixed** by `extend_history_to_present()` in `hydrosheaf/nuclear/tracer_inputs.py`: the
station record is continued past its last observation along the regional post-bomb
decline curve, rescaled to join the record continuously at the splice year, converging to
modern background rather than freezing a bomb-era value. Verified at the four worst
sites — modern (2010) forcing changes from a clamped 14–245 TU to a physical 1.1–12.2 TU.

The full M3 benchmark was re-run under the corrected forcing (author decision).

### Z-1 (both, HIGH). Result artifacts were untracked, and M4's committed posterior columns had silently drifted.

`.gitignore:163` excluded `**/results/` (and `*.csv` at :172), so **none** of the M3/M4
result files backing the tables and figures were version-controlled. Consequence, found
by re-running: `posterior_n_edges_mean` for the M4 Hodge scenario reproduced as **183.03**
today versus the committed **171.86**, and the projected-gradient value as **151.43**
versus **63.77** — because `hydrosheaf/inference/topology_posterior.py` was modified on
2026-07-19 for M7 work, *after* the M4 results were produced. The reported F1 metrics were
unaffected only because `probability_threshold = 0.0` made them independent of the
posterior (B-5), which is why the drift went unnoticed.

**Fixed**: `.gitignore` now re-includes `M3/m3_age_benchmark/{results,tables}` and
`M4/m4_topology_benchmark/{results,tables}` (bulky `results/public_archives/` still
excluded). 157 artifacts are now tracked.

### Z-2 (M3, LOW). A fabricated threshold was introduced and then removed.

While restoring the missing §2.7 symbols I initially wrote "elevation drop of at most
25 m". No such rule exists: `hydrosheaf/graph/build.py:487` only rejects candidates whose
elevation proxy is `>=` the source elevation, and the "25" in the corrupted source text
belongs to the **25 km** distance filter (`run_m3_real_usgs_graph_benchmark.py:232`).
Corrected in the manuscript before saving. Flagged here because it is the kind of error
this exercise exists to prevent.

## Z-3. M4 — code

| File | Change |
|---|---|
| `phase2b_independent_validation.py` | Added `sink_aware_baseline` scenario (informed structural floor) and `negative_misrouted_sink` control; made the two-hop shortcut control record `control_applicable=False` when it emits zero edges, with an explicit "not evidence of rejection" claim string; added `posterior_operating_curve()` and wrote `posterior_operating_curve.csv`; made shortcut construction deterministic (`sorted(set(...))`) and its reference lookup a set; corrected the module docstring, which claimed Hodge pruning "modestly improves precision" |
| `make_publication_tables.py` | Table 2 now includes the two new rows; added `EVIDENCE_LADDER` entries for them; shortcut row now reads INAPPLICABLE with the reason; Supp Table S5 gained a `processing_status` column so the unprocessed Tier 3 stub cannot read as a result |
| `make_publication_figures.py` | Fig 4a no longer clamps zeros to 1 on a log axis (FP = 0 / FN = 0 now render as annotated zeros); Fig 2a replaced (was a redundant duplicate of 2b) with the posterior operating curve; Fig 2b gained the sink-aware floor line and an n/a marker on the empty control; Fig 2c CIs clipped to [0,1] and the x-label/annotation collision fixed |

**Verification**: every pre-existing scenario reproduced bit-identically after the
changes (`tp/fp/fn/tn/precision/recall/f1/n_inferred_edges` unchanged for all nine
original rows). New results: sink-aware baseline **P 0.382 / R 1.000 / F1 0.552**;
misrouted-sink control **F1 0.138**; posterior operating curve shows precision never
exceeds **0.485** at any threshold while recall falls from 0.845 to <0.02.

## Z-4. M4 — manuscript

Backups: `Manuscript-AWR_BACKUP_preAccuracyFix.docx`,
`Supplementary-AWR_BACKUP_preAccuracyFix.docx`.

1. **Abstract** — shortcut claim withdrawn; fan-in structure and the 0.552 sink-aware
   baseline disclosed; the head-gradient advantage restated as ≈0.066 F1 over an informed
   baseline rather than 0.618 over zero.
2. **§2.2** — the 174 reference transitions are now stated to run from 150 source cells
   into three receptor cells, with a note that this determines benchmark difficulty.
3. **§3.2** — shortcut control reported as inapplicable with the structural reason;
   misrouted-sink substitute (F1 0.138) reported; the identical Hodge/projected-gradient
   metrics attributed to `probability_threshold = 0.0` rather than presented as a result;
   operating-curve outcome reported; sink-aware comparison added, including that 301 of
   302 inferred edges already terminate at a reference receptor.
4. **§3.4** — Tier 3 Long Island disclosed as a non-ingested stub contributing no
   evidence.
5. **Conclusion** — control claims corrected; the sink-aware narrowing and the negative
   posterior-thresholding result stated.
6. **Declarations** — "available from the corresponding author upon reasonable request"
   replaced (it contradicted the Data and Code Availability section).
7. **Supplement, Fig. S2 caption (j)** — shortcut control no longer described as having
   "failed as expected".

## Z-5. M3 — code and manuscript

| File | Change |
|---|---|
| `hydrosheaf/nuclear/tracer_inputs.py` | **`extend_history_to_present()`** — the B-6 fix |
| `hydrosheaf/nuclear/input_history.py` | Fixed the synthetic curve: the rise `5 + 100·exp((y−1952)/2)` reached ≈19,000 TU at 1962.5 then dropped discontinuously to the nominal 3,000 TU peak. Now a continuous log-linear rise to exactly 3,000 TU at 1963, and a decay branch continuous with the peak |
| `make_publication_figures.py` | Fig 1a now plots the history the benchmark actually used (nearest-WISER via the same code path) instead of the synthetic curve the caption misattributed to GNIP; annotation matching made robust to irregular station sampling; Fig 6a converted from an uninterpretable dual-axis pair to percentage change from baseline, annotating **both** the depth-constrained and the randomised control; new Fig 6d fills the previously blank quadrant with the RMSE-vs-median-absolute-error decomposition that exposes the tail-driven effect (B-2) |
| `make_publication_tables.py` | `N=43` → `N=272` (6 occurrences); supplementary table outputs renumbered S1/S2/S3 → S2/S3/S4 to match the published supplement |

Manuscript backup: `Manucript_3_BACKUP_preAccuracyFix.docx`. Restored the dropped
symbols and values in §2.7 (elevation-proxy rule; the missing **25 km** distance
threshold), §2.10 (`f = 2` and `f = 10`; the R² symbol, twice), fixed **₅⁵Kr → ⁸⁵Kr** in
§4.3, and disambiguated the `(Jurgens et al., 2022)` citation to `2022b`.

**M3 results text and figures are regenerated from the corrected forcing; see the
re-run summary appended below.**

## Z-6. M3 — effect of the corrected forcing on the reported results

Full design matrix re-run: 14 scenarios × 1,272 rows = 17,808 fits, `--full --age-steps 90`
(matching the age-grid resolution recorded in the original parity manifests).

### Scenario metrics

| Scenario | MdAE before → after | RMSE before → after | R² before → after | wf2 before → after |
|---|---|---|---|---|
| `tracerlpm_parity_agefractions` | 0.169 → **0.166** | 0.736 → **0.799** | 0.684 → **0.628** | 0.612 → **0.630** |
| `tracerlpm_strict_parity` | 0.208 → **0.206** | 0.823 → **0.895** | 0.608 → **0.535** | 0.576 → **0.573** |
| `parity_reported_corrected` | 0.210 → **0.221** | 0.824 → **0.898** | 0.607 → **0.531** | 0.574 → **0.568** |
| `hydrosheaf_selection_corrected` | 0.487 → **0.441** | 0.980 → **0.952** | 0.231 → **0.272** | 0.385 → 0.385 |
| `tracer_young_only` | 1.239 → **1.297** | 1.611 → **1.682** | −0.500 → **−0.632** | 0.193 → **0.210** |

The pattern is consistent and interpretable: **typical-case accuracy improves slightly**
(lower median error, higher within-factor-2) while **tail behaviour worsens** (higher
RMSE, lower R²). Removing the clamp exposes genuine bomb-peak double-valuedness in
young-water tritium inference — the re-run emitted thousands of
`Ambiguous age inference for 3H=…: found 2 peaks` diagnostics that the clamped forcing had
suppressed. This is the honest behaviour of the physics, not a regression.

### Graph benchmark — the central negative result is now reversed

Single-well baseline RMSE improves **1.065 → 0.898** under corrected forcing, and for the
first time a candidate family beats it:

| Family (weak prior) | ΔRMSE before | ΔRMSE after | violations after |
|---|---|---|---|
| `parameter_smooth_aicc` | +0.003 | **−0.006 (improves)** | 87 → 68 |
| `hydraulic_proxy_constrained` | +0.019 | +0.008 | 87 → 88 |
| `depth_constrained` | +0.052 | +0.055 | 338 → 312 |
| `wrong_direction_negative_control` | +0.041 | **+0.069** | 125 → 141 |
| `randomized_negative_control` | +0.279 | **+0.354** (strong: +0.446) | 572 → 590 |

`improved_vs_single` is now **True** for weak parameter smoothing. The manuscript's
headline — *"no tested graph family universally improved log-age RMSE"* — is no longer
what the data show.

**How this should be reported.** The improvement (−0.006 log₁₀) is below the manuscript's
own pre-specified negligibility threshold (§2.9: changes below 0.01 log₁₀ are negligible
unless supported by tracer-withholding evidence and negative-control separation). The
defensible statement is therefore: *a weak, physically constrained prior produced a small
RMSE improvement and reduced age-ordering violations from 87 to 68, while wrong-direction
and randomised controls degraded RMSE by an order of magnitude more (+0.069 and +0.354).*
That is a **stronger** result for the paper's thesis than the previous flat negative,
because the separation between candidate families and controls is now much wider.

## Z-7. Prototype: topology as a bomb-peak branch selector (item 1) — INCONCLUSIVE on this dataset

`M3/m3_age_benchmark/scripts/run_m3_branch_selection_benchmark.py` implements the
recommended redesign: instead of smoothing the age field, the graph restricts each node's
**discrete tritium alias set** to solutions no younger than its selected upstream
neighbour, then takes the best tritium match among the feasible aliases. This is the rule
already demonstrated synthetically in `run_m3_network_dating_demo.py` (within-factor-2
0.63 → 0.84), applied for the first time to real USGS data with a randomised graph as the
negative control.

**Result (352 eligible sites, 184 ambiguous):**

| subset | mode | MdAE | RMSE | wf2 |
|---|---|---|---|---|
| ambiguous | single-node | 0.521 | 1.165 | 0.353 |
| ambiguous | **graph branch selection** | 0.510 | 1.155 | 0.337 |
| ambiguous | randomised control | 0.515 | 1.158 | 0.337 |

Graph selection is indistinguishable from the randomised control. **But the test is not
informative**, and the reason is the finding worth keeping:

- The depth-constrained rule yields **52 edges over 352 sites**, and the graph changes the
  selected alias on only **13 sites (3.7%)**. The remaining 96% are unconstrained and
  identical under all three rules, so the comparison is almost entirely noise.
- Cause: **median nearest-neighbour well spacing is 40–95 km** in 14 of 16 study units,
  against a 25 km edge criterion. USGS public-supply wells are production wells scattered
  across regional aquifers, not monitoring transects. The synthetic demo worked because it
  constructed chains in which every ambiguous node was bracketed by design.
- In the single dense unit (Biscayne, median spacing 2.4 km, 23 edges / 31 sites) the
  mechanism behaves as intended: MdAE 0.750 → 0.685, RMSE 0.814 → 0.719, wf2 0.032 →
  0.065, while the randomised control moves nothing at all (0.750 / 0.818 / 0.032). The
  direction and the control separation are both correct — but the graph acts on only 4
  sites, and Biscayne is a shallow young aquifer where tritium-only PFM dating against LPM
  reference ages is poor in absolute terms (wf2 ≈ 0.03–0.07).

**Verdict.** The mechanism is sound and correctly implemented, and it is the right target,
but it **cannot be validated on the USGS public-supply benchmark**. Do not report it as a
positive M3 result. The defensible statement is a field-design constraint:

> Flow-ordering priors can only resolve tracer ambiguity where well spacing is at or below
> the correlation length of the flow system. In the public-supply network evaluated here,
> spacing exceeds the usable edge criterion by a factor of 2–4, so fewer than 4% of
> ambiguous nodes are reachable by any upstream constraint.

## Z-8. MRVA re-test — the mechanism needs two conditions no public dataset satisfies

### Z-8a. M3 §2.2 declared three data sources; the benchmark used one

| Declared in §2.2 | Present in repo | Loaded by any M3 script |
|---|---|---|
| National public-supply (Jurgens et al., 2022a) | yes | **yes** — 1,272 rows, 20 study units |
| Western Principal Aquifers (Faulkner & Jurgens, 2019, `10.5066/P9U9ZSBN`) | **no** | **no** |
| MRVA alluvial (Gratzer et al., 2025a, `10.5066/P14DPCXE`) | **no** | **no** |

`load_usgs_national_dataset()` is the only loader in the M3 pipeline and reads only
`DataForNationalGroundwaterAge_1_1`. Neither missing DOI appears anywhere in M3's scripts
or configs. The likely origin of the error: the national release contains a study unit
named *"Mississippi embayment–Texas coastal uplands aquifer system (METX)"* (86 rows,
sampled 2014.3–2015.7), which is **not** the separate 2018–20 MRVA alluvial-aquifer
release. Both missing releases have since been downloaded (see Z-8b).

### Z-8b. Datasets retrieved

Downloaded from USGS ScienceBase into
`M2/m2_benchmark/external/usgs_age/input/`:

- **MRVA** (`MRVA_GroundwaterAge_2018_20/`, 0.92 MB) — Tables 2–7 and 9 plus the
  `Table1_Wells` shapefile. The shapefile is attached as a ScienceBase *facet*, not a
  listed file, so it is invisible to the normal file listing; its DBF carries `ALT_VA`,
  `WELL_DEPTH`, `OTOP_dep`/`OBTM_dep` and exact coordinates for 511 wells, and was parsed
  to `Table1_Wells_attributes.csv`. `Table8_LPMs.csv` returns HTTP 404 from ScienceBase
  (server-side; retried four times) and is the one file not obtained.
- **Western Principal Aquifers** (`WesternPrincipalAquifers_2004_2018/`, 1.23 MB) —
  `WPAsTables.zip` extracted to four tables including age interpretations (1.4 MB) and
  computed tracer concentrations (1.2 MB).

### Z-8c. MRVA result: density solved, ambiguity regime lost

MRVA is far denser, as expected — median nearest-neighbour spacing **15.5 km** against
40–95 km nationally, with 74% of sites having a neighbour inside the 25 km criterion, and
genuine nested piezometers (identical coordinates at 5, 8 and 10 m depth). But the test
still fails, for the opposite reason:

| | National public-supply | MRVA alluvial |
|---|---|---|
| sites | 352 | 63 |
| ambiguous (≥2 aliases) | 52% | 48% |
| live ³H (≥ 0.5 TU) | 89% | **75%** |
| reference age ≤ 90 yr (tritium-datable) | 64% | **44%** |
| graph-reachable (alias actually changed) | **3.7%** | 6.3% |
| **reachable AND datable** | **13 sites** | **2 sites** |

MRVA is an **old-water** aquifer: median reference age 100 yr, 75th percentile 730 yr,
maximum 22,000 yr, and 16 of 63 sites have tritium below 0.5 TU. The bomb-peak alias
structure the mechanism exploits simply is not present in most of it. On the
tritium-datable subset (n = 28) the graph changed 2 aliases and the randomised control
outperformed it (MdAE 0.256 against 0.276).

### Z-8d. Conclusion — a two-condition requirement

Graph-based branch selection requires **both** conditions simultaneously:

1. **Well spacing at or below the flow-system correlation length**, so ordering constraints
   reach ambiguous nodes. The national public-supply network fails this (3.7% reachable).
2. **Water inside the bomb-peak window with live tritium**, so a discrete alias ambiguity
   exists to resolve. MRVA fails this (44% datable, 25% dead tritium).

Neither public release satisfies both; the best usable samples are 13 and 2 sites. This is
why the synthetic demo succeeded — it constructed both conditions by design.

**This is the reportable finding.** It is a concrete, quantified field-design constraint
rather than a failed method, and it is more useful than the original inconclusive result:

> Flow-ordering priors can only resolve tracer ambiguity where a dense monitoring network
> coincides with young, tritium-live water. In the two public releases evaluated,
> connectivity and ambiguity are anti-correlated: the network with young water is too
> sparse to constrain it, and the network dense enough to constrain it is too old to be
> ambiguous.

Testing the mechanism therefore requires a purpose-built nested-piezometer transect in a
young-water setting, not another public release. Relaxing the distance criterion to
manufacture edges is not advisable — M4 establishes that false connectivity is the
dominant failure mode.

## Z-9. Three-release integration — attempted, measured, and scoped back

All three releases were harmonised (`M3/m3_age_benchmark/scripts/m3_data_sources.py`) and
the full design matrix was re-run on the union: 1,645 rows × 14 scenarios = 23,030 fits.
The measurement then showed that pooling is not defensible, and the integration was scoped
back to the roles each release can actually support.

### The decisive diagnostic: do scenarios discriminate?

| release | rows | distinct estimates per site across the 14 scenarios |
|---|---|---|
| National public-supply | 1,272 | **5.33** |
| Western Principal Aquifers | 290 | 3.73 |
| MRVA alluvial | 83 | **1.00** |

**1.00 means every scenario returns an identical age.** MRVA carries only tritium and no
reported-LPM table, so strict parity, age-fraction parity, young-only and every ablation
collapse to the same fit. Those rows add sample size and *zero* comparative information —
which is precisely what the design matrix exists to measure. Per-release R² was also
degenerate (Western −9.2, MRVA −2.6), reflecting rows fitted with default assumptions
against reference ages produced by a different workflow.

A defect in the harmonisation contributed to Western's 3.73: the scenario runner reads
`LPM_Name` and `LPM_TracersMod`, but the loader wrote `reported_model_name`, so Western's
parity scenarios ran without the reported-LPM information the release does contain. This
is fixable, but it does not change the conclusion for MRVA.

### Composition, not breadth

Western's net contribution after removing the 332 samples it shares with the national
release is **290 rows, all Central Valley**. Every non-CVAL reference-age sample was
already present nationally. Pooling would raise one Californian aquifer from ~6% to ~22%
of the benchmark and move aggregate metrics for reasons unrelated to method performance.

### What was kept

- **Scenario matrix (Tables 2, 3, 5, 6): national public-supply only.** Restored; the
  numbers in Part Z-6 stand unchanged.
- **Graph benchmark: national-only as primary (Table 4), plus a cross-aquifer replication
  with MRVA (new Table 4b).** This is where the new data earned its place.
- **Cross-validation: unaffected.** Eligibility requires a reported fitted tracer set and,
  for graph families, coordinates, so the CV support is inherently the national 794/260/1103
  rows. Re-running under the combined loader reproduced them exactly.

### The replication result — and it is a negative one

| family (weak prior) | national only (1,272 nodes) | + MRVA (1,346 nodes) |
|---|---|---|
| `parameter_smooth_aicc` | **−0.006 (improves)** | **+0.0003 (neutral)** |
| `hydraulic_proxy_constrained` | +0.008 | +0.008 |
| `wrong_direction_negative_control` | +0.069 | +0.051 |
| `randomized_negative_control` | +0.354 | +0.270 |

**The single positive graph result does not replicate.** Adding 74 wells from a
structurally different aquifer removes the margin entirely, while the control ordering is
preserved. M3 now states the graph result as neutral rather than beneficial, and identifies
the reproducible signal as *rejection of implausible topology*, not gain from plausible
topology. That is a weaker headline and a better-supported one.

Section 2.2 was rewritten to describe this design explicitly, including that Western
publishes no well coordinates, that its net contribution is one aquifer, and that MRVA's
`Table8_LPMs` file is unavailable.

### Z-9a. A reproducibility sensitivity found during the scope-back

Re-running the tracer-withholding CV under the combined loader reproduced the single-well
and candidate-graph values **exactly** (³H 2.161 / 2.428, SF₆ 4.380 / 4.309, ¹⁴C 26.003)
but shifted the **randomised negative control** for ³H from 3.011 to 2.970. The
CV-eligible rows were identical in both runs; the difference is that
`build_random_graph` draws from the *whole loaded node pool*, so adding MRVA's rows and a
21st study unit changed the seeded draw even though none of those rows entered the CV.

The randomised control is therefore sensitive to dataset composition in a way the
candidate families are not. This does not affect any conclusion — the control still
degrades prediction far more than any candidate graph — but it means the control value is
only reproducible against a stated dataset scope. All reported CV numbers are the
national-public-supply run, consistent with the scenario matrix.

---

## Z-10. M7 cross-check — unaffected by the input-history defect

The tritium clamping defect (B-6) affects any module that forces ³H through
`build_site_tracer_histories`. M7 is the thesis's strongest chapter (pre-registered
confirmatory design, protocol committed at `d336e87` before the 5301-series test cases
were generated), so it was checked explicitly.

**Result: M7 is unaffected.** Three independent confirmations:

1. M7 never calls `build_site_tracer_histories`, `build_default_tritium_input`, or any
   WISER loader. Its only tritium source is a constant boundary constructed in
   `strong_inference.py:269`: `InputHistory([1850.0, 2035.0], [6.2, 6.2])`.
2. Because that history is **constant**, the clamping behaviour is inert. Probing recharge
   years from 1700 to 2100 — far outside the declared range, and covering M7's full
   `max_age_years = 200` window from a 2025.5 sample date — returns exactly 6.2 TU in
   every case. Clamping a constant returns the constant.
3. Of the two package files modified during this audit, `input_history.py` was changed
   only inside `build_default_tritium_input` (verified by diff hunk), and
   `tracer_inputs.py` is not imported by M7. The `InputHistory` class that M7 does import
   is unchanged.

No M7 re-run is required.

### Z-10a. Why M7 and M3 disagree about graph priors — a mechanism, not a contradiction

M7's constant 6.2 TU boundary means predicted tritium is `6.2 · exp(-λτ)`, which is
**strictly monotonic in age**. There is no bomb peak, hence no double-valuedness and no
alias ambiguity.

| | tracer–age relation | graph prior effect |
|---|---|---|
| M7 (constant synthetic input) | monotonic, uniquely invertible | **improves** age MAE by 0.062 y; reversed graph harms by 0.282 y |
| M3 (real GNIP/WISER bomb peak) | double-valued across the peak | **neutral to harmful**: +12.3% withheld-³H RMSE |

The graph helps in M7 because age is already identifiable and topology merely refines it.
It fails in M3 because age is *discretely* ambiguous, and an age-smoothing prior cannot
resolve a discrete ambiguity — it can only average across the alternatives.

This sharpens H1 into a conditional law with a named mechanism:

> Graph constraints improve age inference where the tracer–age relation is monotonic and
> topology is correctly specified. Under a bomb-peak input the relation is not monotonic,
> and a smoothing prior is neutral at best regardless of topology quality.

It also retro-justifies the Z-7 branch-selection design: the correct response to a
discrete ambiguity is a discrete *selector*, not a smoother.

---

## Z-11. M5 and M6 audited to the same depth — both sound

### M5 (inverse-reaction benchmark) — verified, no defects

All eleven headline values reproduce exactly from
`results/analysis_summary.json` (regenerated 2026-07-27T01:01Z): 240 PHREEQC
scenarios, 21,600 factorial fits, 160 field pairs, guarded phase F1 0.563,
class F1 0.609, FDR 0.361, Core class F1 0.606 / FDR 0.383, legacy elastic-net
FDR 0.465, Core evidence score 0.622.

**M5 self-reports its own weakest result**, which is exactly the behaviour M3 and
M4 lacked:

> "…mean reaction-evidence score of 0.622…, but only **23.0%** of reaction rows it
> flagged as high evidence corresponded to a truly active reaction…, so a high
> evidence score from this gate should be read as a moderate prior in favour of a
> reaction, not as a reliable positive indicator on its own."

*(Note: an initial pass flagged this as unreported. That was an artifact of the
pandoc plain-text export dropping digits. Verified directly from the docx XML —
the manuscript does report it. Manuscript checks must use python-docx, not
pandoc text.)*

**Baseline check** — the test that exposed M4's inflated headline. Recomputed
from `data_tier_reaction_evidence.csv` (core tier, n=3,840, base rate of truly
active reactions 0.177):

| strategy | F1 | precision | recall |
|---|---|---|---|
| **Hydrosheaf-Core** | **0.588** | 0.568 | 0.611 |
| Informed: archetype-modal reactions | 0.484 | 0.517 | 0.455 |
| Trivial: all reactions active | 0.301 | 0.177 | 1.000 |

Unlike M4, the margin is real (+0.104 F1 over an informed baseline). M5's own
comparator (legacy thermodynamic elastic-net, phase F1 0.551) is harder still,
and against it the honest claim is *comparable F1 with a materially lower
false-discovery rate* (0.361 vs 0.465) — which is what the manuscript says.

**Only recommendation:** the 23.0% figure would be easier to judge if the base
rate (0.177) were stated alongside it — it is a 1.30x enrichment, rising to
2.40x at evidence >= 0.90. The evidence score's discrimination is weak but real
(AUC 0.631). Adding one clause would close this.

### M6 (field-transfer benchmark) — independently reproduced

M6's foundational data-readiness claim was recomputed from the raw workbook
(`data/FieldData/NorthenGhana/NorthernGhana.xlsx`, sheets `Dry` + `Wet`) using an
independent charge-balance calculation:

| quantity | M6 reports | recomputed |
|---|---|---|
| records | 320 | **320** (160 + 160) |
| median abs charge-balance error | 1.5% | **1.49%** |
| quantitative / screening / exploratory | 294 / 19 / 7 | **294 / 19 / 7** |

Talensi and Lower Anayari also reproduce within the difference expected from
ion-set choices (Talensi median |CBE| 30.0% recomputed vs 29.85% reported, 0
quantitative in both; Lower Anayari 3.07% vs 3.13%). Every other headline value
checked — 160 wells, MRS 70.6-71.2, bootstrap stability 0.965-0.976, silicate
weathering 114/160 — is present in the manuscript.

No defects found in M6.

### Audit coverage after this round

| manuscript | audited to depth | headline defects found |
|---|---|---|
| M2 | earlier round | F1 = 0.86 with no computational source |
| M3 | this audit | forcing bug inverted the flagship result |
| M4 | this audit | headline inflated ~10x against the right baseline |
| M5 | this audit | **none** |
| M6 | this audit | **none** |
| M7 | Z-10 | **none** (unaffected by the forcing bug) |

## Z-12. Automated guard against artifact rot

`scripts/check_manuscript_artifacts.py` re-runs the M3 and M4 table builders and
diffs their output against the tracked tables, failing with a non-zero exit
status on any mismatch. It targets the exact failure mode found three times in
this audit: a committed table that no longer corresponds to the code claiming to
produce it.

Verified working in both directions — it passes on the current repository, and a
negative test (injecting `log10_rmse = 9.999` into Table 5) is correctly caught
and localised to the column and row.

**Incident during development, recorded because it matters:** the first version
restored the tracked directory with `shutil.rmtree()` + `copytree()`. On Windows
the `rmdir` failed *after* the contents had been deleted, destroying all ten M3
tables. They were deterministic outputs of committed results and were fully
regenerated and re-verified within minutes, with no loss. The script now restores
file-by-file and never deletes the tracked directory. A tool whose purpose is to
protect artifacts must not be able to destroy them.

---

## Z-13. Framework-level audit — what the manuscripts never exercise

Prompted by the question "are the calibrations working, and did we check the full
power of Hydrosheaf?", the knowledge graph was queried for the framework's
subsystems and each was checked against manuscript usage.

### Z-13a. A fundamental defect: asymmetric input-history bounds (FIXED)

`hydrosheaf/nuclear/lpm.py` guarded the **lower** bound of the tracer input
history and passed an explicit `left=` to `numpy.interp`, but had **no check on
the upper bound and no `right=`**:

```python
if t_recharge < input_years[0]:          # checked AND warned
    warnings.warn(...)
c_in = np.interp(t_recharge, input_years, input_values,
                 left=float(input_values[0]))   # right= absent
```

Out-of-range behaviour had clearly been considered — in one direction only. The
unguarded direction is the one that matters: station records end *before* the
sample date, not after. That asymmetry is the root cause of B-6 and it stayed
silent for 99.8% of M3 samples.

**Fixed**: a symmetric upper-bound warning and an explicit `right=`. Numerics are
unchanged (clamping is preserved, so no result moves), but a truncated input
history can no longer be consumed silently. Verified: the warning fires on a
record ending 1986.8 sampled in 2010, and 66 related tests still pass.

### Z-13b. Test-suite state

Full suite: **613 passed, 5 failed, 13 errors, 2 skipped**.

- **Calibration specifically: 77 passed, 0 failed.** The calibration subsystem is
  tested and green.
- 4 failures + all 13 errors are in `tests/synthetic_data_tests/` and trace to a
  missing data fixture (`test_paths_exist` fails first) — an environment gap, not
  a code defect.
- 1 genuine failure: `test_network_aging::test_bayesian_network_inference`,
  `age_r_hat_max = 1.135 > 1.1`. This is a marginal MCMC convergence miss in a
  deliberately fast 50-sample configuration. Notably it is the *same* R-hat rule
  M7 enforces, and M7 reports every case met it — so the rule works; the unit
  test is simply under-sampled. Worth raising the test's sample count rather than
  relaxing the threshold.

### Z-13c. Roughly half the framework is never benchmarked

| subsystem | lines | benchmark scripts using it |
|---|---|---|
| **calibration** | **6,249** | **0** |
| nuclear | 6,023 | 10 |
| sheaf | 2,660 | 1 |
| **temporal** | **2,579** | **0** |
| inference | 2,481 | 4 |
| models | 2,673 | 3 |
| **vadose** | **2,003** | **0** |
| **graph3d** | **1,730** | **0** |
| **transport** | **1,541** | **0** |
| graph | 1,531 | 8 |
| **reactive_transport** | **1,454** | **0** |
| uncertainty | 1,445 | 1 |
| **tuning** | **527** | **0** |
| **causal** | **459** | **0** |

**About 16,500 lines — including the single largest subsystem — are exercised by
no manuscript benchmark.** The calibration package contains a PEST-style
Gauss-Levenberg-Marquardt engine (`PESTGLM`), an active-learning module (821
lines), bootstrap confidence intervals, topology and isotope-mixing adapters, and
a validation workflow. It is not even referenced from the public `api.py`.

**Implication, two-sided.** Risk: a large body of unvalidated code ships with a
framework whose stated contribution is validation discipline, and a reviewer who
notices will ask why. Opportunity: a tested, unbenchmarked calibration engine is
an entire unclaimed contribution — the obvious next manuscript, and the one that
would answer "does Hydrosheaf calibrate?" rather than only "does Hydrosheaf
falsify?".
