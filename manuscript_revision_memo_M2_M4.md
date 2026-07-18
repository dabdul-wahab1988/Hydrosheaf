# M2–M4 Manuscript Revision Memo

> **STATUS (2026-07-10): M2 topology F1 correction COMPLETE and reproduced.**
> Provenance audit confirmed the submitted M2's no-prior F1 = 0.86 (215/168/47/6) has
> **no computational source** anywhere in the repo or its git history and exceeds the
> benchmark's 147-TP information ceiling. The correct value is **F1 = 0.618**
> (302 cand / 147 TP / 155 FP / 27 FN; P 0.487, R 0.845), reproduced by
> `M2/m2_benchmark/scripts/run_m2_modpath_noprior_topology.py` (matches the M4 benchmark
> exactly). Applied to `M2/M2_ready/Manuscript-CG.docx` (abstract, §3.4, §4.2.2 disc.,
> §4.3 disc., conclusions, Table 4, Table 2). Backup: `Manuscript-CG_BACKUP_preF1fix.docx`.
> The M3-demo 0.63→0.84 and M4 0.618 were verified as REAL.
>
> **ALSO COMPLETE (2026-07-10):**
> - **M2 Figure 2 corrected** — was showing only the prior-assisted 174/0/0 (F1=1.00)
>   mislabelled as topology validation; now a two-panel figure showing independent
>   no-prior (147/155/27, F1=0.618, overconnection) beside the prior-assisted fidelity
>   check. Regenerated from `run_m2_modpath_noprior_topology.py`; figure script updated;
>   old PNG backed up. (Only Fig 2 used topology metrics; other M2 figures unaffected.)
> - **M2 network-age reframed** (¶44, 127, 138, 174, 188) — the 1/√(1+in-degree)
>   σ-reduction is now stated as an uncertainty-partitioning model result *by
>   construction*, not a measured accuracy gain; a limitations caveat cites the companion
>   age study (consistency/violation-reduction, not universal RMSE gain; old-water
>   non-identifiability).
> - **M4 author block fixed** — phantom authors EEK/RW/ZN/ABK/CKB removed; rewritten as
>   the 5 real authors with DA-W lead (matches series). **Verify exact CRediT role split
>   with co-authors.** Backup: `..._BACKUP_preAuthorFix.docx`.
> - **M4 duplicate references removed** (Moracchini, Pang, Pollock — one copy each).
>
> Still open (Medium, revision round): M2 field-R²=0.94 wording → "chemistry-closure,
> no process truth" + class-level identifiability note.



Cross-checked the three companion manuscripts (M2 umbrella, M3 age benchmark, M4
topology benchmark) against the honesty findings established in M6/M7 and in the
network-dating investigations. The three papers are at very different maturity levels
of honesty. Ranked by severity.

---

## CRITICAL — M2 topology F1 contradicts M4 (same benchmark, same authors)

**The single most important fix.** M2 and M4 both benchmark the *identical* USGS Savage
MODPATH reference (174 directed edges), but report different independent-topology skill:

| | Prior-assisted | Independent / no-prior |
|---|---|---|
| **M4 manuscript** (dedicated topology paper, submitted to AWR) | F1 = 1.00 (labelled Level-6, *not* independent) | **F1 = 0.618** (P = 0.487, R = 0.845; 147 TP / 155 FP / 27 FN) |
| **M2 manuscript** (Manuscript-CG.docx §3.4, abstract, conclusions) | F1 = 1.00 (headline) | **F1 = 0.86** (P = 0.78, R = 0.96; 168 TP / 47 FP / 6 FN) |
| **M2 benchmark outputs** (m2_results_summary.md) | F1 = 1.00 | *not present* |

Two problems:
1. **M2's no-prior F1 = 0.86 is not supported by M2's own current benchmark data.** The
   M2 results summary carries only the prior-assisted `TP=174; FP=0; FN=0; F1=1.00`
   row. The 0.86 / 168-TP-47-FP figures appear to be a stale pre-honesty number.
2. **It directly contradicts M4's honest F1 = 0.618 on the same Savage reference.** Any
   reviewer reading both companion papers will catch this immediately. M4's 0.618 is the
   correct, honest independent number (the low precision — FP 155 > TP 147 — is the real
   behaviour of the elevation-as-head steepest-descent rule; see the "degenerate flat-z"
   analysis).

**Action (M2):**
- Abstract: change "topology F1-scores of 1.00 and 0.86 under prior-assisted and
  heuristic-only conditions respectively" → report the independent number as
  **F1 = 0.618 (precision 0.49, recall 0.85)** and explicitly state that the F1 = 1.00
  prior-assisted result is *confirmation that MODPATH priors are transmitted faithfully,
  not independent topology recovery.*
- §3.4 and §4.2.2 and Conclusions: replace the 215/168/47/6 → P=0.78/R=0.96/F1=0.86
  passage with M4's 147 TP / 155 FP / 27 FN → P=0.487 / R=0.845 / F1=0.618, and add the
  overconnection caveat (best independent graph still produces more false positives than
  true positives; an inferred edge is a *candidate hypothesis*, not validated
  connectivity).
- Rerun/point Table 4 to the M4 independent numbers so M2 and M4 are numerically
  consistent.

---

## HIGH — M2 network-Bayesian age claim overstates a benefit M3 shows is conditional

M2 presents the `1/√(1+in-degree)` posterior-uncertainty reduction as a validated benefit
("σ_net consistently lower than single-node prior", "reducing posterior uncertainty
relative to single-node LPM estimates"). But:
- **M3 (dedicated)** found the graph-age RMSE benchmark was *negative*: **no graph family
  universally improved log-age RMSE** on real USGS data; graph priors deliver
  *consistency (violation reduction), not accuracy*, except in specific tracer–topology
  regimes (depth-constrained ³H, 35.9% RMSE gain).
- **The `infer_network_ages_bayesian` MCMC path** is *not identifiable for old water*
  (tritium dead beyond ~60–70 yr → velocity/input-scale/age posterior degenerate). An
  explicit caveat now lives in `network_aging.py`.

M2's claim isn't false (uncertainty *does* shrink by construction), but it conflates
*uncertainty-interval narrowing* with *accuracy improvement*, which M3 explicitly refuses
to claim.

**Action (M2):**
- §3.3 / §4.1 / §4.3: reframe the network-age contribution as **age-order consistency and
  uncertainty-partitioning, not accuracy**. Add one sentence: "Consistent with the
  dedicated age benchmark (companion paper), graph priors reduce age-ordering violations
  and narrow posteriors but do not universally reduce reference-age RMSE; accuracy gains
  are tracer- and topology-specific (e.g. depth-constrained tritium)."
- Note the identifiability limit: network Bayesian dating is reliable in the
  tracer-informative / ambiguity-resolution regime; old-water accuracy needs an
  independent velocity constraint or a longer-lived tracer (¹⁴C).
- Optionally cross-cite M3 for the falsifiable-prior framing.

---

## MEDIUM — M2 reaction-recovery and field R² framing

- Abstract "reaction-recovery R² = 0.86" and "Field deployment … overall R² = 0.94" are
  legitimate (synthetic recovery; downstream-chemistry mass-balance closure) but the field
  R² = 0.94 reads as *validation accuracy* when it is chemistry-reconstruction closure with
  **no independent process truth** (the M2 summary itself flags this). M6/M7 established
  that reaction identifiability is **class-level (equivalence classes)**, not point-level.
- **Action (M2):** in the abstract/discussion, qualify the field R² as
  "chemistry-reconstruction closure (no independent process-truth labels)" and add a clause
  that reaction attributions are resolved to identifiability *classes*, with point-level
  extents non-unique under sparse ions — aligning with the MRS/identifiability language
  used in M5/M6. The Limitations section already says most of this; the fix is to stop the
  abstract headline from implying accuracy.

---

## LOW — M4 author-contribution block is wrong (copy-paste error)

M4 (submitted to AWR) is otherwise the **gold-standard honest paper** — evidence ladder,
F1 = 0.618 independent, precision guardrail (FP > TP), prior-assisted separation, negative
controls all correct. No scientific revision needed. But:

- **Author contributions** read: "EAA conceived the study, led the investigation and
  drafted the manuscript. EAA, DA-W, EEK, RW, ZN, ABK, and CKB reviewed…". The authors
  **EEK, RW, ZN, ABK, CKB do not exist in the M4 author list** (Abdul-Wahab, Asare,
  Adomako, Abass, Ganyaglo). This is a copy-paste block from a different manuscript.
- It also assigns lead/conception to EAA (Asare), inconsistent with M2/M3 where
  Abdul-Wahab is lead and corresponding author.
- **Action (M4):** rewrite the contributions to the actual five authors with DA-W as
  lead/corresponding, matching M2/M3.
- Minor: the reference list has each of Moracchini, Pang, Pollock duplicated — dedupe.

---

## M3 — essentially no revision needed (reference standard)

M3 (Applied Geochemistry) already embodies every honesty finding: graph topology as a
**falsifiable, rejectable prior**; the main graph-age RMSE benchmark reported as
**negative**; negative/wrong-direction controls as the core contribution; age-fraction
parity flagged as **circular (reference-workflow emulation)**; tracer-withholding as the
**only independent** test; explicit non-transfer to Ghana. It is the model the other two
should be aligned to.

- **Optional strengthening only:** add the controlled bomb-peak ambiguity-resolution
  demonstration (`run_m3_network_dating_demo.py`; within-factor-2 0.63 → 0.84 on ambiguous
  nodes) as the constructive "*when* does topology genuinely improve accuracy" case — it
  complements the real-data negative result with a mechanism-level positive. Not required
  for correctness; it would pre-empt the "so when is the graph ever useful for accuracy?"
  reviewer question.

---

## Summary of what to change

| Paper | Severity | Change |
|---|---|---|
| **M2** | Critical | Replace no-prior topology F1 = 0.86 with M4's honest **F1 = 0.618** (P 0.49 / R 0.85); label F1 = 1.00 as prior-assisted, not independent. Reconcile Table 4. |
| **M2** | High | Reframe network-Bayesian age as consistency/uncertainty-partitioning, **not accuracy**; add identifiability + conditional-benefit caveat cross-citing M3. |
| **M2** | Medium | Qualify field R² = 0.94 as chemistry-closure (no process truth); state reaction identifiability is class-level. |
| **M4** | Low | Fix wrong author-contribution names (EEK/RW/ZN/ABK/CKB → actual 5 authors, DA-W lead); dedupe references. |
| **M3** | None / optional | Already honest. Optionally add bomb-peak ambiguity-resolution demo. |
