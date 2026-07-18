# M2–M4 Figure & Table Accuracy Audit (2026-07-10)

## FINAL REPRODUCIBILITY VERIFICATION (all three manuscripts)

Every headline number was recomputed from the committed result files. Status:

**M4 — 100% reproducible.** head_gradient 147/155/27 → F1 0.618 (P 0.487, R 0.845);
head_depth 0.616; hydrostratigraphic 0.614; spatial 0.000; sparse 0.445. All match
`independent_graph_vs_modpath.csv` exactly.

**M3 — 100% reproducible.** age-fraction parity median|log10|=0.167, RMSE=0.740,
R²=0.681, wf2=0.614; strict parity 0.221; network-dating demo ambiguous 0.63→0.84;
depth-constrained withheld-tritium RMSE 12.42→7.97 (35.9%). All match M3 result CSVs.

**M2 — reproducible after this round's fixes.** A second pass found more stale
synthetic-tier numbers (from the pre-upgrade run) and corrected them to reproducible
values:
| Item | Was (stale) | Now (reproducible) | Source |
|---|---|---|---|
| Age R² | 0.86 | **0.98** (MAE 183.3 unchanged; per-class 0.82–0.97) | age_inference_validation.csv |
| Isotopic R² | 0.962 | **δ¹⁸O 0.52, δ²H 0.96** | res_isotopes.csv |
| Reaction identification | "91% of realisations" | **selection precision 0.30, recall 0.45** (class-level; consistent with equifinality) | reaction_recovery.csv |
| Table 18 R² column | 0.95/0.82/0.71/0.86 | **0.97/0.85/0.82/0.98** | age_inference_validation.csv |

Already-verified M2 (unchanged): chemistry 0.955, topology 0.618/1.00, PHREEQC feasible
0.89, field 0.711 (LA 0.60 / Talensi 0.77), PSI (Talensi pyrite 0.49, gypsum 0.00; LA
pyrite 0.00), reaction extent error 0.050, age MAE 183.3.

**All figures regenerated from current data; Figure 7 rebuilt from the geologically-correct
per-site PSI (Talensi excludes gypsum). No stale numbers remain in text, tables, or
figure captions (full scan clean).**

---


## UPDATE — M2 regenerated with current upgraded code (both Ghana sites)

Per author decision, M2 was re-run with the current (null-analysis) Hydrosheaf and the
manuscript updated to the reproducible numbers:

| Tier | Old | Current code (APPLIED to docx + Table 2 + figures) |
|---|---|---|
| Ghana field R² | 0.94 | **0.711** (Lower Anayari 0.60 / Talensi 0.77; both sites, 82+126 edges) |
| Synthetic reaction | R²=0.86, RMSE=1.38 | **chemistry-reconstruction R²=0.955**; extents class-level (extent R²≈0) |
| PHREEQC forward | NSE=0.58, RMSE=1.39 | **mass-balance proxy, feasible fraction 0.89**; live PHREEQC pending |
| Topology no-prior | 0.86 | **0.618** (done earlier) |
| Residence-time | R²=0.86, MAE=183.3 | unchanged — verified stable |

Figures 2,3,4,5,6 regenerated from current data. Field process *selection* still matches
the narrative (LA = exchange+calcite; Talensi = denitrification+pyrite+exchange).

### UNRESOLVED — PSI tier (§3.8, Table 6, Figure 7): cross-benchmark contradiction
- **Figure 7 sources PSI from an M6 file** (`M6/.../m6_phase_stability_index.csv`) and shows
  **Talensi gypsum PSI≈0.99, pyrite low** — geochemically wrong for a crystalline mining
  area and OPPOSITE to the §3.8 text.
- **§3.8 text** (from the old M2 per-reaction PSI run) says **Talensi pyrite 1.00,
  denitrification 1.00**.
- **M2's own current PSI generator** (`run_edge_psi.py` → `top_edges_psi.csv`) outputs only
  *family-level* PSI and, in Monte Carlo, flips SO₄ attribution between gypsum and pyrite
  (the known equifinality). The crystalline gypsum-suppression gate is not disciplining this.
- Figure 7 was **reverted** to the submitted version (byte-identical) to avoid worsening the
  text↔figure state. **This inconsistency is PRE-EXISTING in the submitted manuscript**, not
  introduced by the rerun.
- **Decision needed:** which PSI source is authoritative (M2 vs M6), and fix the
  gypsum-vs-pyrite equifinality (enforce crystalline gypsum suppression for Talensi) before
  §3.8/Table 6/Fig 7 can be made self-consistent. NOT auto-edited.

---


Cross-checked every numeric figure and table in the three manuscripts against the
current committed result files.

## M4 — CLEAN ✅
- **Table 2 (evidence ladder)** and **Table 3 (topology summary)**: every value matches
  `independent_graph_vs_modpath.csv` exactly (head-gradient 0.487/0.845/0.618/155 FP/27 FN;
  head-depth 0.616/156; hydrostratigraphic 0.614/158; negatives 0.006/0/0; sparse 0.445).
- **Figs 2–3** render the same numbers. No action needed.

## M3 — CLEAN ✅
- **Tables 10–14** are internally consistent with the manuscript text and abstract
  (age-fraction parity 0.167/0.740/0.681/0.614; strict parity 0.221; Hydrosheaf selection
  0.487; negative controls degrade as reported). Honest reference paper; no action.

## M2 — MIXED: topology fixed, but a deeper sync gap found ⚠️

### Verified correct (manuscript matches current data)
- **Table 3 / Table 18 (residence-time)**: age MAEs match `age_inference_validation.csv`
  **exactly** (overall 183.31; young 1.25; intermediate 22.82; old 333.26); class R² in
  log-space (young 0.94, etc.). ✅
- **Table 4 / Table 19 + Figure 2 (topology)**: CORRECTED to the reproduced
  147/155/27 → F1 0.618, and the figure rebuilt to show both modes honestly. ✅
- Transport-parameter error (0.061) and reaction-extent error (0.050) match
  `transport_recovery.csv` / `reaction_recovery.csv`. ✅

### NOT reproducible from current artifacts — earlier-run provenance ❌
The current M2 package (regenerated in commit 3102680, *after* the manuscript) produces
materially different numbers for three tiers. The manuscript's values appear in **no**
current result file, table, or script — the same situation as the topology 0.86 was:

| Tier | Manuscript (docx) | Current pipeline / CSVs |
|---|---|---|
| Synthetic **reaction R² = 0.86, RMSE = 1.38** | abstract, Table 2, §3.2 | `edge_fit chemistry_r2` median = **0.999**; reaction error 0.05 |
| PHREEQC **NSE = 0.58, RMSE = 1.39** | abstract, Table 2, §3.6, conclusions | `phreeqc_forward` NSE = **0.999**, RMSE = **0.046** |
| Ghana field **R² = 0.94** (LA 0.88 / Talensi 0.91) | abstract, Table 2, §3.7, §3.6 | `field_discovery` R² = **0.991**, **Lower-Anayari only — no Talensi edges present** |

Consequences:
1. **Manuscript text ↔ its own figures are inconsistent** for these tiers: Figs 3/5/S2/S3
   are generated from the current CSVs (≈0.999) while the text/tables report the earlier
   run (0.86 / 0.58 / 0.94).
2. The current numbers are *cleaner* (0.999), so — unlike topology — the manuscript's more
   conservative values are **not** obviously wrong; they may be the better science and the
   re-run may be over-idealised.
3. **Talensi field results are absent** from the current package, yet the manuscript builds
   a two-site (Lower Anayari + Talensi) narrative (PSI tables, Fig 4, Fig 7).

### Recommendation (needs author decision — NOT auto-edited)
Do **not** blindly replace 0.86/0.58/0.94 with the current 0.999 — that would degrade the
paper and remove Talensi. Instead:
- **Locate or re-run the benchmark that produced the manuscript numbers** (the pre-3102680
  run, including the two-site Ghana field analysis), regenerate figures + tables from it,
  and lock the manuscript to that single reproducible run; **or**
- if the current run is to be authoritative, re-run the **full** two-site field analysis,
  regenerate all figures/tables, and update every text number together — accepting the
  0.999 synthetic recovery must then be justified (it looks too perfect for a defensible
  synthetic benchmark).

Only the **topology** tier could be pinned to a verifiable ground truth (via M4) and was
therefore the only one corrected in place. The reaction/PHREEQC/field mismatch is a
reproducibility-provenance decision for the author.
