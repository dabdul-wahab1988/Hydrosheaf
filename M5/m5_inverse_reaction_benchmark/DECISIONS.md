# M5 manuscript decisions (r2m Track B)

## 2026-07-26 — Field-data source correction

The user reported suspected errors in the field data used for the M5-M7
analyses and asked that all three use data strictly from
`data/FieldData/` (three canonical datasets: `LowerAnayari/manu.csv`,
`Talensi_MiningArea/talensi.csv`, `NorthenGhana/NorthernGhana.xlsx`). M5 was
fixed first, module by module, per explicit user instruction.

**What was actually wrong.** Two distinct problems, not one:

1. **Broken paths.** `TALENSI_CSV`, `LOWER_ANAYARI_CSV`, and the Northern
   Ghana workbook constant in `run_m5_inverse_reaction_benchmark.py` pointed
   at `data/Talensi_MiningArea/`, `data/LowerAnayari/`, and
   `data/NorthenGhana/` (no `FieldData/` prefix). None of those paths exist
   on disk any more; only `data/FieldData/...` does. Every field-dependent
   stage of the pipeline (`run_ghana_field_demonstration`,
   `_legacy_northern_ghana_external_edges`, the Talensi/Lower Anayari ELRI
   edges) would raise `FileNotFoundError` if rerun today. Talensi and Lower
   Anayari were a pure prefix fix: the canonical CSVs have exactly the
   columns the code already expected (verified by reading both files).

2. **A non-canonical Northern Ghana source.** `run_ghana_field_demonstration`
   / `load_ghana_pairs` read `data/NorthenGhana/Aquifers_Dataset_Mendeley.xlsx`
   (sheets `Hydrochemistry_Seasonal` + `Wells_Nodes`), not the plain
   `NorthernGhana.xlsx` (sheets `Dry`/`Wet`) that `_legacy_northern_ghana_
   external_edges` already used. That Mendeley-schema file does not exist in
   `data/FieldData/` or anywhere else in the repository. Reading its
   accompanying `data/FieldData/NorthenGhana/SI.pdf` (17 pages) established
   what it actually was: the Supplementary Information for a **different,
   antecedent study** ("Graph-inverted Ghanaian aquifers under
   aridification"), covering the same 160 boreholes but adding that study's
   own aquifer-type/geology-group/lithology/land-use classification,
   precomputed saturation indices, 35 years of separate climate/aridification
   data, and — critically — **397 already-computed graph edges**. Reusing it
   would have imported another study's classifications and its own
   connectivity graph into what is supposed to be an independent Hydrosheaf
   benchmark; `M6/m6_field_transfer_benchmark/DECISIONS.md` had already
   flagged unresolved overlap with this same antecedent project as a
   reviewer concern without tracing it to this specific dependency.

**User decision (asked explicitly, given the scope change):** rebuild the
Northern Ghana loader from the raw `Dry`/`Wet` sheets only, dropping every
column that only the Mendeley file supplied, rather than sourcing or
reconstructing that workbook. Confirmed via `AskUserQuestion` before any
code was rewritten.

**What changed in code** (`run_m5_inverse_reaction_benchmark.py`):

- `TALENSI_CSV`, `LOWER_ANAYARI_CSV`, and the Northern Ghana workbook
  constant (renamed `NORTHERN_GHANA_WORKBOOK`, dropping the separate
  `GHANA_WORKBOOK`/`LEGACY_NORTHERN_GHANA_WORKBOOK` split) now all resolve
  under `data/FieldData/`.
- `load_ghana_pairs()` rebuilt to concatenate the raw `Dry`/`Wet` sheets
  (`Season` tagged directly) instead of reading the Mendeley sheets. It now
  returns one DataFrame instead of `(hydro, wells)`, since the `wells` half
  had no other caller.
- `run_ghana_field_demonstration()`: removed the `Data_Class == "Quantitative
  inverse modelling"` filter (that column does not exist in the raw file —
  every well with a complete-enough panel in both seasons is used instead,
  which is the same >=6-measured-ion gate the rest of the function already
  applied); removed `Sampling_Date`-based year validation and the
  `wet_date`/`dry_date` output columns (no date column in the raw sheets;
  the 2025 March/August campaign dates are independently confirmed from
  `SI.pdf` Supplementary Table 19, not invented); replaced every
  `Aquifer_Type`/`Geology_Group`/`Lithology` output field with `Region`/
  `District`, which the raw sheets do carry.
- **Saturation indices**: the function previously built `upstream_si` from
  `Calcite_SI`/`Dolomite_SI`/`Gypsum_SI`/`Halite_SI` columns that only
  existed in the Mendeley file. No PHREEQC installation is available in this
  environment to recompute them (`find_phreeqc()` fails; verified), so
  `upstream_si` is now an explicit empty dict — the same treatment
  `_external_field_edge_outputs` already used for the Talensi/Lower
  Anayari/legacy-Ghana ELRI edges. This is not a no-op with side effects
  hidden elsewhere: `thermodynamic_bounds()` and
  `count_thermodynamic_violations()` in `m5_common.py` were read directly
  and both `continue` (no bound tightened, zero violations counted) whenever
  a phase's SI is non-finite, so the thermodynamic-bound and SI
  evidence-gate components are genuinely inactive for this section, not
  silently biased. This is disclosed in Methods and Discussion rather than
  left implicit.

**Consumers updated to match** (all "aquifer"/"lithology" facets replaced
with "region"): `make_m5_publication_figures.py` (Figure 6 and Figures
S13/S14/S16), `make_m5_publication_tables.py` (Table S10, Table S13),
`r_figures/plot_m5_publication_figures.R` (`aquifer_names`/`clean_aquifer`
-> `region_names`/`clean_region`, Figure 6 panel d). The DuckDB export/figure
scripts (`export_m5_results_database.py`, `make_m5_database_figures.py`) are
schema-agnostic and needed no changes.

**Numeric effect** (rerun with `--reuse-synthetic`; synthetic-tier PHREEQC
scenarios are unaffected by field data and were confirmed byte-identical to
a synthetic-only rerun of the pre-fix code — see below): Northern Ghana
wet-to-dry pairs 138 -> 160 (all wells now used, not a 294-record
Data_Class-filtered subset that no longer exists); median Hydrosheaf-Core
evidence score 0.681 -> 0.664; median TDS consistency score 0.941 -> 0.942;
all pairs remain classified `partially_identifiable` (0% fully identifiable)
in both the old and new run. The external field ELRI transfer numbers
(NorthernGhana.xlsx 160 edges/0.072 median ELRI, Talensi 85/0.072, Lower
Anayari 49/0.010) are **unchanged**, because `_legacy_northern_ghana_
external_edges` already read the raw Dry/Wet sheets before this fix — only
the path prefix was wrong there.

**Isolating the fix from pre-existing staleness.** A second, unrelated issue
surfaced during verification: `mechanism_resolution_scores.csv`,
`mrs_calibration_model.json`, `data_tier_experiment.csv`, and
`identifiability_diagnostics.csv` (all synthetic-tier, no field-data
dependency) differed from the committed git baseline. Stashing this
session's script changes and rerunning with `--reuse-synthetic --skip-field`
(old, pre-fix code, field stage skipped entirely) reproduced those same
differences and was byte-identical to this session's fixed-code rerun of the
same files — proving those files were already stale in the repository
before this session touched anything, unrelated to the field-data fix, and
simply got refreshed as a side effect of rerunning the pipeline. Not
otherwise investigated; flagged here rather than silently taken credit for.

**Manuscript propagation.** Edited `manuscript/sections/{00-abstract,
03-methods,04-results,05-discussion}/section.md` and `Outline.md` for the
138->160 count, the 0.681/0.941->0.664/0.942 scores, replaced every
"aquifer classes and lithologies" framing with "regions", and added an
explicit Methods/Discussion disclosure that thermodynamic/SI evidence gating
does not apply to this section. Updated `manuscript/methods/methods_
evidence.csv` (METH-15, METH-19, METH-21) and `manuscript/artifact_
registry.json` (FIG-1, FIG-6) to match. Re-ran `methods_gate_check.py`
(plan+draft), `citation_audit.py` (passes with `--bib-files manuscript/
LITERATURE.bib`; the default bib-file glob does not pick up that filename,
unrelated to this fix), `abstract_gate_check.py`, and
`assemble_manuscript.py` — all pass; `Manuscript-Final.md` was regenerated
and scanned clean of stale numbers. `manuscript/Manuscript.md` (the
pre-revision snapshot, superseded by `Manuscript-Final.md`) and
`manuscript/review/M5_manuscript_review.md` (a dated historical peer-review
report) were deliberately left untouched as frozen records; the latter's
"Numerical/traceability gap" comment about the 294-vs-138 aggregation rule
is now moot rather than retroactively edited away.

`manuscript/methods/human_review.json` did not exist before this session
(the manifest declared its path but it had never been created), which was
blocking the B4-draft gate independently of this fix. Created it as a
self-review attestation (no independent human reviewer available in this
session), following the M6/M7 precedent.

**Known pre-existing gaps, not addressed here** (out of scope for a
data-correction pass): `M5/manuscript/analysis_plan.json` does not exist and
`artifact_validation.json` does not match the B2/B3 gate's expected shape,
so `track_b_gate_check.py --stage B2/B3` fail. This mirrors M7's own
documented pattern (`M7/.../DECISIONS.md`, "registering already-locked
results, not new computation") — M5's manuscript was evidently drafted the
same retrospective way, gate artifacts and all, before this session began.
Backfilling the full B1-B3 artifact chain was not requested and was not
attempted.

**M2 has the same broken-path bug** (`data/LowerAnayari/manu.csv`,
`data/Talensi_MiningArea/talensi.csv` in `scripts/analysis/
run_m2_field_benchmarks.py` and elsewhere) but M2 was explicitly out of the
user's stated M5-M7 scope and was not touched.

## 2026-07-26 — Superseded: real PHREEQC-derived SI, not empty

The user pushed back on leaving `upstream_si` empty and asked whether
saturation indices could have been computed instead. They could, and now
are; this section supersedes the "Saturation indices" bullet above and the
corresponding "thermodynamic screening was not applied" language that was
briefly in the manuscript.

**PHREEQC was not installed in this environment** (confirmed: `find_phreeqc()`
raised `FileNotFoundError`, and the hardcoded candidate path
`C:\Program Files\USGS\phreeqc-3.7.3-15968-x64\...` did not exist). Given
explicit user approval, downloaded the official USGS PHREEQC 3.7.3-15968
Windows x64 MSI from `github.com/phreeqc-dev/phreeqc3` (an apparent rename
of `usgs-coupled/phreeqc3`: same GitHub repository ID, same release
tag/asset, confirmed via the GitHub API before downloading, and separately
confirmed by the user). Installing it failed with MSI error 1603
(`LaunchConditions` return value 3) because this session has no
administrator rights (`C:\Program Files\` is not writable; confirmed not a
member of `S-1-5-32-544`) and could not self-elevate. The user then
installed PHREEQC themselves, but got the current official release,
**3.8.6-17100**, not the exact 3.7.3-15968 build the code hardcodes (that
older build's own installer needs elevation this session doesn't have
either, so the version mismatch is a direct consequence of the earlier
elevation failure, not a substitution of convenience).

- Added the 3.8.6-17100 paths as additional fallback candidates in
  `m5_common.py`'s `find_phreeqc()` (kept the original 3.7.3-15968
  candidates too, so this remains portable to a machine with either
  version installed).
- Added `GHANA_SI_PHASES = ["Calcite", "Dolomite", "Gypsum", "Halite",
  "Fluorite"]` and `compute_ghana_field_saturation_indices()` to
  `run_m5_inverse_reaction_benchmark.py`: builds one PHREEQC `SOLUTION`
  block per wet-season borehole directly from its own measured pH,
  temperature, and major-ion panel (mirroring the existing
  `_solution_lines`/`build_phreeqc_input` pattern used for the synthetic
  benchmark), runs it through `run_phreeqc()`, and parses the
  `SELECTED_OUTPUT`/`-saturation_indices` table (verified column-naming
  convention `si_<Phase>` and row-order-matches-solution-order behaviour
  with a standalone two-solution test before wiring it into the real
  160-well batch). Cached to `results/ghana_field_saturation_indices.csv`
  and reused on subsequent `--reuse-synthetic` runs rather than
  recomputed every time.
- **Phase selection was deliberately narrowed, not copied from the retired
  Mendeley 4-phase set.** `M5/m5_common.py`'s `REACTIONS` dictionary
  actually declares 9 `si_phase` values (Albite, Anorthite, Calcite,
  Dolomite, Fluorapatite, Fluorite, Gypsum, Halite, K-feldspar), but the
  raw field panel has no Al (rules out the three aluminosilicates) and no
  PO4 (rules out fluorapatite) — confirmed by checking the actual
  `NorthernGhana.xlsx` column list. Requesting SI for phases the measured
  chemistry cannot constrain would produce meaningless values, so only the
  5 phases the panel actually supports were requested (one more,
  Fluorite, than the retired Mendeley file ever provided, since Ca and F
  are both measured).
- `run_ghana_field_demonstration()` now takes an optional
  `saturation_indices` DataFrame parameter and looks up each well's own
  computed SI by `well_id`, falling back to an empty dict only if none was
  supplied (so the function still degrades gracefully if PHREEQC is ever
  unavailable again).

**Numeric effect of switching from empty to real SI**: median
Hydrosheaf-Core evidence score 0.664 -> 0.690 (the SI-dependent scoring
terms in `hydrosheaf_core_evidence` are no longer uniformly skipped);
TDS-consistency score unchanged at 0.942 (does not depend on SI);
`thermodynamic_violations` summed to 0 across all 160 wells post-fix
(the SI bounds constrain the fit before violations can occur, not a
post-hoc coincidence — these are dilute, mostly undersaturated waters with
respect to calcite/dolomite/gypsum/halite/fluorite, consistent with the
raw EC/TDS statistics). Fluoride's held-out measurement-value ranking is
unaffected (still first by a wide margin, 2.12 vs 1.02 for the runner-up).

**Manuscript re-propagation**: reworded the Methods/Discussion sentences
that had described SI as unavailable/inactive to instead describe the
live-PHREEQC computation and its Al/PO4 limitation; updated the Results
score to 0.690; added `METH-24` to `methods_evidence.csv` and registered
`parkhurst2013phreeqc` against it in `literature_search.json`'s
`supports_method_ids` (the plan-stage gate failed once with "unmapped
citation" until this was added); rewrote Supplementary Methods S7's first
paragraph, which still described the retired Mendeley-derived aquifer/
lithology/294-record/graph-edge narrative in detail even after the
Methods-section fix earlier in this session — this was missed in the
first pass and only caught when re-auditing after the SI change. Left S7's
second paragraph (the separate "1,019 candidate Hydrosheaf-generated
edges" / hydraulic-residence-time exploratory analysis, sourced from
`m5_data_readiness_audit.md`) untouched: it is a disclosed limitation
already excluded from the reported results, does not appear to be produced
by `run_m5_inverse_reaction_benchmark.py` at all (no code path generating
graph edges was found there), and auditing whatever script did produce it
was not attempted in this session — flagged here as a possible follow-up,
not resolved. Re-ran the full Methods/citation/abstract gate suite and
`assemble_manuscript.py` after every change in this section; all pass.
Regenerated all Python, DuckDB, and R figures/tables a second time from the
SI-corrected results.

Software-version note: `tableS9_software_environment.csv`'s `PHREEQC`
row was hardcoded to the literal string `"3.7.3-15968"` in
`write_summary()`, independent of the `executable`/`database` paths
actually resolved and reported in the same table — meaning it silently
went stale (still claiming 3.7.3-15968 while the executable/database paths
below it pointed at 3.8.6-17100) until this was caught and fixed to derive
the version string from the resolved executable path by regex. The
already-locked synthetic benchmark (240 PHREEQC scenarios) was generated
under 3.7.3-15968 and was not rerun or revalidated under 3.8.6-17100 in
this session; only the newly added Northern Ghana SI computation used the
newer version. This is disclosed in `manuscript/methods/human_review.json`
rather than left implicit.

## 2026-07-27 — Saturation indices computed for Talensi and Lower Anayari

The user asked whether the manuscript claim that Talensi/Lower Anayari
lacked saturation indices "because it was needed, not because it was not
available" held up. Checked directly: `run_external_field_transfer()`'s
`_external_field_edge_outputs()` passed `upstream_si={}` unconditionally,
for every dataset including `NorthernGhana.xlsx`, in this specific
external-transfer/ELRI pathway (a separate code path from
`run_ghana_field_demonstration`, which already had real Ghana SI from the
2026-07-26 fix). `LEGACY_PRIMARY_METHOD = "thermo_elastic_net"` was already
requesting thermodynamic bounds, but with `upstream_si={}`,
`thermodynamic_bounds()` found every phase's SI missing and fell back to
unconstrained bounds — the "thermo" method name was a no-op in practice.
Talensi and Lower Anayari both measure complete pH, temperature and
major-ion panels (the only inputs SI computation needs; Sr/SiO2 are not
required), so nothing about them justified leaving SI uncomputed.

- Added `compute_external_field_saturation_indices()`, mirroring the
  existing `compute_ghana_field_saturation_indices()` PHREEQC input-deck
  pattern, for Talensi (63 samples) and Lower Anayari (41 samples)
  combined. Found and fixed two real bugs while building it: (1) one field
  value in the source CSVs is recorded as a below-detection-limit string
  (`"<0.001"`), which crashed a naive `float()` cast — fixed by parsing it
  as half the stated detection limit, matching the convention already used
  in `scripts/analysis/run_m2_field_benchmarks.py` for the same datasets;
  (2) PHREEQC's raw executable reports a `-999(.999)` sentinel for a phase
  it cannot evaluate (Talensi measures no fluoride, so fluorite SI is
  undefined for all 63 samples) — this was not being filtered and would
  have been read downstream as a genuine, extremely undersaturated SI value
  rather than "not computed"; fixed by treating any value below -900 as
  `NaN`. Verified Northern Ghana's existing SI cache has zero such sentinel
  rows (its ion panel is complete for every sample), so this bug did not
  affect any already-reported Ghana number.
- `_external_field_edge_outputs()` gained an `upstream_si` parameter (both
  its `hydrosheaf_core_evidence()` and `fit_inverse()` calls previously
  hardcoded `{}`); `run_external_field_transfer()` now builds a
  `(dataset, sample_id) -> SI dict` lookup from both the existing Ghana SI
  and the new external SI, keyed at each edge's upstream sample
  (consistent with the Ghana field demonstration's own convention of
  computing SI from the wet-season/upstream chemistry), and passes the
  real value through instead of `{}`.
- **Numeric effect**: median ELRI for `NorthernGhana.xlsx` and Lower
  Anayari both 0.072 (Lower Anayari up from 0.010); Talensi 0.147 (up from
  0.072). Edge counts (160/85/49) are unchanged. Propagated to Results
  (04-results/section.md external field-transfer paragraph) and
  Supplementary Methods (S7's external-transfer paragraph).
- All figures and tables regenerated (`make_m5_publication_tables.py`,
  `make_m5_publication_figures.py`, `r_figures/plot_m5_publication_figures.R`).
  All gates (methods draft, citation audit, abstract, assemble) and
  `pytest M5/m5_inverse_reaction_benchmark/tests/test_m5_common.py` (11/11)
  re-verified green.
- M6 had the identical issue (Talensi/Lower Anayari's SI hardcoded absent
  in `m6_common.py`, gated behind Sr/SiO2 by the tier ladder's own nested
  construction rather than by genuine data absence) and was fixed in the
  same session; see `M6/m6_field_transfer_benchmark/DECISIONS.md`.

## Missing References list fixed (2026-07-27)

The user reported that M5's manuscript had no reference list. `manuscript/
Manuscript-Final.md` never renders one by design — Pandoc citeproc only
resolves `[@key]` citations and generates the bibliography at DOCX
conversion time (verified this is also true of M6/M7, which both already
had a `Manuscript-Final.docx`). The actual gap: M5 had no
`Manuscript-Final.docx` at all, because two prerequisites were missing:

- `M5/References.bib` (the citeproc-input file the docx converter reads)
  did not exist because `citation_audit.py` had only been run without
  `--write-references-bib` this session. Re-ran with that flag.
- `manuscript/methods/apa.csl` did not exist, even though `methods_manifest.
  json` declares `"csl_path": "manuscript/methods/apa.csl"` (matching
  M6/M7's convention exactly) — a genuine setup gap specific to M5. Copied
  M6's copy of the same standard, project-independent APA 7 CSL file.
- Running `convert_to_docx.py --project-root M5 --csl manuscript/methods/
  apa.csl` then produced a DOCX with a correctly formatted bibliography but
  **no "References" heading above it**. Root cause: the r2m skill's
  `convert_to_docx.py` never passes Pandoc's `--metadata reference-section-
  title=References`, and this flag is not exposed as a script argument at
  all — confirmed by testing the equivalent raw `pandoc` invocation with
  and without that metadata flag. M6/M7's existing `Manuscript-Final.docx`
  files (both already present before this session) must have been produced
  by some other invocation that included this flag. Rather than patch the
  shared, externally-provided skill script, ran `pandoc` directly with the
  same bibliography/CSL inputs plus `--metadata reference-section-
  title=References` to produce the corrected `Manuscript-Final.docx`,
  matching M6/M7's actual output exactly (verified: 16 `ref-` bookmarks,
  full APA-formatted list under a "References" heading).
- `methods_gate_check.py --stage docx` re-verified PASS.
