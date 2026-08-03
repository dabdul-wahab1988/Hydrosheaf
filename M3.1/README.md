# M3.1 manuscript package

Version upgrade of `M3` (accuracy-locked 2026-07-28, commit `47fbf7c`,
submission package in `M3/M3_geochemistry/`). M3 is retained unchanged as the
prior locked record; M3.1 supersedes it once regeneration and review are
complete. See `DECISIONS.md` for the full audit trail of what changed in the
codebase since the lock and why an upgrade was warranted.

## Working title

Unchanged from M3: **Conditional Value of Graph Priors in Multi-Tracer
Groundwater-Age Inference: A Hydrosheaf Benchmark of Public Aquifer Datasets**

## What is different from M3

1. Every reportable statistic is regenerated from the current codebase (not
   inherited from the M3 lock) and diffed against it; see `DECISIONS.md` D2 for
   the outcome.
2. A new site map (FIG-1) shows the geographic distribution of the USGS
   public-supply sites used — M3 had no map figure.
3. A compact new subsection reports an independent set-valued compatibility
   audit (agrees with the main graph-regularisation finding by a different
   method) and a 27.85% tracer-level local-infeasibility rate that seven
   tested explanations failed to account for. Both are explicitly scoped as
   exploratory/development-stage, never as confirmatory or field-validated.
4. The main text is rewritten in reader-friendly language at a 6,000-word
   budget (see `Outline.md`); full mathematics, derivations and the
   seven-explanation exclusion protocol move to the Supplementary Information,
   which carries no length cap.
5. All figures and maps are produced in R from regenerated CSV exports,
   following the `M2.3/analysis/r/` convention in this repository
   (`_theme.R`, `_map.R`, one script per figure), rather than the mixed
   Python/matplotlib pipeline M3 used.

## Layout

```
M3.1/
  DECISIONS.md               locked decisions and audit trail
  Outline.md                 positioning, argument spine, word budget
  proposal.normalized.json   objectives, inputs, risks
  m3_age_benchmark/          computational authority (copied from M3, rerun
                              against current HEAD; see its own README.md)
  analysis/
    r/                       figure/map authority; consumes CSV exports only
      _theme.R
      _map.R                 US study-area map for Figure 1
      fetch_boundaries.R     one-time boundary fetch (only network step)
      fig01..fig07
      make_all_figures.R
  manuscript/
    Manuscript-M3.1.md        assembled main manuscript
    sections/                 authoritative section sources
    supplementary/            Supplementary Information (methods, math, figs)
    artifacts/
      data/                   tidy CSV exports (the evidence interface)
      figures/                PDF and PNG at 600 dpi
      TAB-*.md                generated tables
    methods/apa.csl            APA 7 CSL (shared with M2.3)
    LITERATURE.bib
```

## Reproduction

From the repository root:

```powershell
.venv\Scripts\python.exe M3.1\m3_age_benchmark\scripts\run_m3_manuscript_analysis.py --full --age-steps 90
Rscript M3.1/analysis/r/fetch_boundaries.R        # once; needs network
Rscript M3.1/analysis/r/make_all_figures.R
```

Python owns computation inside `m3_age_benchmark/` and writes the locked,
hash-manifested results. R consumes exported CSVs under
`manuscript/artifacts/data/` and owns every figure and map; no figure
recomputes a reported statistic.

## Status

Analysis regeneration in progress. See `DECISIONS.md` D2 for the numerical
comparison against the M3 lock once the rerun completes.
