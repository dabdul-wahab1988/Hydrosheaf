# M2.3 manuscript package

Merged, submission-track manuscript for **Computers & Geosciences**. M2.3
supersedes the legacy `M2` Word package and the `M2_1` markdown package, both of
which are retained unchanged as historical records and neither of which should be
submitted.

## Working title

**HydroSheaf: a claim-gated evidence-integration framework for groundwater
connectivity, residence-time and reaction inference in data-limited aquifers**

## What is different from M2 and M2_1

M2 and M2_1 described the same software and reported the same Ghanaian result,
and could not both have been submitted. M2.3 merges them and rebuilds the
evidence from primary records. Three things changed materially.

Every reported value was recomputed from primary run records and primary data
files. Eight quantities disagreed with earlier internal reporting and are
disclosed in Table 4, including a reaction-recovery coefficient of determination
that earlier reporting gave as 0.57 and that recomputes to approximately zero.

The Northern Ghana workbook was found to carry a reconstructed rather than
measured seasonal split: static water level is identical in all 160 wells across
the Dry and Wet sheets. The author confirmed the chemistry is measured and the
seasonal attribute constructed. The measured unit is now the well, the Wet panel
is excluded from inference, and no seasonal or repeated-head result is reported.

The paper's central claim is now structural rather than empirical. The
stoichiometry matrix of the reaction dictionary has rank 8 against 14 unknown
terms, leaving a six-dimensional null space, so point reaction extents are
non-identifiable from major-ion chemistry regardless of solver or penalty.

## Layout

```
M2.3/
  DECISIONS.md              locked decisions and claim ledger
  Outline.md                positioning, argument spine, word budget
  proposal.normalized.json  objectives, inputs, risks
  analysis/
    python/                 computation authority; writes read-only exports
      derive_field_evidence.py
      derive_model_evidence.py
      derive_identifiability.py
      build_manuscript.py
    r/                      figure authority; consumes exports only
      _theme.R              shared theme and validated palette
      _map.R                study-area map for Figure 2(a)
      fetch_boundaries.R    one-time boundary fetch (only network step)
      check_coordinates.R   point-in-polygon gate
      fig01..fig07          one script per figure
      make_all_figures.R
      data/                 cached boundary GeoPackage
  manuscript/
    Manuscript-M2.3.md      assembled manuscript
    sections/               authoritative section sources
    supplementary/          Supplementary Methods
    artifacts/
      data/                 tidy CSV exports (the evidence interface)
      figures/              PDF and PNG at 600 dpi
      TAB-*.md              generated tables
    review/                 reviewer report and resolution ledger
    LITERATURE.bib          28 entries, all cited, none uncited
```

## Reproduction

From the repository root:

```
Rscript M2.3/analysis/r/fetch_boundaries.R        # once; needs network
.venv/Scripts/python.exe M2.3/analysis/python/derive_field_evidence.py
.venv/Scripts/python.exe M2.3/analysis/python/derive_model_evidence.py
.venv/Scripts/python.exe M2.3/analysis/python/derive_identifiability.py
Rscript M2.3/analysis/r/check_coordinates.R      # gate: point-in-polygon
Rscript M2.3/analysis/r/make_all_figures.R
.venv/Scripts/python.exe M2.3/analysis/python/build_manuscript.py
```

`check_coordinates.R` exits non-zero if any sample falls more than 5 km outside
the Ghana national boundary. It exists because a rendered map is not a
coordinate check: the Talensi longitudes were stored without their negative sign
and plotted outside the country without that being caught by eye.

`fetch_boundaries.R` is the only step that touches the network. It pulls Ghana's
national and first-level administrative boundaries plus three neighbouring
countries from the geoBoundaries database (William & Mary, CC BY 4.0) and writes
`analysis/r/data/study_area_boundaries.gpkg`. Once that file exists, every other
step runs offline. R packages `sf`, `ggspatial` and `geobounds` come from the
project-local `.r-lib` library, which `_theme.R` adds to the library path.

Python owns computation and writes the CSV exports under
`manuscript/artifacts/data/`. R consumes those exports and recomputes no
reported statistic. No number is typed by hand into the manuscript; tables are
generated and figure values are read from the same exports the prose cites.

## Review status

An adversarial review is in `manuscript/review/manuscript_reviewer_report.md`
(score 68/100, return for moderate revisions). Four issues were classified
Critical; three were resolved with data in hand and are documented in
`manuscript/review/reviewer_issue_resolution.md`. The resolution of the third
strengthened the paper's central claim rather than weakening it.

## Outstanding before submission

The repository DOI and URL in the availability statement are placeholders.
Computers & Geosciences requires a live, documented, public repository at
submission and treats a private or incomplete one as grounds for rejection; this
is a blocker, not a formality.

Also outstanding: confirmation of the Talensi longitude sign and charge balance
against the source publication; clean-tree regeneration of the locked programme
run and a hash comparison; named CRediT author roles; expansion of the bibliography toward post-2022 work; a
statement of what was predicted before the recovery benchmark was run; and
confirmation of redistribution terms for the two published field datasets.
