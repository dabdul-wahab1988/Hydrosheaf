# M2.3 decisions and claim ledger

M2.3 supersedes the legacy `M2` Word package and the `M2_1` markdown package as
the submission-track manuscript. Neither earlier package is deleted; both are
retained as historical records.

## Locked decisions

### D1. Package status

M2.3 is a single merged manuscript targeting Computers & Geosciences. It replaces
the plan to submit `M2` and `M2_1` as two papers, which would have duplicated the
software description and reported the same Ghanaian result twice.

### D2. Northern Ghana data provenance (user decision, 2026-08-02)

Integrity checks on `data/FieldData/NorthenGhana/NorthernGhana.xlsx` established
that the static water level is identical in 160 of 160 wells between the Dry and
Wet sheets, that 158-159 of 160 wells decrease in every conservative solute, that
no sampling dates, DOI, laboratory or source publication is embedded, and that the
referenced `SI.pdf` provenance document does not exist in the repository.

The author confirmed that **the underlying chemistry is real but the dry/wet
seasonal separation was reconstructed rather than independently sampled**.

Consequences, all binding:

1. The measured unit is **160 wells**, not 320 seasonal samples. The workbook is
   described as 160 wells presented as two reconstructed seasonal panels.
2. The **Dry** sheet is used as the primary measured panel for the reader-facing
   data description. The Wet sheet is reported as a reconstructed counterpart and
   is excluded from every inferential result.
3. **No seasonal-change finding may be reported.** Wet-dry contrasts are
   properties of the reconstruction, not of the aquifer. The seasonal
   hold-forward prediction diagnostic carried in the earlier packages is
   withdrawn.
4. **No repeated hydraulic-head claim.** Static water level is not season
   resolved.
5. The citation `@anku2008ghana` is **not** a source for this workbook. It may be
   cited for regional hydrogeological context only. The workbook itself is cited
   as an author-compiled dataset with a reconstructed seasonal attribute.
6. The reconstruction is stated in the Data section, in the limitations, and in
   the data-availability statement. It is never implied to be measured
   seasonality.

### D3. Real-data application rests on Lower Anayari and Talensi

These two datasets trace to published field studies, show the analytical
heterogeneity expected of measured data, and carry no reconstructed attribute.
They carry the field-application claim. Northern Ghana contributes a compiled
regional chemistry panel.

### D4. USGS public-age benchmark reporting (user decision, 2026-08-02)

Three different summaries of this benchmark were in circulation across earlier
packages. M2.3 reports the **held-out uncalibrated strict-parity** comparison as
the primary result, and reports the calibrated result beside it explicitly labelled
as emulation of the reference outputs rather than independent agreement. The
earlier abstract values that mixed the two are not reproduced.

### D5. Evidence recomputation

No number is copied from `M2`, `M2_1` or any earlier summary document. Every
reported value is recomputed from primary run records and primary data files by
scripts under `M2.3/analysis/python/`. Where a recomputed value disagrees with an
earlier manuscript, the recomputed value stands and the disagreement is recorded.

### D6. Computation and figure authority

Python owns computation and writes read-only tidy CSV exports. R consumes those
exports and owns all figures. No figure recomputes a reported statistic.

### D7. Terminology

"Sheaf-inspired" or "sheaf-style" is used throughout. No formal sheaf theorem is
claimed. RAPM is identified as project-specific notation.

### D8. Claims excluded by construction

Field validation of groundwater age, directed connectivity or unique reaction
mechanism; universal superiority over MODFLOW, MODPATH, TracerLPM or PHREEQC;
measured seasonality in the Northern Ghana workbook; and formal sheaf theory.

## Claim ledger

| ID | Claim | Evidence | Status |
|---|---|---|---|
| C1 | The package implements topology, age, reaction, kinetic, calibration, discrepancy, averaging, abstention and decision-utility components. | Package source and exports | PASS |
| C2 | Age inference passes a bounded locked synthetic gate with calibrated intervals and lower specialist error than a competence-matched baseline. | Locked run performance gate | PASS, bounded |
| C3 | Reaction-family inference passes a bounded calibrated log-loss, coverage and selective-risk gate. | Locked run performance gate | PASS, bounded |
| C4 | Kinetic prediction passes; parameter identification is conditional on independent surface-area information. | Locked run performance gate | PASS, conditional |
| C5 | Integrated decision selection passes a bounded prospective cost-adjusted utility rule against random. | Locked run performance gate | PASS, bounded |
| C6 | No-prior topology inference against a particle-tracking reference is screening-level, with recall far exceeding precision. | Independent topology comparison | PASS, screening only |
| C7 | Agreement with published lumped-parameter age outputs is screening-level on held-out sites. | Held-out uncalibrated parity | PASS, screening only |
| C8 | Lower Anayari and Talensi support candidate screening and failure-mode diagnostics under sparse and analytically heterogeneous input. | Field application | PASS, field scope |
| C9 | Northern Ghana supplies a compiled regional chemistry panel; its seasonal attribute is reconstructed. | Data audit, author statement | PASS, with stated reconstruction |
| C10 | Field effectiveness or universal superiority is established. | No independent field truth | ABSTAIN |
| C11 | Measured seasonal change is established for Northern Ghana. | Seasonal split is reconstructed | ABSTAIN |
