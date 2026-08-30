# M2_1 manuscript package

This is the new, evidence-locked M2_1 software and validation manuscript package.
It is deliberately separate from the legacy M2 Word manuscript, which reports an
earlier benchmark and should not be silently overwritten.

## Working title

**HydroSheaf: an auditable evidence-integration framework for groundwater age,
connectivity and reaction inference with controlled-synthetic validation**

## Scope

M2_1 explains how the HydroSheaf software was assembled and what its current
implementation can support. It reports the current bounded controlled-synthetic
programme run `RUN-INTEGRATION-FULL-20260802-15` and a conservative application to
the real data under `data/FieldData`. The synthetic result is not field validation,
and the field application is not treated as truth for age, directed connectivity,
or unique reaction mechanisms.

## Word contract

- Main manuscript: 6,000 prose words including the abstract and statements, and
  excluding headings, references, tables, figure captions and this README.
- Supplementary Methods: 3,000 prose words, excluding headings, references,
  tables and equations.
- UK English, past tense, non-personal scientific voice, APA 7 citation keys.

The authoritative section files are under `manuscript/sections/`. The assembled
manuscript is `manuscript/Manuscript-Final.md`; the supplementary methods are
`manuscript/supplementary/Supplementary-Methods.md`.

## Evidence boundary

The primary software claim is bounded: HydroSheaf passed the predeclared age,
reaction, kinetics and integrated decision gates on the locked synthetic test
cases in the current programme run. The result supports competence and decision
utility under the tested generators and observation regimes. It does not establish
universal superiority, field effectiveness, or validity outside the tested
synthetic domain.

The field data support data-readiness, seasonal chemistry prediction,
reaction-family plausibility, measurement-value and edge-sensitivity diagnostics.
They do not contain environmental age tracers, screen-resolved intervals,
repeated hydraulic-head series, independently verified directed connectivity, or
field reaction-mechanism truth.

## Review and reproducibility

The package follows the Track-B r2m contract. The final review report is kept under
`manuscript/review/` and is produced with the manuscript-reviewer skill after the
first complete draft. Every issue classified as critical or major must be resolved
or explicitly retained as a claim limitation before the package is finalised.

The synthetic source snapshot is identified by SHA-256 hashes in the locked run
manifest. The manifest records revision `8718d669363bc4955cba91ab729108c7f604e610`
and a dirty-worktree flag. The current checkout is commit
`2d4b8af4b31cb0673ad5067c6169cba14f9a53bd`; therefore the run is traceable, but it
was not generated from the current clean commit and the manuscript does not imply
clean-commit regeneration.
