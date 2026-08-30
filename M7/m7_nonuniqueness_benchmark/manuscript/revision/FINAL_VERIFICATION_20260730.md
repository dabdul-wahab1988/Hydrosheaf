# M7 final reviewer-revision verification, 30 July 2026

## Outcome

All 36 reviewer comments have a verbatim point-by-point response and a routed
change record. Thirty-three are fully addressed. Three are technically or
structurally corrected but remain partial because the final versioned public
release, persistent identifiers and author declarations do not yet exist.

No locked confirmatory test was rerun during the revision.

## Scientific and statistical checks

- Exact local-first/global-fallback estimator-diagnostic runner SHA-256:
  `a0ef13bde5391af62698927211cb4e701123affebb108d331795ce8596e2e191`.
  This matches the confirmatory lock and content-addressed source archive.
- Post-review audit status: PASS; 10,000 shared case-block resamples, complete
  120-contrast representation-benchmark and 560-contrast estimator-diagnostic
  families, and 20,000 planning simulations. The audit does not regenerate
  locked tests.
- Focused benchmark tests: 14 passed in 2.59 seconds.
- Main manuscript: seven figures and four tables.
- Supplement: thirteen tables.
- Citation audit: PASS.
- Abstract gate: PASS at 329 whitespace-delimited words (maximum 350).
- Methods plan, draft and DOCX gates: PASS.
- Track B artifact gates B1-B3: PASS.
- Track D revision gates D1-D6: PASS.
- `git diff --check` for the M7 package: PASS, with line-ending warnings only.

## Document interoperability and visual QA

- Microsoft Word: complete 33-page main PDF and 22-page supplement PDF.
- LibreOffice: complete 31-page main PDF and 21-page supplement PDF.
- The former LibreOffice `libpng error: Write Error` did not recur.
- All 55 Word-rendered manuscript and supplement pages were visually inspected.
- The 24-page second-round review was also rendered in Word and LibreOffice and
  visually inspected.

## Reader-facing terminology verification

- Internal labels `M7.2` through `M7.5` were removed from the manuscript,
  supplement, publication tables, figure titles and captions.
- They were replaced by self-describing scientific terms: independent
  supporting-validation experiment, process-based integration benchmark,
  competence-matched sheaf-versus-weighted-graph representation benchmark, and
  local-first/global-fallback estimator diagnostic.
- Exact immutable run identifiers retain `M7` only where required for
  reproducibility and provenance.
- A plain-text audit of both final DOCX files found no `M7.2`-`M7.5` labels.
- The rebuilt LibreOffice render contains 31 main-manuscript pages and 21
  supplementary pages; every page of the latest render was visually checked.
- The dated and canonical main-manuscript DOCX files are byte-identical at
  SHA-256 `3900e499f1abb20380f254a99af421fef19d5094ad341ae7d0879555d7f2bcbe`;
  the corresponding supplementary files are byte-identical at SHA-256
  `43c77d784ee4f7047a98f8236045a28eb1ae3c01e02d1adf4eab2752600c0229`.

Artifact hashes and exact QA paths are recorded in
`manuscript/revision/docx_interoperability_qa.json`.

## Remaining submission actions

1. Obtain final CRediT author roles, funding and competing-interest statements.
2. Commit and tag the complete representation-benchmark and
   estimator-diagnostic technical package.
3. Deposit the release and data package, mint persistent software and dataset
   identifiers, and replace the explicit pending statements in the manuscript.
4. Verify the public release from a clean checkout before journal submission.

These actions require author declarations or external repository publication;
they were not fabricated or performed as part of this technical revision.
