# M7 post-refactor scientific audit

Date: 2026-07-29

## Decision

**PASS for a bounded scientific-workflow claim; FAIL for general sheaf
superiority.** The revised manuscript now answers what the sheaf layer adds
beyond an ordinary weighted graph with a prospectively locked,
competence-matched comparison. It keeps system execution, representation
value, predictive performance, and field validation as separate claims.

## Claim-bearing evidence

- The public pipeline completed every stage on six cases, but full-sheaf
  PR-AUC was lower than hydraulic-only PR-AUC by 0.0197 (95% CI -0.0355 to
  -0.0039), so execution acceptance is not described as overall superiority.
- In the M7.4 identity limit, native affine-sheaf and identity-graph residuals
  and predictions were exactly equal.
- Across 64 locked cases, affine-sheaf PR-AUC exceeded the identity graph
  Laplacian by 0.0854 (95% CI 0.0666-0.1050) and exceeded the permuted-map
  control by 0.0909 (95% CI 0.0705-0.1117).
- Against the strong edge-local graph, overall PR-AUC differed by 0.0097 (95%
  CI -0.0054-0.0248), Brier score by 0.0005 (95% CI -0.0033-0.0044), and log
  loss was worse by 0.0117 (95% CI 0.0008-0.0232).
- In planted incompatible-cycle cases, the sheaf improved conflict-
  localisation PR-AUC over the edge-local graph by 0.0689 (95% CI
  0.0466-0.0914).

## Defensible thesis statement

Under static scalar controlled-affine conditions, non-identity restriction
maps and a global section add explicit relation composition and conflict
localisation beyond an identity weighted graph. They do not establish general
predictive superiority over a strong edge-local weighted graph. Temporal,
three-dimensional, vadose-zone and independent field-validity claims remain
outside scope.

## Manuscript coverage checked

The title, Abstract, Introduction, Methods, Results, Discussion, Conclusion,
Open Research statement, Supplementary Methods, Supplementary Tables,
literature ledger, evidence ledger, outline, analysis plan, artifact registry,
artifact validation and generation entrypoint were revised to the same claim
boundary. Figure 6, Table 8, Supplementary Tables S7-S8, the protocol lock,
immutable first claim-bearing run, and replay commands provide the associated
audit trail. The generated 29-page main manuscript and 14-page supplement were
rendered through LibreOffice and inspected page by page. Complete supplementary
CSVs remain authoritative; the printable supplement uses compact views so that
per-case and per-edge audit tables remain legible without weakening provenance.
