# Reader-facing terminology audit

## Purpose

The manuscript and supplementary methods were checked as if they were being read without access to the project workspace, Codex, milestone folders or internal run logs. Internal implementation labels were removed from the reader-facing text. They remain available in the software repository and machine-readable provenance records where they are useful for reproduction.

## Findings and resolution

Internal milestone labels such as M2 through M8 and the manuscript working label M2_1 were not meaningful to an external reader. They were replaced with terms such as “earlier analyses,” “historical analysis records,” “this study” and “the present study.”

Repository-specific paths and names, including `.codex_work`, package-directory paths, benchmark-script paths, the field-data root and the assembly-script command, were removed from the manuscript and supplement. The text now describes reusable package code, control scripts, machine-readable evidence, and the field-data archive in ordinary scientific language.

The internal run identifier, raw revision hashes and the machine flag `git_worktree_dirty` were removed from the prose. The manuscript now states that the provenance manifest recorded the generating revision, cryptographic hashes and whether source changes had been committed. The exact values remain in the provenance artefacts.

Machine status strings such as `MODEL_DISAGREEMENT`, `FIT_NOT_CONVERGED`, `NON_IDENTIFIABLE` and `ABSTAIN_NO_SUPERIORITY_GATE` were removed. PASS, WEAK and ABSTAIN are defined in the Introduction, and the integrated result is reported as a bounded PASS with the general-superiority claim withheld.

Technical shorthand was checked. RAPM is identified as project-specific regularised adjusted plus-minus notation; KKT is expanded as the Karush–Kuhn–Tucker condition; R² is defined as the coefficient of determination; RMSE, MAE and expected calibration error are written out in the main performance table; and k×A is explained as rate-constant and surface-area confounding.

## Verification

The final assembled manuscript and supplementary methods contain no matches for internal milestone labels, Codex paths, raw run identifiers, machine status strings or repository-specific source paths. The final word counts are 5,974 for the main manuscript and 2,490 for the supplementary methods, both within their declared limits.
