# M2_1 decisions and claim ledger

## Locked decisions

1. M2_1 is a new manuscript package. The legacy `M2/M2_ready/Manuscript-CG_updated.docx`
   remains unchanged and is not silently treated as the current manuscript.
2. The main-text limit is 6,000 prose words including the abstract and statements.
   Supplementary Methods is capped at 3,000 prose words.
3. The current synthetic evidence is identified by
   `RUN-INTEGRATION-FULL-20260802-15`, with its source SHA-256 hashes, recorded
   revision `8718d669363bc4955cba91ab729108c7f604e610` and dirty-worktree flag.
   The current checkout is `2d4b8af4b31cb0673ad5067c6169cba14f9a53bd`; a clean
   commit-level regeneration remains a final reproducibility action.
4. The current field source is `data/FieldData`; the actual Northern Ghana folder
   is `data/FieldData/NorthenGhana` and the spelling is retained for provenance.
5. Field results are application and scope diagnostics, never independent truth
   for groundwater age, directed connectivity, or unique reaction mechanism.
6. “Sheaf-style” or “sheaf-inspired” is used unless a formal sheaf construction is
   demonstrated in a separate mathematical manuscript.
7. Claims of superiority are excluded. The synthetic programme supports bounded
   competence, calibration and prospective utility under its tested domain.
8. Existing synthetic metrics are not silently pooled with the older M2 100-run
   benchmark. The run ID and protocol are stated beside every current result.

## Current claim ledger

| ID | Safe claim | Evidence | Status |
|---|---|---|---|
| C1 | HydroSheaf contains modular topology, age, reaction, kinetic, calibration, discrepancy, averaging, abstention and decision-utility components. | Package source and exports | PASS |
| C2 | Age output passes the bounded locked synthetic gate with calibrated intervals and lower specialist MAE than its competence-matched baseline. | Performance gate, run 20260802-15 | PASS, bounded |
| C3 | Reaction-family output passes a bounded calibrated log-loss/coverage/selective-risk gate. | Performance gate, run 20260802-15 | PASS, bounded |
| C4 | Kinetic prediction passes the bounded predictive gate, while parameter identification remains conditional on independent surface-area information. | Performance gate, run 20260802-15 | PASS, conditional |
| C5 | Integrated decision selection passes the bounded prospective synthetic utility gate against random and is non-inferior to the strongest specialist comparator under the locked rule. | Performance gate, run 20260802-15 | PASS, bounded |
| C6 | Ghana FieldData support chemistry/readiness/identifiability and measurement-value diagnostics. | FieldData audit and M6 synthesis | PASS, field scope |
| C7 | Current evidence establishes field effectiveness or universal superiority. | No independent field truth or universal benchmark | ABSTAIN |
