# Track B B3 gate

- Status: **FAIL**
- Purpose: Register final artifacts

| Check | Result | Detail |
|---|---|---|
| artifact_registry.json contains records | PASS | count=7 |
| Registry matches the locked plan | FAIL | missing=[] extra=['FIG-1', 'FIG-2', 'FIG-3', 'FIG-4', 'FIG-5', 'FIG-6', 'TAB-1'] |
| Registry descriptions satisfy contract | FAIL | caption + production + interpretation |
| Registered artifact paths are contained and exist | PASS | ok |
