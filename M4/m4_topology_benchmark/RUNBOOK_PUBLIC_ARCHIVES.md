# Runbook: Public Archive Validation Pipeline

Follow these commands to run and verify the validation pipeline:

## Step 1: Inspect Archive Integrity

To scan all configurations and log readiness:
```powershell
python M4/m4_topology_benchmark/scripts/inspect_public_archive.py
```

To scan a single configuration:
```powershell
python M4/m4_topology_benchmark/scripts/inspect_public_archive.py --config M4/m4_topology_benchmark/configs/savage.yaml
```

---

## Step 2: Run the Ingestion & Validation Pipeline

To run the validation pipeline on a single config (e.g. Great Miami):
```powershell
python M4/m4_topology_benchmark/scripts/run_public_archive_pipeline.py --config M4/m4_topology_benchmark/configs/great_miami.yaml
```

To run the full orchestrated sequence (Savage, Great Miami, and Long Island) and merge summaries:
```powershell
python M4/m4_topology_benchmark/scripts/run_m4_external_modpath_archive_validation.py
```

---

## Step 3: Regenerate Manuscript Figures & Tables

To run the entire M4 analysis suite, which runs all validation stages and regenerates all figures and tables:
```powershell
python M4/m4_topology_benchmark/scripts/run_m4_manuscript_analysis.py
```

Check outputs in:
- `M4/m4_topology_benchmark/tables/Manuscript_Ready/`
- `M4/m4_topology_benchmark/figures/Manuscript_Ready/`

---

## Step 4: Run the Test Suites

To execute the unit tests verifying pipeline logic, reference graph analysis, and scientific claim boundaries:
```powershell
pytest tests/test_m4_public_archive_pipeline.py tests/test_modpath_archive_reference_graph.py tests/test_m4_claim_guardrails.py
```
