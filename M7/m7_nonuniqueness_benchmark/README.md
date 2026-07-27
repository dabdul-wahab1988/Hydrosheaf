# M7 — Conditional integration and non-uniqueness

This is the single maintained M7 code package for the Q1 manuscript. The
superseded original M7/M7.1 implementation has been removed. The independent
generator and confirmatory validation formerly developed as M7.2 are retained
here only because they are direct dependencies of the final benchmark and its
supporting evidence.

M7.3 tests the thesis claim that integration is useful when independent
evidence is complementary, but can be redundant or harmful when an evidence
stream is weak or misspecified. It does not require every additional stream to
increase point accuracy.

The benchmark has four linked audits:

1. Fresh independent MODFLOW 6/MODPATH 7 cases compare hydraulic, age and
   chemistry evidence alone, in pairs and together.
2. Correct, partial and reversed flow graphs are imposed on the same local
   tracer likelihoods to quantify how connectivity assumptions change age
   uncertainty and accuracy.
3. Reaction-family recovery is bootstrapped under core and enhanced chemistry
   panels, with special reporting for carbonate processes.
4. The Northern Ghana workbook is audited as a data-limited application:
   supported component diagnostics are separated from unavailable field truth.

The final locked protocol is in
[docs/m7_3_protocol.md](docs/m7_3_protocol.md). Results must be written only
after the protocol and runner have been committed.

The locked interpretation and decision table are in
[docs/m7_3_results.md](docs/m7_3_results.md).

## Package layout

- `scripts/run_m7_3_nonuniqueness.py`: final conditional-integration benchmark.
- `scripts/independent_modflow_generator.py`: truth generator used by the final
  benchmark; it intentionally imports no HydroSheaf code.
- `scripts/strong_inference.py`: blind HydroSheaf inference adapter.
- `scripts/run_supporting_validation.py`: confirmatory supporting experiment,
  retained for reproducibility rather than as a separate M7 version.
- `results/m7_3_locked`: final locked M7 results.
- `results/supporting_validation`: confirmatory evidence used by M7 and M6.

The supporting experiment's
[protocol](docs/supporting_validation_protocol.md),
[fresh-seed amendment](docs/supporting_validation_amendment.md), and
[results](docs/supporting_validation_results.md) remain as the audit trail for
the negative incremental-age finding.

## Replay the final benchmark

Place `mf6.exe` and `mp7.exe` in `.codex_work/modflow-bin`, then run:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_m7_3_nonuniqueness.py
```

A non-claim-bearing smoke run is available with `--quick`.

## Replay the supporting validation

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_supporting_validation.py --confirmatory
```

To regenerate Supplementary Figure S1 for representative locked realization
4101:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' M7\m7_nonuniqueness_benchmark\scripts\make_supporting_model_domain_map.R
```

This figure uses synthetic model-space metres. It is not georeferenced and
must not be presented on a Ghana or other geographic basemap.
