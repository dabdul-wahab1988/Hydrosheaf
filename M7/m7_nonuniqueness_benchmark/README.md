# M7.3 — Conditional integration and non-uniqueness

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

The locked protocol is in
[docs/m7_3_protocol.md](docs/m7_3_protocol.md). Results must be written only
after the protocol and runner have been committed.

The locked interpretation and decision table are in
[docs/m7_3_results.md](docs/m7_3_results.md).

## Replay

Place `mf6.exe` and `mp7.exe` in `.codex_work/modflow-bin`, then run:

```powershell
.venv\Scripts\python.exe M7\m7_nonuniqueness_benchmark\scripts\run_m7_3_nonuniqueness.py
```

A non-claim-bearing smoke run is available with `--quick`.
