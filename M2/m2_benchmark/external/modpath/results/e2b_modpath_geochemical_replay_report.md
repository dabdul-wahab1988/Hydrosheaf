# E2b MODPATH-Conditioned Geochemical Replay Report

Run timestamp UTC: 2026-05-05T10:27:00.959587+00:00

Physical topology source: USGS Savage MODPATH archive, DOI `10.5066/F7J102FK`.

## Scope

This is a semi-synthetic integration validation. The graph topology and travel-time priors are real MODPATH outputs from E2. The geochemical observations are controlled known-truth injections placed on those real edges because the Savage archive does not provide paired chemistry/tracer samples for the same graph nodes.

## Outputs

- Replay results: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\modpath\results\modpath_geochemical_replay.csv`
- Replay summary: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\modpath\results\modpath_geochemical_replay_summary.csv`
- Figure S2B: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\figures\figure_s2b_modpath_conditioned_geochemical_replay.png`

## Summary

|   n_edges |   gamma_mae |   mean_extent_mae |   max_extent_mae |   reaction_tp |   reaction_fp |   reaction_fn |   reaction_precision |   reaction_recall |   reaction_f1 |   median_objective_score |
|----------:|------------:|------------------:|-----------------:|--------------:|--------------:|--------------:|---------------------:|------------------:|--------------:|-------------------------:|
|        80 |  0.00458189 |        0.00246646 |       0.00627001 |           240 |             0 |             0 |                    1 |                 1 |             1 |              3.88846e-05 |

## Interpretation

E2b closes the integration gap left by E2: it tests whether Hydrosheaf can recover known geochemical transport and reaction signals when conditioned on a real MODPATH-derived physical graph. It does not replace a future fully paired physical-geochemical field dataset, but it prevents the M2 paper from implying that the MODPATH archive alone validates chemistry.
