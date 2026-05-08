# E5 DGMETA Dissolved-Gas Validation

Run timestamp UTC: 2026-05-06T04:27:38.886179+00:00

Source: official USGS DGMETA example workbooks from the USGS GitLab repository.

Run scope: all parsed rows.

## Scope

Hydrosheaf reads cached DGMETA `*_ModOut` outputs and compares recharge temperature/excess-air solutions against Hydrosheaf's dissolved-gas fitter using the same example dissolved-gas observations. Row-specific DGMETA solubility values are passed into Hydrosheaf so this validation isolates fitting/correction mechanics from Hydrosheaf's compact fallback solubility approximations.

## Summary

| group   |   n_rows |   temperature_mae_c |   temperature_rmse_c |   excess_air_mae_cm3kg |   median_hydrosheaf_rmse_standardized |   large_residual_fraction |
|:--------|---------:|--------------------:|---------------------:|-----------------------:|--------------------------------------:|--------------------------:|
| all     |      409 |            0.508409 |              1.235   |                2.4403  |                              0.36964  |                 0.0537897 |
| CE      |      114 |            0.530702 |              1.17447 |                2.50202 |                              0.201412 |                 0.0175439 |
| UA      |      295 |            0.499794 |              1.25761 |                2.31237 |                              0.41377  |                 0.0677966 |

## Interpretation Guardrail

Agreement in this validation supports Hydrosheaf's dissolved-gas workflow mechanics. Exact DGMETA parity still requires direct verification against DGMETA's full coefficient tables, macros, and all workbook options.

Hydrosheaf now uses the standard closed-system-equilibration gas-volume equation for CE rows. Exact DGMETA parity still depends on the full DGMETA coefficient tables, workbook option handling, and independent verification of every macro pathway.
