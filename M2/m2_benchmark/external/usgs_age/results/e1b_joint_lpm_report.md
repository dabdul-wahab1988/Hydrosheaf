# E1b Joint-LPM Validation Report

Run timestamp: 2026-05-05T11:56:46Z

## Metrics

| Metric | Value |
| :--- | :--- |
| Samples Evaluated | 69 |
| Log10 Age RMSE | 0.9562 |
| Median Log10 Bias | nan |
| Model Family Match Rate | 10.14% |

## Accuracy by Model Family

| ref_model   |   count |       mean |        std |
|:------------|--------:|-----------:|-----------:|
| BMM-DM-DM   |      12 |  0.349781  |   0.825567 |
| BMM-EMM-DM  |       1 | -2.21304   | nan        |
| BMM-PEM-PEM |       1 |  1.43378   | nan        |
| DM          |      49 |  0.0248317 |   0.946728 |
| EPM         |       3 |  0.343583  |   0.716933 |
| PEM         |       2 | -0.77501   |   0.790571 |

## Interpretation
The Joint-LPM validation tests Hydrosheaf's ability to jointly fit multiple tracers using standard lumped parameter model families, compared to official USGS TracerLPM fits.
