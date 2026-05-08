# E4b Public Groundwater-Chemistry Pilot Report

Run timestamp: 2026-05-06T04:24:17Z

Source: USGS NAWQA groundwater-quality public data release, DOI `10.5066/P9J7I9DH`.
Selected network: `vpdcpas1`.
Samples fitted: 30.
Candidate edge fits: 150.

## Top 20 Inferred Connectivity

| edge_id                  | transport_model   |   objective_score |   chemistry_r2 |
|:-------------------------|:------------------|------------------:|---------------:|
| VPDCPAS1-61->VPDCPAS1-30 | evap              |          0.501491 |       0.99981  |
| VPDCPAS1-45->VPDCPAS1-46 | evap              |          0.505006 |       0.999222 |
| VPDCPAS1-45->VPDCPAS1-15 | evap              |          0.507973 |       0.999304 |
| VPDCPAS1-59->VPDCPAS1-30 | evap              |          0.508217 |       0.99897  |
| VPDCPAS1-52->VPDCPAS1-53 | evap              |          0.508381 |       0.997557 |
| VPDCPAS1-61->VPDCPAS1-59 | evap              |          0.510884 |       0.998566 |
| VPDCPAS1-16->VPDCPAS1-14 | evap              |          0.517115 |       0.998637 |
| VPDCPAS1-43->VPDCPAS1-14 | evap              |          0.518663 |       0.998514 |
| VPDCPAS1-30->VPDCPAS1-13 | evap              |          0.520144 |       0.998662 |
| VPDCPAS1-60->VPDCPAS1-61 | evap              |          0.523967 |       0.995358 |
| VPDCPAS1-23->VPDCPAS1-20 | evap              |          0.531211 |       0.999131 |
| VPDCPAS1-60->VPDCPAS1-30 | evap              |          0.532364 |       0.995956 |
| VPDCPAS1-43->VPDCPAS1-30 | evap              |          0.540017 |       0.995001 |
| VPDCPAS1-16->VPDCPAS1-30 | evap              |          0.540738 |       0.994911 |
| VPDCPAS1-52->VPDCPAS1-39 | evap              |          0.56892  |       0.992994 |
| VPDCPAS1-60->VPDCPAS1-59 | evap              |          0.572496 |       0.99048  |
| VPDCPAS1-61->VPDCPAS1-60 | evap              |          0.576695 |       0.982936 |
| VPDCPAS1-16->VPDCPAS1-59 | evap              |          0.577257 |       0.989856 |
| VPDCPAS1-43->VPDCPAS1-40 | evap              |          0.580402 |       0.948968 |
| VPDCPAS1-40->VPDCPAS1-43 | evap              |          0.582927 |       0.946553 |

## Reviewer Interpretation

This run replaces the local E4 pilot with a fully public groundwater-chemistry dataset. It validates Hydrosheaf's public-data ingestion, unit harmonization, sparse candidate graph construction, and chemistry-fit behavior under real field observations. It does not provide an independently known reaction-path truth graph, so it should be described as public-field demonstration evidence rather than full process-truth validation.
