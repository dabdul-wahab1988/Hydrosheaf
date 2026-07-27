## Predeclared tuning, sample size, and software provenance

Regularisation used a Cartesian grid of nine $\lambda_1$ values logarithmically
spaced from 0.0003 to 0.08, three $\lambda_2$ values (0.005, 0.02, and 0.08),
and six support thresholds from 0.015 to 0.040 mmol/L. Selection used the
one-standard-error rule stated above. Sixty independent scenarios per
archetype were fixed before analysis as a compromise between factorial
coverage and live-PHREEQC cost. Thus, each predeclared noise-panel-method
comparison contains 240 independent scenarios; the 21,600 fits are not
treated as independent replicates.

Two software environments generated different, explicitly scoped artifacts.
The locked synthetic benchmark and its Table 1 statistics used Python 3.10.17,
NumPy 2.2.6, pandas 2.3.3, and PHREEQC 3.7.3-15968. The later field-transfer
saturation-index computations used Python 3.14.6, NumPy 2.4.6, pandas 2.3.3,
and PHREEQC 3.8.6-17100. This version difference does not mix synthetic and
field results: each reported artifact retains its generating scope, and Table
S9 records both environments.

Reviewer-requested post-hoc analyses read locked outputs without changing the
primary estimands. They quantify bounded-least-squares performance by
convergence status, the majority-class MRS comparator, the exact denominator
for thermodynamic-screening failure, and realised next-measurement value. The
complete machine-readable results are archived in
`results/m5_review_sensitivity.csv`.
