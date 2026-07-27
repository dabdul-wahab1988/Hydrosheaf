# M8 confirmatory claim decision

Run: `RUN-M8-CONFIRM-20260728-01`
Status: PASS
Scope: controlled synthetic evaluation through production HydroSheaf adapters;
not field validation.

## Decision

The proposed headline that scalar identifiability metrics conceal opposing
parameter-specific optimal sampling times is **not confirmed**. After using a
noise-whitened Fisher matrix in log-parameter space, removing random
start-centred priors, and applying the locked decision rule, every
parameter-specific and joint local criterion selected the same 50 d front
observation from the three-point late-time base.

This is not a null result. It gives two supported findings:

1. **A front observation repaired the late-only transport design for both
   parameters.** Under local starts, the common 50 d choice reduced median
   dispersivity absolute log10 error from 0.2154 to 0.0173. The paired median
   difference was -0.1952 (95% bootstrap interval -0.2874 to -0.1656). Decay
   error decreased from 0.0182 to 0.0146, with a paired difference of -0.0040
   (-0.0071 to -0.00003). All calibrations succeeded, and the distant-start
   sensitivity produced essentially the same median errors.
2. **The kinetic adapter demonstrated a structural limit of sampling design.**
   One and four residence-time experiments both produced rank-one information
   for calcite `k` and `A`. Predictions along the constant-`k*A` ridge differed
   by at most 1.33e-6 declared observation standard deviations, whereas doubling
   the product changed predictions by 6.42 standard deviations. An independent
   surface-area observation increased the information rank to two.

The fixed four-point transport sweep still showed parameter-specific placement:
`c90_s45` had the lowest dispersivity error (0.0172), whereas `very_late` had the
lowest decay error (0.0140). That supports reporting per-parameter resolution,
but it does not rescue the stronger claim of opposed conditioning responses or
different optimal next observations. Replicate-wise Spearman correlations with
the corrected condition number were positive for both parameters (median 0.587
for dispersivity and 0.490 for decay); the 95% replicate interval excluded zero
only for dispersivity.

## Manuscript consequence

Withdraw these claims:

- scalar condition, D-, A- and E-optimality were actively misleading in the
  transport confirmation;
- dispersivity- and decay-targeted criteria selected different next samples;
- the earlier negative decay-versus-conditioning correlation was a robust
  physical result.

Retain this bounded claim:

> In a controlled 1-D transport benchmark, a noise-whitened log-parameter
> Fisher analysis selected one front observation that improved recovery of both
> dispersivity and decay from a late-only base. In the PHREEQC kinetic adapter,
> residence-time sampling could not separate a rate constant from reactive
> surface area because the forward law identified only their product; an
> independent surface-area constraint restored local identifiability.

The analytic toy benchmark may remain a positive control, but its 60-fold result
must not be presented as the real-adapter effect. The next optional extension is
independent-model robustness—generate transport truth with the numerical solver
and calibrate with the analytical adapter—rather than adding more hand-selected
sampling designs to the present matched-model experiment.
