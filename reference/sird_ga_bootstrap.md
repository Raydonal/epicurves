# Parametric bootstrap confidence intervals for the SIRD-GA MLE

Two bootstrap schemes are supported. `method = "gillespie"` simulates
`B` independent epidemic trajectories and surveillance records from the
fitted parameters via the exact jump-process simulator
([`sird_ga_simulate_obs()`](https://raydonal.github.io/epicurves/reference/sird_ga_simulate_obs.md))
and re-estimates the model on each replicate — appropriate when
extinction risk is low. `method = "residual"` resamples the standardised
one-step-ahead Kalman residuals of the fitted model with replacement,
reconstructs bootstrap observation series around the fitted predictive
mean, and re-estimates the model on each — appropriate for real-data
settings near the epidemic threshold, where re-simulating the jump
process can produce early extinction.

## Usage

``` r
sird_ga_bootstrap(
  fit,
  y,
  N,
  Delta = 1,
  method = c("gillespie", "residual"),
  B = 200L,
  level = 0.95,
  S0 = NULL,
  I0 = NULL,
  R0 = NULL
)
```

## Arguments

- fit:

  Object returned by
  [`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md).

- y:

  Original observed series used to produce `fit`.

- N, Delta:

  Passed through to
  [`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md)
  / the simulator.

- method:

  `"gillespie"` or `"residual"`.

- B:

  Number of bootstrap replicates.

- level:

  Confidence level.

- S0, I0, R0:

  Initial head counts, required for `method = "gillespie"` (derived from
  `fit$theta` and `N` if omitted).

## Value

A list with `boot_theta` (`B x p` matrix of bootstrap estimates) and
`ci` (matrix of percentile confidence limits, one row per parameter).
