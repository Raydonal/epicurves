# Maximum (approximate) likelihood estimation for the SIRD-GA model

Estimates `theta = (lambda, gamma, alpha, i0, p, tau2)` (and optionally
`r0`) by Nelder-Mead maximisation of the Kalman-filter approximate
log-likelihood
([`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md)),
following the estimation strategy described in the SIRD-GA manuscript.
Positivity and unit-interval constraints are enforced through a
log/logit reparametrisation, so the optimiser itself runs unconstrained.

## Usage

``` r
sird_ga_fit(
  y,
  N,
  Delta = 1,
  start,
  r0_fixed = 0,
  estimate_r0 = FALSE,
  p_max = 0.999,
  n_sub = 10L,
  Q0_diag = 1e-08,
  control = list(maxit = 3000, reltol = 1e-10),
  hessian = TRUE
)
```

## Arguments

- y:

  Numeric vector of observed infectious counts, `y_0, ..., y_n`.

- N:

  Population size.

- Delta:

  Time step between observations.

- start:

  Named numeric vector/list with starting values for `lambda`, `gamma`,
  `alpha`, `i0`, `p`, `tau2` (and `r0` if `estimate_r0 = TRUE`).

- r0_fixed:

  Value at which `r0` is held fixed when `estimate_r0 = FALSE` (default
  `0`, appropriate for data starting at or near the onset of the
  epidemic).

- estimate_r0:

  Logical; estimate `r0` jointly instead of fixing it.

- p_max:

  Upper bound enforced on the reporting probability `p` (default
  `0.999`). Real daily case series are rarely reported with near-100%
  completeness; without an upper cap strictly below 1, the optimiser can
  drive `p -> 1` and `tau2` to a degenerate extreme (perfect reporting,
  all residual variance absorbed by the measurement-noise term), a
  boundary artefact rather than a genuine optimum. Setting `p_max`
  around 0.95-0.98 is recommended for real data; simulation studies
  where the reporting rate is itself the target of recovery can use
  `p_max = 0.999`.

- n_sub, Q0_diag:

  Passed to
  [`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md).

- control:

  List of control parameters passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- hessian:

  Logical; compute the numerical Hessian in the unconstrained
  parametrisation (used only for a quick asymptotic diagnostic — see
  [`sird_ga_bootstrap()`](https://raydonal.github.io/epicurves/reference/sird_ga_bootstrap.md)
  for the inferential confidence intervals reported in the manuscript).

## Value

A list with elements `theta` (named vector of MLEs in natural scale),
`loglik`, `convergence`, `optim_fit` (raw
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) output in the
unconstrained parametrisation), and `free_se` (delta-free standard
errors on the unconstrained scale, or `NULL`).
