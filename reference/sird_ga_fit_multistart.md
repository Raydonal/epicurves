# Multi-start SIRD-GA fitting, with a practical-identifiability summary

Fits the SIRD-GA model from several starting points spanning a range of
assumed mean infectious periods, and reports all converged optima sorted
by log-likelihood. This is useful because, on real surveillance data,
the SIRD-GA log-likelihood can have a near-flat ridge along which
`lambda` trades off against `gamma + alpha` while `R0` and the shape of
the fitted curve stay essentially unchanged — a known
weak-identifiability feature of SIR/SIRD models fit to incidence-only
series (the removal rate is identified mainly through its sum with the
transmission rate and the curve's inflection timing, not independently).
Reporting the single best-fitting optimum without checking for this
ridge can be misleading; this function makes the ridge (if present)
visible.

## Usage

``` r
sird_ga_fit_multistart(
  y,
  N,
  Delta = 1,
  r0_guess = 1.1,
  infectious_period_days = c(2, 3, 4, 5, 7, 10),
  alpha_frac = 0.1,
  p0 = 0.9,
  tau2_0 = 1,
  i0_0 = 1e-04,
  ...
)
```

## Arguments

- y:

  Numeric vector of observed infectious counts, `y_0, ..., y_n`.

- N:

  Population size.

- Delta:

  Time step between observations.

- r0_guess:

  Starting guess for \\R_0\\ (used to scale `lambda` across starts,
  together with each candidate infectious period).

- infectious_period_days:

  Numeric vector of candidate mean infectious periods
  \\1/(\gamma+\alpha)\\ (in days) to start from (default spans 2-10
  days, covering typical respiratory-infection priors).

- alpha_frac:

  Fraction of the total removal rate `gamma + alpha` assigned to `alpha`
  at each starting point (default 0.1).

- p0, tau2_0, i0_0:

  Starting values shared across all starts for `p`, `tau2` and `i0`.

- ...:

  Further arguments passed to
  [`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md).

## Value

A list with `results` (a list of
[`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md)
outputs, one per start, in the same order as `infectious_period_days`),
`summary` (a data frame with the starting infectious period, final
log-likelihood, `R0` and key parameter estimates for each start, sorted
best-first) and `best` (the single best-fitting result).
