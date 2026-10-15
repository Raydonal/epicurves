# SIRD-GA with piecewise-constant transmission, fitted to cases and deaths

Maximum likelihood for the SIRD-GA model in which the transmission rate
is constant within blocks of observation intervals,
\\\lambda(t)=\lambda_b\\ for \\t\\ in block \\b\\, while the removal
rates \\\gamma\\, \\\alpha\\, the initial infectious fraction, the case
reporting probability and the two noise scales are shared. Parameters
listed in `fixed` are held at the given values (e.g. a profile over
`p`).

## Usage

``` r
sird_ga_fit_tv(
  cases = NULL,
  deaths = NULL,
  N,
  blocks,
  use = c("both", "cases", "deaths"),
  start,
  fixed = list(),
  q = 1,
  p_max = 0.999,
  Delta = 1,
  nsub = 10L,
  maxit = 3000
)
```

## Arguments

- cases, deaths, N, use, Delta, nsub:

  See
  [`sird_ga_loglik_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik_inc.md).

- blocks:

  Integer vector, one entry per observation interval, giving the block
  of each interval (`1, ..., B`).

- start:

  Named list of starting values for `gamma`, `alpha`, `i0`, `p`, `tau2`,
  `tau2_d` and a vector `lambda` of length `B`.

- fixed:

  Named list of parameters held fixed (any of `gamma`, `alpha`, `i0`,
  `p`, `tau2`, `tau2_d`); `q` (default 1) is always fixed.

- q:

  Death reporting probability (fixed).

- p_max:

  Upper bound for `p`.

- maxit:

  Iteration limit.

## Value

A list with `lambda` (per block), `theta`, `loglik`, `par`, `hessian`,
`R_t` (per interval, \\\lambda(t)s(t)/(\gamma+\alpha)\\) and the full
Kalman output `kf`.
