# Profile log-likelihood of one SIRD-GA parameter

Maximises the log-likelihood over the remaining free parameters on a
grid of fixed values of `param`, warm-starting each grid point from the
previous optimum.

## Usage

``` r
sird_ga_profile(fit, cases = NULL, deaths = NULL, param, grid, maxit = 3000)
```

## Arguments

- fit:

  Output of
  [`sird_ga_fit_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit_inc.md).

- cases, deaths:

  Data used for `fit`.

- param:

  Name of the profiled parameter (one of `fit$free`).

- grid:

  Values at which to profile.

- maxit:

  Nelder-Mead iteration limit per grid point.

## Value

Data frame with `value`, `loglik` and the re-optimised parameters.
