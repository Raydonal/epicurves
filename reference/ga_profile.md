# Profile log-likelihood for one parameter of a `ga_fit`

Profile log-likelihood for one parameter of a `ga_fit`

## Usage

``` r
ga_profile(fit, param, grid, maxit = 3000)
```

## Arguments

- fit:

  A
  [`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
  result.

- param:

  Name of a free, scalar parameter.

- grid:

  Values at which the profile is computed. Each point is warm-started
  from the previous one, so an ordered grid that starts near the
  estimate works best.

- maxit:

  Nelder-Mead iteration limit per point.

## Value

A data frame with the grid value, the profile log-likelihood and the
other estimates.
