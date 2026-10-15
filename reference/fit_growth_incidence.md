# Fit one growth curve to daily counts by maximum likelihood

Poisson or negative-binomial likelihood for the daily counts `x` of a
segment, with mean given by
[`growth_incidence()`](https://raydonal.github.io/epicurves/reference/growth_incidence.md).
Parameters are optimised on the scale \\(\log y_0, \log(K - y_0),
\log\gamma\[, \log\phi\])\\ from a small grid of starting values.

## Usage

``` r
fit_growth_incidence(
  x,
  model = c("logistic", "gompertz", "richards"),
  family = c("poisson", "negbin"),
  beta = 1.5,
  start = NULL
)
```

## Arguments

- x:

  Daily counts of the segment (local times `1, ..., length(x)`).

- model:

  Growth model.

- family:

  `"poisson"` or `"negbin"`.

- beta:

  Richards shape (held fixed).

- start:

  Optional warm start on the unconstrained scale (e.g. `par` of a
  previous fit), tried in addition to the default grid.

## Value

A list with `coef`, `loglik`, `aic`, `bic`, `k` (number of estimated
parameters), `fitted`, `model`, `family`, `beta`, `n`.
