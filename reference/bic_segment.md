# Offline segmentation of daily counts by penalised likelihood

Binary segmentation of the daily-count series: a segment is split at the
candidate day maximising the summed ensemble log-likelihood of the two
halves, and the split is kept only if it lowers the Bayesian information
criterion, with the split charged for the parameters of one additional
curve and the change location. The parameter count of an ensemble is the
weight-averaged count of its members (three per curve with Poisson
counts, four with negative binomial counts): the members are alternative
descriptions of a single mean, not additive components. Each half is
fitted by
[`fit_ensemble_incidence()`](https://raydonal.github.io/epicurves/reference/fit_ensemble_incidence.md).

## Usage

``` r
bic_segment(
  x,
  d_min = 20L,
  family = c("poisson", "negbin"),
  weights = c("akaike", "inverse_mse"),
  beta = 1.5,
  step = 2L,
  max_cuts = 10L
)
```

## Arguments

- x:

  Daily counts.

- d_min:

  Minimum segment length.

- family, weights, beta:

  Passed to
  [`fit_ensemble_incidence()`](https://raydonal.github.io/epicurves/reference/fit_ensemble_incidence.md).

- step:

  Spacing of the coarse candidate grid (days); the optimum is then
  refined at unit spacing within `step` days.

- max_cuts:

  Maximum number of change points.

## Value

A list with `cuts` and `segments`.
