# Bootstrap confidence intervals for each growth model's parameters

Resamples the ensemble fit (see
[`ensemble_bootstrap_ci()`](https://raydonal.github.io/epicurves/reference/ensemble_bootstrap_ci.md)
for the choice between `"poisson"` and `"residual"` resampling) and,
instead of (or in addition to) the ensemble curve, collects each of the
three growth models' own parameter estimates (\\Y_0\\, \\\gamma\\,
\\K\\) and weight on every replicate, returning percentile confidence
limits for each. This is what the manuscript's per-segment parameter
tables report.

## Usage

``` r
ensemble_bootstrap_params(
  fit,
  method = c("poisson", "residual"),
  B = 200L,
  level = 0.95,
  beta_fixed = 1.5
)
```

## Arguments

- fit:

  An `ensemble_fit` object.

- method:

  `"poisson"` (default) or `"residual"`.

- B:

  Number of bootstrap replicates.

- level:

  Confidence level.

- beta_fixed:

  Passed to
  [`fit_ensemble()`](https://raydonal.github.io/epicurves/reference/fit_ensemble.md)
  for each replicate.

## Value

A data frame with one row per model x parameter, giving the point
estimate, lower/upper confidence limits, and number of successful
replicates (`B` minus failed refits).
