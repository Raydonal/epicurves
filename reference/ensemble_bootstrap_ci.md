# Bootstrap confidence band for an ensemble fit

Two resampling schemes are supported. `method = "poisson"` implements
the parametric bootstrap of the manuscript's simulation study: bootstrap
pseudo-observations are drawn as \\y_k^{(b)} \sim
\mathrm{Poisson}(f\_{ens}(t_k))\\, appropriate when the data-generating
process is itself (approximately) Poisson, as in the simulated growth
curves used there. `method = "residual"` instead resamples the
standardised fitted residuals with replacement and adds them back onto
the fitted curve, which is appropriate for real surveillance counts: the
Poisson scheme assumes \\\mathrm{Var}(y_k)=f\_{ens}(t_k)\\, which is
routinely violated by over-dispersed daily death/case counts (reporting
batches, weekday effects), producing artificially narrow bands (see the
real-data application in the manuscript, where this is checked
directly). In both cases the full ensemble (all three growth models,
re-weighted) is refit to each replicate, and percentile confidence
limits are taken across the `B` replicates at every requested time
point.

## Usage

``` r
ensemble_bootstrap_ci(
  fit,
  newdata = NULL,
  method = c("poisson", "residual"),
  B = 500L,
  level = 0.95,
  beta_fixed = 1.5
)
```

## Arguments

- fit:

  An `ensemble_fit` object.

- newdata:

  Numeric vector of times at which to build the band (defaults to
  `fit$t`; include out-of-window times for a forecast band).

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

A list with `boot_curves` (`B x length(newdata)` matrix) and `ci` (data
frame with `time`, `lower`, `upper`, `fit`).
