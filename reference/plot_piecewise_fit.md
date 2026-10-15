# Plot a fitted piecewise ensemble curve, with confirmed partition points

Plot a fitted piecewise ensemble curve, with confirmed partition points

## Usage

``` r
plot_piecewise_fit(fit, ci = NULL, forecast_from = NULL)
```

## Arguments

- fit:

  A `piecewise_fit` object.

- ci:

  Optional data frame with columns `time`, `lower`, `upper` (e.g. built
  by calling
  [`ensemble_bootstrap_ci()`](https://raydonal.github.io/epicurves/reference/ensemble_bootstrap_ci.md)
  per segment and row-binding the results) to draw as a shaded
  confidence band.

- forecast_from:

  Optional time beyond which a dashed vertical line marks the start of
  an out-of-sample forecast window.
