# Fit the full ensemble of growth-curve models to one data window

Fits Logistic, Gompertz and Richards curves
([`fit_growth_model()`](https://raydonal.github.io/epicurves/reference/fit_growth_model.md))
to a single window of cumulative counts, computes MSE- and AIC-based
weights, and returns the ensemble (weighted-average) curve together with
the individual fits.

## Usage

``` r
fit_ensemble(t, y, criterion = c("mse", "aic"), beta_fixed = 1.5)
```

## Arguments

- t, y:

  As in
  [`fit_growth_model()`](https://raydonal.github.io/epicurves/reference/fit_growth_model.md).

- criterion:

  `"mse"` (default, matching the main text) or `"aic"` (matching the
  manuscript's comparative analysis).

- beta_fixed:

  Shape parameter passed to the Richards fit (see
  [`fit_growth_model()`](https://raydonal.github.io/epicurves/reference/fit_growth_model.md)).

## Value

An object of class `ensemble_fit`: a list with `models` (named list of
three `growth_fit` objects), `weights`, `criterion`, `fitted` (ensemble
fitted values at `t`), `mse` (ensemble MSE) and `t`, `y`.
