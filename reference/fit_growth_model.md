# Fit a single growth-curve model by nonlinear least squares

Fits
[`logistic_curve()`](https://raydonal.github.io/epicurves/reference/logistic_curve.md),
[`gompertz_curve()`](https://raydonal.github.io/epicurves/reference/gompertz_curve.md)
or
[`richards_curve()`](https://raydonal.github.io/epicurves/reference/richards_curve.md)
to a cumulative-incidence series by nonlinear least squares
([`minpack.lm::nlsLM()`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html)),
which is markedly more robust to the identification issues of these
S-shaped curves than base
[`stats::nls()`](https://rdrr.io/r/stats/nls.html).

## Usage

``` r
fit_growth_model(
  t,
  y,
  model = c("logistic", "gompertz", "richards"),
  beta_fixed = 1.5,
  start = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- t:

  Numeric vector of time indices (`1, 2, ..., n`, relative to the start
  of the fitting window).

- y:

  Numeric vector of cumulative counts, same length as `t`.

- model:

  One of `"logistic"`, `"gompertz"`, `"richards"`.

- beta_fixed:

  If not `NULL` (Richards only), the shape parameter is held fixed at
  this value instead of estimated (improves convergence on
  short/fast-saturating phases, as used throughout the simulation
  study).

- start:

  Optional named list overriding the default starting values `y0`, `K`,
  `gamma` (and `beta` for Richards).

- lower, upper:

  Optional named lists overriding the default box constraints.

## Value

An object of class `growth_fit`: a list with `model`, the fitted `nls`
object, `coef` (named coefficient vector), `fitted` (fitted values at
`t`), `resid`, `mse` (mean squared error, Equation for EQM in the
manuscript), `aic`, and `n_par` (number of estimated parameters, used by
the AIC-based ensemble weights).
