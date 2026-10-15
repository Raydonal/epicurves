# Ensemble of growth curves fitted to daily counts

Fits the Logistic, Gompertz and Richards incidence models by maximum
likelihood and combines them with Akaike weights \\w_i \propto
\exp(-\Delta\mathrm{AIC}\_i/2)\\ (`weights = "akaike"`) or with
inverse-MSE weights on the count scale (`weights = "inverse_mse"`).

## Usage

``` r
fit_ensemble_incidence(
  x,
  family = c("poisson", "negbin"),
  weights = c("akaike", "inverse_mse"),
  beta = 1.5,
  models = c("gompertz", "richards", "logistic"),
  warm = NULL
)
```

## Arguments

- x:

  Daily counts.

- family:

  `"poisson"` or `"negbin"`.

- weights:

  Weighting scheme.

- beta:

  Richards shape (held fixed).

- models:

  Models included.

- warm:

  Optional previous `incidence_ensemble` whose estimates are used as
  warm starts.

## Value

An object of class `incidence_ensemble`.
