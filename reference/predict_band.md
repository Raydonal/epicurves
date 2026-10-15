# Forecast band of a growth-incidence fit

Delta-method band \\\hat\mu_u \pm z\\\widehat{\mathrm{se}}(\hat\mu_u)\\,
with the standard error propagated from the inverse observed information
of the unconstrained parameters by a numerical gradient.

## Usage

``` r
predict_band(object, u, z = 2)
```

## Arguments

- object:

  A `growth_incidence_fit`.

- u:

  Local times.

- z:

  Band half-width multiplier.

## Value

Matrix with columns `mu`, `lower`, `upper`.
