# Jacobian of the SIRD drift vector

Analytic Jacobian \\\nabla_x b(x)\\ of
[`sird_drift()`](https://raydonal.github.io/epicurves/reference/sird_drift.md),
used to build the state-transition (resolvent) matrix \\\Phi(t,s)\\ of
the SIRD-GA model.

## Usage

``` r
sird_jacobian(state, parms)
```

## Arguments

- state:

  Numeric vector `c(s, i, r)`.

- parms:

  Named numeric vector/list with elements `lambda`, `gamma`, `alpha`.

## Value

A 3x3 numeric matrix.
