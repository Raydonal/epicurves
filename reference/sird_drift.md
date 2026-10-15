# Deterministic SIRD drift vector

Right-hand side of the normalised SIRD ODE system \\\dot s = -\lambda s
i\\, \\\dot i = \lambda s i - (\gamma+\alpha) i\\, \\\dot r = \gamma
i\\.

## Usage

``` r
sird_drift(state, parms)
```

## Arguments

- state:

  Numeric vector `c(s, i, r)`.

- parms:

  Named numeric vector/list with elements `lambda`, `gamma`, `alpha`.

## Value

Numeric vector of length 3, `b(x)`.
