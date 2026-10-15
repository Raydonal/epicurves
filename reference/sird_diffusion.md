# Diffusion covariance matrix of the SIRD jump process

Symmetric matrix \\\Sigma(x)\\ from the second-order moments of the jump
rates, obtained from the three elementary events (infection, recovery,
disease-induced death).

## Usage

``` r
sird_diffusion(state, parms)
```

## Arguments

- state:

  Numeric vector `c(s, i, r)`.

- parms:

  Named numeric vector/list with elements `lambda`, `gamma`, `alpha`.

## Value

A 3x3 symmetric numeric matrix.
