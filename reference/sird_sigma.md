# Diffusion "square-root" factor of the SIRD jump process

Lower-triangular-like matrix \\\sigma(x)\\ such that
\\\sigma(x)\sigma(x)^\top = \Sigma(x)\\, driven by a 3-dimensional
Brownian motion with one independent noise channel per elementary event.

## Usage

``` r
sird_sigma(state, parms)
```

## Arguments

- state:

  Numeric vector `c(s, i, r)`.

- parms:

  Named numeric vector/list with elements `lambda`, `gamma`, `alpha`.

## Value

A 3x3 numeric matrix.
