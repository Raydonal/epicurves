# Inverse-criterion ensemble weights

Computes weights proportional to the inverse of each model's
goodness-of- fit criterion (mean squared error or AIC), so that
better-fitting models receive more weight, following Chowell et al.
(2021) as adapted in the manuscript.

## Usage

``` r
ensemble_weights(crit)
```

## Arguments

- crit:

  Named numeric vector of per-model criterion values (MSE or AIC; lower
  is better for both).

## Value

Named numeric vector of weights summing to 1.
