# Mean squared error of a split at a single candidate partition point

Fits the growth-curve ensemble separately to the two sides of a
candidate split `d` and returns the pooled MSE, as used by the grid
search in
[`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md).

## Usage

``` r
split_mse(t, y, d, criterion = "mse", beta_fixed = 1.5)
```

## Arguments

- t, y:

  Time index and cumulative counts of the (open) window being searched,
  both starting at 1.

- d:

  Candidate partition point, `1 <= d < length(y)`.

- criterion, beta_fixed:

  Passed to
  [`fit_ensemble()`](https://raydonal.github.io/epicurves/reference/fit_ensemble.md).

## Value

A list with `mse` (pooled MSE, `NA` if either side fails to converge)
and `fit1`, `fit2` (the two `ensemble_fit` objects, or `NULL` on
failure).
