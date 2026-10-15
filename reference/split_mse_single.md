# Best single-model split MSE at a candidate partition point (fast proxy)

A single-growth-model analogue of
[`split_mse()`](https://raydonal.github.io/epicurves/reference/split_mse.md),
roughly three times cheaper (one model fit per side instead of three),
used internally by
[`permutation_test_partition()`](https://raydonal.github.io/epicurves/reference/permutation_test_partition.md)
to generate the null distribution of the grid-maximum MSE-reduction
statistic at a cost that keeps a full permutation test tractable. The
model used is `null_model` (default `"gompertz"`), fixed across the
whole null simulation, not re-selected per replicate.

## Usage

``` r
split_mse_single(t, y, d, null_model = "gompertz", beta_fixed = 1.5)
```

## Arguments

- t, y:

  Time index and cumulative counts of the (open) window being searched,
  both starting at 1.

- d:

  Candidate partition point, `1 <= d < length(y)`.

- null_model:

  One of `"logistic"`, `"gompertz"`, `"richards"`.
