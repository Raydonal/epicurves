# Permutation test for a single candidate structural change point

Tests whether the observed reduction in pooled MSE from the best
candidate split in
[`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md)
is larger than what a single homogeneous (unsplit) growth process would
produce by chance, properly accounting for the multiple comparisons
implicit in scanning the whole candidate grid: the null distribution is
built from the *grid-maximum* MSE reduction in each of `B`
parametric-bootstrap replicates simulated from the unsplit null fit,
following the same logic as the permutation test of
[`ecp::e.divisive`](https://rdrr.io/pkg/ecp/man/e.divisive.html)
(Matteson & James, 2013) but specialised to parametric growth curves
rather than a general multivariate distributional distance. For
computational tractability, the null replicates' grid search uses a
single growth model
([`split_mse_single()`](https://raydonal.github.io/epicurves/reference/split_mse_single.md))
rather than the full three-model ensemble used for the reported fit;
this is an approximation, disclosed here rather than left implicit.

## Usage

``` r
permutation_test_partition(
  t,
  y,
  min_obs = 12L,
  B = 99L,
  null_model = NULL,
  beta_fixed = 1.5,
  use_ensemble_null = TRUE
)
```

## Arguments

- t, y:

  Time index and cumulative counts of the open window.

- min_obs:

  Minimum segment length (as in
  [`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md)).

- B:

  Number of parametric-bootstrap null replicates.

- null_model:

  Growth model used for the (cheaper) null-replicate grid searches;
  defaults to whichever of the three models best fits the unsplit
  window.

- beta_fixed:

  Passed to the growth-curve fits.

- use_ensemble_null:

  Logical. If `TRUE`, the null replicates are generated from, and re-fit
  with, the full three-model ensemble (matching the observed statistic
  exactly, at roughly three times the cost); if `FALSE` (default), the
  cheaper single-model proxy of
  [`split_mse_single()`](https://raydonal.github.io/epicurves/reference/split_mse_single.md)
  is used. We found the single-model proxy under-powered in simulation
  (it can fail to confirm genuine change points detected by the full
  ensemble); `use_ensemble_null = TRUE` is recommended whenever the
  extra cost is affordable.

## Value

A list with `p_value`, `observed_reduction` (MSE_nosplit minus MSE at
the best split), `null_reductions` (the `B` simulated null maxima) and
`d_hat` (the best split, as in
[`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md)).
