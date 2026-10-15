# Grid search for a single structural change (partition) point

Scans a grid of candidate partition points with increments of size 1, as
described in the manuscript, and selects the point minimising the pooled
mean squared error, subject to a minimum number of observations on each
side.

## Usage

``` r
find_partition_point(t, y, min_obs = 12L, criterion = "mse", beta_fixed = 1.5)
```

## Arguments

- t, y:

  Time index and cumulative counts of the open window.

- min_obs:

  Minimum number of observations required on each side of the candidate
  split (default 12, as in Simulation Scenarios 1 and 3 of the
  manuscript; use 15 to match Scenario 2).

- criterion, beta_fixed:

  Passed to
  [`fit_ensemble()`](https://raydonal.github.io/epicurves/reference/fit_ensemble.md).

## Value

A list with `grid` (candidate points searched), `mse_profile` (MSE at
each grid point), `d_hat` (selected partition point), and `mse_split`
(MSE at `d_hat`).
