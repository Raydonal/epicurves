# Plot an MSE (or AIC) profile from a partition-point grid search

Reproduces the manuscript's EQM-profile figures.

## Usage

``` r
plot_partition_profile(search, true_d = NULL)
```

## Arguments

- search:

  Object returned by
  [`find_partition_point()`](https://raydonal.github.io/epicurves/reference/find_partition_point.md).

- true_d:

  Optional true partition point (vertical reference line).
