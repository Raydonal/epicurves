# Plot observed vs. fitted SIRD-GA incidence, with an optional bootstrap band

Plot observed vs. fitted SIRD-GA incidence, with an optional bootstrap
band

## Usage

``` r
plot_sird_ga_fit(y, kf, ci = NULL, cumulative = FALSE)
```

## Arguments

- y:

  Observed counts.

- kf:

  Output of
  [`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md)
  called with `keep_states = TRUE`.

- ci:

  Optional data frame (`time`, `lower`, `upper`) on the head-count
  scale, e.g. built from
  [`sird_ga_bootstrap()`](https://raydonal.github.io/epicurves/reference/sird_ga_bootstrap.md)
  replicate trajectories.

- cumulative:

  Logical; plot cumulative (default `FALSE`, i.e. daily) incidence.
