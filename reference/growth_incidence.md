# Daily incidence implied by a growth curve

For a cumulative growth curve \\C(u)\\ with \\C(0)=y_0\\, returns the
expected count in each unit interval, \\\mu_u = C(u) - C(u-1)\\, \\u =
1, 2, \ldots\\.

## Usage

``` r
growth_incidence(u, model, y0, K, gamma, beta = 1.5)
```

## Arguments

- u:

  Integer local times (1 = first day of the segment).

- model:

  One of `"logistic"`, `"gompertz"`, `"richards"`.

- y0, K, gamma:

  Curve parameters (`0 < y0 < K`, `gamma > 0`).

- beta:

  Richards shape parameter.
