# Bootstrap confidence intervals for a `ga_boot`

Percentile intervals use the quantiles of the replicates; basic
intervals reflect them about the original estimate. Replicates with a
missing value of a quantity are dropped for that quantity.

## Usage

``` r
# S3 method for class 'ga_boot'
confint(object, parm, level = 0.95, type = c("percentile", "basic"), ...)
```

## Arguments

- object:

  A
  [`ga_bootstrap()`](https://raydonal.github.io/epicurves/reference/ga_bootstrap.md)
  result.

- parm:

  Names of the quantities; all by default.

- level:

  Confidence level.

- type:

  `"percentile"` or `"basic"`.

- ...:

  Unused.
