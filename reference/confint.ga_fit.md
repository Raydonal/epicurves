# Wald confidence intervals for a `ga_fit`

Intervals are computed on the unconstrained scale of the optimisation
and transformed back, so they respect the parameter bounds.

## Usage

``` r
# S3 method for class 'ga_fit'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A
  [`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
  result.

- parm:

  Ignored (all parameters are returned).

- level:

  Confidence level.

- ...:

  Unused.
