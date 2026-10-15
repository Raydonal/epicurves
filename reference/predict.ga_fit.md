# One-step-ahead predictions from a `ga_fit`

One-step-ahead predictions from a `ga_fit`

## Usage

``` r
# S3 method for class 'ga_fit'
predict(object, ...)
```

## Arguments

- object:

  A
  [`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
  result.

- ...:

  Unused.

## Value

A list with `fitted` and `sd` (matrices of predicted means and standard
deviations of the counts, one column per stream) and `x` (path of the
mean, as fractions).
