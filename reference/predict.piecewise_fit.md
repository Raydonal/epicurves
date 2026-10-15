# Predict from a piecewise ensemble fit

Predict from a piecewise ensemble fit

## Usage

``` r
# S3 method for class 'piecewise_fit'
predict(object, newdata, ...)
```

## Arguments

- object:

  A `piecewise_fit` object.

- newdata:

  Numeric vector of times (original time scale); times beyond the last
  confirmed segment's final observation are forecast using that
  segment's ensemble fit.

- ...:

  Unused; present for S3 consistency with
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html).
