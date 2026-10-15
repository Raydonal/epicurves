# Predict the ensemble curve at arbitrary times

Predict the ensemble curve at arbitrary times

## Usage

``` r
# S3 method for class 'ensemble_fit'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An `ensemble_fit` object.

- newdata:

  Numeric vector of times (defaults to `object$t`).

- ...:

  Unused; present for S3 consistency with
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html).
