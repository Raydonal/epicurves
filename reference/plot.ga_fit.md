# Plot the one-step-ahead predictions of a `ga_fit`

One panel per observation stream: observed counts (points), predicted
means (line) and predicted means plus and minus two standard deviations
(band), all from the Kalman filter at the estimates.

## Usage

``` r
# S3 method for class 'ga_fit'
plot(x, labels = NULL, ...)
```

## Arguments

- x:

  A
  [`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
  result.

- labels:

  Optional panel labels, one per stream.

- ...:

  Unused.

## Value

A `ggplot` object.
