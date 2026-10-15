# Centered moving-average smoother

A simple centered rolling mean, used to attenuate day-of-week reporting
artefacts (lower weekend counts, Monday catch-up spikes) in daily
COVID-19 surveillance series before SIRD-GA estimation. Edge windows are
truncated (not padded), so the smoothed series has the same length as
the input.

## Usage

``` r
smooth_series(x, window = 7L)
```

## Arguments

- x:

  Numeric vector.

- window:

  Window width in days (default 7). Must be odd for a symmetric window;
  even values are incremented by 1 with a warning.

## Value

A numeric vector of the same length as `x`.
