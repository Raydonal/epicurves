# Online segmentation of an epidemic curve with e-detectors

Sequentially monitors daily counts against a forecast band of a
growth-curve ensemble fitted to the current segment, and opens a new
segment when an anytime-valid e-detector (package edetect) raises an
alarm. Within a segment the ensemble is refitted every `refit_every`
days on data observed before the day being monitored, so each band is
predictable with respect to the past.

## Usage

``` r
edetect_segment(
  x,
  d_min = 20L,
  arl = 1000,
  refit_every = 7L,
  chart = c("poisson", "subgaussian"),
  z = 2,
  w_min = 0.05,
  weights = c("akaike", "inverse_mse"),
  type = c("sr", "cusum"),
  beta = 1.5
)
```

## Arguments

- x:

  Daily counts.

- d_min:

  Minimum segment length before monitoring starts.

- arl:

  Target average run length to a false alarm.

- refit_every:

  Refit interval (days) within a segment.

- chart:

  `"poisson"` or `"subgaussian"`.

- z:

  Half-width multiplier of each model's band (`0` monitors the point
  forecast).

- w_min:

  Minimum ensemble weight for a model to enter the envelope.

- weights:

  Ensemble weighting scheme.

- type:

  e-detector type, `"sr"` or `"cusum"`.

- beta:

  Richards shape.

## Value

A list with `cuts` (last day of each closed segment), `alarms` (alarm
times), `delay`, `band` (matrix of monitored forecast, lower and upper
band, `NA` on unmonitored days), `segments` and the final ensemble of
every segment (`fits`).

## Details

The band \\\[\ell_t, u_t\]\\ is the envelope, over the models carrying
at least `w_min` ensemble weight, of their delta-method bands
([`predict_band()`](https://raydonal.github.io/epicurves/reference/predict_band.md)).
Two one-sided detectors run in parallel, each with target average run
length `2 * arl`: an upward detector tests the composite null
"conditional mean at most \\u_t\\" and a downward detector the null "at
least \\\ell_t\\". If the conditional mean lies in the band, the
probability of any false alarm by time \\t\\ is at most
\\t/\mathrm{arl}\\ (union bound over the two detectors). With
`chart = "poisson"` the detectors are Poisson u charts with exposure
\\u_t\\ (resp. \\\ell_t\\); with `chart = "subgaussian"` they monitor
the negative-binomial Pearson residuals \\(x_t-u_t)/s(u_t)\\ and
\\(x_t-\ell_t)/s(\ell_t)\\, \\s(\mu)=\sqrt{\mu(1+\mu/\hat\phi)}\\,
declared sub-Gaussian with unit scale (an approximation).
