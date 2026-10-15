# Richards growth curve

Richards growth curve

## Usage

``` r
richards_curve(t, y0, K, gamma, beta)
```

## Arguments

- t:

  Time.

- y0:

  Cumulative count at the start of the fitting window.

- K:

  Final epidemic size (carrying capacity).

- gamma:

  Growth rate.

- beta:

  Shape parameter.

## Details

Parametrised so that `richards_curve(0, y0, K, gamma, beta) == y0` and
`richards_curve(t, ...) -> K` as `t -> Inf`, consistent with the
Logistic and Gompertz curves in this package and with the standard
generalised-logistic form (Richards, 1959). Note this uses
`exp(-beta * gamma * t)` (decaying exponential); the thesis chapter this
package accompanies transcribes the equation with
`exp(+beta * gamma * t)`, which drives the curve to 0 rather than to `K`
as `t` grows — a sign slip relative to the standard form, corrected here
(and in the revised manuscript).
