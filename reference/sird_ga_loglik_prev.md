# Fast Kalman log-likelihood for prevalence observations

Compiled counterpart of
[`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md)
(same model and observation equation, RK4 transition quantities).

## Usage

``` r
sird_ga_loglik_prev(theta, y, N, Delta = 1, nsub = 10L)
```

## Arguments

- theta:

  Named numeric vector with elements `lambda`, `gamma`, `alpha`, `s0`,
  `i0`, `r0`, `p`, `tau2` (all strictly positive except `r0`, which may
  be zero; `p` in `(0,1)`).

- y:

  Numeric vector of observed infectious **head counts** \\y_0, y_1,
  \ldots, y_n\\ (not normalised; normalisation by `N` is performed
  internally, matching the \\O_N(t_k)=O(t_k)/N\\ observation model of
  the manuscript).

- N:

  Population size.

- Delta:

  Time step between consecutive observations (default `1`).

- nsub:

  RK4 steps per interval.
