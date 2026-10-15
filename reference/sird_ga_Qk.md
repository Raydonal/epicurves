# State-noise covariance \\Q_k\\ of interval `k`

Computes \\Q_k = N^{-1}\int\_{t\_{k-1}}^{t_k}
\Phi(t_k,s)\\\Sigma(x(s))\\ \Phi(t_k,s)^\top ds\\ by Simpson's rule,
with \\\Phi(t_k,s) = \Phi(t_k,t\_{k-1})\\\Phi(s,t\_{k-1})^{-1}\\, using
the locally-anchored fundamental matrix stored in `tr$segments[[k]]`.

## Usage

``` r
sird_ga_Qk(tr, k, N)
```

## Arguments

- tr:

  Object returned by
  [`sird_ga_transition()`](https://raydonal.github.io/epicurves/reference/sird_ga_transition.md).

- k:

  Integer interval index, `1, ..., n`.

- N:

  Population size.
