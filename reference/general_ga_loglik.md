# Kalman-filter log-likelihood for a general jump model

Dimension-generic analogue of
[`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md).
A single compartment is observed (binomial under-reporting plus Gaussian
measurement error, as in the SIRD-GA observation model), identified by
name in `obs_compartment`.

## Usage

``` r
general_ga_loglik(
  theta,
  y,
  N,
  model,
  obs_compartment,
  Delta = 1,
  n_sub = 10L,
  Q0_diag = 1e-08,
  keep_states = FALSE
)
```

## Arguments

- theta:

  Named numeric vector: the model's rate parameters, its initial
  condition (one named entry per `model$compartments`), `p` and `tau2`.

- y, N, Delta, n_sub, Q0_diag, keep_states:

  As in
  [`sird_ga_loglik()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik.md).

- model:

  A `jump_model`.

- obs_compartment:

  Name of the observed compartment (must be one of
  `model$compartments`).
