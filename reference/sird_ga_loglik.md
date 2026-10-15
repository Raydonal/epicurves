# Kalman-filter log-likelihood of the SIRD-GA state-space model

Runs the predict/update Kalman recursion of the linear-Gaussian
state-space representation of the SIRD-GA model and returns the
approximate log-likelihood of the observed (partially reported, noisy)
infectious counts.

## Usage

``` r
sird_ga_loglik(
  theta,
  y,
  N,
  Delta = 1,
  n_sub = 10L,
  Q0_diag = 1e-08,
  keep_states = FALSE
)
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

- n_sub:

  Internal quadrature resolution, passed to
  [`sird_ga_transition()`](https://raydonal.github.io/epicurves/reference/sird_ga_transition.md).

- Q0_diag:

  Variance assigned to each component of the initial-state distribution
  (small, reflects near-certainty about `x0`).

- keep_states:

  Logical; if `TRUE`, also return the filtered state trajectory (in
  normalised units) and the predictive mean on the original (head-count)
  scale.

## Value

If `keep_states` is `FALSE` (default), the scalar approximate
log-likelihood. Otherwise a list with elements `loglik`, `filtered`
(matrix of filtered means \\\bar X_k\\, normalised units), `predicted`
(matrix of one-step-ahead predicted means \\\hat X_k\\, normalised
units) and `fitted` (predictive mean of `y`, on the original head-count
scale, i.e. `N * p * predicted[, "i"]`).

## Details

The state-space system is \$\$X_k = F_k(\eta,\Delta) +
A\_{k-1}(\eta,\Delta) X\_{k-1} + U_k, \quad U_k \sim N_3(0,
Q_k(\eta,\Delta)),\$\$ \$\$Y_k/N = p \\ X\_{k,2} + V_k, \quad V_k \sim
N(0, R_k(\theta)),\$\$ with \\X_k=(S_k,I_k,R_k)/N\\, observation
operator `B = (0, p, 0)` and \\R_k(\theta) = N^{-1}(p(1-p)+\tau^2) \\
i(\eta,t_k)\\, following the observation model of Narci et al. (2021)
extended to the three-event SIRD jump process.
