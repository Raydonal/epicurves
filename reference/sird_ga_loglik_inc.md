# SIRD-GA log-likelihood for daily reported cases and/or deaths

Kalman-filter approximate log-likelihood of the SIRD-GA state-space
model when the observations are daily *increments* (new reported cases,
new reported deaths), as recorded by routine surveillance systems. Daily
increments are linear in the augmented state \\Z_k=(X_k, X\_{k-1})\\, so
the linear-Gaussian structure is preserved. Transition quantities are
computed in compiled code (RK4 on the joint system for the ODE solution,
the resolvent and the Lyapunov equation of the state-noise covariance).

## Usage

``` r
sird_ga_loglik_inc(
  theta,
  cases = NULL,
  deaths = NULL,
  N,
  use = c("both", "cases", "deaths"),
  Delta = 1,
  nsub = 10L,
  full = FALSE,
  lambda_t = NULL
)
```

## Arguments

- theta:

  Named vector with `lambda`, `gamma`, `alpha`, `s0`, `i0`, `r0`, `p`,
  `tau2`, `q`, `tau2_d`.

- cases, deaths:

  Daily counts for intervals `1, ..., n` (`NA` allowed). Either may be
  `NULL` when `use` excludes it.

- N:

  Population size.

- use:

  Which streams enter the likelihood: `"both"`, `"cases"` or `"deaths"`.

- Delta:

  Observation interval.

- nsub:

  RK4 steps per interval.

- full:

  If `TRUE`, return the full list (fitted one-step-ahead means and
  standard deviations, ODE path) instead of the scalar.

- lambda_t:

  Optional vector of transmission rates, one per observation interval;
  when given it overrides `theta[["lambda"]]` (piecewise-constant
  transmission).
