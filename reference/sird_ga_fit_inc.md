# Maximum likelihood for SIRD-GA from daily cases and/or deaths

Nelder-Mead maximisation of
[`sird_ga_loglik_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik_inc.md)
on an unconstrained scale (log rates, logit probabilities and
fractions), from one or several starting points; the best optimum is
polished by a final BFGS run. Parameters not listed in `free` are held
at their values in `fixed`.

## Usage

``` r
sird_ga_fit_inc(
  cases = NULL,
  deaths = NULL,
  N,
  use = c("both", "cases", "deaths"),
  starts,
  free = c("lambda", "gamma", "alpha", "i0", "p", "tau2", "tau2_d"),
  fixed = list(s0 = NULL, r0 = 0, q = 1, p = 1, tau2_d = 0),
  p_max = 0.999,
  Delta = 1,
  nsub = 10L,
  maxit = 4000,
  hessian = TRUE
)
```

## Arguments

- cases, deaths, N, use, Delta, nsub:

  See
  [`sird_ga_loglik_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_loglik_inc.md).

- starts:

  A list of named lists of starting values (one per start).

- free:

  Character vector of parameters to estimate, a subset of
  `c("lambda", "gamma", "alpha", "i0", "s0", "p", "tau2", "q", "tau2_d")`.

- fixed:

  Named list with values for the remaining parameters (`r0` is always
  derived as `1 - s0 - i0` when `s0` is free, otherwise taken from
  `fixed`).

- p_max:

  Upper bound for `p`.

- maxit:

  Iteration limit for each Nelder-Mead run.

- hessian:

  Logical; return the numerical Hessian (unconstrained scale) at the
  optimum.

## Value

A list with `theta`, `loglik`, `starts` (data frame of the optimum
reached from each start), `hessian`, `free` and the transformation
functions.
