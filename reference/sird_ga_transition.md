# Solve the SIRD-GA state-space transition quantities

Jointly solves the deterministic SIRD trajectory \\x(t)\\ and the
fundamental (sensitivity) matrix of the linearised dynamics in a
**single** [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html)
call spanning the whole observation window, using a *reset event* that
sets the fundamental matrix back to the identity at the start of every
observation interval. This gives the same locally-anchored resolvent
\\\Phi(s, t\_{k-1})\\ as restarting the solver once per interval (and
the same numerical stability — the matrix never has to span more than
one interval, which matters for long series), while paying the fixed
per-call overhead of
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) only once
instead of `n` times. This is the single biggest cost in repeated
likelihood evaluation (MLE optimisation, bootstrap), so the single-call
form is the default implementation.

## Usage

``` r
sird_ga_transition(times, x0, parms, n_sub = 10L)
```

## Arguments

- times:

  Numeric vector of observation times, `t_0 < t_1 < ... < t_n`
  (typically `0:n`).

- x0:

  Initial state `c(s0, i0, r0)` at `times[1]`.

- parms:

  Named list/vector with `lambda`, `gamma`, `alpha`.

- n_sub:

  Number of internal subdivisions per observation interval, used for the
  [`sird_ga_Qk()`](https://raydonal.github.io/epicurves/reference/sird_ga_Qk.md)
  quadrature. Defaults to 10 (even, required by Simpson's rule).

## Value

A list of class `sird_ga_transition` with components `x_obs` (matrix of
states at `t_0, ..., t_n`), `segments` (one list per interval, with
fine-grid `t`, `x` and `P`, `P[,,1]` always the identity), and `parms`.
