# Solve the Gaussian-approximation transition quantities for a general jump model

Dimension-generic analogue of
[`sird_ga_transition()`](https://raydonal.github.io/epicurves/reference/sird_ga_transition.md):
identical event-reset construction (a single
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) call per
model, with the fundamental matrix reset to the identity at the start of
each observation interval), but driven by a
[`new_jump_model()`](https://raydonal.github.io/epicurves/reference/new_jump_model.md)
object of any dimension `d = length(model$compartments)` instead of the
hard-coded `d = 3` SIRD case.
[`sird_ga_transition()`](https://raydonal.github.io/epicurves/reference/sird_ga_transition.md)
remains the validated, hand-derived implementation used for the main
SIRD-GA results; this generic engine is what
[`general_ga_fit()`](https://raydonal.github.io/epicurves/reference/general_ga_fit.md)
uses for extensions (SEIRD, SIRDV, ...).

## Usage

``` r
general_ga_transition(times, x0, model, parms, n_sub = 10L)

general_ga_Ak(tr, k)

general_ga_Fk(tr, k)

general_ga_Qk(tr, k, N)
```

## Arguments

- times, x0, n_sub:

  As in
  [`sird_ga_transition()`](https://raydonal.github.io/epicurves/reference/sird_ga_transition.md).

- model:

  A `jump_model` (see
  [`new_jump_model()`](https://raydonal.github.io/epicurves/reference/new_jump_model.md)).

- parms:

  Named list/vector of the model's rate parameters.

- tr:

  Object returned by `general_ga_transition()`.

- k:

  Integer interval index.

- N:

  Population size.

## Functions

- `general_ga_Ak()`: Resolvent matrix of interval `k`.

- `general_ga_Fk()`: Drift increment of interval `k`.

- `general_ga_Qk()`: State-noise covariance of interval `k` (Simpson's
  rule).
