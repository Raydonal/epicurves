# Solve the deterministic normalised SIRD ODE system

Solve the deterministic normalised SIRD ODE system

## Usage

``` r
sird_ode_solve(times, x0, parms)
```

## Arguments

- times:

  Numeric vector of output times (must include `t0`).

- x0:

  Initial state `c(s0, i0, r0)`.

- parms:

  Named list/vector with `lambda`, `gamma`, `alpha`.

## Value

A data frame with columns `time`, `s`, `i`, `r`, `d` (data frame as
returned by
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html), with the
implicit death fraction added).
