# Numerical Jacobian of a jump model's drift

Central-difference Jacobian of
[`jm_drift()`](https://raydonal.github.io/epicurves/reference/jm_drift.md).
Used in place of an analytic Jacobian (as hand-derived for SIRD in
[`sird_jacobian()`](https://raydonal.github.io/epicurves/reference/sird_jacobian.md))
so that
[`new_jump_model()`](https://raydonal.github.io/epicurves/reference/new_jump_model.md)
generalises to arbitrary compartmental structures without per-model
algebra. The cost is a small, well- controlled numerical error (step
`h`), validated in the package tests against the analytic SIRD Jacobian.

## Usage

``` r
jm_jacobian(model, state, parms, h = 1e-06)
```

## Arguments

- model:

  A `jump_model`.

- state, parms:

  As in
  [`sird_drift()`](https://raydonal.github.io/epicurves/reference/sird_drift.md).

- h:

  Finite-difference step.
