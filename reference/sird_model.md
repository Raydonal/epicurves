# SIRD compartmental model, as a `jump_model`

The same SIRD model as
[`sird_drift()`](https://raydonal.github.io/epicurves/reference/sird_drift.md)/[`sird_jacobian()`](https://raydonal.github.io/epicurves/reference/sird_jacobian.md)/
[`sird_diffusion()`](https://raydonal.github.io/epicurves/reference/sird_diffusion.md),
expressed as a
[`new_jump_model()`](https://raydonal.github.io/epicurves/reference/new_jump_model.md),
purely so that
[`general_ga_transition()`](https://raydonal.github.io/epicurves/reference/general_ga_transition.md)/[`general_ga_loglik()`](https://raydonal.github.io/epicurves/reference/general_ga_loglik.md)
can be validated against the hand-derived, analytic-Jacobian SIRD-GA
engine (which remains the implementation used for the main SIRD-GA
results). See `test-general-engine.R` for the numerical cross-check.

## Usage

``` r
sird_model()
```

## Value

A `jump_model` with compartments `c("s","i","r")` and parameters
`lambda`, `gamma`, `alpha`.
