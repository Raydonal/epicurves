# Define a density-dependent Markov jump compartmental model

A general representation of a compartmental epidemic model as a set of
elementary events, each with an integer (here, real-valued after
normalisation) jump vector and a rate function of the current state and
parameters.
[`sird_drift()`](https://raydonal.github.io/epicurves/reference/sird_drift.md),
[`sird_jacobian()`](https://raydonal.github.io/epicurves/reference/sird_jacobian.md)
and
[`sird_diffusion()`](https://raydonal.github.io/epicurves/reference/sird_diffusion.md)
are the SIRD-specific, hand-derived analogues of
[`jm_drift()`](https://raydonal.github.io/epicurves/reference/jm_drift.md),
[`jm_jacobian()`](https://raydonal.github.io/epicurves/reference/jm_jacobian.md)
and
[`jm_diffusion()`](https://raydonal.github.io/epicurves/reference/jm_diffusion.md)
below; a `jump_model` makes the same Gaussian-diffusion-approximation
machinery available for **any** compartmental structure (SEIR, SEIRD,
SIRDV, age- or risk-stratified models, ...), at the cost of a numerical
rather than analytic Jacobian.
[`general_ga_fit()`](https://raydonal.github.io/epicurves/reference/general_ga_fit.md)
and friends use a `jump_model` exactly as
[`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md)
uses the hand-derived SIRD functions.

## Usage

``` r
new_jump_model(compartments, events)
```

## Arguments

- compartments:

  Character vector naming the explicit state compartments (the model's
  state space has one *implicit* compartment, recovered as
  `1 - sum(explicit compartments)`, exactly as `d` is implicit in
  [`sird_drift()`](https://raydonal.github.io/epicurves/reference/sird_drift.md)).

- events:

  A list of elementary events, each a list with `jump` (a numeric vector
  the same length as `compartments`, giving the change to each explicit
  compartment when the event fires) and `rate` (a function of
  `(state, parms)` returning the event's instantaneous rate).

## Value

An object of class `jump_model`.
