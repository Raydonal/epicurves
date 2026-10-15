# Exact stochastic simulation of a general jump model (Gillespie algorithm)

Dimension-generic analogue of
[`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md),
driven by a
[`new_jump_model()`](https://raydonal.github.io/epicurves/reference/new_jump_model.md)
of any dimension.

## Usage

``` r
general_gillespie(
  model,
  N,
  x0_counts,
  parms,
  Tmax,
  Delta = 1,
  max_events = 50L * N
)
```

## Arguments

- model:

  A `jump_model`.

- N:

  Population size.

- x0_counts:

  Named integer vector of initial head counts, one per
  `model$compartments`.

- parms:

  Named list of rate parameters.

- Tmax, Delta:

  As in
  [`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md).

- max_events:

  Safety cap on the number of simulated events.

## Value

A data frame with column `time` and one column per compartment in
`model$compartments` (head counts), recorded on
`seq(0, Tmax, by = Delta)`.
