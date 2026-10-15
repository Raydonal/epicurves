# Simulate partially-observed surveillance data from a general jump model

Simulate partially-observed surveillance data from a general jump model

## Usage

``` r
general_ga_simulate_obs(
  model,
  N,
  x0_counts,
  parms,
  obs_compartment,
  p,
  tau2,
  Tmax,
  Delta = 1
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

- obs_compartment:

  Compartment observed (with under-reporting/noise).

- p, tau2:

  As in
  [`sird_ga_simulate_obs()`](https://raydonal.github.io/epicurves/reference/sird_ga_simulate_obs.md).

- Tmax, Delta:

  As in
  [`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md).
