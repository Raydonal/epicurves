# Maximum-likelihood estimation for a general jump model

Dimension-generic analogue of
[`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md).
All explicit compartments except one (`free_init`) are held fixed at
`start`'s values for the initial condition; by default only the observed
compartment's initial fraction is estimated, with all others fixed
(mirroring
[`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md)'s
treatment of `r0`).

## Usage

``` r
general_ga_fit(
  y,
  N,
  model,
  obs_compartment,
  Delta = 1,
  start,
  free_init = obs_compartment,
  p_max = 0.999,
  n_sub = 10L,
  Q0_diag = 1e-08,
  control = list(maxit = 3000, reltol = 1e-10),
  hessian = FALSE
)
```

## Arguments

- y, N, Delta:

  As in
  [`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md).

- model:

  A `jump_model`.

- obs_compartment:

  Name of the observed compartment.

- start:

  Named list with one entry per rate parameter, one entry per
  compartment in `model$compartments` (the initial condition), `p` and
  `tau2`.

- free_init:

  Character vector of compartment names whose initial fraction is
  estimated (default: just `obs_compartment`); all other compartments
  are held fixed at their `start` value, and the "susceptible-like"
  first compartment is always treated as
  `1 - sum( other explicit compartments)` so the initial condition stays
  on the simplex.

- p_max, n_sub, Q0_diag, control, hessian:

  As in
  [`sird_ga_fit()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit.md).
