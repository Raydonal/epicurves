# Simulate daily reported cases and deaths from the SIRD jump process

Runs the exact Gillespie simulator
[`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md)
and turns the latent daily increments into surveillance counts: reported
cases are `Binomial(new infections, p)` plus
`N(0, tau2 * new infections)` measurement error, reported deaths are
`Binomial(new deaths, q)` plus `N(0, tau2_d * new deaths)`. Counts are
rounded and truncated at zero.

## Usage

``` r
sird_ga_simulate_inc(
  N,
  S0,
  I0,
  R0,
  lambda,
  gamma,
  alpha,
  p,
  tau2,
  q = 1,
  tau2_d = 0,
  Tmax,
  Delta = 1
)
```

## Arguments

- N:

  Population size (integer head count).

- S0, I0, R0:

  Initial compartment counts (integers, `S0+I0+R0 <= N`).

- lambda, gamma, alpha:

  Transmission, recovery and disease-induced mortality rates.

- p, tau2:

  Case reporting probability and case measurement-error scale.

- q, tau2_d:

  Death reporting probability and death measurement-error scale.

- Tmax:

  Final observation time.

- Delta:

  Spacing of the recorded observation grid (default `1`).

## Value

A list with `latent` (the Gillespie output), `cases` and `deaths`
(numeric vectors of length `Tmax / Delta`, one count per interval
`(t_{k-1}, t_k]`), `new_inf` and `new_deaths` (latent increments).
