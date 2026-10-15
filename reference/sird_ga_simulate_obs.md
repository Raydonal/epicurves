# Simulate partially-observed SIRD-GA surveillance data

Wraps
[`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md)
with the two-part observation-error model (binomial under-reporting plus
Gaussian measurement error) used throughout the SIRD-GA manuscript.

## Usage

``` r
sird_ga_simulate_obs(
  N,
  S0,
  I0,
  R0,
  lambda,
  gamma,
  alpha,
  p,
  tau2,
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

- p:

  Reporting probability in `(0,1)`.

- tau2:

  Measurement-error variance parameter.

- Tmax:

  Final observation time.

- Delta:

  Spacing of the recorded observation grid (default `1`).

## Value

A list with `latent` (output of
[`sird_gillespie()`](https://raydonal.github.io/epicurves/reference/sird_gillespie.md))
and `y` (numeric vector of simulated observations, one per grid time).
