# Exact stochastic simulation of the SIRD jump process (Gillespie algorithm)

Simulates one exact trajectory of the continuous-time, density-dependent
Markov jump process underlying the SIRD-GA model, using the direct
Gillespie method, and returns the compartment counts recorded at an
equally-spaced observation grid.

## Usage

``` r
sird_gillespie(
  N,
  S0,
  I0,
  R0,
  lambda,
  gamma,
  alpha,
  Tmax,
  Delta = 1,
  max_events = 50L * N
)
```

## Arguments

- N:

  Population size (integer head count).

- S0, I0, R0:

  Initial compartment counts (integers, `S0+I0+R0 <= N`).

- lambda, gamma, alpha:

  Transmission, recovery and disease-induced mortality rates.

- Tmax:

  Final observation time.

- Delta:

  Spacing of the recorded observation grid (default `1`).

- max_events:

  Safety cap on the number of simulated events.

## Value

A data frame with columns `time`, `S`, `I`, `R`, `D`, recorded on the
grid `seq(0, Tmax, by = Delta)` (state held constant between events, as
is standard for jump processes).
