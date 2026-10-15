# Gaussian-approximation log-likelihood of a compartmental model

Kalman-filter log-likelihood of counts observed at times
`Delta, 2 Delta, ..., n Delta` under the Gaussian (linear noise)
approximation of the density-dependent Markov jump process defined by
`model`. Over each interval the ordinary differential equation of the
mean, its resolvent and the Lyapunov equation of the state-noise
covariance are integrated jointly by the classical Runge-Kutta method in
compiled code. Incidence observations are handled with event counters
and an augmented state that contains the state at the start and the end
of each interval.

## Usage

``` r
ga_loglik(
  model,
  theta,
  init,
  data,
  obs,
  N,
  Delta = 1,
  nsub = 10L,
  varying = NULL,
  full = FALSE
)
```

## Arguments

- model:

  A
  [`ga_model()`](https://raydonal.github.io/epicurves/reference/ga_model.md).

- theta:

  Named list or vector with the rate parameters of `model` and any
  parameters named in `obs`.

- init:

  Named vector of initial fractions, one per compartment (normalised to
  sum to one).

- data:

  Matrix or data frame of counts, one column per element of `obs` (in
  the same order), one row per interval. `NA` marks a missing count.

- obs:

  List of
  [`ga_incidence()`](https://raydonal.github.io/epicurves/reference/ga_incidence.md)
  /
  [`ga_prevalence()`](https://raydonal.github.io/epicurves/reference/ga_incidence.md)
  streams.

- N:

  Population size.

- Delta:

  Length of the observation interval.

- nsub:

  Runge-Kutta steps per interval.

- varying:

  Optional named list of time-varying rate parameters: each element is a
  vector with one value per interval, overriding the constant value in
  `theta`.

- full:

  If `TRUE`, return the list with one-step-ahead predicted means
  (`fitted`) and standard deviations (`sd`) of the counts and the path
  of the mean (`x`, fractions) instead of the scalar.

## Value

The log-likelihood, or a list if `full = TRUE`.
