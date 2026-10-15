# Exact simulation of a compartmental model and its surveillance counts

Gillespie simulation of the Markov jump process defined by `model`, with
counts aggregated over intervals of length `Delta` and reported counts
drawn as binomial thinning plus Gaussian measurement error (rounded and
truncated at zero), the data-generating process assumed by
[`ga_loglik()`](https://raydonal.github.io/epicurves/reference/ga_loglik.md).

## Usage

``` r
ga_simulate(model, theta, init_counts, obs, n, Delta = 1, varying = NULL)
```

## Arguments

- model:

  A
  [`ga_model()`](https://raydonal.github.io/epicurves/reference/ga_model.md).

- theta:

  Named list or vector with the rate parameters of `model` and any
  parameters named in `obs`.

- init_counts:

  Named integer vector of initial counts per compartment.

- obs:

  List of
  [`ga_incidence()`](https://raydonal.github.io/epicurves/reference/ga_incidence.md)
  /
  [`ga_prevalence()`](https://raydonal.github.io/epicurves/reference/ga_incidence.md)
  streams.

- n:

  Number of intervals.

- Delta:

  Length of the observation interval.

- varying:

  Optional named list of rate parameters that change between intervals,
  each a vector of length `n` (as in
  [`ga_loglik()`](https://raydonal.github.io/epicurves/reference/ga_loglik.md));
  a rate listed here may be absent from `theta`.

## Value

A list with `data` (matrix of reported counts, one column per stream),
`latent` (matrix of latent counts per stream) and `states` (compartment
counts at the end of each interval, first row the initial state).
