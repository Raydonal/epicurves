# Maximum likelihood for a compartmental model under the Gaussian approximation

Maximises
[`ga_loglik()`](https://raydonal.github.io/epicurves/reference/ga_loglik.md)
over the parameters in `start`, on an unconstrained scale (logarithms of
positive parameters, logits of those in `unit`), by a Nelder-Mead search
(relative tolerance `1e-8`) from each starting point followed by BFGS
refinement of the best optimum. Parameters in `fixed` are held at their
values.

## Usage

``` r
ga_fit(
  model,
  data,
  obs,
  N,
  init,
  start,
  fixed = list(),
  unit = NULL,
  blocks = NULL,
  Delta = 1,
  nsub = 10L,
  hessian = TRUE,
  maxit = 1500
)
```

## Arguments

- model:

  A
  [`ga_model()`](https://raydonal.github.io/epicurves/reference/ga_model.md).

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

- init:

  A function of the parameter list returning the named vector of initial
  fractions, e.g. `function(p) c(S = 1 - p$i0, I = p$i0, R = 0, D = 0)`.

- start:

  Named list of starting values of the free parameters, or a list of
  such lists (several starts). A parameter listed in `blocks` takes a
  vector with one value per block (a scalar is recycled).

- fixed:

  Named list of parameters held fixed.

- unit:

  Names of free parameters constrained to `(0, 1)`. By default, the
  reporting probabilities named in `obs` and any free parameter whose
  name ends in `0` (initial fractions such as `i0`).

- blocks:

  Optional named list giving, for a time-varying rate parameter, an
  integer vector with the block (`1, ..., B`) of each observation
  interval; the parameter then takes one value per block.

- Delta:

  Length of the observation interval.

- nsub:

  Runge-Kutta steps per interval.

- hessian:

  Logical; compute the numerical Hessian at the optimum.

- maxit:

  Iteration limit of each Nelder-Mead run.

## Value

An object of class `ga_fit` with elements `coef` (estimates on the
natural scale, block-specific values named `name[b]`), `loglik`,
`convergence` (the `optim` code of the final run), `par` and `hessian`
(unconstrained scale), `vcov` (on that scale), `starts` (log-likelihood
reached from each start) and the inputs needed by
[`ga_profile()`](https://raydonal.github.io/epicurves/reference/ga_profile.md),
[`confint.ga_fit()`](https://raydonal.github.io/epicurves/reference/confint.ga_fit.md)
and
[`predict.ga_fit()`](https://raydonal.github.io/epicurves/reference/predict.ga_fit.md).
