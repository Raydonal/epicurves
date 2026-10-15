# Parametric bootstrap for a `ga_fit`

Simulates `B` data sets from the fitted model with
[`ga_simulate()`](https://raydonal.github.io/epicurves/reference/ga_simulate.md),
that is, from the Markov jump process and the observation model assumed
by
[`ga_loglik()`](https://raydonal.github.io/epicurves/reference/ga_loglik.md),
refits each by
[`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
started at the estimates, and records a statistic of every refit.
Missing values in the original data are kept missing in the simulated
data sets, and rates that vary between blocks are simulated with their
block values.

## Usage

``` r
ga_bootstrap(
  fit,
  B = 200,
  statistic = stats::coef,
  init_counts = NULL,
  cores = 1L,
  seed = NULL,
  maxit = 1500
)
```

## Arguments

- fit:

  A
  [`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md)
  result.

- B:

  Number of bootstrap replicates.

- statistic:

  Function of a `ga_fit` returning a named numeric vector, applied to
  `fit` and to every refit. The default returns the estimates.

- init_counts:

  Named vector of initial counts per compartment, summing to the
  population size of the fit. By default the initial fractions at the
  estimates times the population size, rounded, with the rounding
  remainder added to the largest compartment.

- cores:

  Number of processes for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  which forks; use 1 on Windows.

- seed:

  Optional seed set before the replicate seeds are drawn.

- maxit:

  Nelder-Mead iteration limit of each refit.

## Value

An object of class `ga_boot` with elements `t0` (statistic at `fit`),
`t` (matrix with one row per replicate, `NA` where the simulation or the
refit failed), `loglik` and `convergence` (the `optim` code of each
refit) and `init_counts`.

## Details

Each replicate draws its own seed from the current random-number stream,
so the result is reproducible under
[`set.seed()`](https://rdrr.io/r/base/Random.html) (or `seed`) and does
not depend on `cores`.

## See also

[`confint.ga_boot()`](https://raydonal.github.io/epicurves/reference/confint.ga_boot.md),
[`ga_profile()`](https://raydonal.github.io/epicurves/reference/ga_profile.md).

## Examples

``` r
# \donttest{
sir <- ga_model(c("S", "I", "R"), list(
  infection = ga_event("S", "I", "lambda", by = "I"),
  removal   = ga_event("I", "R", "gamma")))
obs <- list(cases = ga_incidence("infection", prob = "p", noise = "tau2"))
set.seed(1)
sim <- ga_simulate(sir, list(lambda = 0.5, gamma = 0.2, p = 0.5, tau2 = 0.1),
                   c(S = 4990, I = 10, R = 0), obs, n = 60)
init <- function(p) c(S = 1 - p$i0, I = p$i0, R = 0)
fit <- ga_fit(sir, sim$data, obs, 5000, init,
              start = list(lambda = 0.4, gamma = 0.25, i0 = 0.005, p = 0.4, tau2 = 0.2))
bt <- ga_bootstrap(fit, B = 20, seed = 2)
confint(bt)
#>               2.5 %      97.5 %
#> lambda 3.113220e-01 0.534183195
#> gamma  5.521845e-02 0.246993159
#> i0     1.283805e-03 0.005556414
#> p      4.402647e-01 0.530667869
#> tau2   6.413786e-07 0.423937863
# }
```
