# SIRD-GA from daily cases and deaths

``` r
library(epicurves)
```

This vignette shows the functions behind the article on what reported
cases and deaths can identify in a stochastic SIRD model (Ospina,
Oliveira, Ferraz and Gondim). The heavier computations of the article
are shown but not run.

## Simulated cases and deaths

[`sird_ga_simulate_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_simulate_inc.md)
simulates the SIRD jump process exactly and turns its daily increments
into reported cases (binomial thinning with probability `p` plus
measurement error with variance `tau2` times the latent count) and
reported deaths (probability `q`, variance scale `tau2_d`).

``` r
set.seed(1)
N <- 20000
si <- sird_ga_simulate_inc(N, S0 = N - 10, I0 = 10, R0 = 0, lambda = 0.6, gamma = 0.3,
                           alpha = 0.1, p = 0.5, tau2 = 0.16, q = 1, tau2_d = 0.1, Tmax = 100)
c(cases = sum(si$cases), deaths = sum(si$deaths))
#>  cases deaths 
#>   5826   2988
```

## Fitting to cases and deaths

[`sird_ga_fit_inc()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit_inc.md)
maximises the Kalman-filter log-likelihood of the Gaussian
approximation, with the death reporting probability held at `q = 1` (the
product of `q` and the death rate is what the death counts identify).

``` r
fit <- sird_ga_fit_inc(si$cases, si$deaths, N, use = "both",
                       starts = list(list(lambda = 0.5, gamma = 0.25, alpha = 0.15, i0 = 5e-4,
                                          p = 0.6, tau2 = 0.3, tau2_d = 0.3)),
                       nsub = 5L)
round(fit$theta[c("lambda", "gamma", "alpha", "p")], 3)
#> lambda  gamma  alpha      p 
#>  0.662  0.323  0.110  0.491
fit$loglik
#> [1] 1536.861
```

## The death rate is not identifiable from cases alone

With cases only, recovery and death enter the likelihood through
`gamma + alpha`. The profile log-likelihood of `alpha` is flat when
fitted to cases and has a maximum when deaths are added:

``` r
fit_c <- sird_ga_fit_inc(si$cases, NULL, N, use = "cases",
                         starts = list(list(lambda = 0.5, gamma = 0.25, alpha = 0.15,
                                            i0 = 5e-4, p = 0.6, tau2 = 0.3)),
                         free = c("lambda", "gamma", "alpha", "i0", "p", "tau2"))
grid <- c(0.01, 0.05, 0.1, 0.2, 0.3)
sird_ga_profile(fit_c, si$cases, NULL, "alpha", grid)
sird_ga_profile(fit, si$cases, si$deaths, "alpha", grid)
```

## Paraíba: weekly data and transmission changing every four weeks

The application aggregates the daily counts of 2020 into complete
epidemiological weeks (daily counts carry a strong weekly reporting
cycle) and lets the transmission rate change every four weeks with
[`sird_ga_fit_tv()`](https://raydonal.github.io/epicurves/reference/sird_ga_fit_tv.md).
Rates are then per week.

``` r
data(paraiba_covid19)
data(paraiba_population)
d <- paraiba_covid19
wk <- d$epi_week
keep <- wk %in% names(which(table(wk) == 7))
cases <- as.numeric(tapply(d$new_cases[keep], wk[keep], sum))
deaths <- as.numeric(tapply(d$new_deaths[keep], wk[keep], sum))
blk <- ceiling(seq_along(cases) / 4)
c(weeks = length(cases), blocks = max(blk))
#>  weeks blocks 
#>     40     10
```

``` r
rem <- 7 / 7  # removal rate per week for a mean infectious period of 7 days
st <- list(lambda = rep(1.2 * rem, max(blk)), gamma = 0.98 * rem, alpha = 0.02 * rem,
           i0 = 5e-6, p = 0.2, tau2 = 20, tau2_d = 2)
f <- sird_ga_fit_tv(cases, deaths, paraiba_population, blk, "both", st, nsub = 20L)
7 / (f$theta[["gamma"]] + f$theta[["alpha"]])   # mean infectious period, days
tapply(f$R_t, blk, mean)                         # reproduction number by block
```

The same model can be written with the general engine
([`vignette("compartmental-models")`](https://raydonal.github.io/epicurves/articles/compartmental-models.md)):
[`ga_model()`](https://raydonal.github.io/epicurves/reference/ga_model.md)
with infection, recovery and death events,
[`ga_incidence()`](https://raydonal.github.io/epicurves/reference/ga_incidence.md)
streams for infection and death, and `blocks = list(lambda = blk)` in
[`ga_fit()`](https://raydonal.github.io/epicurves/reference/ga_fit.md).
