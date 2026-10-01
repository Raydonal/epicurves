# epicurves

Reference implementation accompanying two companion manuscripts derived from
Oliveira's doctoral thesis on COVID-19 epidemic curve modelling (Ospina,
Oliveira & Ferraz):

- **SIRD-GA**: a stochastic SIRD model estimated by Gaussian-diffusion
  approximation of the underlying density-dependent Markov jump process,
  with maximum-likelihood estimation via a Kalman filter under partial
  observation and under-reporting.
- **Piecewise ensemble growth model**: Logistic, Gompertz and Richards
  growth curves combined by MSE- or AIC-based weights, with structural
  change points ("partition points") identified by sequential grid search.

The package reproduces every simulation study, table and figure reported in
the two manuscripts, and includes a `paraiba_covid19` dataset (the real
COVID-19 surveillance series for Paraíba, Brazil, used in both
applications).

## Installation

```r
# from a local clone:
devtools::install("path/to/epicurves")
```

## Usage

```r
library(epicurves)

# SIRD-GA: simulate, then estimate
sim <- sird_ga_simulate_obs(N = 20000, S0 = 19980, I0 = 20, R0 = 0,
                             lambda = 0.4, gamma = 0.15, alpha = 0.08,
                             p = 0.6, tau2 = 0.2, Tmax = 100)
fit <- sird_ga_fit(sim$y, N = 20000,
                    start = list(lambda = 0.3, gamma = 0.1, alpha = 0.1,
                                 i0 = 0.002, p = 0.5, tau2 = 0.3))
sird_R0(fit$theta[["lambda"]], fit$theta[["gamma"]], fit$theta[["alpha"]])

# Piecewise ensemble: identify structural change points
data(paraiba_covid19)
pw <- fit_piecewise_ensemble(paraiba_covid19$deaths[paraiba_covid19$deaths > 0 |
                                                       seq_along(paraiba_covid19$deaths) > 1])
pw$cuts
```

## Status

Research code accompanying manuscripts in preparation. See `NEWS.md` for
notable methodological findings made while validating this package against
the thesis (including a sign correction to the Richards growth-curve
equation and a changepoint-confirmation robustness analysis).
