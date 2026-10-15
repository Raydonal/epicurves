# Worked example of the R Journal article: a SEIRD epidemic simulated with ga_simulate(),
# fitted with ga_fit() to cases and deaths and to cases only, with profile likelihoods of alpha.
suppressMessages(library(epicurves))
seird <- ga_model(c("S", "E", "I", "R", "D"), list(
  infection   = ga_event("S", "E", "lambda", by = "I"),
  progression = ga_event("E", "I", "sigma"),
  recovery    = ga_event("I", "R", "gamma"),
  death       = ga_event("I", "D", "alpha")))
N <- 50000
truth <- list(lambda = 0.5, sigma = 0.25, gamma = 0.2, alpha = 0.05, p = 0.4, tau2 = 0.2, tau2_d = 0.1)
obs <- list(cases = ga_incidence("progression", prob = "p", noise = "tau2"),
            deaths = ga_incidence("death", prob = 1, noise = "tau2_d"))
set.seed(2026)
sim <- ga_simulate(seird, truth, c(S = N - 20, E = 0, I = 20, R = 0, D = 0), obs, n = 120)
init <- function(p) c(S = 1 - p$i0, E = 0, I = p$i0, R = 0, D = 0)
start <- list(lambda = 0.4, sigma = 0.3, gamma = 0.15, alpha = 0.08, i0 = 1e-3, p = 0.5, tau2 = 0.3, tau2_d = 0.3)
t0 <- proc.time()[3]
fit_both <- ga_fit(seird, sim$data, obs, N, init, start = start)
secs_both <- proc.time()[3] - t0
start_c <- start[setdiff(names(start), "tau2_d")]
cases_only <- sim$data; cases_only[, 2] <- NA
fit_cases <- ga_fit(seird, cases_only, obs, N, init, start = start_c, fixed = list(tau2_d = 0.1))
grid <- c(0.05, 0.04, 0.03, 0.02, 0.01, 0.06, 0.08, 0.1, 0.12, 0.15)
prof_both <- ga_profile(fit_both, "alpha", grid)
prof_cases <- ga_profile(fit_cases, "alpha", grid)
saveRDS(list(sim = sim, truth = truth, N = N, fit_both = fit_both, fit_cases = fit_cases,
             secs_both = secs_both, prof_both = prof_both, prof_cases = prof_cases),
        "results/example_seird.rds")
print(fit_both); print(fit_cases)
cat("seconds (both):", secs_both, "\n")
print(prof_both[order(prof_both$value), c("value", "loglik", "gamma")])
print(prof_cases[order(prof_cases$value), c("value", "loglik", "gamma")])
