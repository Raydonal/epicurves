# Profile log-likelihood of alpha: cases only vs cases + deaths (Scenario 1, N = 20000).
suppressMessages(library(epicurves))
set.seed(2024)
N <- 20000; I0 <- 10
repeat { si <- sird_ga_simulate_inc(N, N - I0, I0, 0, 0.6, 0.3, 0.1, p = 0.5, tau2 = 0.16,
                                    q = 1, tau2_d = 0.1, Tmax = 100)
         if (sum(si$new_inf) >= 0.05 * N) break }
st <- list(lambda = 0.6, gamma = 0.3, alpha = 0.1, i0 = I0 / N, p = 0.5, tau2 = 0.16, tau2_d = 0.1)
fb <- sird_ga_fit_inc(si$cases, si$deaths, N, "both", starts = list(st), nsub = 10L)
fc <- sird_ga_fit_inc(si$cases, NULL, N, "cases", starts = list(st),
                      free = c("lambda", "gamma", "alpha", "i0", "p", "tau2"), nsub = 10L)
grid <- c(0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3)
pb <- sird_ga_profile(fb, si$cases, si$deaths, "alpha", grid)
pc <- sird_ga_profile(fc, si$cases, NULL, "alpha", grid)
out <- rbind(data.frame(data = "cases and deaths", pb[, c("value", "loglik", "gamma", "lambda")], max = fb$loglik),
             data.frame(data = "cases only", pc[, c("value", "loglik", "gamma", "lambda")], max = fc$loglik))
print(out, digits = 6)
saveRDS(list(profile = out, fit_both = fb$theta, fit_cases = fc$theta, sim = si),
        "results/profile_alpha_sim.rds")
