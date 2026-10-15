# Parametric bootstrap of the SEIRD example fit (cases and deaths) with ga_bootstrap():
# B = 200 replicates, statistic = estimates plus infection fatality ratio and basic reproduction
# number; compared with the Wald intervals of the fit. Output: results/example_bootstrap.{rds,log}.
suppressMessages(library(epicurves))
ex <- readRDS("results/example_seird.rds")
fit <- ex$fit_both
stat <- function(f) { th <- f$theta
  c(coef(f), ifr = th$alpha / (th$gamma + th$alpha), R0 = th$lambda / (th$gamma + th$alpha)) }
t0 <- proc.time()[3]
bt <- ga_bootstrap(fit, B = 200, statistic = stat, cores = 4, seed = 2026)
secs <- proc.time()[3] - t0
truth <- with(ex$truth, c(lambda = lambda, sigma = sigma, gamma = gamma, alpha = alpha, i0 = 20 / ex$N,
                          p = p, tau2 = tau2, tau2_d = tau2_d, ifr = alpha / (gamma + alpha),
                          R0 = lambda / (gamma + alpha)))
wald <- confint(fit)
saveRDS(list(boot = bt, wald = wald, truth = truth, secs = secs), "results/example_bootstrap.rds")
print(bt)
cat("seconds:", secs, "\n")
pc <- confint(bt)
print(data.frame(truth = signif(truth, 4), est = signif(bt$t0, 4), boot_l = signif(pc[, 1], 4), boot_u = signif(pc[, 2], 4),
                 wald_l = signif(c(wald[, 1], NA, NA), 4), wald_u = signif(c(wald[, 2], NA, NA), 4)))
