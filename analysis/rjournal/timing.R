# Cost of one likelihood evaluation: hand-coded SIRD engine vs generic engine, for several
# population sizes and model sizes (median of repeated evaluations, single core).
suppressMessages(library(epicurves))
sird <- ga_model(c("S", "I", "R", "D"), list(
  infection = ga_event("S", "I", "lambda", by = "I"),
  recovery  = ga_event("I", "R", "gamma"),
  death     = ga_event("I", "D", "alpha")))
seird <- ga_model(c("S", "E", "I", "R", "D"), list(
  infection   = ga_event("S", "E", "lambda", by = "I"),
  progression = ga_event("E", "I", "sigma"),
  recovery    = ga_event("I", "R", "gamma"),
  death       = ga_event("I", "D", "alpha")))
th <- list(lambda = 0.5, sigma = 0.25, gamma = 0.2, alpha = 0.05, p = 0.4, tau2 = 0.2, tau2_d = 0.1,
           q = 1, s0 = 0.999, i0 = 0.001, r0 = 0)
out <- list()
for (N in c(1e4, 1e5, 1e6)) for (n in c(100, 200)) {
  set.seed(1)
  dat <- cbind(rpois(n, 50), rpois(n, 2))
  b <- bench::mark(
    sird_hand = sird_ga_loglik_inc(th, dat[, 1], dat[, 2], N, "both"),
    sird_generic = ga_loglik(sird, th, c(S = 0.999, I = 0.001, R = 0, D = 0), dat,
                             list(ga_incidence("infection", "p", "tau2"), ga_incidence("death", 1, "tau2_d")), N),
    seird_generic = ga_loglik(seird, th, c(S = 0.999, E = 0, I = 0.001, R = 0, D = 0), dat,
                              list(ga_incidence("progression", "p", "tau2"), ga_incidence("death", 1, "tau2_d")), N),
    check = FALSE, min_iterations = 30)
  out[[length(out) + 1]] <- data.frame(N = N, n = n, engine = as.character(b$expression),
                                       ms = as.numeric(b$median) * 1000)
}
res <- do.call(rbind, out)
saveRDS(res, "results/timing.rds")
print(res)
