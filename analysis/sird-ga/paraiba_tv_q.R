# Sensitivity of the time-varying Paraíba fit to the death reporting probability q.
suppressMessages(library(epicurves))
r <- readRDS("results/paraiba_tv.rds")
best <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
st <- c(list(lambda = best$lambda), as.list(best$theta))
out <- do.call(rbind, lapply(c(1, 0.9, 0.8, 0.6), function(q) {
  f <- sird_ga_fit_tv(r$cases, r$deaths, r$N, r$blk, "both", st, q = q, nsub = 20L)
  data.frame(q = q, loglik = f$loglik, inf_period = 7 / (f$theta[["gamma"]] + f$theta[["alpha"]]),
             ifr = f$theta[["alpha"]] / (f$theta[["gamma"]] + f$theta[["alpha"]]), p = f$theta[["p"]],
             R1 = mean(f$R_t[r$blk == 1]), Rlast = mean(f$R_t[r$blk == max(r$blk)]),
             attack = 1 - f$kf$x[nrow(f$kf$x), 1])
}))
print(out, digits = 4)
saveRDS(out, "results/paraiba_tv_q.rds")
