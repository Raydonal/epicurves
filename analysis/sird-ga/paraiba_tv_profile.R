# Profile likelihood over the case reporting probability p (time-varying SIRD-GA, weekly data).
suppressMessages(library(epicurves))
r <- readRDS("results/paraiba_tv.rds")
best <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
grid <- c(0.02, 0.03, 0.045, 0.066, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8)
prof <- parallel::mclapply(grid, function(pp) {
  st <- c(list(lambda = best$lambda), as.list(best$theta))
  f <- tryCatch(sird_ga_fit_tv(r$cases, r$deaths, r$N, r$blk, "both", st,
                               fixed = list(p = pp), nsub = 20L),
                error = function(e) { message("p=", pp, ": ", conditionMessage(e)); NULL })
  if (is.null(f)) return(data.frame(p = pp, loglik = NA, inf_period = NA, ifr = NA, R1 = NA,
                                    Rlast = NA, attack = NA))
  data.frame(p = pp, loglik = f$loglik, inf_period = 7 / (f$theta[["gamma"]] + f$theta[["alpha"]]),
             ifr = f$theta[["alpha"]] / (f$theta[["gamma"]] + f$theta[["alpha"]]),
             R1 = mean(f$R_t[r$blk == 1]), Rlast = mean(f$R_t[r$blk == max(r$blk)]),
             attack = 1 - f$kf$x[nrow(f$kf$x), 1])
}, mc.cores = 2)
prof <- do.call(rbind, prof)
print(prof, digits = 4)
saveRDS(prof, "results/paraiba_tv_profile_p.rds")
