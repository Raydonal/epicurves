# Profile log-likelihood of p for the daily day-of-week fit (M1 of paraiba_daily_dow.R). The M1 fits
# reach two families of optima within 0.2 log-likelihood units (p near 0.10 and near 0.26), with an
# indefinite Hessian. At each p of the grid the model is refitted with p fixed from two warm starts,
# the best M1 fit of each family with p overwritten; the higher fit is kept.
# Each grid point is saved to results/paraiba_daily_prof_<p>.rds (skipped if present); the summary
# goes to stdout (-> results/paraiba_daily_profile_p.log).
# Usage: Rscript paraiba_daily_profile_p.R [cores]
suppressMessages(library(epicurves))
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
options(width = 200)

r <- readRDS("results/paraiba_daily_dow.rds")
fits <- lapply(list.files("results", "^paraiba_daily_M1_.*\\.rds$", full.names = TRUE), readRDS)
fits <- Filter(function(f) is.null(f$error) && is.finite(f$loglik), fits)
low <- Filter(function(f) f$theta[["p"]] < 0.18, fits)
high <- Filter(function(f) f$theta[["p"]] >= 0.18 & f$theta[["p"]] < 0.5, fits)
pick <- function(l) l[[which.max(sapply(l, `[[`, "loglik"))]]
fam <- list(low = pick(low), high = pick(high))

as_start <- function(f) c(list(lambda = f$lambda), as.list(f$theta[c("gamma", "alpha", "i0", "tau2", "tau2_d")]),
                          list(dow_cases = f$dow_weights$cases, dow_deaths = f$dow_weights$deaths))
grid <- c(0.04, 0.055, 0.07, 0.085, 0.10, 0.12, 0.15, 0.18, 0.22, 0.26, 0.32, 0.40)
jobs <- expand.grid(p = grid, fam = names(fam), stringsAsFactors = FALSE)

one <- function(i) {
  pv <- jobs$p[i]; fm <- jobs$fam[i]
  out <- sprintf("results/paraiba_daily_prof_%.3f_%s.rds", pv, fm)
  if (file.exists(out)) return(readRDS(out))
  t0 <- proc.time()[3]
  f <- tryCatch(sird_ga_fit_tv(r$cases, r$deaths, r$N, r$blk, "both", as_start(fam[[fm]]),
                               fixed = list(p = pv), Delta = 1, nsub = 5L, dow = r$dow),
                error = function(e) NULL)
  res <- if (is.null(f)) data.frame(p = pv, fam = fm, loglik = NA) else {
    g <- f$theta[["gamma"]] + f$theta[["alpha"]]
    data.frame(p = pv, fam = fm, loglik = f$loglik, conv = f$convergence, inf_period_days = 1 / g,
               ifr_pct = 100 * f$theta[["alpha"]] / g, attack_pct = 100 * (1 - f$kf$x[nrow(f$kf$x), 1]),
               R_block1 = mean(f$R_t[f$blocks == 1]), tau2 = f$theta[["tau2"]], i0N = f$theta[["i0"]] * r$N,
               minutes = round((proc.time()[3] - t0) / 60, 1))
  }
  saveRDS(res, out)
  res
}
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), one, mc.cores = cores, mc.preschedule = FALSE))
best <- do.call(rbind, lapply(split(res, res$p), function(x) x[which.max(x$loglik), ]))
mx <- max(c(best$loglik, sapply(fits, `[[`, "loglik")), na.rm = TRUE)
best$rel <- best$loglik - mx
cat("Profile of p, daily data with day-of-week weights (best of the two warm starts)\n")
cat(sprintf("Unrestricted maximum (best M1 fit): %.2f\n\n", max(sapply(fits, `[[`, "loglik"))))
print(format(best[, c("p", "fam", "rel", "loglik", "conv", "inf_period_days", "ifr_pct", "attack_pct",
                      "R_block1", "tau2", "i0N")], digits = 4), row.names = FALSE)
inside <- best$p[best$rel >= -1.92]
cat(sprintf("\nGrid values within 1.92 of the maximum: %s\n", paste(inside, collapse = ", ")))
cat("\nBoth starts at each p (loglik):\n")
print(reshape(res[, c("p", "fam", "loglik")], idvar = "p", timevar = "fam", direction = "wide"), row.names = FALSE)
saveRDS(list(all = res, best = best), "results/paraiba_daily_profile_p.rds")
cat("DAILYPROFILEDONE\n")
