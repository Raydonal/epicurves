# Paraíba: sensitivity of the 4-week time-varying fit to the lower bound on the expected
# increments in the observation variance R_k. The compiled filter replaces the expected weekly
# increment (as a fraction of N) by max(increment, var_floor / N), i.e. var_floor is a bound in
# persons; the package default is var_floor = 1e-3. It is not an argument of the exported R
# functions, so the fit is reproduced in regime_helpers.R (sird_ga_fit_tv_floor) calling the
# internal epicurves:::sird_ga_loglik_inc_cpp with an explicit var_floor; no package change.
# Floors: 1e-4, 1e-3 (default), 1e-2 (requested), plus 1e-1 and 1 to show where the bound starts to
# bind. Each floor: three generic starts (paraiba_tv.R), best kept, Table 3 quantities with the
# intervals of paraiba_report.R (4000 draws). Also counts the weeks in which the bound binds at
# the optimum. Output: results/paraiba_floor_<floor>.rds (skipped if present), results/paraiba_floor.rds.
# Usage: Rscript paraiba_floor.R [cores]
source("regime_helpers.R")
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
dat <- paraiba_weekly(); n <- length(dat$cases); N <- dat$N
blk <- ceiling(seq_len(n) / 4)
floors <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

# check: with the default floor the local copy reproduces the package log-likelihood
r <- readRDS("results/paraiba_tv.rds")
f0 <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
th0 <- c(as.list(f0$theta), list(lambda = f0$lambda[1], s0 = 1 - f0$theta[["i0"]], r0 = 0))
ll_pkg <- sird_ga_loglik_inc(th0, dat$cases, dat$deaths, N, "both", 1, 20L, lambda_t = f0$lambda[blk])
ll_loc <- loglik_floor(th0, dat$cases, dat$deaths, N, 1, 20L, f0$lambda[blk], 1e-3)
cat(sprintf("check at paraiba_tv.rds optimum: package ll = %.6f, local (floor 1e-3) = %.6f\n", ll_pkg, ll_loc))
for (vf in floors) cat(sprintf("  ll at paraiba_tv.rds optimum with floor %g: %.6f\n", vf,
                               loglik_floor(th0, dat$cases, dat$deaths, N, 1, 20L, f0$lambda[blk], vf)))

for (vf in floors) {
  out_file <- sprintf("results/paraiba_floor_%g.rds", vf)
  if (file.exists(out_file)) next
  t0 <- proc.time()[3]
  fits <- parallel::mclapply(generic_starts_tv(max(blk)), function(st)
    tryCatch(sird_ga_fit_tv_floor(dat$cases, dat$deaths, N, blk, st, var_floor = vf),
             error = function(e) NULL), mc.cores = cores)
  ll <- sapply(fits, function(f) if (is.null(f)) NA else f$loglik)
  f <- fits[[which.max(ll)]]
  w <- wald_derived(f, dat$cases, dat$deaths, N, blk)
  x <- f$kf$x; inc_c <- -diff(x[, 1]) * N; inc_d <- diff(1 - rowSums(x)) * N
  saveRDS(list(var_floor = vf, logliks = ll, fit = f, wald = w,
               binds_cases = sum(inc_c < vf), binds_deaths = sum(inc_d < vf),
               min_inc_cases = min(inc_c), min_inc_deaths = min(inc_d)), out_file)
  cat(format(Sys.time()), sprintf("floor %g: logliks %s (%.0fs)\n", vf,
      paste(sprintf("%.3f", ll), collapse = " "), proc.time()[3] - t0))
}

options(width = 200)
res <- lapply(floors, function(vf) readRDS(sprintf("results/paraiba_floor_%g.rds", vf)))
summ <- do.call(rbind, lapply(res, function(z) data.frame(
  var_floor = z$var_floor, loglik = z$fit$loglik, ll_spread_starts = diff(range(z$logliks, na.rm = TRUE)),
  conv = z$fit$convergence, weeks_bound_cases = z$binds_cases, weeks_bound_deaths = z$binds_deaths,
  min_exp_weekly_infections = z$min_inc_cases, min_exp_weekly_deaths = z$min_inc_deaths)))
cat("\nFloor sensitivity: log-likelihood and binding\n"); print(format(summ, digits = 6), row.names = FALSE)
for (z in res) {
  cat(sprintf("\nvar_floor = %g: Table 3 quantities (estimate, 2.5%%, 97.5%% of 4000 draws)\n", z$var_floor))
  print(signif(data.frame(estimate = z$wald$est, lower = z$wald$lower, upper = z$wald$upper), 4))
}
saveRDS(list(summary = summ, details = res), "results/paraiba_floor.rds")
cat("FLOORDONE\n")
