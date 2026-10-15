# Monte Carlo comparison of segmentation rules for the piecewise ensemble article.
# Usage: Rscript mc_study.R <scenario: 1|2|3|1c|2c|3c|nc_log|nc_gom> <J> <cores>
# 1-3: thesis scenarios (one-day spike at each change point); 1c-3c: continuous versions.
suppressMessages(library(epicurves))
source("scenarios.R")
args <- commandArgs(trailingOnly = TRUE)
sc <- args[1]; J <- as.integer(args[2]); cores <- as.integer(args[3])

truth <- if (sc %in% c("1", "2", "3")) scenario_curve(as.integer(sc)) else
  if (sc %in% c("1c", "2c", "3c")) scenario_curve(as.integer(substr(sc, 1, 1)), continuous = TRUE) else {
  C <- if (sc == "nc_log") logistic_curve(1:200, 10, 3000, 0.05) else gompertz_curve(1:200, 10, 3000, 0.03)
  list(C = C, mu = pmax(diff(c(10, C)), 0.01), cuts = integer(), y0 = 10)
}

timed <- function(expr) { t0 <- proc.time()[3]; v <- tryCatch(expr, error = function(e) NULL); list(v = v, s = proc.time()[3] - t0) }

one <- function(j) {
  set.seed(5000L + 100L * j + match(sc, c("1", "2", "3", "nc_log", "nc_gom", "1c", "2c", "3c")))
  x <- stats::rpois(200, truth$mu)
  y <- truth$y0 + cumsum(x)
  out <- list()
  r <- timed(fit_piecewise_ensemble(y, window = 90, update_every = 20, min_obs = 12,
                                    criterion = "mse", threshold = 0.5,
                                    confirm_criterion = "mse_threshold"))
  out$thesis_flat <- list(cuts = r$v$cuts, secs = r$s)
  r <- timed(fit_piecewise_ensemble(y, window = 90, update_every = 20, min_obs = 12,
                                    criterion = "mse", confirm_criterion = "aic"))
  out$cum_bic <- list(cuts = r$v$cuts, secs = r$s)
  r <- timed(bic_segment(x, d_min = 12, family = "poisson", step = 5))
  out$inc_bic <- list(cuts = r$v$cuts, secs = r$s)
  for (z in c(0, 1)) {
    r <- timed(edetect_segment(x, d_min = 20, arl = 1000, refit_every = 7, z = z))
    out[[paste0("edetect_z", z)]] <- list(cuts = r$v$cuts, alarms = r$v$alarms, secs = r$s)
  }
  out$x <- x
  out
}
res <- parallel::mclapply(seq_len(J), one, mc.cores = cores, mc.preschedule = FALSE)
saveRDS(list(res = res, truth = truth, scenario = sc),
        sprintf("results/mc_%s.rds", sc))
cat("done", sc, "\n")
