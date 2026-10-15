# Summaries of the Monte Carlo study (results/mc_*.rds) for the piecewise ensemble article.
# A true change point counts as detected if some estimated cut lies within `tol` days of it;
# an estimated cut is false if no true change point lies within `tol` days.
# Usage: Rscript summarise_mc.R [tol]
args <- commandArgs(trailingOnly = TRUE)
tol <- if (length(args)) as.numeric(args[1]) else 7
rules <- c("thesis_flat", "cum_bic", "inc_bic", "edetect_z0", "edetect_z1")
scs <- c("1c", "2c", "3c", "nc_log", "nc_gom", "1", "2", "3")

summ_rule <- function(res, truth, rule) {
  tc <- truth$cuts
  per <- lapply(res, function(r) {
    o <- r[[rule]]; failed <- is.null(o$cuts)  # the rule threw an error
    ec <- if (failed) numeric() else o$cuts
    det <- vapply(tc, function(c) any(abs(ec - c) <= tol), logical(1))
    err <- vapply(tc, function(c) if (length(ec)) min(abs(ec - c)) else NA_real_, numeric(1))
    false <- sum(vapply(ec, function(e) !any(abs(tc - e) <= tol), logical(1)))
    delay <- if (!is.null(o$alarms) && length(tc)) vapply(tc, function(c) {
      k <- which(abs(ec - c) <= tol); if (length(k)) o$alarms[k[1]] - c else NA_real_ }, numeric(1)) else NA_real_
    list(failed = failed, n = length(ec), det = det, err = err[det], false = false, delay = delay, secs = unname(o$secs))
  })
  n <- vapply(per, `[[`, numeric(1), "n")
  data.frame(rule = rule, failed = sum(vapply(per, `[[`, logical(1), "failed")), true_cuts = length(tc), mean_cuts = mean(n),
             exact = mean(n == length(tc)),
             detect = if (length(tc)) mean(unlist(lapply(per, `[[`, "det"))) else NA,
             mean_abs_err = if (length(tc)) mean(unlist(lapply(per, `[[`, "err")), na.rm = TRUE) else NA,
             false_per_rep = mean(vapply(per, `[[`, numeric(1), "false")),
             median_delay = suppressWarnings(median(unlist(lapply(per, `[[`, "delay")), na.rm = TRUE)),
             median_secs = median(vapply(per, `[[`, numeric(1), "secs")))
}

all <- list()
for (sc in scs) {
  f <- sprintf("results/mc_%s.rds", sc)
  if (!file.exists(f)) next
  m <- readRDS(f)
  ok <- !vapply(m$res, inherits, logical(1), "try-error")
  s <- do.call(rbind, lapply(rules, summ_rule, res = m$res[ok], truth = m$truth))
  s <- cbind(scenario = sc, J = sum(ok), s)
  all[[sc]] <- s
  cat(sprintf("\n== scenario %s  J = %d  true cuts: %s  (tol = %g days)\n", sc, sum(ok),
              paste(m$truth$cuts, collapse = " "), tol))
  print(s[, -(1:2)], digits = 3, row.names = FALSE)
  # per-change-point detection rates
  if (length(m$truth$cuts)) {
    pc <- sapply(rules, function(rl) vapply(m$truth$cuts, function(c)
      mean(vapply(m$res[ok], function(r) any(abs(r[[rl]]$cuts - c) <= tol), logical(1))), numeric(1)))
    pc <- matrix(pc, nrow = length(m$truth$cuts), dimnames = list(paste0("cut", m$truth$cuts), rules))
    # size of each change: relative jump of the daily mean from the last day of a phase to the
    # first day of the next (the thesis scenarios add a one-day spike there; we report day c+2)
    mu <- m$truth$mu; cc <- m$truth$cuts
    jump <- if (sc %in% c("1", "2", "3")) (mu[cc + 2] - mu[cc]) / mu[cc] else (mu[cc + 1] - mu[cc]) / mu[cc]
    print(round(cbind(rel_jump = jump, pc), 2))
  }
}
saveRDS(do.call(rbind, all), sprintf("results/mc_summary_tol%g.rds", tol))
