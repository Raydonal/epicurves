# Paraíba: the segmentation rules of the Monte Carlo study applied to daily cases and deaths.
# Usage: Rscript paraiba_rules.R <cases|deaths>
suppressMessages(library(epicurves))
series <- commandArgs(trailingOnly = TRUE)[1]
data(paraiba_covid19)
d <- paraiba_covid19
x <- if (series == "cases") d$new_cases else d$new_deaths
y <- cumsum(x)
cat(series, ": n =", length(x), "days,", format(min(d$date)), "to", format(max(d$date)), "\n")

timed <- function(expr) { t0 <- proc.time()[3]; v <- tryCatch(expr, error = function(e) { message(conditionMessage(e)); NULL }); list(v = v, s = proc.time()[3] - t0) }
out <- list()
out$thesis_flat <- timed(fit_piecewise_ensemble(y, window = 90, update_every = 20, min_obs = 12,
                                                criterion = "mse", threshold = 0.5,
                                                confirm_criterion = "mse_threshold"))
out$cum_bic <- timed(fit_piecewise_ensemble(y, window = 90, update_every = 20, min_obs = 12,
                                            criterion = "mse", confirm_criterion = "aic"))
for (fam in c("poisson", "negbin"))
  out[[paste0("inc_bic_", fam)]] <- timed(bic_segment(x, d_min = 12, family = fam, step = 5))
for (ch in c("poisson", "subgaussian")) for (z in c(0, 1))
  out[[paste0("edetect_", ch, "_z", z)]] <- timed(edetect_segment(x, d_min = 20, arl = 1000,
                                                                  refit_every = 7, chart = ch, z = z))
for (nm in names(out)) {
  v <- out[[nm]]$v
  cat(sprintf("%-22s cuts: %-40s (%.0fs)\n", nm,
              if (is.null(v)) "FAILED" else paste(v$cuts, collapse = " "), out[[nm]]$s))
  if (!is.null(v$cuts) && length(v$cuts)) cat(sprintf("%-22s dates: %s\n", "", paste(format(d$date[v$cuts]), collapse = " ")))
}
# per-segment negative binomial ensembles for the offline penalised-likelihood segmentation
cuts <- out$inc_bic_negbin$v$cuts
bounds <- cbind(c(1, cuts + 1), c(cuts, length(x)))
seg <- lapply(seq_len(nrow(bounds)), function(j) {
  e <- fit_ensemble_incidence(x[bounds[j, 1]:bounds[j, 2]], family = "negbin")
  data.frame(start = format(d$date[bounds[j, 1]]), end = format(d$date[bounds[j, 2]]),
             days = bounds[j, 2] - bounds[j, 1] + 1, t(round(e$weights, 3)), loglik = e$loglik)
})
print(do.call(rbind, seg))
saveRDS(list(out = out, x = x, dates = d$date, series = series),
        sprintf("results/paraiba_rules_%s.rds", series))
