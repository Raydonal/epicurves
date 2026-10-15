# Monte Carlo under overdispersion: series without change point (Logistic or Gompertz mean of the
# no-change scenarios of mc_study.R) with negative binomial counts, variance mu (1 + mu / phi).
# Rules: penalised likelihood (Poisson and negative binomial) and e-detectors (Poisson z = 1,
# sub-Gaussian z = 0 and z = 1), with the settings of mc_study.R and paraiba_rules.R.
# Usage: Rscript mc_negbin.R <log|gom> <phi> <J> <cores>
suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
curve <- args[1]; phi <- as.numeric(args[2]); J <- as.integer(args[3]); cores <- as.integer(args[4])
C <- if (curve == "log") logistic_curve(1:200, 10, 3000, 0.05) else gompertz_curve(1:200, 10, 3000, 0.03)
mu <- pmax(diff(c(10, C)), 0.01)

timed <- function(expr) { t0 <- proc.time()[3]; v <- tryCatch(expr, error = function(e) NULL); list(v = v, s = proc.time()[3] - t0) }
one <- function(j) {
  set.seed(9000L + 100L * j + 10L * (curve == "gom") + round(phi))
  x <- stats::rnbinom(200, mu = mu, size = phi)
  out <- list(x = x)
  for (fam in c("poisson", "negbin")) {
    r <- timed(bic_segment(x, d_min = 12, family = fam, step = 5))
    out[[paste0("pl_", fam)]] <- list(cuts = r$v$cuts, ok = !is.null(r$v), secs = r$s)
  }
  for (cfg in list(c("poisson", 1), c("subgaussian", 0), c("subgaussian", 1))) {
    r <- timed(edetect_segment(x, d_min = 20, arl = 1000, refit_every = 7, chart = cfg[1], z = as.numeric(cfg[2])))
    out[[sprintf("ed_%s_z%s", cfg[1], cfg[2])]] <- list(cuts = r$v$cuts, alarms = r$v$alarms, ok = !is.null(r$v), secs = r$s)
  }
  out
}
res <- parallel::mclapply(seq_len(J), one, mc.cores = cores, mc.preschedule = FALSE)
f <- sprintf("results/mc_negbin_%s_phi%g.rds", curve, phi)
saveRDS(list(res = res, mu = mu, curve = curve, phi = phi), f)
rules <- setdiff(names(res[[1]]), "x")
cat(sprintf("== %s, phi = %g, J = %d\n", curve, phi, J))
for (r in rules) {
  ok <- vapply(res, function(z) isTRUE(z[[r]]$ok), logical(1))
  k <- vapply(res, function(z) length(z[[r]]$cuts), numeric(1))
  cat(sprintf("%-18s failed %2d  mean cuts %.2f  share with any cut %.2f  median secs %.1f\n", r, sum(!ok),
              mean(k[ok]), mean(k[ok] > 0), median(vapply(res, function(z) z[[r]]$secs, numeric(1)))))
}
