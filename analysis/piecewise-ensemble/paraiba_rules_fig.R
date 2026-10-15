# Figure and segment table for the Paraíba application (results/paraiba_rules_*.rds).
suppressMessages(library(epicurves))
figdir <- if (dir.exists("../SMMR_submission/figures")) "../SMMR_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
pdf(file.path(figdir, "paraiba_segments.pdf"), width = 10, height = 3.8)
op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 1.5, 0.5), las = 1)
for (s in c("cases", "deaths")) {
  r <- readRDS(sprintf("results/paraiba_rules_%s.rds", s)); x <- r$x; dt <- r$dates
  cuts <- r$out$inc_bic_negbin$v$cuts
  bounds <- cbind(c(1, cuts + 1), c(cuts, length(x)))
  plot(dt, x, pch = 19, cex = 0.35, col = "grey45", xlab = "", ylab = paste("Daily", s),
       main = if (s == "cases") "Cases" else "Deaths")
  cat("\n==", s, "\n")
  for (j in seq_len(nrow(bounds))) {
    idx <- bounds[j, 1]:bounds[j, 2]
    e <- fit_ensemble_incidence(x[idx], family = "negbin")
    lines(dt[idx], e$fitted, lwd = 2)
    cat(sprintf("segment %d: %s to %s (%d days)  w(G,R,L)=%s  size=%.1f  mean count=%.1f\n", j,
                format(dt[idx[1]]), format(dt[max(idx)]), length(idx),
                paste(round(e$weights[c("gompertz", "richards", "logistic")], 3), collapse = "/"),
                e$size, mean(x[idx])))
  }
  abline(v = dt[cuts], lty = 2)
  al <- r$out$edetect_subgaussian_z1$v$alarms
  rug(dt[al], side = 3, lwd = 2, col = "black", ticksize = 0.05)
  cat("sub-Gaussian z=1 alarms:", format(dt[al]), " delays:", r$out$edetect_subgaussian_z1$v$delay, "\n")
}
par(op); dev.off()
