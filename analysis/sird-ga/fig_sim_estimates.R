# Figure: sampling distribution of SIRD-GA estimates relative to the truth (estimate / true value, log scale),
# by scenario and population size, from results/sim_s<sc>_N<N>.rds. All replicates are used, as in summarise_sim.R.
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
Ns <- c(5000, 20000, 100000); scs <- 1:5
cols <- c("#C6DBEF", "#6BAED6", "#08519C")            # one hue, light to dark = increasing N (greyscale-safe)
pars <- c(lambda = expression(lambda), gamma = expression(gamma), alpha = expression(alpha),
          p = expression(p), R0 = expression(R[0]), ifr = "IFR")
get <- function(sc, N) {
  x <- readRDS(sprintf("results/sim_s%d_N%d.rds", sc, N)); t <- attr(x, "truth")
  d <- list(lambda = x$lambda / t[["lambda"]], gamma = x$gamma / t[["gamma"]], alpha = x$alpha / t[["alpha"]],
            p = x$p / t[["p"]],
            R0 = (x$lambda / (x$gamma + x$alpha)) / (t[["lambda"]] / (t[["gamma"]] + t[["alpha"]])),
            ifr = (x$alpha / (x$gamma + x$alpha)) / (t[["alpha"]] / (t[["gamma"]] + t[["alpha"]])))
  d
}
dat <- lapply(scs, function(sc) lapply(Ns, function(N) get(sc, N)))
pdf(file.path(figdir, "sim_estimates.pdf"), width = 7.5, height = 8)
par(mfrow = c(3, 2), mar = c(3.2, 4.4, 1.6, 0.6), oma = c(0, 0, 0, 0), las = 1, cex.axis = 0.85, mgp = c(2.4, 0.7, 0))
ylab <- "Estimate / true value"
for (k in names(pars)) {
  allv <- unlist(lapply(dat, function(a) lapply(a, function(b) b[[k]])))
  yl <- range(allv); yl <- c(max(yl[1], 0.05), min(yl[2], 20))
  plot(NA, xlim = c(0.4, 5.6), ylim = yl, log = "y", xaxt = "n", xlab = "", ylab = ylab, bty = "l")
  abline(h = 1, col = "grey30", lty = 2)
  for (sc in scs) for (j in 1:3) {
    v <- dat[[sc]][[j]][[k]]
    boxplot(v, at = sc + (j - 2) * 0.27, add = TRUE, boxwex = 0.24, col = cols[j], border = "grey20",
            axes = FALSE, outline = TRUE, pch = 20, cex = 0.35, lwd = 0.7, staplewex = 0, log = "y")
  }
  axis(1, at = scs, labels = scs, tick = TRUE); mtext("Scenario", 1, 2.1, cex = 0.7)
  mtext(pars[[k]], 3, 0.2, adj = 0, cex = 1.0)
  if (k == "lambda") legend("topright", paste("N =", formatC(Ns, format = "d", big.mark = ",")), fill = cols, border = "grey20",
                           bty = "n", cex = 0.75)
}
dev.off()
for (sc in scs) for (j in 1:3) cat(sc, Ns[j], "median est/true:", paste(names(pars), round(sapply(dat[[sc]][[j]], median), 2), collapse = " "), "\n")
