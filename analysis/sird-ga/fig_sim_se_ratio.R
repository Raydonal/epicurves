# Figure: median Wald SE / Monte Carlo SD (optimisation scale) and Wald 95% coverage, per parameter, scenario, N.
# Same definitions as sim_se_diagnostic.R and summarise_sim.R; coverage uses replicates with finite SE.
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
tf <- list(lambda = log, gamma = log, alpha = log, p = function(v) qlogis(v / 0.999))
pn <- c("lambda", "gamma", "alpha", "p"); Ns <- c(5000, 20000, 100000); scs <- 1:5
cols <- c("#C6DBEF", "#6BAED6", "#08519C"); pchs <- c(21, 22, 24)
res <- list()
for (sc in scs) for (j in 1:3) {
  x <- readRDS(sprintf("results/sim_s%d_N%d.rds", sc, Ns[j])); tr <- attr(x, "truth")
  for (pp in pn) {
    u <- x[[paste0("u_", pp)]]; se <- x[[paste0("se_", pp)]]; ok <- is.finite(se)
    res[[length(res) + 1]] <- data.frame(sc = sc, j = j, par = pp, ratio = median(se[ok]) / sd(u[ok]),
      cover = mean(abs(u[ok] - tf[[pp]](tr[[pp]])) <= 1.96 * se[ok]))
  }
}
res <- do.call(rbind, res)
lab <- list(lambda = expression(lambda), gamma = expression(gamma), alpha = expression(alpha), p = expression(p))
pdf(file.path(figdir, "sim_se_ratio.pdf"), width = 7.5, height = 3.6)
par(mfrow = c(1, 2), mar = c(4, 4.4, 1.6, 0.6), las = 1, cex.axis = 0.85, mgp = c(2.4, 0.7, 0))
panel <- function(col, ylab, ref, yl, main) {
  plot(NA, xlim = c(0.5, 4.5), ylim = yl, xaxt = "n", xlab = "", ylab = ylab, bty = "l")
  abline(h = ref, col = "grey30", lty = 2)
  abline(v = 1:3 + 0.5, col = "grey85", lwd = 0.8)
  for (i in seq_along(pn)) for (sc in scs) for (j in 1:3) {
    r <- res[res$par == pn[i] & res$sc == sc & res$j == j, ]
    xx <- i + (sc - 3) * 0.15
    points(xx + (j - 2) * 0.035, r[[col]], pch = pchs[j], bg = cols[j], col = "grey20", cex = 0.7, lwd = 0.6)
  }
  axis(1, at = 1:4, labels = FALSE)
  for (i in 1:4) mtext(lab[[i]], 1, 1.1, at = i, cex = 0.9)
  mtext("Parameter (scenarios 1-5 left to right within each)", 1, 2.5, cex = 0.7)
  mtext(main, 3, 0.2, adj = 0, cex = 0.9)
}
panel("ratio", "Median Wald SE / Monte Carlo SD", 1, c(0, 1.1), "(a) Standard error ratio")
panel("cover", "Coverage of Wald 95% interval", 0.95, c(0, 1), "(b) Coverage")
legend("bottom", paste("N =", formatC(Ns, format = "d", big.mark = ",")), pch = pchs, pt.bg = cols, col = "grey20", bty = "n", cex = 0.75, horiz = TRUE)
dev.off()
print(aggregate(cbind(ratio, cover) ~ j + par, res, range))
