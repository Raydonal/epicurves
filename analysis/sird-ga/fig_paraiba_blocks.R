# Figure: sensitivity of the time-varying Paraiba fit to the block length L (weeks), from
# results/paraiba_blocks.rds (generic starts) and results/paraiba_blocks_warm.rds (best of generic and
# warm starts); best log-likelihood per L across both. AIC = -2 ll + 2k, BIC = -2 ll + k log 80.
a <- readRDS("results/paraiba_blocks.rds"); w <- readRDS("results/paraiba_blocks_warm.rds")$table
cols <- c("L", "k", "loglik", "p", "inf_period_days", "attack")
d <- do.call(rbind, lapply(sort(unique(a$L)), function(L) {
  x <- rbind(a[a$L == L, cols], w[w$L == L, cols]); x[which.max(x$loglik), ]
}))
d$AIC <- -2 * d$loglik + 2 * d$k; d$BIC <- -2 * d$loglik + d$k * log(80)
d$dAIC <- d$AIC - min(d$AIC); d$dBIC <- d$BIC - min(d$BIC)
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
pdf(file.path(figdir, "paraiba_blocks.pdf"), width = 9, height = 3.0)
par(mfrow = c(1, 4), mar = c(4, 4.4, 1.6, 0.6), las = 1, cex.axis = 0.85, mgp = c(2.4, 0.7, 0))
blue <- "#0072B2"; verm <- "#D55E00"; shade <- adjustcolor(verm, 0.13)
pan <- function(ylab, ylim, tag) {
  plot(d$L, d$L, type = "n", ylim = ylim, xlab = "Block length L (weeks)", ylab = ylab, bty = "l", xaxt = "n")
  axis(1, d$L); u <- par("usr"); rect(3.5, u[3], 4.5, u[4], col = shade, border = NA); mtext(tag, 3, 0.2, adj = 0)
}
pan(expression(Delta*"AIC, "*Delta*"BIC"), c(0, max(d$dAIC, d$dBIC) * 1.05), "(a)")
lines(d$L, d$dAIC, type = "b", pch = 19, col = blue); lines(d$L, d$dBIC, type = "b", pch = 17, lty = 2, col = "grey20")
legend("topright", c("AIC", "BIC"), pch = c(19, 17), lty = c(1, 2), col = c(blue, "grey20"), bty = "n", cex = 0.8)
pan(expression(hat(p)), c(0, max(d$p) * 1.1), "(b)"); lines(d$L, d$p, type = "b", pch = 19, col = blue)
pan("Attack rate by week 52", c(0, 1), "(c)"); lines(d$L, d$attack, type = "b", pch = 19, col = blue)
pan("Infectious period (days)", c(0, max(d$inf_period_days) * 1.1), "(d)"); lines(d$L, d$inf_period_days, type = "b", pch = 19, col = blue)
dev.off()
options(width = 160); print(round(d, 4), row.names = FALSE)
