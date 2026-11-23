# Figure: refined profile log-likelihood of p (Paraiba, time-varying fit, 4-week blocks) from
# results/paraiba_profile_p_refined.rds: best of the starts (generic, forward and backward warm sweeps)
# at each fixed p; infectious period and attack rate at that best fit. Shaded: profile interval.
r <- readRDS("results/paraiba_profile_p_refined.rds")
f <- readRDS("results/paraiba_tv.rds")$res; best <- f[[which.max(sapply(f, `[[`, "loglik"))]]
pm <- unname(best$theta["p"]); lm <- best$loglik
pr <- r$profile; pr <- pr[order(pr$p), ]
llmax <- max(c(pr$loglik, lm)); rel <- pr$loglik - llmax
# interval limits recomputed by linear interpolation of the profile at llmax - 1.92
cut <- llmax - 1.92; inside <- pr$loglik >= cut; i1 <- min(which(inside)); i2 <- max(which(inside))
ip <- function(a, b) pr$p[a] + (cut - pr$loglik[a]) * (pr$p[b] - pr$p[a]) / (pr$loglik[b] - pr$loglik[a])
lo <- ip(i1 - 1, i1); hi <- ip(i2 + 1, i2)
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
pdf(file.path(figdir, "paraiba_profile_p.pdf"), width = 9, height = 3.4)
par(mfrow = c(1, 3), mar = c(4, 4.4, 1.6, 0.6), las = 1, cex.axis = 0.85, mgp = c(2.4, 0.7, 0))
blue <- "#0072B2"; verm <- "#D55E00"; shade <- adjustcolor(verm, 0.13)
band <- function() { u <- par("usr"); rect(lo, 10^u[3], hi, 10^u[4], col = shade, border = NA) }
frame3 <- function(y, ylab, ylim, tag, ...) {
  plot(pr$p, y, log = "x", type = "n", ylim = ylim, xlab = expression(p~"(fixed, log scale)"), ylab = ylab, bty = "l")
  u <- par("usr"); rect(lo, u[3], hi, u[4], col = shade, border = NA)
  lines(pr$p, y, type = "b", pch = 19, col = blue); abline(v = pm, col = verm, lty = 2)
  mtext(tag, 3, 0.2, adj = 0)
}
frame3(rel, "Profile log-likelihood minus maximum", c(-45, 2), "(a)")
abline(h = -1.92, col = "grey40", lty = 3)
text(0.0205, -1.92, "-1.92", pos = 3, cex = 0.75, col = "grey30")
text(pm * 0.93, -44, sprintf("MLE %.3f", pm), adj = 1, cex = 0.8, col = verm)
text(sqrt(lo * hi), -30, sprintf("[%.3f, %.3f]", lo, hi), cex = 0.75, col = verm, srt = 90)
frame3(pr$inf_period_days, "Infectious period (days)", c(0, max(pr$inf_period_days) * 1.05), "(b)")
frame3(pr$attack, "Attack rate by week 52", c(0, 1), "(c)")
dev.off()
print(data.frame(p = pr$p, rel = round(rel, 3), inf = round(pr$inf_period_days, 2), attack = round(pr$attack, 3)))
cat("MLE p =", pm, " loglik =", lm, " llmax =", llmax, " interval = [", lo, hi, "]  stored:", r$ci_refined[1:2], "\n")
