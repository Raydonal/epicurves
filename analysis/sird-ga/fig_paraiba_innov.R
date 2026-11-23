# Figure and log: standardised one-step-ahead (Kalman) innovations of the best time-varying Paraiba fit.
# Standardised innovation = (observed - fitted one-step-ahead mean) / predictive SD, per series (kf$fitted_*, kf$sd_*).
r <- readRDS("results/paraiba_tv.rds"); f <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
kf <- f$kf
ser <- list(Cases = (r$cases - kf$fitted_cases) / kf$sd_cases, Deaths = (r$deaths - kf$fitted_deaths) / kf$sd_deaths)
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
blue <- "#0072B2"; verm <- "#D55E00"
pdf(file.path(figdir, "paraiba_innovations.pdf"), width = 9, height = 5.6)
par(mfrow = c(2, 3), mar = c(4, 4.4, 1.6, 0.6), las = 1, cex.axis = 0.85, mgp = c(2.4, 0.7, 0))
for (nm in names(ser)) {
  z <- ser[[nm]]
  plot(r$weeks, z, type = "h", col = blue, lwd = 1.5, xlab = "Epidemiological week", ylab = "Standardised innovation", bty = "l",
       ylim = range(c(z, -2, 2)))
  points(r$weeks, z, pch = 19, cex = 0.5, col = blue); abline(h = c(-2, 0, 2), lty = c(3, 1, 3), col = "grey50")
  mtext(paste0(nm, ": innovations"), 3, 0.2, adj = 0, cex = 0.85)
  a <- acf(z, lag.max = 12, plot = FALSE)
  plot(a$lag[-1], a$acf[-1], type = "h", lwd = 2, col = blue, ylim = c(-0.6, 0.6), xlab = "Lag (weeks)", ylab = "Autocorrelation", bty = "l")
  abline(h = 0); abline(h = c(-1, 1) * 1.96 / sqrt(length(z)), lty = 2, col = verm)
  mtext(paste0(nm, ": ACF"), 3, 0.2, adj = 0, cex = 0.85)
  qq <- qqnorm(z, plot.it = FALSE)
  plot(qq$x, qq$y, pch = 19, cex = 0.6, col = blue, xlab = "Normal quantile", ylab = "Sample quantile", bty = "l")
  qqline(z, col = "grey30", lty = 2); mtext(paste0(nm, ": normal QQ"), 3, 0.2, adj = 0, cex = 0.85)
}
dev.off()
sink("results/paraiba_innovations.log")
cat("Standardised one-step-ahead innovations, best fit of results/paraiba_tv.rds (loglik", round(f$loglik, 2), "), n =", length(r$cases), "weeks\n")
cat("z = (y - fitted_one_step_mean) / predictive_sd, per series (cases and deaths marginally)\n")
for (nm in names(ser)) {
  z <- ser[[nm]]
  cat(sprintf("\n%s: mean %.3f, sd %.3f, min %.2f, max %.2f\n", nm, mean(z), sd(z), min(z), max(z)))
  for (L in c(4, 8)) { b <- Box.test(z, lag = L, type = "Ljung-Box"); cat(sprintf("  Ljung-Box lag %d: Q = %.3f, p = %.4f\n", L, b$statistic, b$p.value)) }
  s <- shapiro.test(z); cat(sprintf("  Shapiro-Wilk: W = %.4f, p = %.4f\n", s$statistic, s$p.value))
  cat("  lag-1 autocorrelation:", round(acf(z, plot = FALSE)$acf[2], 3), "\n")
}
sink(); cat(readLines("results/paraiba_innovations.log"), sep = "\n")
