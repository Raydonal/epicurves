# Tables and figure for the Paraíba application (time-varying transmission, weekly data).
# Uncertainty of derived quantities: 4000 draws from the asymptotic normal distribution of the
# estimates on the unconstrained scale (inverse numerical Hessian), mapped through the model.
suppressMessages(library(epicurves))
r <- readRDS("results/paraiba_tv.rds")
f <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
B <- max(r$blk); p_max <- 0.999
V <- solve(f$hessian)
cat("Hessian eigenvalues (min, max):", signif(range(eigen(f$hessian, only.values = TRUE)$values), 3), "\n")
unpack <- function(v) {
  th <- list(q = 1, r0 = 0, lambda = exp(v[1]))
  bw <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = plogis,
             p = function(z) p_max * plogis(z))
  for (j in seq_along(f$sc_names)) th[[f$sc_names[j]]] <- bw[[f$sc_names[j]]](unname(v[B + j]))
  th$s0 <- 1 - th$i0 - th$r0
  list(th = th, lam = exp(v[seq_len(B)]))
}
derived <- function(v) {
  u <- unpack(v); th <- u$th
  kf <- sird_ga_loglik_inc(th, r$cases, r$deaths, r$N, "both", 1, 20L, full = TRUE,
                           lambda_t = u$lam[r$blk])
  s <- kf$x[-1, 1]
  Rt <- tapply(u$lam[r$blk] * s / (th$gamma + th$alpha), r$blk, mean)
  c(inf_period_days = 7 / (th$gamma + th$alpha), ifr = th$alpha / (th$gamma + th$alpha),
    p = th$p, i0N = th$i0 * r$N, tau2 = th$tau2, tau2_d = th$tau2_d,
    attack = 1 - s[length(s)], setNames(Rt, paste0("R_block", seq_len(B))))
}
est <- derived(f$par)
set.seed(1)
L <- chol(V)
draws <- t(replicate(4000, derived(f$par + drop(rnorm(length(f$par)) %*% L))))
tab <- data.frame(estimate = est, lower = apply(draws, 2, quantile, 0.025),
                  upper = apply(draws, 2, quantile, 0.975))
print(signif(tab, 3))
blk_weeks <- tapply(r$weeks, r$blk, function(w) paste(range(w), collapse = "-"))
cat("block epidemiological weeks:", blk_weeks, "\n")
cat("weekly lambda:", round(f$lambda, 3), "\n")
cat("loglik:", f$loglik, " convergence:", f$convergence, "\n")

# constant-transmission comparison (weekly, s0 fixed at 1 - i0)
w <- readRDS("results/paraiba_weekly_s0fixed_q1.00.rds")$fits$both
cat(sprintf("constant lambda: loglik=%.2f lambda=%.3f gamma=%.3f alpha=%.4f R0=%.4f p=%.3f tau2=%.0f\n",
            w$loglik, w$theta["lambda"], w$theta["gamma"], w$theta["alpha"],
            w$theta["lambda"] / (w$theta["gamma"] + w$theta["alpha"]), w$theta["p"], w$theta["tau2"]))
cat("parameters: time-varying", length(f$par), " constant", length(w$par), "\n")
saveRDS(list(table = tab, draws = draws, blk_weeks = blk_weeks), "results/paraiba_report.rds")

# Figure: one-step-ahead predictions and R_t by block
figdir <- if (dir.exists("../CSDA_submission/figures")) "../CSDA_submission/figures" else "figures"
dir.create(figdir, showWarnings = FALSE)
pdf(file.path(figdir, "paraiba_tv_fit.pdf"), width = 10, height = 3.4)
op <- par(mfrow = c(1, 3), mar = c(4, 4.2, 1.5, 0.5), las = 1, cex.axis = 0.85)
kf <- f$kf
plot(r$weeks, r$cases / 1000, type = "n", xlab = "Epidemiological week (2020)",
     ylab = "Reported cases (thousands)", main = "Cases", ylim = c(0, max(r$cases) * 1.15 / 1000))
polygon(c(r$weeks, rev(r$weeks)), c(kf$fitted_cases - 2 * kf$sd_cases, rev(kf$fitted_cases + 2 * kf$sd_cases)) / 1000,
        col = "grey85", border = NA)
lines(r$weeks, kf$fitted_cases / 1000, lwd = 2); points(r$weeks, r$cases / 1000, pch = 19, cex = 0.6)
plot(r$weeks, r$deaths, type = "n", xlab = "Epidemiological week (2020)", ylab = "Reported deaths",
     main = "Deaths", ylim = c(0, max(r$deaths) * 1.15))
polygon(c(r$weeks, rev(r$weeks)), c(kf$fitted_deaths - 2 * kf$sd_deaths, rev(kf$fitted_deaths + 2 * kf$sd_deaths)),
        col = "grey85", border = NA)
lines(r$weeks, kf$fitted_deaths, lwd = 2); points(r$weeks, r$deaths, pch = 19, cex = 0.6)
Rcols <- grep("^R_block", colnames(draws))
mid <- tapply(r$weeks, r$blk, mean)
plot(mid, est[Rcols], ylim = range(tab[Rcols, c("lower", "upper")], 1), pch = 19,
     xlab = "Epidemiological week (2020)", ylab = expression(R[t]~"(block mean)"), main = "Reproduction number")
arrows(mid, tab$lower[Rcols], mid, tab$upper[Rcols], angle = 90, code = 3, length = 0.03)
abline(h = 1, lty = 2)
par(op); dev.off()
