# Summaries of the SIRD-GA simulation: mean, empirical 90% interval, bias, RMSE,
# coverage of Wald 95% intervals (on the unconstrained scale) for each parameter.
summ <- function(file) {
  x <- readRDS(file); tr <- attr(x, "truth")
  pars <- c("lambda", "gamma", "alpha", "p", "tau2", "tau2_d")
  tf <- list(lambda = log, gamma = log, alpha = log, tau2 = log, tau2_d = log,
             p = function(v) qlogis(v / 0.999))
  out <- do.call(rbind, lapply(pars, function(pp) {
    est <- x[[pp]]; u <- x[[paste0("u_", pp)]]; se <- x[[paste0("se_", pp)]]
    ok <- is.finite(se)
    cover <- mean(abs(u[ok] - tf[[pp]](tr[[pp]])) <= 1.96 * se[ok])
    data.frame(param = pp, true = tr[[pp]], mean = mean(est), median = median(est),
               q05 = quantile(est, 0.05), q95 = quantile(est, 0.95),
               bias = mean(est) - tr[[pp]], rmse = sqrt(mean((est - tr[[pp]])^2)),
               cover95 = cover, n_se = sum(ok))
  }))
  R0 <- x$lambda / (x$gamma + x$alpha); R0t <- tr[["lambda"]] / (tr[["gamma"]] + tr[["alpha"]])
  ifr <- x$alpha / (x$gamma + x$alpha); ifrt <- tr[["alpha"]] / (tr[["gamma"]] + tr[["alpha"]])
  out <- rbind(out,
    data.frame(param = "R0", true = R0t, mean = mean(R0), median = median(R0), q05 = quantile(R0, .05),
               q95 = quantile(R0, .95), bias = mean(R0) - R0t, rmse = sqrt(mean((R0 - R0t)^2)), cover95 = NA, n_se = NA),
    data.frame(param = "ifr", true = ifrt, mean = mean(ifr), median = median(ifr), q05 = quantile(ifr, .05),
               q95 = quantile(ifr, .95), bias = mean(ifr) - ifrt, rmse = sqrt(mean((ifr - ifrt)^2)), cover95 = NA, n_se = NA))
  rownames(out) <- NULL
  attr(out, "J") <- nrow(x); attr(out, "tries") <- mean(x$tries); attr(out, "secs") <- median(x$secs)
  out
}
if (sys.nframe() == 0) {
  for (f in sort(list.files("results", "^sim_s.*rds$", full.names = TRUE))) {
    s <- summ(f); cat("\n==", f, " J =", attr(s, "J"), " mean tries =", round(attr(s, "tries"), 2),
                       " median secs =", round(attr(s, "secs")), "\n")
    print(s, digits = 3)
  }
}
