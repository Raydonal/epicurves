# Paraíba: SIRD-GA with piecewise-constant transmission (4-week blocks), weekly data.
suppressMessages(library(epicurves))
data(paraiba_covid19); data(paraiba_population)
N <- paraiba_population; d <- paraiba_covid19
wk <- d$epi_week; keep <- wk %in% names(which(table(wk) == 7))
cases <- as.numeric(tapply(d$new_cases[keep], wk[keep], sum))
deaths <- as.numeric(tapply(d$new_deaths[keep], wk[keep], sum))
weeks <- as.numeric(names(table(wk[keep])))
n <- length(cases); blk <- ceiling(seq_len(n) / 4)
cat("weeks", range(weeks), "n =", n, "blocks =", max(blk), "\n")
res <- list()
for (D in c(5, 7, 10)) {            # starting mean infectious period (days)
  rem <- 7 / D                      # per-week removal rate
  st <- list(lambda = rep(1.2 * rem, max(blk)), gamma = 0.98 * rem, alpha = 0.02 * rem,
             i0 = 5e-6, p = 0.2, tau2 = 20, tau2_d = 2)
  t0 <- proc.time()[3]
  f <- tryCatch(sird_ga_fit_tv(cases, deaths, N, blk, "both", st, Delta = 1, nsub = 20L),
                error = function(e) { message(e$message); NULL })
  if (is.null(f)) next
  cat(sprintf("start D=%d: ll=%.2f  1/(g+a)=%.2f days  alpha/(g+a)=%.4f  p=%.3f  i0N=%.0f tau2=%.3g tau2d=%.3g (%.0fs)\n",
              D, f$loglik, 7 / (f$theta["gamma"] + f$theta["alpha"]),
              f$theta["alpha"] / (f$theta["gamma"] + f$theta["alpha"]), f$theta["p"],
              f$theta["i0"] * N, f$theta["tau2"], f$theta["tau2_d"], proc.time()[3] - t0))
  cat("  R_t by block:", round(tapply(f$R_t, blk, mean), 2), "\n")
  res[[as.character(D)]] <- f
}
saveRDS(list(res = res, cases = cases, deaths = deaths, weeks = weeks, blk = blk, N = N),
        "results/paraiba_tv.rds")
