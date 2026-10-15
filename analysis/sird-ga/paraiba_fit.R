# Paraíba first wave: SIRD-GA fitted to reported cases, deaths, and both.
# Usage: Rscript paraiba_fit.R <daily|weekly> <s0fixed|s0free> [q]
suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
agg <- args[1]; s0mode <- args[2]
q <- if (length(args) >= 3) as.numeric(args[3]) else 1
data(paraiba_covid19); data(paraiba_population)
N <- paraiba_population
d <- paraiba_covid19
if (agg == "daily") {
  cases <- d$new_cases; deaths <- d$new_deaths; Delta <- 1
} else {
  # complete epidemiological weeks only
  wk <- d$epi_week
  keep <- wk %in% names(which(table(wk) == 7))
  cases <- as.numeric(tapply(d$new_cases[keep], wk[keep], sum))
  deaths <- as.numeric(tapply(d$new_deaths[keep], wk[keep], sum))
  Delta <- 7
}

periods <- c(3, 5, 7, 10)
mk_starts <- function(R0g = 1.3, s0 = FALSE) lapply(periods, function(D) {
  rem <- 1 / D
  st <- list(lambda = R0g * rem, gamma = 0.97 * rem, alpha = 0.03 * rem,
             i0 = 2e-5, p = 0.3, tau2 = 50, tau2_d = 1)
  if (s0) st$s0 <- 0.2
  st
})

free_base <- c("lambda", "gamma", "alpha", "i0", "p", "tau2", "tau2_d")
fits <- list()
for (use in c("cases", "deaths", "both")) {
  free <- switch(use,
    cases = c("lambda", "gamma", "alpha", "i0", "p", "tau2"),
    deaths = c("lambda", "gamma", "alpha", "i0", "tau2_d"),
    both = free_base)
  if (s0mode == "s0free") free <- c(free, "s0")
  fixed <- list(r0 = 0, q = q, p = 1, tau2 = 1, tau2_d = 1)
  if (s0mode == "s0free") fixed$s0 <- NULL
  st <- mk_starts(s0 = s0mode == "s0free")
  t0 <- proc.time()[3]
  f <- sird_ga_fit_inc(if (use != "deaths") cases, if (use != "cases") deaths, N,
                       use = use, starts = st, free = free, fixed = fixed,
                       p_max = 0.999, Delta = Delta, nsub = 10L, maxit = 6000)
  f$secs <- proc.time()[3] - t0
  f$full <- sird_ga_loglik_inc(f$theta, if (use != "deaths") cases, if (use != "cases") deaths,
                               N, use, Delta, 10L, full = TRUE)
  fits[[use]] <- f
  cat(sprintf("[%s %s q=%.2f] %-6s ll=%.2f  lambda=%.4f gamma=%.4f alpha=%.4f R0=%.3f s0=%.4f i0N=%.1f p=%.3f tau2=%.3g tau2d=%.3g  (%.0fs)\n",
              agg, s0mode, q, use, f$loglik, f$theta["lambda"], f$theta["gamma"], f$theta["alpha"],
              f$theta["lambda"] / (f$theta["gamma"] + f$theta["alpha"]), f$theta["s0"],
              f$theta["i0"] * N, f$theta["p"], f$theta["tau2"], f$theta["tau2_d"], f$secs))
  print(f$starts[, c("start", "loglik", intersect(c("lambda", "gamma", "alpha"), names(f$starts)))], digits = 5)
}
saveRDS(list(fits = fits, cases = cases, deaths = deaths, Delta = Delta, N = N, q = q),
        sprintf("results/paraiba_%s_%s_q%.2f.rds", agg, s0mode, q))
