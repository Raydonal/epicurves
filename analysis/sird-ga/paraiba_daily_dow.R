# Daily instead of weekly data for the Paraíba fit, with and without day-of-week reporting weights.
# Same period as paraiba_tv.R (epidemiological weeks 13-52 of 2020, 280 days), ten blocks of 28 days
# for lambda(t), Delta = 1 day (rates per day), 5 RK4 steps per day.
#   M0: daily counts, no day-of-week effect.
#   M1: daily counts, seven reporting weights with mean one for cases and seven for deaths
#       (sird_ga_fit_tv(dow =)), which multiply the mean and variance of that day's counts.
# Starts: the three generic starts of paraiba_tv.R on the daily scale (mean infectious period 5, 7,
# 10 days) and a warm start from the weekly fit (rates divided by 7); M1 also from the best M0 fit
# with empirical weekday ratios. Each fit is saved to results/paraiba_daily_<model>_<start>.rds
# (skipped if present); the summary goes to stdout (-> results/paraiba_daily_dow.log).
# Usage: Rscript paraiba_daily_dow.R [cores]
suppressMessages(library(epicurves))
args <- commandArgs(TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 3L
options(width = 200)

data(paraiba_covid19, package = "epicurves"); data(paraiba_population, package = "epicurves")
d <- paraiba_covid19
d <- d[d$epi_week >= 13 & d$epi_week <= 52, ]
stopifnot(nrow(d) == 280, all(table(d$epi_week) == 7))
cases <- d$new_cases; deaths <- d$new_deaths; N <- paraiba_population
dow <- as.integer(format(d$date, "%u"))                 # 1 = Monday, ..., 7 = Sunday
blk <- rep(1:10, each = 28); B <- 10L
nsub <- 5L

wk <- readRDS("results/paraiba_tv.rds")
fw <- wk$res[[which.max(sapply(wk$res, `[[`, "loglik"))]]

generic <- lapply(c(5, 7, 10), function(D) {
  rem <- 1 / D
  list(lambda = rep(1.2 * rem, B), gamma = 0.98 * rem, alpha = 0.02 * rem,
       i0 = 5e-6, p = 0.2, tau2 = 20, tau2_d = 2)
})
warm <- list(lambda = fw$lambda / 7, gamma = fw$theta[["gamma"]] / 7, alpha = fw$theta[["alpha"]] / 7,
             i0 = fw$theta[["i0"]], p = fw$theta[["p"]], tau2 = fw$theta[["tau2"]],
             tau2_d = fw$theta[["tau2_d"]])
starts <- c(setNames(generic, c("gen5", "gen7", "gen10")), list(weekly = warm))

# empirical weekday ratios (from week 16 on, when counts are large), mean one
late <- d$epi_week >= 16
emp <- function(x) { m <- tapply(x[late], dow[late], mean); m / mean(m) }

fit_one <- function(model, sname, st) {
  out <- sprintf("results/paraiba_daily_%s_%s.rds", model, sname)
  if (file.exists(out)) return(readRDS(out))
  t0 <- proc.time()[3]
  f <- tryCatch(sird_ga_fit_tv(cases, deaths, N, blk, "both", st, Delta = 1, nsub = nsub,
                               dow = if (model == "M1") dow else NULL),
                error = function(e) structure(list(error = conditionMessage(e)), class = "fit_error"))
  f$secs <- unname(proc.time()[3] - t0); f$model <- model; f$start_name <- sname
  saveRDS(f, out)
  f
}

jobs <- c(lapply(names(starts), function(s) list("M0", s, starts[[s]])),
          lapply(names(starts), function(s) list("M1", s, starts[[s]])))
res <- parallel::mclapply(jobs, function(j) fit_one(j[[1]], j[[2]], j[[3]]),
                          mc.cores = cores, mc.preschedule = FALSE)
ok <- function(f) is.list(f) && is.null(f$error) && is.finite(f$loglik)
best_of <- function(m) { r <- Filter(function(f) ok(f) && f$model == m, res); r[[which.max(sapply(r, `[[`, "loglik"))]] }
m0 <- best_of("M0")
st <- c(list(lambda = m0$lambda), as.list(m0$theta[c("gamma", "alpha", "i0", "p", "tau2", "tau2_d")]),
        list(dow_cases = as.numeric(emp(cases)), dow_deaths = as.numeric(emp(deaths))))
res <- c(res, list(fit_one("M1", "fromM0", st)))
m1 <- best_of("M1")

cat("Daily Paraíba fits, weeks 13-52 of 2020 (280 days), 10 blocks of 28 days\n\n")
tab <- do.call(rbind, lapply(Filter(ok, res), function(f) data.frame(model = f$model, start = f$start_name,
  loglik = f$loglik, conv = f$convergence, p = f$theta[["p"]], minutes = round(f$secs / 60, 1))))
print(tab[order(tab$model, -tab$loglik), ], row.names = FALSE)
bad <- Filter(function(f) !ok(f), res)
if (length(bad)) cat("failed:", paste(sapply(bad, function(f) paste(f$model, f$start_name)), collapse = ", "), "\n")

k0 <- length(m0$par); k1 <- length(m1$par); n <- 2 * 280
cat(sprintf("\nM0: loglik %.2f, k = %d, AIC %.1f, BIC %.1f\n", m0$loglik, k0, -2 * m0$loglik + 2 * k0, -2 * m0$loglik + log(n) * k0))
cat(sprintf("M1: loglik %.2f, k = %d, AIC %.1f, BIC %.1f\n", m1$loglik, k1, -2 * m1$loglik + 2 * k1, -2 * m1$loglik + log(n) * k1))
cat(sprintf("Likelihood ratio M1 vs M0: %.1f on %d df\n", 2 * (m1$loglik - m0$loglik), k1 - k0))

cat("\nDay-of-week weights (M1; Monday ... Sunday) and empirical ratios from week 16\n")
print(round(rbind(cases_fit = m1$dow_weights$cases, cases_emp = emp(cases),
                  deaths_fit = m1$dow_weights$deaths, deaths_emp = emp(deaths)), 3))

# estimates and Wald SEs (delta method on the unconstrained scale) next to the weekly fit
summ <- function(f, days_per_unit) {
  g <- f$theta[["gamma"]] + f$theta[["alpha"]]
  s_end <- f$kf$x[nrow(f$kf$x), 1]
  R_blk <- tapply(f$R_t, f$blocks, mean)
  est <- c(inf_period_days = days_per_unit / g, ifr_pct = 100 * f$theta[["alpha"]] / g, p = f$theta[["p"]],
           i0N = f$theta[["i0"]] * N, attack_pct = 100 * (1 - s_end), tau2 = f$theta[["tau2"]],
           tau2_d = f$theta[["tau2_d"]], setNames(R_blk, paste0("R_block", seq_along(R_blk))))
  se_p <- NA
  H <- f$hessian
  if (!is.null(H) && all(is.finite(H))) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    j <- length(f$lambda) + match("p", f$sc_names)
    if (!is.null(V) && V[j, j] > 0) {
      eta <- f$par[j]; pr <- plogis(eta)
      se_p <- 0.999 * pr * (1 - pr) * sqrt(V[j, j])
    }
  }
  list(est = est, se_p = se_p, min_eig = if (is.null(H)) NA else
       min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values))
}
sw <- summ(fw, 7); s0 <- summ(m0, 1); s1 <- summ(m1, 1)
cmp <- round(cbind(weekly = sw$est, daily_M0 = s0$est, daily_M1 = s1$est), 4)
print(cmp)
cat(sprintf("\nWald SE of p: weekly %.4f, M0 %.4f, M1 %.4f; min Hessian eigenvalue: %.3g, %.3g, %.3g\n",
            sw$se_p, s0$se_p, s1$se_p, sw$min_eig, s0$min_eig, s1$min_eig))

# standardised one-step-ahead innovations
innov <- function(f) list(cases = (cases - f$kf$fitted_cases) / f$kf$sd_cases,
                          deaths = (deaths - f$kf$fitted_deaths) / f$kf$sd_deaths)
cat("\nInnovations: autocorrelation at lags 1 and 7, Ljung-Box p-values at lags 7 and 14\n")
for (nm in c("M0", "M1")) {
  e <- innov(if (nm == "M0") m0 else m1)
  for (s in c("cases", "deaths")) {
    a <- acf(e[[s]], lag.max = 14, plot = FALSE)$acf
    cat(sprintf("  %s %-6s acf1 %.2f acf7 %.2f  LB7 %.3f LB14 %.3f  sd %.2f\n", nm, s, a[2], a[8],
                Box.test(e[[s]], 7, "Ljung-Box")$p.value, Box.test(e[[s]], 14, "Ljung-Box")$p.value, sd(e[[s]])))
  }
}
saveRDS(list(m0 = m0, m1 = m1, table = tab, compare = cmp, dow = dow, blk = blk, cases = cases,
             deaths = deaths, N = N, dates = d$date), "results/paraiba_daily_dow.rds")
cat("DAILYDONE\n")
