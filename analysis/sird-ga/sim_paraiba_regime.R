# Monte Carlo study in the regime of the Paraíba application: 40 weekly observations, 10 blocks of
# 4 weeks for lambda(t), N = 4,060,000, truth = best fit in results/paraiba_tv.rds.
# Scenarios:  A  truth as fitted (p = 0.066)                                   default B = 200
#             B  same, but p = 0.15                                            default B = 100
#             C  deaths reported with q = 0.8, fitted with q = 1               default B = 50
# Simulation: chain-binomial SIRD with step 0.1 day (as paraiba_bootstrap.R); reported cases
# Bin(new infections, p) + N(0, tau2 * new infections), deaths Bin(new deaths, q) +
# N(0, tau2_d * new deaths), rounded and truncated at 0.
# Fit: sird_ga_fit_tv from the three generic starts of paraiba_tv.R (never from the truth), best
# kept; all three logliks/estimates and convergence codes recorded. Intervals for the best fit:
# Wald for the unconstrained parameters, and for the Table 3 quantities the 4000-draw intervals of
# paraiba_report.R (regime_helpers.R::wald_derived), with Hessian positive-definiteness.
# Truths stored: unconstrained parameters, derived quantities on the deterministic path, attack
# rate and block R_t from the simulated path, and for C the pseudo-truth of a q = 1 fit
# (IFR scaled by q, tau2_d -> (q(1-q) + tau2_d)/q, other quantities unchanged).
# Output: results/regime_<A|B|C>_<first rep>.rds chunks (skipped if present).
# Usage: Rscript sim_paraiba_regime.R <A|B|C> [Brep] [cores] [chunk] [first] [last]
source("regime_helpers.R")
args <- commandArgs(TRUE)
sc <- args[1]
Brep <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else c(A = 200L, B = 100L, C = 50L)[[sc]]
cores <- if (length(args) >= 3) as.integer(args[3]) else 3L
chunk <- if (length(args) >= 4) as.integer(args[4]) else 6L
first <- if (length(args) >= 5) as.integer(args[5]) else 1L
last <- if (length(args) >= 6) as.integer(args[6]) else Brep
outdir <- Sys.getenv("REGIME_OUTDIR", "results")

r <- readRDS("results/paraiba_tv.rds")
f0 <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
N <- r$N; blk <- r$blk; n <- length(r$cases); B <- max(blk)
th <- as.list(f0$theta); lam <- f0$lambda
p_true <- if (sc == "B") 0.15 else th$p
q_true <- if (sc == "C") 0.8 else 1
p_max <- 0.999
sc_names <- c("gamma", "alpha", "i0", "p", "tau2", "tau2_d")
g_a <- th$gamma + th$alpha

# truth on the unconstrained scale (the fitted model has q = 1)
alpha_t <- th$alpha; gamma_t <- th$gamma; tau2d_t <- th$tau2_d
if (sc == "C") {                               # pseudo-truth of the q = 1 fit
  ifr_ps <- q_true * th$alpha / g_a
  alpha_t <- ifr_ps * g_a; gamma_t <- g_a - alpha_t
  tau2d_t <- (q_true * (1 - q_true) + th$tau2_d) / q_true
}
par_true <- c(log(lam), log(gamma_t), log(alpha_t), qlogis(th$i0), qlogis(p_true / p_max),
              log(th$tau2), log(tau2d_t))
par_nominal <- c(log(lam), log(th$gamma), log(th$alpha), qlogis(th$i0), qlogis(p_true / p_max),
                 log(th$tau2), log(th$tau2_d))
names(par_true) <- names(par_nominal) <- c(paste0("log_lambda", seq_len(B)), "log_gamma", "log_alpha",
                                           "logit_i0", "logit_p", "log_tau2", "log_tau2_d")
der_true <- derived_from_par(par_true, sc_names, r$cases, r$deaths, N, blk)
der_nominal <- derived_from_par(par_nominal, sc_names, r$cases, r$deaths, N, blk)

simulate_once <- function() {
  h <- 1 / 70                                  # 0.1 day, in weeks
  I <- round(th$i0 * N); S <- N - I
  inf <- dth <- Swk <- numeric(n)
  pr <- 1 - exp(-(th$gamma + th$alpha) * h); pd <- th$alpha / (th$gamma + th$alpha)
  for (k in seq_len(n)) {
    for (j in 1:70) {
      ni <- rbinom(1, S, 1 - exp(-lam[blk[k]] * I / N * h))
      nr <- rbinom(1, I, pr); nd <- rbinom(1, nr, pd)
      S <- S - ni; I <- I + ni - nr
      inf[k] <- inf[k] + ni; dth[k] <- dth[k] + nd
    }
    Swk[k] <- S
  }
  obs <- function(x, pp, t2) pmax(round(rbinom(n, x, pp) + rnorm(n, 0, sqrt(t2 * x))), 0)
  list(cases = obs(inf, p_true, th$tau2), deaths = obs(dth, q_true, th$tau2_d),
       infections = inf, true_deaths = dth, S = Swk)
}

one_rep <- function(b) {
  set.seed(20261121L + 10000L * match(sc, c("A", "B", "C")) + b)
  sim <- simulate_once()
  t0 <- proc.time()[3]
  fits <- lapply(generic_starts_tv(B), function(st)
    tryCatch(sird_ga_fit_tv(sim$cases, sim$deaths, N, blk, "both", st, Delta = 1, nsub = 20L),
             error = function(e) NULL))
  ll <- sapply(fits, function(f) if (is.null(f)) NA_real_ else f$loglik)
  starts <- do.call(rbind, lapply(seq_along(fits), function(s) {
    f <- fits[[s]]
    if (is.null(f)) return(data.frame(start = s, loglik = NA))
    g <- f$theta[["gamma"]] + f$theta[["alpha"]]
    data.frame(start = s, loglik = f$loglik, conv = f$convergence, p = f$theta[["p"]],
               ifr = f$theta[["alpha"]] / g, inf_period_days = 7 / g,
               attack = 1 - f$kf$x[nrow(f$kf$x), 1], t(setNames(f$par, names(par_true))))
  }))
  if (all(is.na(ll))) return(list(b = b, ok = FALSE, starts = starts))
  best <- fits[[which.max(ll)]]
  w <- wald_derived(best, sim$cases, sim$deaths, N, blk, ndraw = 4000, seed = b)
  R_sim <- tapply(lam[blk] * sim$S / N / g_a, blk, mean)
  list(b = b, ok = TRUE, starts = starts, best_start = which.max(ll), conv = best$convergence,
       pd = w$pd, min_eig = w$min_eig, par = setNames(best$par, names(par_true)), se_par = w$se_par,
       est = w$est, lower = w$lower, upper = w$upper, sd = w$sd,
       attack_sim = sum(sim$infections) / N, R_block_sim = setNames(R_sim, paste0("R_block", seq_len(B))),
       total_cases = sum(sim$cases), total_deaths = sum(sim$deaths),
       secs = unname(proc.time()[3] - t0))
}

cat(format(Sys.time()), "scenario", sc, "reps", first, "-", last, "cores", cores, "\n")
meta_file <- file.path(outdir, sprintf("regime_%s_truth.rds", sc))
if (!file.exists(meta_file))
  saveRDS(list(scenario = sc, p = p_true, q = q_true, par_true = par_true, par_nominal = par_nominal,
               derived_true = der_true, derived_nominal = der_nominal, N = N, blk = blk), meta_file)
for (c0 in seq(first, last, by = chunk)) {
  out_file <- file.path(outdir, sprintf("regime_%s_%03d.rds", sc, c0))
  if (file.exists(out_file)) next
  ids <- c0:min(c0 + chunk - 1, last)
  res <- parallel::mclapply(ids, one_rep, mc.cores = cores)
  saveRDS(res, out_file)
  cat(format(Sys.time()), "scenario", sc, "chunk", c0, "done; secs per rep:",
      paste(round(sapply(res, function(z) if (is.list(z) && !is.null(z$secs)) z$secs else NA)), collapse = " "), "\n")
}
cat("REGIMEDONE", sc, "\n")
