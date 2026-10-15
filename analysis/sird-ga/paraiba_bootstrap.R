# Paraíba: parametric bootstrap of the time-varying-transmission fit (weekly data, 4-week blocks).
# Each replicate simulates the SIRD jump process from the fitted parameters by a chain-binomial
# scheme with step 0.1 day, draws reported cases and deaths from the observation model (binomial
# thinning plus Gaussian error), refits the model from the fitted values, and records the derived
# quantities of Table 3. Replicates are saved in chunks; finished chunks are skipped on rerun.
# Usage: Rscript paraiba_bootstrap.R <B> [cores] [chunk]
suppressMessages(library(epicurves))
args <- commandArgs(TRUE)
Btot <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
chunk <- if (length(args) >= 3) as.integer(args[3]) else 20L

r <- readRDS("results/paraiba_tv.rds")
f <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
N <- r$N; blk <- r$blk; n <- length(r$cases); B <- max(blk)
th <- as.list(f$theta); lam <- f$lambda
start <- list(lambda = lam, gamma = th$gamma, alpha = th$alpha, i0 = th$i0, p = th$p,
              tau2 = th$tau2, tau2_d = th$tau2_d)

derived <- function(fit) {
  s <- fit$kf$x[-1, 1]
  g <- fit$theta[["gamma"]] + fit$theta[["alpha"]]
  c(inf_period_days = 7 / g, ifr = fit$theta[["alpha"]] / g, p = fit$theta[["p"]],
    i0N = fit$theta[["i0"]] * N, tau2 = fit$theta[["tau2"]], tau2_d = fit$theta[["tau2_d"]],
    attack = 1 - s[length(s)],
    setNames(tapply(fit$R_t, fit$blocks, mean), paste0("R_block", seq_len(B))),
    loglik = fit$loglik)
}

simulate_once <- function() {
  h <- 1 / 70                                  # 0.1 day, in weeks
  I <- round(th$i0 * N); S <- N - I
  inf <- dth <- numeric(n)
  pr <- 1 - exp(-(th$gamma + th$alpha) * h); pd <- th$alpha / (th$gamma + th$alpha)
  for (k in seq_len(n)) {
    for (j in 1:70) {
      ni <- rbinom(1, S, 1 - exp(-lam[blk[k]] * I / N * h))
      nr <- rbinom(1, I, pr); nd <- rbinom(1, nr, pd)
      S <- S - ni; I <- I + ni - nr
      inf[k] <- inf[k] + ni; dth[k] <- dth[k] + nd
    }
  }
  obs <- function(x, p, t2) pmax(round(rbinom(n, x, p) + rnorm(n, 0, sqrt(t2 * x))), 0)
  list(cases = obs(inf, th$p, th$tau2), deaths = obs(dth, 1, th$tau2_d), infections = inf)
}

one_rep <- function(b) {
  set.seed(20261015 + b)
  sim <- simulate_once()
  t0 <- proc.time()[3]
  fit <- tryCatch(sird_ga_fit_tv(sim$cases, sim$deaths, N, blk, "both", start, Delta = 1, nsub = 20L),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(b = b, ok = 0))
  c(b = b, ok = 1, derived(fit), conv = fit$convergence, secs = unname(proc.time()[3] - t0),
    attack_true = sum(sim$infections) / N)
}

cat(format(Sys.time()), "B =", Btot, "cores =", cores, "\n")
for (c0 in seq(1, Btot, by = chunk)) {
  out_file <- sprintf("results/paraiba_boot_%03d.rds", c0)
  if (file.exists(out_file)) next
  ids <- c0:min(c0 + chunk - 1, Btot)
  res <- parallel::mclapply(ids, one_rep, mc.cores = cores)
  saveRDS(res, out_file)
  cat(format(Sys.time()), "chunk", c0, "done\n")
}
cat("BOOTDONE\n")
