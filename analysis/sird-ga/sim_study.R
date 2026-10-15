# Simulation study for the SIRD-GA article (joint daily cases + deaths).
# Usage: Rscript sim_study.R <scenario> <N> [J] [cores]
# Writes results/sim_s<scenario>_N<N>.rds, one row per accepted replicate.

suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
sc <- as.integer(args[1]); N <- as.integer(args[2])
J <- if (length(args) >= 3) as.integer(args[3]) else 100L
cores <- if (length(args) >= 4) as.integer(args[4]) else 4L

scen <- list(
  `1` = list(lambda = 0.6, gamma = 0.3,  alpha = 0.1, p = 0.5, tau2 = 0.16,   Tmax = 100),
  `2` = list(lambda = 0.6, gamma = 0.3,  alpha = 0.1, p = 0.7, tau2 = 0.16,   Tmax = 100),
  `3` = list(lambda = 0.6, gamma = 0.3,  alpha = 0.1, p = 0.9, tau2 = 0.16,   Tmax = 100),
  `4` = list(lambda = 0.3, gamma = 0.15, alpha = 0.1, p = 0.7, tau2 = 0.0625, Tmax = 300),
  `5` = list(lambda = 0.7, gamma = 0.2,  alpha = 0.2, p = 0.8, tau2 = 0.16,   Tmax = 100)
)[[as.character(sc)]]
tau2_d <- 0.1
I0 <- 10L
truth <- c(lambda = scen$lambda, gamma = scen$gamma, alpha = scen$alpha,
           i0 = I0 / N, p = scen$p, tau2 = scen$tau2, tau2_d = tau2_d)

one_rep <- function(j) {
  set.seed(100000L * sc + 10L * j + N %/% 1000L)
  tries <- 0L
  repeat {
    tries <- tries + 1L
    si <- sird_ga_simulate_inc(N, N - I0, I0, 0, scen$lambda, scen$gamma, scen$alpha,
                               p = scen$p, tau2 = scen$tau2, q = 1, tau2_d = tau2_d,
                               Tmax = scen$Tmax)
    if (sum(si$new_inf) >= 0.05 * N) break   # discard minor outbreaks
  }
  # one start: truth perturbed multiplicatively by exp(N(0, 0.25^2))
  st <- as.list(truth * exp(stats::rnorm(length(truth), 0, 0.25)))
  st$p <- min(st$p, 0.95)
  t0 <- proc.time()[3]
  f <- tryCatch(sird_ga_fit_inc(si$cases, si$deaths, N, use = "both", starts = list(st),
                                nsub = 5L, p_max = 0.999),
                error = function(e) NULL)
  el <- proc.time()[3] - t0
  if (is.null(f)) return(NULL)
  se <- rep(NA_real_, length(f$free))
  if (!is.null(f$hessian)) {
    V <- tryCatch(solve(f$hessian), error = function(e) NULL)
    if (!is.null(V) && all(diag(V) > 0)) se <- sqrt(diag(V))
  }
  names(se) <- paste0("se_", f$free)
  data.frame(rep = j, tries = tries, loglik = f$loglik, secs = el,
             t(f$theta[c("lambda", "gamma", "alpha", "i0", "p", "tau2", "tau2_d")]),
             t(setNames(f$par, paste0("u_", f$free))), t(se), check.names = FALSE)
}

res <- parallel::mclapply(seq_len(J), one_rep, mc.cores = cores)
res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
attr(res, "truth") <- truth
attr(res, "scenario") <- scen
out <- sprintf("results/sim_s%d_N%d.rds", sc, N)
saveRDS(res, out)
cat("done", out, nrow(res), "reps\n")
