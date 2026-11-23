# SIRD-GA simulation study (sim_study.R) refitted from generic starting values instead of the
# perturbed truth. Same scenario definitions, same data generation and seeds as sim_study.R, so
# replicate j has the same simulated data as replicate j of results/sim_s<sc>_N<N>.rds; this is
# checked for every replicate by evaluating the log-likelihood at the original estimate on the
# regenerated data and comparing it with the stored maximum (column ll_check_diff).
# Starts (the same for every scenario and replicate; no true value is used):
#   mean infectious period 1/(gamma + alpha) = 2, 4, 8 days (brackets the plausible range for a
#   respiratory infection on a daily scale), R0 = 1.3 (lambda = 1.3 (gamma + alpha)), fatality
#   ratio alpha/(gamma + alpha) = 0.1, p = 0.4, tau2 = 0.1, tau2_d = 0.5, and
#   i0 = max(1, mean(cases[1:3]) / (p lambda)) / N (early growth: reported daily cases are about
#   p lambda I0), a data-based value. p = 0.5 was avoided because it is the truth of Scenario 1,
#   tau2_d = 0.1 because it is the truth of all scenarios.
# Fit: sird_ga_fit_inc with the three starts (Nelder-Mead from each, BFGS polish of the best),
# nsub = 5, p_max = 0.999, as in sim_study.R. Chunked and resumable.
# Usage: Rscript sim_generic_starts.R <scenario> <N> [M] [cores] [chunk]
#   -> results/simgs_s<sc>_N<N>_<first>.rds chunks, results/sim_generic_s<sc>_N<N>.rds combined.
suppressMessages(library(epicurves))
args <- commandArgs(trailingOnly = TRUE)
sc <- as.integer(args[1]); N <- as.integer(args[2])
J <- if (length(args) >= 3) as.integer(args[3]) else 100L
cores <- if (length(args) >= 4) as.integer(args[4]) else 3L
chunk <- if (length(args) >= 5) as.integer(args[5]) else 12L
outdir <- Sys.getenv("SIMGS_OUTDIR", "results")

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
orig <- readRDS(sprintf("results/sim_s%d_N%d.rds", sc, N))

generic_starts <- function(cases) {
  lapply(c(2, 4, 8), function(D) {
    ga <- 1 / D; lam <- 1.3 * ga; p0 <- 0.4
    list(lambda = lam, gamma = 0.9 * ga, alpha = 0.1 * ga,
         i0 = max(1, mean(cases[1:3]) / (p0 * lam)) / N, p = p0, tau2 = 0.1, tau2_d = 0.5)
  })
}

one_rep <- function(j) {
  set.seed(100000L * sc + 10L * j + N %/% 1000L)        # as sim_study.R
  tries <- 0L
  repeat {
    tries <- tries + 1L
    si <- sird_ga_simulate_inc(N, N - I0, I0, 0, scen$lambda, scen$gamma, scen$alpha,
                               p = scen$p, tau2 = scen$tau2, q = 1, tau2_d = tau2_d,
                               Tmax = scen$Tmax)
    if (sum(si$new_inf) >= 0.05 * N) break
  }
  # same-data check against the original fit of replicate j
  o <- orig[orig$rep == j, , drop = FALSE]
  ll_check <- NA_real_
  if (nrow(o) == 1) {
    tho <- list(lambda = o$lambda, gamma = o$gamma, alpha = o$alpha, i0 = o$i0, s0 = 1 - o$i0, r0 = 0,
                p = o$p, tau2 = o$tau2, q = 1, tau2_d = o$tau2_d)
    ll_check <- sird_ga_loglik_inc(tho, si$cases, si$deaths, N, "both", 1, 5L) - o$loglik
  }
  sts <- generic_starts(si$cases)
  t0 <- proc.time()[3]
  f <- tryCatch(sird_ga_fit_inc(si$cases, si$deaths, N, use = "both", starts = sts,
                                nsub = 5L, p_max = 0.999),
                error = function(e) NULL)
  el <- proc.time()[3] - t0
  if (is.null(f)) return(data.frame(rep = j, tries = tries, ll_check_diff = ll_check, failed = TRUE))
  se <- rep(NA_real_, length(f$free))
  pd <- FALSE
  if (!is.null(f$hessian)) {
    pd <- isTRUE(min(eigen((f$hessian + t(f$hessian)) / 2, TRUE, only.values = TRUE)$values) > 0)
    V <- tryCatch(solve(f$hessian), error = function(e) NULL)
    if (!is.null(V) && all(diag(V) > 0)) se <- sqrt(diag(V))
  }
  names(se) <- paste0("se_", f$free)
  stab <- f$starts
  st_ll <- setNames(rep(NA_real_, 3), paste0("ll_start", 1:3)); st_ll[stab$start] <- stab$loglik
  st_p <- setNames(rep(NA_real_, 3), paste0("p_start", 1:3)); st_p[stab$start] <- stab$p
  data.frame(rep = j, tries = tries, loglik = f$loglik, secs = el, failed = FALSE,
             ll_check_diff = ll_check, hessian_pd = pd,
             t(f$theta[c("lambda", "gamma", "alpha", "i0", "p", "tau2", "tau2_d")]),
             t(setNames(f$par, paste0("u_", f$free))), t(se), t(st_ll), t(st_p),
             i0_start = sts[[1]]$i0, check.names = FALSE)
}

cat(format(Sys.time()), "scenario", sc, "N", N, "M", J, "cores", cores, "\n")
for (c0 in seq(1, J, by = chunk)) {
  out_file <- file.path(outdir, sprintf("simgs_s%d_N%d_%03d.rds", sc, N, c0))
  if (file.exists(out_file)) next
  ids <- c0:min(c0 + chunk - 1, J)
  res <- parallel::mclapply(ids, one_rep, mc.cores = cores)
  saveRDS(res, out_file)
  cat(format(Sys.time()), "chunk", c0, "done\n")
}
files <- list.files(outdir, sprintf("^simgs_s%d_N%d_[0-9]+\\.rds$", sc, N), full.names = TRUE)
res <- unlist(lapply(files, readRDS), recursive = FALSE)
res <- Filter(is.data.frame, res)
all_cols <- unique(unlist(lapply(res, names)))
res <- do.call(rbind, lapply(res, function(d) { d[setdiff(all_cols, names(d))] <- NA; d[all_cols] }))
res <- res[order(res$rep), ]
attr(res, "truth") <- truth
attr(res, "scenario") <- scen
saveRDS(res, file.path(outdir, sprintf("sim_generic_s%d_N%d.rds", sc, N)))
cat("GENERICDONE", sc, N, nrow(res), "reps\n")
