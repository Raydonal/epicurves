# Optimiser failure or competing maximum? For every replicate of sim_paraiba_regime.R, regenerate
# the data from its seed, check that the stored best log-likelihood is reproduced, and evaluate the
# log-likelihood at the true parameters (for C, the q = 1 pseudo-truth). If the best of the three
# generic starts is below the value at the truth, the optimiser missed a better point; if it is
# above, the degenerate solution is a maximum the data prefer.
# With refit = 1, degenerate replicates (attack rate < 5%) and any whose best fit is more than 0.5
# below the truth are refitted from the
# truth (the start used by the parametric bootstrap of paraiba_bootstrap.R).
# Output: results/regime_truth_check.rds; summary to stdout (-> results/regime_truth_check.log).
# Usage: Rscript regime_truth_check.R [refit 0|1] [cores]
source("regime_helpers.R")
args <- commandArgs(TRUE)
refit <- length(args) >= 1 && args[1] == "1"
cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
options(width = 200)

r <- readRDS("results/paraiba_tv.rds")
f0 <- r$res[[which.max(sapply(r$res, `[[`, "loglik"))]]
N <- r$N; blk <- r$blk; n <- length(r$cases); B <- max(blk)
th <- as.list(f0$theta); lam <- f0$lambda
sc_names <- c("gamma", "alpha", "i0", "p", "tau2", "tau2_d")
p_max <- 0.999

# simulate_once() of sim_paraiba_regime.R, with the scenario's p and q passed in
simulate_once <- function(p_true, q_true) {
  h <- 1 / 70
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
  list(cases = obs(inf, p_true, th$tau2), deaths = obs(dth, q_true, th$tau2_d))
}

ll_par <- function(v, cases, deaths) {
  t <- list(q = 1, r0 = 0, lambda = exp(v[1]))
  bw <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = plogis,
             p = function(z) p_max * plogis(z))
  for (j in seq_along(sc_names)) t[[sc_names[j]]] <- bw[[sc_names[j]]](unname(v[B + j]))
  t$s0 <- 1 - t$i0 - t$r0
  sird_ga_loglik_inc(t, cases, deaths, N, "both", 1, 20L, full = TRUE,
                     lambda_t = exp(v[seq_len(B)])[blk])$loglik
}

start_from_par <- function(v) {
  bw <- list(gamma = exp, alpha = exp, tau2 = exp, tau2_d = exp, i0 = plogis,
             p = function(z) p_max * plogis(z))
  s <- list(lambda = exp(v[seq_len(B)]))
  for (j in seq_along(sc_names)) s[[sc_names[j]]] <- bw[[sc_names[j]]](unname(v[B + j]))
  s
}

rows <- list()
for (sc in c("A", "B", "C")) {
  tr <- readRDS(sprintf("results/regime_%s_truth.rds", sc))
  files <- list.files("results", sprintf("^regime_%s_[0-9]+\\.rds$", sc), full.names = TRUE)
  reps <- Filter(function(z) is.list(z) && isTRUE(z$ok), unlist(lapply(files, readRDS), recursive = FALSE))
  one <- function(z) {
    set.seed(20261121L + 10000L * match(sc, c("A", "B", "C")) + z$b)
    sim <- simulate_once(tr$p, tr$q)
    best_ll <- max(z$starts$loglik, na.rm = TRUE)
    data.frame(scenario = sc, b = z$b, same_data = sum(sim$cases) == z$total_cases &&
                 sum(sim$deaths) == z$total_deaths,
               ll_best = best_ll, ll_best_recomputed = ll_par(z$par, sim$cases, sim$deaths),
               ll_truth = ll_par(tr$par_true, sim$cases, sim$deaths),
               p_best = z$est[["p"]], attack_best = z$est[["attack"]], conv = z$conv)
  }
  rows[[sc]] <- do.call(rbind, parallel::mclapply(reps, one, mc.cores = cores))
}
d <- do.call(rbind, rows)
d$gap <- d$ll_best - d$ll_truth                 # > 0: best fit beats the truth
d$degenerate <- d$attack_best < 0.05

cat("Data regenerated identically:", sum(d$same_data), "of", nrow(d),
    "; max |stored - recomputed best loglik| =", signif(max(abs(d$ll_best - d$ll_best_recomputed)), 3), "\n\n")
summ <- do.call(rbind, lapply(split(d, d$scenario), function(x) data.frame(
  scenario = x$scenario[1], n = nrow(x),
  best_below_truth_gt_0.5 = sum(x$gap < -0.5), best_above_truth_gt_0.5 = sum(x$gap > 0.5),
  median_gap = median(x$gap), degenerate = sum(x$degenerate),
  degenerate_below_truth = sum(x$degenerate & x$gap < -0.5),
  degenerate_above_truth = sum(x$degenerate & x$gap > 0.5))))
cat("Best of three generic starts against the true parameters (gap = loglik best - loglik truth)\n")
print(format(summ, digits = 3), row.names = FALSE)
cat("\nDegenerate best fits (attack rate < 5%):\n")
print(format(d[d$degenerate, c("scenario", "b", "ll_best", "ll_truth", "gap", "p_best", "attack_best", "conv")],
             digits = 4), row.names = FALSE)

if (refit) {
  todo <- d[d$degenerate | d$gap < -0.5, ]
  cat("\nRefitting", nrow(todo), "degenerate replicates from the truth\n")
  rf <- parallel::mclapply(seq_len(nrow(todo)), function(i) {
    sc <- todo$scenario[i]; b <- todo$b[i]
    tr <- readRDS(sprintf("results/regime_%s_truth.rds", sc))
    set.seed(20261121L + 10000L * match(sc, c("A", "B", "C")) + b)
    sim <- simulate_once(tr$p, tr$q)
    f <- tryCatch(sird_ga_fit_tv(sim$cases, sim$deaths, N, blk, "both", start_from_par(tr$par_true),
                                 Delta = 1, nsub = 20L), error = function(e) NULL)
    if (is.null(f)) return(data.frame(scenario = sc, b = b, ll_refit = NA, p_refit = NA, attack_refit = NA))
    data.frame(scenario = sc, b = b, ll_refit = f$loglik, p_refit = f$theta[["p"]],
               attack_refit = 1 - f$kf$x[nrow(f$kf$x), 1])
  }, mc.cores = cores)
  rf <- do.call(rbind, rf)
  d <- merge(d, rf, all.x = TRUE)
  x <- d[!is.na(d$ll_refit), ]
  cat("Refit from truth minus best generic fit: median", signif(median(x$ll_refit - x$ll_best), 3),
      "; refit better by > 0.5 in", sum(x$ll_refit - x$ll_best > 0.5), "of", nrow(x), "\n")
  print(format(x[, c("scenario", "b", "ll_best", "ll_truth", "ll_refit", "p_best", "p_refit",
                     "attack_best", "attack_refit")], digits = 4), row.names = FALSE)
}
saveRDS(d, "results/regime_truth_check.rds")
cat("TRUTHCHECKDONE\n")
